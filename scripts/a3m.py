#!/usr/bin/env python


class A3MFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)


class A3M_Container:
    RESIDUES = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    VALID_MATCH_STATES = set(RESIDUES)
    VALID_INSERTION_STATES = set(RESIDUES.lower())
    VALID_GAP_STATES = set("-.")
    VALID_SS_CONF_STATES = set("0123456789")
    VALID_SS_STATES = set("ECH")
    VALID_DSSP_STATES = set("CHBEGITS-")

    def __init__(self):
        self.header = None
        self.annotations = dict()
        self.consensus = None
        self.sequences = []
        self.nr_match_states = None

    @property
    def number_sequences(self):
        """get the current number of protein sequences"""
        return len(self.sequences)

    def check_and_add_sequence(self, header, sequence):
        try:
            if (not self.check_and_add_annotation(header, sequence) and
                    not self.check_and_add_consensus(header, sequence)):
                self.check_sequence(sequence)
                self.sequences.append((header, sequence))
        except A3MFormatError as e:
            raise e

    def check_and_add_consensus(self, header, sequence):
        header_name = header[1:].split()[0]
        if header_name.endswith("_consensus"):
            if self.consensus:
                raise A3MFormatError("Multiple definitions of consensus!")
            else:
                self.check_sequence(sequence)
                self.consensus = (header, sequence)
                return True
        else:
            return False

    def check_and_add_annotation(self, header, sequence):
        annotation_classes = [
            ("ss_conf", self.check_ss_conf),
            ("ss_pred", self.check_ss_pred),
            ("ss_dssp", self.check_dssp)
        ]

        for (annotation_name, check) in annotation_classes:
            if(header[1:].startswith(annotation_name)):
                if(annotation_name in self.annotations):
                    raise A3MFormatError(
                        "Multiple definitions of {}!".format(annotation_name)
                    )
                elif check(sequence):
                    self.annotations[annotation_name] = sequence
                    return True
        return False

    def check_match_states(self, match_states):
        if not self.nr_match_states:
            self.nr_match_states = match_states

        if match_states == 0:
            raise A3MFormatError("Sequence with zero match states!")
        elif match_states != self.nr_match_states:
            raise A3MFormatError(
                ("Sequence with diverging number "
                 "of match states ({} vs. {})!").format(
                    match_states,
                    self.nr_match_states
                )
            )


    def check_ss_conf(self, sequence):
        count_match_states = sum((c in self.VALID_SS_CONF_STATES
                                 or c in self.VALID_GAP_STATES)
                                 for c in sequence)
        self.check_match_states(count_match_states)

        invalid_states = set(sequence) - self.VALID_SS_CONF_STATES
        invalid_states -= self.VALID_GAP_STATES

        if len(invalid_states):
            raise A3MFormatError(
                ("Undefined character(s) '{}' in predicted "
                 "secondary structure confidence!").format(invalid_states))
        else:
            return True

    def check_ss_pred(self, sequence):
        count_match_states = sum((c in self.VALID_SS_STATES
                                 or c in self.VALID_GAP_STATES)
                                 for c in sequence)
        self.check_match_states(count_match_states)

        invalid_states = set(sequence) - self.VALID_SS_STATES
        invalid_states -= self.VALID_GAP_STATES

        if len(invalid_states):
            raise A3MFormatError(
               ("Undefined character(s) '{}' in predicted "
                "secondary structure!").format(invalid_states))
        else:
            return True

    def check_dssp(self, sequence):
        count_match_states = sum(
                (c in self.VALID_DSSP_STATES) for c in sequence)
        self.check_match_states(count_match_states)

        invalid_states = set(sequence) - self.VALID_DSSP_STATES

        if len(invalid_states):
            raise A3MFormatError(
                ("Undefined character(s) '{}' in "
                 "dssp annotation!").format(invalid_states))
        else:
            return True

    def check_sequence(self, sequence):
        count_match_states = sum((c in self.VALID_MATCH_STATES
                                 or c in self.VALID_GAP_STATES)
                                 for c in sequence)
        self.check_match_states(count_match_states)


        invalid_states = set(sequence) - self.VALID_MATCH_STATES
        invalid_states -= self.VALID_GAP_STATES
        invalid_states -= self.VALID_INSERTION_STATES

        if len(invalid_states):
            raise A3MFormatError(
               ("Undefined character(s) '{}' in "
                "protein sequence!").format(invalid_states))
        else:
            return True

    def get_sub_sequence(self, sequence, limits):
        sub_sequence = []

        for (start, end) in limits:
            start_pos = 0
            pos = -1
            for i in range(len(sequence)):
                if (sequence[i] in self.VALID_MATCH_STATES or
                        sequence[i] in self.VALID_GAP_STATES):
                    pos += 1

                    if pos + 1 == start:
                        start_pos = i
                        break

            end_pos = 0
            pos = -1
            for i in range(len(sequence)):
                if (sequence[i] in self.VALID_MATCH_STATES or
                        sequence[i] in self.VALID_GAP_STATES):
                    pos += 1
                    if pos + 1 == end:
                        end_pos = i
                        break
            sub_sequence.append(sequence[start_pos:end_pos+1])

        return "".join(sub_sequence)

    def __str__(self):
        content = []

        if self.header:
            content.append(self.header)

        if self.consensus:
            content.append(self.consensus[0])
            content.append(self.consensus[1])

        for (header, sequence) in self.sequences:
            content.append(header)
            content.append(sequence)

        return "\n".join(content)

    def split_a3m(self, limits):
        new_a3m = A3M_Container()

        if self.consensus:
            new_consensus_sequence = self.get_sub_sequence(self.consensus[1],
                                                           limits)
            new_a3m.consensus = (self.consensus[0], new_consensus_sequence)

        for (header, sequence) in self.sequences:
            new_sequence = self.get_sub_sequence(sequence, limits)
            new_a3m.sequences.append((header, new_sequence))

        return new_a3m

    def read_a3m(self, fh):
        lines = fh.readlines()
        self.read_a3m_from_lines(lines)
        fh.close()

    def read_a3m_from_lines(self, lines):
        sequence_header = None
        sequence = []

        is_first_line = True

        for line in lines:
            line = line.strip()
            if len(line) == 0:
                continue
            elif line[0] == "#":
                if is_first_line:
                    self.header = line
            elif line[0] == ">":
                if sequence_header:
                    self.check_and_add_sequence(sequence_header,
                                                "".join(sequence))
                    sequence = []
                sequence_header = line.rstrip()
            else:
                sequence.append(line.strip().strip("\x00"))

            is_first_line = False

        if sequence_header:
            self.check_and_add_sequence(sequence_header, "".join(sequence))
