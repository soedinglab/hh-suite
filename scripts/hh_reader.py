#!/usr/bin/env python

"""
Parser for hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
"""


import sys
from collections import namedtuple


__author__ = 'Markus Meier (markus.meier@mpibpc.mpg.de)'
__version__ = '1.0'
__license__ = "GPL-3"


hhr_alignment = namedtuple('hhr_alignment', ['query_id', 'query_length', 'query_neff',
                                             'template_id', 'template_length', 'template_info',
                                             'template_neff', 'query_ali', 'template_ali', 
                                             'start', 'end', 'probability', 'evalue', 'score', 
                                             'aligned_cols', 'identity', 'similarity', 'sum_probs'])


class HHRFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)


def get_sequence_name(header):
    name = header.replace(">", "").split()[0]
    return name


def parse_result(lines):
    results = []

    query_id = None
    query_length = None
    query_neff = None
    query_seq = []
    template_id = None
    template_length = None
    template_seq = []
    template_info = None
    query_start = None
    query_end = None
    template_start = None
    template_end = None
    probability = None
    evalue = None
    score = None
    identity = None
    similarity = None
    template_neff = None
    sum_probs = None
    aligned_cols = None

    skipped_ali_tags = ["ss_dssp", "ss_pred", "Consensus"]

    is_alignment_section = False

    for line in lines:
        if(line.startswith("Query")):
            query_id = line.split()[1]
        elif(line.startswith("Match_columns")):
            query_length = int(line.split()[1])
        elif(line.startswith("Neff")):
            query_neff = float(line.split()[1])
        elif(is_alignment_section and (line.startswith("No") or line.startswith("Done!"))):
            if query_start is not None:
                result = hhr_alignment(query_id, query_length, query_neff,
                                       template_id, template_length, template_info, template_neff,
                                       query_seq, template_seq, (query_start, template_start),
                                       (query_end, template_end), probability, evalue, score,
                                       aligned_cols, identity, similarity, sum_probs)
                results.append(result)
            template_id = None
            template_info = None
            query_seq = []
            template_seq = []

            query_start = None
            query_end = None
            template_start = None
            template_end = None
        elif(line.startswith("Probab")):
            tokens = line.split()
            probability = float(tokens[0].split("=")[1])
            evalue = float(tokens[1].split("=")[1])
            score = float(tokens[2].split("=")[1])
            aligned_cols = int(tokens[3].split("=")[1])
            identity = float(tokens[4].split("=")[1].replace("%", "")) / 100.0
            similarity = float(tokens[5].split("=")[1])
            sum_probs = float(tokens[6].split("=")[1])
            if(len(tokens) > 7):
                template_neff = float(tokens[7].split("=")[1])
            continue
        elif(line.startswith(">")):
            is_alignment_section = True
            template_id = line[1:].split()[0]
            template_info = line
        elif(line.startswith("Q")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of query alignment").format(tokens[2]))

            if query_start is None:
                query_start = token_2
            query_start = min(query_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of query alignment").format(tokens[4]))

            if query_end is None:
                query_end = token_4
            query_end = max(query_end, token_4)
            query_seq.append(tokens[3])
        elif(line.startswith("T")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue
            template_seq.append(tokens[3])

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of template alignment").format(tokens[2]))

            if template_start is None:
                template_start = token_2
            template_start = min(template_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of template alignment").format(tokens[4]))

            if template_end is None:
                template_end = token_4
            template_end = max(template_end, token_4)

            try:
                token_5 = tokens[4].replace("(", "").replace(")", "")
                token_5 = int(token_5)
            except:
                raise HHRFormatError(("Converting failure of template length ({}) "
                                      "in template alignment").format(tokens[5]))
            template_length = token_5


    if(template_id is not None and query_start is not None):
        result = hhr_alignment(query_id, query_length, query_neff,
                               template_id, template_length, template_info, template_neff,
                               query_seq, template_seq, (query_start, template_start),
                               (query_end, template_end), probability, evalue, score,
                               aligned_cols, identity, similarity, sum_probs)
        results.append(result)

    return results


def read_result(input_file):
    with open(input_file) as fh:
        lines = fh.readlines()
        return parse_result(lines)


def main():
    counter = 0
    for result in read_result(sys.argv[1]):
        print("Alignment " + str(counter) + "\t evalue: " + str(result.evalue) +
              "\t probability: " + str(result.probability))

        print(result.query_id + "\t" + str(result.start[0]) + "\t" +
              result.query_ali + "\t" +
              str(result.end[0]))

        print(result.template_id + "\t" + str(result.start[1]) + "\t" +
              result.template_ali + "\t" +
              str(result.end[1]))

        counter += 1


if __name__ == "__main__":
    main()
