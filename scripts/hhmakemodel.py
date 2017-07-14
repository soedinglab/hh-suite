#!/usr/bin/env python

from hh_reader import read_result
from copy import deepcopy
from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
import re, os, sys, tempfile, glob

from operator import itemgetter  # hzhu
from itertools import groupby    # hzhu

EMPTY = '*'
GAP = '-'
DEBUG_MODE = False

class Gap: 
    """ A gap is a continuous stretch of indels.
        It is defined by a opening position and a size/length
    """
    def __init__(self, open_pos, size):
        self.open_pos = open_pos   # gap opening position
        self.size = size           # num of indels in the gap

    def __repr__(self):
        return 'Gap opening pos = %d, size = %d' % (self.open_pos, self.size)
        
class Grid:
    """
    Implementation of 2D grid of cells
    Includes boundary handling
    """
    
    def __init__(self, grid_height, grid_width):
        """
        Initializes grid to be empty, take height and width of grid as parameters
        Indexed by rows (left to right), then by columns (top to bottom)
        """
        
        self._grid_height = grid_height
        self._grid_width = grid_width
        self._cells = [ [ EMPTY for dummy_col in range(self._grid_width) ]
                       for dummy_row in range(self._grid_height)]
                
    def __str__(self):
        """ Return multi-line string represenation for grid """
        
        ans = ''
        for row in range(self._grid_height):
            ans += ''.join(self._cells[row])
            ans += '\n'
        return ans

    def clear(self):
        """ Clears grid to be empty """
        
        self._cells = [[EMPTY for dummy_col in range(self._grid_width)]
                       for dummy_row in range(self._grid_height)]
    
    def get_grid_height(self):
        """ Return the height of the grid """
        
        return self._grid_height

    def get_grid_width(self):
        """ Return the width of the grid """
        
        return self._grid_width
    
    def get_cell(self, row, col):
        return self._cells[row][col]

    def get_seq_start(self, row):
        """ Returns the start position of the sequence """
        
        index = 0
        for pos in self._cells[row]:
            if pos != EMPTY: 
                return index
            index += 1

        return None

    def get_seq_end(self, row):
        """ Returns the end position of the sequence """
        
        index = 0
        for pos in reversed(self._cells[row]):
            if pos != EMPTY: 
                return self.get_grid_width() - index
            index += 1

        return None

    def get_gaps(self, row): 
        """ Return the position of gaps in a row """
        
        gaps = list()

        index = 0
        for pos in self._cells[row]:
            if pos == GAP: 
                gaps.append(index)
            index += 1

        return gaps

    def get_gaps_ref_gapless(self, row):
        """ Return the pos of gaps in a row.
            The opening positions of the gaps are wrt. the gapless seq
        """
        # get all the indels
        indels = self.get_gaps(row)  
        gaps = []
        # combine continuous indels into a gap
        for k,i in groupby( enumerate(indels), lambda x: x[0]-x[1] ):
            g = list(map(itemgetter(1), i))
            gaps.append( Gap(g[0], len(g)) )

        # offset the gap opening positions
        for i in range(1, len(gaps)):
            # offset by total gap number before
            gaps[i].open_pos -= sum([gaps[j].size for j in range(i)])
            
        return gaps  # a list of Gap instances
    
    def get_seq_indeces(self, row):

        seq = list()
        for pos, res in enumerate(self._cells[row]):
            if res != EMPTY and res != GAP:
                seq.append(pos)

        return seq
            
    ## def get_gap_list(self):  # hzhu commented this out. wrote a new version
    ##     """ Returns a list of list of all gap positions in the sequence grid. """
    ##     gap_pos = set()

    ##     for row in range(self.get_grid_height()):
    ##         for gap in self.get_gaps(row):
    ##             gap_pos.add(gap)

    ##     gap_pos = list(sorted(gap_pos))
    ##     boundaries = [ (x + 1) for x, y in zip(gap_pos, gap_pos[1:]) if y - x != 1 ]
 
    ##     gap_list = list()
    ##     prev = 0
 
    ##     for boundary in boundaries:
    ##         sub_list = [ pos for pos in gap_pos[prev:] if pos < boundary ]
    ##         gap_list.append(sub_list)
    ##         prev += len(sub_list)
 
    ##     gap_list.append([ x for x in gap_pos[prev:]])
    
    ##     return gap_list

    def get_gap_list(self):
        """ Returns a list of Gap instances for all rows in the grid
        """
        gap_dict = dict()  # each position should occur as gap at most once
                           # keys are gap openning positions
                           # values are Gap instances
        gap_list = []
        for row in range(self.get_grid_height()):
            gap_pos = []
            gaps = self.get_gaps_ref_gapless(row)
             
            for g in gaps:
                if g.open_pos in gap_dict: # if there is already gaps at this open pos
                    if g.size > gap_dict[g.open_pos].size: # if new gap is bigger
                        gap_dict[g.open_pos] = g  # keep the larger gap as they overlap
                else:
                    gap_dict[g.open_pos] = g
                    
        gap_list = sorted(list(gap_dict.values()), key=lambda x: x.open_pos) # sort according to start position
        return gap_list  # a list of Gap instances
        
    def set_gap(self, row, col):
        """ Set cell with index (row, col) to be a gap """
        
        self._cells[row][col] = GAP

    def set_empty(self, row, col):
        """ Set cell with index (row, col) to be a gap """
        
        self._cells[row][col] = EMPTY
    
    def set_cell(self, row, col, res):
        """ Set cell with index (row, col) to be full """
        
        self._cells[row][col] = res

    def is_empty(self, row, col):
        """ Checks whether cell with index (row, col) is empty """
        
        return self._cells[row][col] == EMPTY

    def is_gap(self, row, col):
        """ Checks whetehr cell with indxex (row, col) is a gap """
        
        return self._cells[row][col] == GAP

    def insert_gaps(self, cols):
        """ Inserts a gaps into a column of the template grid """
        
        for col in cols:
            for row in range(self._grid_height):
                if col >= self.get_seq_start(row) and col < self.get_seq_end(row):
                    self._cells[row].insert(col, GAP)
                else:
                    self._cells[row].insert(col, EMPTY)
            
            self._grid_width += 1

    def insert_gaps_row(self, cols, row):
        """ Intert gaps into cols only for certain row"""
        for col in cols:
            if col >= self.get_seq_start(row) and col < self.get_seq_end(row):
                self._cells[row].insert(col, GAP)
            else:
                self._cells[row].insert(col, EMPTY)
            # NOTE: grid_with should not be changed after every row is updated.    
            #self._grid_width += 1

    def clean_trail_empty(self):
        """ Remove all trailing EMPTY and pad grid to same width"""
        # first find out the max length (exluding trailing EMPTY)
        max_width = 0
        for row in range(self._grid_height):
            for i in range(len(self._cells[row])-1, -1, -1):
                if self._cells[row][i] != EMPTY:
                    break
            if i+1 > max_width:
                max_width = i+1
                
        # delete excessive EMPTY        
        for row in range(self._grid_height):
            del self._cells[row][max_width:]

        # then pad all rows to the same length
        [self._cells[row].append( EMPTY * (max_width-len(self._cells[row])) ) \
                                 for row in range(self._grid_height) if len(self._cells[row]) < max_width]
        self._grid_width = max_width
        return
    
    def remove_gaps(self, keep_width=True): # hzhu add keep_width option
        """ Removes all gaps from the grid. """

        for row in range(self.get_grid_height()):
            not_gap = list()
            for col in range(self.get_grid_width()):
                if not self.is_gap(row, col):
                    not_gap.append(col)

            self._cells[row] = [ self._cells[row][col] for col in not_gap ]

            if keep_width:  # hzhu only pad to original width if desired
                for del_pos in range(self._grid_width - len(not_gap)):
                    self._cells[row].append(EMPTY)

        if not keep_width: # hzhu if width is not kept, make sure width is consistent
            self.clean_trail_empty()
                    
        return
                

class QueryGrid(Grid):

    def __init__(self, grid_height, grid_width):
        Grid.__init__(self, grid_height, grid_width)
    
    def get_query_start(self, row):
        """ Returns the query start position """
        return self.get_seq_start(row) + 1

    def get_query_end(self, row):
        """ Returns the query end postion """
        
        return self.get_seq_end(row) - len(self.get_gaps(row))

    def get_col_residue(self, col):
        """ Tries to find a the query residue in a given column. Used by derive_global_seq() to
        identify the global query sequence """
        
        for row in range(self.get_grid_height()):
            if not self.is_empty(row, col):
                return self._cells[row][col]

        return GAP

class TemplateGrid(Grid):

    def __init__(self, grid_height, grid_width):
        Grid.__init__(self, grid_height, grid_width)

        self._start = list()
        self._end = list()
        self._pdb_code = list()
        self._chain = list()
        self._organism = list()
        self._resolution = list()
        
    def display(self):
        """ Return multi-line string represenation for grid """
        
        ans = ''
        for row in range(self._grid_height):
            ans += '>P1;{p}\nstructure:{p}:{s}:{c}:{e}:{c}::{o}:{r}:\n{a}*\n'.format(
                p = self._pdb_code[row],
                s = add_white_space_end(self.get_template_start(row), 4),
                e = add_white_space_end(self.get_template_end(row), 4),
                c = self._chain[row],
                o = self._organism[row],
                r = self._resolution[row], 
                a = ''.join(self._cells[row]).replace(EMPTY, GAP).replace('#', GAP))

        return ans

    def debug(self, row):
        """ Return multi-line string represenation for grid, for debugging purposes """

        ans = '{p}\nInternal: {s}, {e} Query: {qs}, {qe} Gaps ({g1}): {g2}\n{seq}\n'.format(
            p = self._pdb_code[row],
            s = self.get_seq_start(row),
            e = self.get_seq_end(row),
            qs = self.get_template_start(row),
            qe = self.get_template_end(row),
            g1 = len(self.get_gaps(row)),
            g2 = ', '.join([str(gap) for gap in self.get_gaps(row)]),
            seq = ''.join(self._cells[row]))

        return ans 

    def set_metadata(self, row, start, end, pdb_code, chain, organism, resolution):
        """ Used by create_template_grid() to setup metadata of pir template """
        
        self._start.append(start)
        self._end.append(end)
        self._pdb_code.append(pdb_code)
        self._chain.append(chain)
        self._organism.append(organism)
        self._resolution.append(resolution)

    def set_map(self, row, start, end):
        
        self._start[row] = start
        self._end[row] = end

    def get_template_start(self, row):
        """ Returns the template start position """
        
        return self._start[row]

    def get_template_end(self, row):
        """ Return sthe template end position """
        
        return self._end[row]

    def del_row(self, row):
        """ Removes a complete template entry from the grid """
        
        del self._cells[row]
        del self._start[row]
        del self._end[row]
        del self._pdb_code[row]
        del self._chain[row]
        del self._organism[row]
        del self._resolution[row]
        self._grid_height -= 1

# Helper functions

def add_white_space_end(string, length):
    """ Adds whitespaces to a string until it has the wished length"""
    
    edited_string = str(string)
    
    if len(edited_string) >= length:
        return string
    else:
        while len(edited_string) != length:
            edited_string += ' '
    
    return edited_string

def convert_aa_code(three_letter, convert):
    """
    Assumes a string that contains a three letter aminoacid code and
    returns the corresponding one letter code.
    """    
    
    aa_code = {
        'CYS': 'C', 
        'ASP': 'D', 
        'SER': 'S', 
        'GLN': 'Q', 
        'LYS': 'K',
        'ILE': 'I', 
        'PRO': 'P', 
        'THR': 'T', 
        'PHE': 'F', 
        'ASN': 'N', 
        'GLY': 'G', 
        'HIS': 'H', 
        'LEU': 'L', 
        'ARG': 'R', 
        'TRP': 'W', 
        'ALA': 'A', 
        'VAL': 'V', 
        'GLU': 'E', 
        'TYR': 'Y', 
        'MET': 'M',
      } 

    non_canonical = {
        'MSE': 1,
        'HYP': 2,
        'MLY': 3,
        'SEP': 4,
        'TPO': 5,
        'CSO': 6,
        'PTR': 7,
        'KCX': 8,
        'CME': 9,
        'CSD': 10,
        'CAS': 11,
        'MLE': 12,
        'DAL': 13,
        'CGU': 14,
        'DLE': 15,
        'FME': 16,
        'DVA': 17,
        'OCS': 18,
        'DPR': 19,
        'MVA': 20,
        'TYS': 21,
        'M3L': 22,
        'SMC': 23,
        'ALY': 24,
        'CSX': 25,
        'DCY': 26,
        'NLE': 27,
        'DGL': 28,
        'DSN': 29,
        'CSS': 30,
        'DLY': 31,
        'MLZ': 32,
        'DPN': 33,
        'DAR': 34,
        'PHI': 35,
        'IAS': 36,
        'DAS': 37,
        'HIC': 38,
        'MP8': 39,
        'DTH': 40,
        'DIL': 41,
        'MEN': 42,
        'DTY': 43,
        'CXM': 44,
        'DGN': 45,
        'DTR': 46,
        'SAC': 47,
        'DSG': 48,
        'MME': 49,
        'MAA': 50,
        'YOF': 51,
        'FP9': 52,
        'FVA': 53,
        'MLU': 54,
        'OMY': 55,
        'FGA': 56,
        'MEA': 57,
        'CMH': 58,
        'DHI': 59,
        'SEC': 60,
        'OMZ': 61,
        'SCY': 62,
        'MHO': 63,
        'MED': 64,
        'CAF': 65,
        'NIY': 66,
        'OAS': 67,
        'SCH': 68,
        'MK8': 69,
        'SME': 70,
        'LYZ': 71
    }

    if three_letter in aa_code.keys():
        return aa_code[three_letter]
    elif convert and (three_letter in non_canonical.keys()):
        return non_canonical[three_letter]
    else:
        return '-'


def get_query_name(hhr_file):

    with open(hhr_file) as fh:
        for line in fh:
            if line.startswith('Query'):
                # match the PDB Code
                m = re.search('(\d[A-Z0-9]{3})_(\S)', line) 

                if m: 
                    pdb_code = m.group(1)
                    chain = m.group(2)
                else: 
                    pdb_code = 'UKNP'
                    chain = 'A'
                    # raise ValueError('Input HHR-File Does not seem to be a PDB-Structure')

                break

    return pdb_code, chain

def get_cif_files(folder):
    """ Gets all cif files located in folder. """
    
    return glob(os.path.join(folder, '*.cif'))

def open_cif(cif_file):
    """ Assumes a mmCif file and returns a data block used for subsequent procedures """
    # The "usual" procedure to open a mmCIF with pdbX/mmCIF
    with open(cif_file) as cif_fh:
        data = []
        reader = PdbxReader(cif_fh)
        reader.read(data)

    block = data[0]
    
    return block

def get_pdb_entry_id(block):
    """ Extracts the PDB entry information of a cif file and returns it as a string """

    entry = block.getObj('entry')
    entry_id = entry.getValue('id')
    
    return entry_id


def template_id_to_pdb(template_id):
    """
    Extracts PDB ID and chain name from the provided template id
    """
    # match PDBID without chain (8fab, 1a01)
    m = re.match(r'/^(\d[A-Za-z0-9]{3})$', template_id)
    if m:
        return m.group(1).upper(), 'A'
    
    # PDB CODE with chain Identifier
    m = re.match(r'^(\d[A-Za-z0-9]{3})_(\S)$', template_id)
    if m:
        return m.group(1).upper(), m.group(2).upper()

    # Match DALI ID
    m = re.match(r'^(\d[A-Za-z0-9]{3})([A-Za-z0-9]?)_\d+$', template_id)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    
    # No PDB code and chain identified
    return None, None


def create_template_grid(hhr_data):
    """ Creates a template grid """

    total_seq = len(hhr_data)
    templ_max = max( [ hhr.start[0] + len(to_seq(hhr.template_ali)) for hhr in hhr_data ] ) - 1


    template_grid = TemplateGrid(total_seq, templ_max)

    for row, template in enumerate(hhr_data):
        seq_start = template.start[0] - 1
        templatealignment = to_seq(template.template_ali)
        seq_end = seq_start + len(templatealignment)

        # Load Meta Data
        start = template.start[1]
        end = template.end[1]

        # Get pdb_code and chain identifier of template
        pdb_code, chain =  template_id_to_pdb(template.template_id)

        m = re.search("(\d+.\d+)A", template.template_info) # try to extract resolution of the structure
        
        if m: 
            resolution = m.group(1)
        else: 
            resolution = ""

        m = re.search("\{(.*)\}", template.template_info) # try to extract the organism
        if m: 
            organism = m.group(1).replace(":", " ") # make sure that no colons are in the organism
        else: 
            organism = ""

        template_grid.set_metadata(row, start, end, pdb_code, chain, organism, resolution)

        # Write sequence into the grid
        for pos, col in enumerate(range(seq_start, seq_end)):
            template_grid.set_cell(row, col, templatealignment[pos])

    return template_grid


def to_seq(ali):
    if isinstance(ali, list):
        return ''.join(ali)
    else:
        return ali
    

def create_query_grid(hhr_data):
    """ Creates a Query Grid """
    
    total_seq = len(hhr_data)
    query_max = max( [ hhr.start[0] + len(to_seq(hhr.query_ali)) for hhr in hhr_data ] ) - 1

    query_grid = QueryGrid(total_seq, query_max)

    for row, query in enumerate(hhr_data):

        queryalignment = to_seq(query.query_ali)
        query_start = query.start[0] - 1
        query_end = query_start + len(queryalignment)

        for pos, col in enumerate(range(query_start, query_end)):
            if queryalignment[pos] not in ['Z', 'U', 'O', 'J', 'X', 'B']: # CAUTION

                query_grid.set_cell(row, col, queryalignment[pos])

    return query_grid

def create_gapless_grid(grid):
    """ Returns a gapless grid """

    gapless = deepcopy(grid)
    gapless.remove_gaps(keep_width=False)  # hzhu:  shrink grid

    return gapless

def process_query_grid(query_grid, gapless_grid):
    """ Processes a query grid sucht that it contains all gaps
    """
    gaplist = query_grid.get_gap_list()
    off_set = 0
    
    for g in gaplist:
        gapless_grid.insert_gaps([ p + off_set for p in range(g.open_pos, g.open_pos+g.size) ])
        off_set += g.size
        
    return gapless_grid

def derive_global_seq(processed_query_grid, query_name, query_chain):

    global_seq = list()

    for col in range(processed_query_grid.get_grid_width()):
        global_seq.append(processed_query_grid.get_col_residue(col))

    # this is the query entry
    header = '>P1;{q}\nsequence:{q}:1    :{c}:{l}  :{c}::::\n'.format(
        q = query_name, 
        l = len(global_seq),
        c = query_chain)

    return header + ''.join(global_seq) + '*'

def process_template_grid(query_grid, template_grid):
    """ Insertes Gaps into the template grid
        Only add gaps from **other** query_grids into template grid (NOT gapless)
    """
    gaplist = query_grid.get_gap_list()  # use this to keep the offset
    
    for row in range(template_grid.get_grid_height()):
        # do NOT consider gaps in current query row         
        gaplist_row = query_grid.get_gaps_ref_gapless(row)
        gapdict_row = dict(zip([g.open_pos for g in gaplist_row],
                               [g.size     for g in gaplist_row]))
        off_set = 0
        for g in gaplist:
            # if there is a gap with same opening position in the current row,
            # only consider g if it is larger than the on in the current row
            if g.open_pos in gapdict_row:
                if g.size > gapdict_row[g.open_pos]: 
                    template_grid.insert_gaps_row([ p + off_set for p in range(g.open_pos,
                                                                              g.open_pos+g.size-gapdict_row[g.open_pos]) ], row)
            else:
                template_grid.insert_gaps_row([ p + off_set for p in range(g.open_pos, g.open_pos+g.size) ], row)
                 
            off_set += g.size  # even if the gaps are not inserted, the offset should be adjusted

    template_grid.clean_trail_empty()  # clean the redundant trailing EMPTY char
    
    return template_grid

def compare_with_cifs(template_grid, folder, output_path, convert, threshold):
    """
    Compare the PIR Alignment with Atomsection of a mmCIF file. To make the ATOM-Section of 
    a mmCIF file compatible with MODELLER, each residue has in the ATOM-Section has to match
    corresponding positions in the PIR-Alignment
    """

    # glob the mmCif files from given directory and map the PDB identifier to the path
    cif_files = glob.glob(os.path.join(folder, '*.cif'))
    cif_paths = { path.split('/')[-1].split('.')[0].upper() : path for path in cif_files }
    cif_edits = dict()
    

    # create the path where renumbered cifs are saved to
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    # if the cif does not contain any residue of the por alignment we delete it
    del_row = list()

    for row in range(template_grid.get_grid_height()):
        # get the pdb code and strand id from the current template
        pdb_code = template_grid._pdb_code[row]
        chain = template_grid._chain[row]  # hhr users pdb chain ID


        # load mmCif file accordingly
        if pdb_code in cif_edits.keys():
            block = cif_edits[pdb_code]
        else:
            try:
                block = open_cif(cif_paths[pdb_code])
            except KeyError:
                del_row.append(row)
                print ('! Did not find the mmCIF file for {pdb}. Removing it from the alignment.'.format(
                    pdb = pdb_code))
                continue

        # Create a mapping of the atom site
        atom_site = block.getObj('atom_site')

        ########################################################################
        ## Get the mapping of the residues in the atom section                ##
        ########################################################################

        cif_seq = dict()
        # For the case that we have to rename a chain
        cif_chains = set([]) 

        # Iterate through the atomsection of the cif file
        for atom_row in range(0, atom_site.getRowCount()):
            
            try:
                if atom_site.getValue('label_comp_id', atom_row) == 'HOH':
                    continue
                cif_chain = atom_site.getValue('label_asym_id', atom_row)
                pdb_chain = atom_site.getValue('auth_asym_id', atom_row)   # use PDB chain ID
            except IndexError:
                pass

            cif_chains.add(cif_chain)

            # We do not care about the residues apart from the chain
            #if cif_chain != chain: # hzhu
            if pdb_chain != chain:  # hhr uses PDB chain, not the cif chain! hzhu
                continue     
            # and update the chain id from pdb_chain to cif_chain
            if atom_site.getValue('group_PDB', atom_row).startswith('ATOM'):  # hzhu in case HETATM ruins ch id
                template_grid._chain[row] = cif_chain
                
            # get the residue and the residue number 
            try:
                res_num = int(atom_site.getValue("label_seq_id", atom_row))
            except ValueError:
                continue

            residue = atom_site.getValue('label_comp_id', atom_row)
            residue = convert_aa_code(residue, convert)         


            if res_num not in cif_seq.keys():
                cif_seq[res_num] = residue
            elif res_num in cif_seq.keys() and cif_seq[res_num] == residue:
                continue
            elif res_num in cif_seq.keys() and cif_seq[res_num] != residue:
                cif_seq[res_num] = '-'
                
                if DEBUG_MODE:
                    print ('! {p} {c}: mmCIF contains a residue position that is assigned {cr} to two residues. Removing it.'.format(
                        p = pdb_code,
                        c = chain,
                        cr = res_num))

        ########################################################################
        ## Rename chain if necessary                                          ##
        ########################################################################

        chain_idx = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        if len(template_grid._chain[row]) != 1:
            i = 0
            new_chain = 0
            
            while i < len(chain_idx):
                if chain_idx[i] in cif_chains:
                    if DEBUG_MODE:
                        print ('! {p} {c}: Chain identifier {i} is already taken.'.format(
                            p = pdb_code,
                            c = chain,
                            i = chain_idx[i]))
                    i += 1
                else:
                    new_chain = chain_idx[i]
                    break

            if new_chain == 0:
                if DEBUG_MODE:
                    print ('! {p} {c}: Could not use {p}. The chain identifier {c} is not compatible with MODELLER (2 letters) and could not be renanmed.'.format(
                        p = pdb_code, 
                        c = chain))
                
                del_row.append(row)
                continue

            if new_chain != 0:
                print ('Selected new chain name {c}'.format(c = new_chain))

            #TODO

        ########################################################################
        ## Compare cif positions with the atom positions                      ##
        ########################################################################

        del_pos = list()
        mod_pos = dict()
        mapping = dict()

        for pos_cif, pos_tem in zip(range(template_grid.get_template_start(row), 
            template_grid.get_template_end(row) + 1), template_grid.get_seq_indeces(row)):

            res_tem = template_grid.get_cell(row, pos_tem)

            try:
                res_cif = cif_seq[pos_cif]
            except KeyError:
                res_cif = -1
                

            match = True if res_tem == res_cif else False

            if not match:
                if res_cif == 1 and res_tem == 'M':
                    mod_pos[pos_cif] = 1
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 2 and res_tem == 'P':

                    mod_pos[pos_cif] = 2
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'P')

                elif res_cif == 3 and res_tem == 'K':

                    mod_pos[pos_cif] = 3
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 4 and res_tem == 'S':

                    mod_pos[pos_cif] = 4
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'S')

                elif res_cif == 5 and res_tem == 'T':

                    mod_pos[pos_cif] = 5
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'T')

                elif res_cif == 6 and res_tem == 'C':

                    mod_pos[pos_cif] = 6
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 7 and res_tem == 'Y':

                    mod_pos[pos_cif] = 7
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 8 and res_tem == 'K':

                    mod_pos[pos_cif] = 8
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 9 and res_tem == 'C':

                    mod_pos[pos_cif] = 9
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 10 and res_tem == 'A':

                    mod_pos[pos_cif] = 10
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'A')

                elif res_cif == 11 and res_tem == 'C':

                    mod_pos[pos_cif] = 11
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 12 and res_tem == 'L':

                    mod_pos[pos_cif] = 12
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'L')

                elif res_cif == 13 and res_tem == 'A':

                    mod_pos[pos_cif] = 13
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'A')

                elif res_cif == 14 and res_tem == 'E':

                    mod_pos[pos_cif] = 14
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'E')

                elif res_cif == 15 and res_tem == 'L':

                    mod_pos[pos_cif] = 15
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'L')

                elif res_cif == 16 and res_tem == 'M':

                    mod_pos[pos_cif] = 16
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 17 and res_tem == 'V':

                    mod_pos[pos_cif] = 17
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'V')

                elif res_cif == 18 and res_tem == 'C':

                    mod_pos[pos_cif] = 18
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 19 and res_tem == 'P':

                    mod_pos[pos_cif] = 19
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'P')

                elif res_cif == 20 and res_tem == 'V':

                    mod_pos[pos_cif] = 20
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'V')

                elif res_cif == 21 and res_tem == 'Y':

                    mod_pos[pos_cif] = 21
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 22 and res_tem == 'K':

                    mod_pos[pos_cif] = 22
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 23 and res_tem == 'C':

                    mod_pos[pos_cif] = 23
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 24 and res_tem == 'K':

                    mod_pos[pos_cif] = 24
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 25 and res_tem == 'C':

                    mod_pos[pos_cif] = 25
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 26 and res_tem == 'C':

                    mod_pos[pos_cif] = 26
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 27 and res_tem == 'L':

                    mod_pos[pos_cif] = 27
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'L')

                elif res_cif == 28 and res_tem == 'E':

                    mod_pos[pos_cif] = 28
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'E')

                elif res_cif == 29 and res_tem == 'S':

                    mod_pos[pos_cif] = 29
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'S')

                elif res_cif == 30 and res_tem == 'C':

                    mod_pos[pos_cif] = 30
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 31 and res_tem == 'K':

                    mod_pos[pos_cif] = 31
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 32 and res_tem == 'K':

                    mod_pos[pos_cif] = 32
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                elif res_cif == 33 and res_tem == 'F':

                    mod_pos[pos_cif] = 33
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'F')

                elif res_cif == 34 and res_tem == 'R':

                    mod_pos[pos_cif] = 34
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'R')

                elif res_cif == 35 and res_tem == 'F':

                    mod_pos[pos_cif] = 35
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'F')

                elif res_cif == 36 and res_tem == 'D':

                    mod_pos[pos_cif] = 36
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'D')

                elif res_cif == 37 and res_tem == 'D':

                    mod_pos[pos_cif] = 37
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'D')

                elif res_cif == 38 and res_tem == 'H':

                    mod_pos[pos_cif] = 38
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'H')

                elif res_cif == 39 and res_tem == 'P':

                    mod_pos[pos_cif] = 39
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'P')

                elif res_cif == 40 and res_tem == 'T':

                    mod_pos[pos_cif] = 40
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'T')

                elif res_cif == 41 and res_tem == 'I':

                    mod_pos[pos_cif] = 41
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'I')

                elif res_cif == 42 and res_tem == 'N':

                    mod_pos[pos_cif] = 42
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'N')

                elif res_cif == 43 and res_tem == 'Y':

                    mod_pos[pos_cif] = 43
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 44 and res_tem == 'M':

                    mod_pos[pos_cif] = 44
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 45 and res_tem == 'G':

                    mod_pos[pos_cif] = 45
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'G')

                elif res_cif == 46 and res_tem == 'W':

                    mod_pos[pos_cif] = 46
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'W')

                elif res_cif == 47 and res_tem == 'S':

                    mod_pos[pos_cif] = 47
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'S')

                elif res_cif == 48 and res_tem == 'N':

                    mod_pos[pos_cif] = 48
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'N')

                elif res_cif == 49 and res_tem == 'M':

                    mod_pos[pos_cif] = 49
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 50 and res_tem == 'A':

                    mod_pos[pos_cif] = 50
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'A')

                elif res_cif == 51 and res_tem == 'Y':

                    mod_pos[pos_cif] = 51
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 52 and res_tem == 'P':

                    mod_pos[pos_cif] = 52
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'P')

                elif res_cif == 53 and res_tem == 'V':

                    mod_pos[pos_cif] = 53
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'V')

                elif res_cif == 54 and res_tem == 'L':

                    mod_pos[pos_cif] = 54
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'L')

                elif res_cif == 55 and res_tem == 'Y':

                    mod_pos[pos_cif] = 55
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 56 and res_tem == 'E':

                    mod_pos[pos_cif] = 56
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'E')
                
                elif res_cif == 57 and res_tem == 'F':

                    mod_pos[pos_cif] = 57
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'F')

                elif res_cif == 58 and res_tem == 'C':

                    mod_pos[pos_cif] = 58
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 59 and res_tem == 'H':

                    mod_pos[pos_cif] = 59
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'H')

                elif res_cif == 60 and res_tem == 'C':

                    mod_pos[pos_cif] = 60
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 61 and res_tem == 'Y':

                    mod_pos[pos_cif] = 61
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 62 and res_tem == 'C':

                    mod_pos[pos_cif] = 62
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 63 and res_tem == 'M':

                    mod_pos[pos_cif] = 63
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 64 and res_tem == 'M':

                    mod_pos[pos_cif] = 64
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 65 and res_tem == 'C':

                    mod_pos[pos_cif] = 65
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 66 and res_tem == 'Y':

                    mod_pos[pos_cif] = 66
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'Y')

                elif res_cif == 67 and res_tem == 'S':

                    mod_pos[pos_cif] = 67
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'S')

                elif res_cif == 68 and res_tem == 'C':

                    mod_pos[pos_cif] = 68
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'C')

                elif res_cif == 69 and res_tem == 'L':

                    mod_pos[pos_cif] = 69
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'L')

                elif res_cif == 70 and res_tem == 'M':

                    mod_pos[pos_cif] = 70
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'M')

                elif res_cif == 71 and res_tem == 'K':

                    mod_pos[pos_cif] = 71
                    mapping[(pos_tem, res_tem)] = (pos_cif, 'K')

                else:
                    # insert a gap
                    template_grid.set_empty(row, pos_tem)
                    mapping[(pos_tem, res_tem)] = (pos_cif, res_cif)
                    
                    if DEBUG_MODE:
                        print ('! {p} {c}: template pos {pt} ({rt}) does not match cif pos {pc} ({rc}). Replacing with gap.'.format( 
                            p = pdb_code, 
                            c = chain,
                            pt = pos_tem, 
                            rt = res_tem, 
                            pc = pos_cif, 
                            rc = res_cif if res_cif != -1 else 'DNE'))
                    
                    if res_cif != -1:
                        del_pos.append(pos_cif)
            else:
                mapping[(pos_tem, res_tem)] = (pos_cif, res_cif)


        # adjust template start and end positions
        correct_mapping = { key:value for key, value in mapping.items() if key[1] == value[1] }

        try:
            tstart = correct_mapping[sorted(correct_mapping.keys())[0]][0]
            tend = correct_mapping[sorted(correct_mapping.keys())[-1]][0]
            template_grid.set_map(row, tstart, tend)
        except IndexError:
            # This exception handles cases in which all residues were deleted
            if DEBUG_MODE:
                print ('! {p} {c}: Removing {p} from alignment. No residues matched the alignment sequence.'.format( 
                    p = pdb_code, 
                    c = chain))

            del_row.append(row)
            continue

        ########################################################################
        ## Delete rows from the PIR Alignment if the residue ratio is to low  ##
        ########################################################################

        if threshold > 0:
            
            gaps = 0
            res = 0

            for col in range(template_grid.get_grid_width()):
                if template_grid.is_empty(row, col):
                    template_grid.set_gap(row, col)

                if template_grid.is_gap(row, col):
                    gaps += 1
                else:
                    res += 1

            ratio = res/float(gaps + res)
        
            if ratio > threshold:
                print ('! Template {p} successfully passed residue ratio ({r:.2f} / {t}).'.format( 
                    p = pdb_code,
                    r = ratio,
                    t = threshold ))
            else:
                print ('! Template {p} did not passed residue ratio ({r:.2f} / {t}). Removing it from pir Alignment.'.format( 
                    p = pdb_code,
                    r = ratio,
                    t = threshold ))            
            
                if row not in del_row:
                    del_row.append(row)
                    continue

        ########################################################################
        ## Edit cif files                                                     ##
        ########################################################################

        rem_row = list() # verbosity: saves information about removed residues
        mod_row = list() # verbosity: saves information about modified residues
        cha_row = list() # verbosity: saves any other changes 

        for atom_row in reversed(range(0, atom_site.getRowCount())):
            
            try:
                cif_chain = atom_site.getValue('label_asym_id', atom_row)
            except IndexError:
                pass

            # We do not care about the residues apart from the chain
            if cif_chain != chain:
                continue    

            # get the residue number 
            try:
                res_num = int(atom_site.getValue("label_seq_id", atom_row))
            except ValueError:
                continue

            # pdb_PDB_model_num has to be set to 1
            try:
                model_num = int(atom_site.getValue('pdbx_PDB_model_num', atom_row))
            except IndexError:
                model_num = 1 # if we cannot extract, assume that it is alright

            try:
                ins_code = atom_site.getValue('pdbx_PDB_ins_code', atom_row)
            except IndexError:
                ins_code = '?' # assume it has no insertion code

            group_PDB = atom_site.getValue('group_PDB', atom_row)
            residue = atom_site.getValue('label_comp_id', atom_row)
            residue = convert_aa_code(residue, convert)         

            # MODELLER accepts only structures if pdbx_PDB_model_num is set to 1   
            if model_num != 1:

                if (res_num, residue, 'model_num') not in cha_row:
                    cha_row.append((res_num, residue, 'model_num'))             
                
                atom_site.setValue(1, "pdbx_PDB_model_num", atom_row)

            if ins_code != '?':
                
                if (res_num, residue, 'ins_code') not in cha_row:
                    cha_row.append((res_num, residue, 'ins_code'))
                
                atom_site.setValue('?', "pdbx_PDB_ins_code", atom_row)

            if group_PDB != 'ATOM':

                if (res_num, residue, 'group_PDB') not in cha_row:
                    cha_row.append((res_num, residue, 'group_PDB'))

                atom_site.setValue('ATOM', 'group_PDB', atom_row)               

            ########################################################################
            ## Delete residues                                                    ##
            ########################################################################

            if res_num in del_pos:
                if (res_num, residue) not in rem_row:
                    rem_row.append((res_num, residue))
                
                atom_site.removeRow(atom_row)

            ########################################################################
            ## Modify residues                                                    ##
            ########################################################################

            if res_num in mod_pos.keys():
                
                # Get the data
                type_symbol = atom_site.getValue('type_symbol', atom_row)
                label_atom_id = atom_site.getValue('label_atom_id', atom_row)
                auth_atom_id = atom_site.getValue('auth_atom_id', atom_row)

                if mod_pos[res_num] == 1: # try to convert MSE to M 
                                        
                    atom_site.setValue('MET', 'label_comp_id', atom_row)
                    
                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if type_symbol == 'SE':
                        atom_site.setValue('S', 'type_symbol', atom_row)
                    if label_atom_id == 'SE':
                        atom_site.setValue('S', 'label_atom_id', atom_row)
                    if auth_atom_id == 'SE':
                        atom_site.setValue('S', 'auth_atom_id', atom_row)

                    if (res_num, residue, 'MSE -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'MSE -> MET'))

                elif mod_pos[res_num] == 2: # try to convert HYP to PRO
                    # apparently it is enough to rename the label_comp_id to PRO to get 
                    # MODELLER working with Hydroxyprolines (HYP)

                    atom_site.setValue('PRO', 'label_comp_id', atom_row)
                    
                    try:
                        atom_site.setValue('PRO', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'HYP -> PRO') not in mod_row:
                        mod_row.append((res_num, residue, 'HYP -> PRO'))

                elif mod_pos[res_num] == 3: # try to convert MLY to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)
                    
                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MLY -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'MLY -> LYS'))

                elif mod_pos[res_num] == 4: # converts Phosphoserine to Serine

                    atom_site.setValue('SER', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('SER', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SEP -> SER') not in mod_row:
                        mod_row.append((res_num, residue, 'SEP -> SER'))

                elif mod_pos[res_num] == 5: # converts Phosphothreonine to Threonine

                    atom_site.setValue('THR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('THR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'TPO -> THR') not in mod_row:
                        mod_row.append((res_num, residue, 'TPO -> THR'))

                elif mod_pos[res_num] == 6: # converts S-HYDROXYCYSTEINE to Cysteine

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CSO -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CSO -> CYS'))

                elif mod_pos[res_num] == 7: # converts O-PHOSPHOTYROSINE to Tyrosine

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'PTR -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'PTR -> TYR'))

                elif mod_pos[res_num] == 8: # converts LYSINE NZ-CARBOXYLIC ACID to Lysine

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'KCX -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'KCX -> LYS'))

                elif mod_pos[res_num] == 9: # converts S,S-(2-HYDROXYETHYL)THIOCYSTEINE to Cysteine

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CME -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CME -> CYS'))

                elif mod_pos[res_num] == 10: # converts 3-SULFINOALANINE to Alanine

                    atom_site.setValue('ALA', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ALA', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CSD -> ALA') not in mod_row:
                        mod_row.append((res_num, residue, 'CSD -> ALA'))

                elif mod_pos[res_num] == 11: # converts S-(DIMETHYLARSENIC)CYSTEINE to Cysteine

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CAS -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CAS -> CYS'))

                elif mod_pos[res_num] == 12: # converts N-METHYLLEUCINE (MLE) to Leucine

                    atom_site.setValue('LEU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LEU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MLE -> LEU') not in mod_row:
                        mod_row.append((res_num, residue, 'MLE -> LEU'))

                elif mod_pos[res_num] == 13: # converts D-ALANINE (DAL) to ALA

                    atom_site.setValue('ALA', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ALA', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DAL -> ALA') not in mod_row:
                        mod_row.append((res_num, residue, 'DAL -> ALA'))

                elif mod_pos[res_num] == 14: # converts GAMMA-CARBOXY-GLUTAMIC ACID (CGU) to GLU

                    atom_site.setValue('GLU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('GLU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CGU -> GLU') not in mod_row:
                        mod_row.append((res_num, residue, 'CGU -> GLU'))

                elif mod_pos[res_num] == 15: # converts D-LEUCINE (DLE) to LEU

                    atom_site.setValue('LEU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LEU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DLE -> LEU') not in mod_row:
                        mod_row.append((res_num, residue, 'DLE -> LEU'))

                elif mod_pos[res_num] == 16: # converts N-FORMYLMETHIONINE (FME) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'FME -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'FME -> MET'))

                elif mod_pos[res_num] == 17: # converts D-VAL (DVA) to VAL

                    atom_site.setValue('VAL', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('VAL', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DVA -> VAL') not in mod_row:
                        mod_row.append((res_num, residue, 'DVA -> VAL'))

                elif mod_pos[res_num] == 18: # converts CYSTEINESULFONIC ACID (OCS) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'OCS -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'OCS -> CYS'))

                elif mod_pos[res_num] == 19: # converts D-PROLINE (DPR) to PRO

                    atom_site.setValue('PRO', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PRO', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DPR -> PRO') not in mod_row:
                        mod_row.append((res_num, residue, 'DPR -> PRO'))

                elif mod_pos[res_num] == 20: # converts N-METHYLVALINE (MVA) to VAL

                    atom_site.setValue('VAL', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('VAL', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MVA -> VAL') not in mod_row:
                        mod_row.append((res_num, residue, 'MVA -> VAL'))

                elif mod_pos[res_num] == 21: # converts O-SULFO-L-TYROSINE (TYS) to VAL

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'TYS -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'TYS -> TYR'))

                elif mod_pos[res_num] == 22: # converts N-TRIMETHYLLYSINE (M3L) to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'M3L -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'M3L -> LYS'))

                elif mod_pos[res_num] == 23: # converts S-METHYLCYSTEINE (SMC) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SMC -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'SMC -> CYS'))

                elif mod_pos[res_num] == 24: # converts N(6)-ACETYLLYSINE (ALY) to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'ALY -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'ALY -> LYS'))

                elif mod_pos[res_num] == 25: # converts S-OXY CYSTEINE (CSX) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CSX -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CSX -> CYS'))

                elif mod_pos[res_num] == 26: # converts D-CYSTEINE (DCY) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DCY -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'DCY -> CYS'))

                elif mod_pos[res_num] == 27: # converts NORLEUCINE (NLE) to LEU

                    atom_site.setValue('LEU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LEU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'NLE -> LEU') not in mod_row:
                        mod_row.append((res_num, residue, 'NLE -> LEU'))

                elif mod_pos[res_num] == 28: # converts D-GLUTAMIC ACID (DGL) to GLU

                    atom_site.setValue('GLU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('GLU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DGL -> GLU') not in mod_row:
                        mod_row.append((res_num, residue, 'DGL -> GLU'))

                elif mod_pos[res_num] == 29: # converts D-SERINE (DSN) to SER

                    atom_site.setValue('SER', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('SER', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DSN -> SER') not in mod_row:
                        mod_row.append((res_num, residue, 'DSN -> SER'))

                elif mod_pos[res_num] == 30: # converts S-MERCAPTOCYSTEINE (CSS) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CSS -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CSS -> CYS'))

                elif mod_pos[res_num] == 31: # converts D-LYSINE (DLY) to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DLY -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'DLY -> LYS'))

                elif mod_pos[res_num] == 32: # converts N-METHYL-LYSINE (MLZ) to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MLZ -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'MLZ -> LYS'))

                elif mod_pos[res_num] == 33: # converts D-PHENYLALANINE (DPN) to PHE

                    atom_site.setValue('PHE', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PHE', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DPN -> PHE') not in mod_row:
                        mod_row.append((res_num, residue, 'DPN -> PHE'))

                elif mod_pos[res_num] == 34: # converts D-ARGININE (DAR) to ARG

                    atom_site.setValue('ARG', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ARG', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DAR -> ARG') not in mod_row:
                        mod_row.append((res_num, residue, 'DAR -> ARG'))

                elif mod_pos[res_num] == 35: # converts IODO-PHENYLALANINE (PHI) to PHE

                    atom_site.setValue('PHE', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PHE', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'PHI -> PHE') not in mod_row:
                        mod_row.append((res_num, residue, 'PHI -> PHE'))

                elif mod_pos[res_num] == 36: # converts BETA-L-ASPARTIC ACID (IAS) to ASP

                    atom_site.setValue('ASP', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ASP', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'IAS -> ASP') not in mod_row:
                        mod_row.append((res_num, residue, 'IAS -> ASP'))

                elif mod_pos[res_num] == 37: # converts D-ASPARTIC ACID (DAS) to ASP

                    atom_site.setValue('ASP', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ASP', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DAS -> ASP') not in mod_row:
                        mod_row.append((res_num, residue, 'DAS -> ASP'))

                elif mod_pos[res_num] == 38: # converts 4-METHYL-HISTIDINE (HIC) to HIS

                    atom_site.setValue('HIS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('HIS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'HIC -> HIS') not in mod_row:
                        mod_row.append((res_num, residue, 'HIC -> HIS'))

                elif mod_pos[res_num] == 39: # converts (4R)-4-methyl-L-proline (MP8) to PRO

                    atom_site.setValue('PRO', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PRO', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MP8 -> PRO') not in mod_row:
                        mod_row.append((res_num, residue, 'MP8 -> PRO'))

                elif mod_pos[res_num] == 40: # converts D-THREONINE (DTH) to THR

                    atom_site.setValue('THR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('THR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DTH -> THR') not in mod_row:
                        mod_row.append((res_num, residue, 'DTH -> THR'))

                elif mod_pos[res_num] == 41: # converts D-ISOLEUCINE (DIL) to ILE

                    atom_site.setValue('ILE', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ILE', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DIL -> ILE') not in mod_row:
                        mod_row.append((res_num, residue, 'DIL -> ILE'))

                elif mod_pos[res_num] == 42: # converts N-METHYL ASPARAGINE (MEN) to ASN

                    atom_site.setValue('ASN', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ASN', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MEN -> ASN') not in mod_row:
                        mod_row.append((res_num, residue, 'MEN -> ASN'))

                elif mod_pos[res_num] == 43: # converts D-TYROSINE (DTY) to TYR

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DTY -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'DTY -> TYR'))

                elif mod_pos[res_num] == 44: # converts N-CARBOXYMETHIONINE (CXM) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CXM -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'CXM -> MET'))

                elif mod_pos[res_num] == 45: # converts D-GLUTAMINE (DGN) to MET

                    atom_site.setValue('GLN', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('GLN', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DGN -> GLN') not in mod_row:
                        mod_row.append((res_num, residue, 'DGN -> GLN'))

                elif mod_pos[res_num] == 46: # converts D-TRYPTOPHAN (DTR) to TRP

                    atom_site.setValue('TRP', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TRP', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DTR -> TRP') not in mod_row:
                        mod_row.append((res_num, residue, 'DTR -> TRP'))

                elif mod_pos[res_num] == 47: # converts N-ACETYL-SERINE (SAC) to SER

                    atom_site.setValue('SER', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('SER', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SAC -> SER') not in mod_row:
                        mod_row.append((res_num, residue, 'SAC -> SER'))

                elif mod_pos[res_num] == 48: # converts D-ASPARAGINE (DSG) to ASN

                    atom_site.setValue('ASN', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ASN', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DSG -> ASN') not in mod_row:
                        mod_row.append((res_num, residue, 'DSG -> ASN'))

                elif mod_pos[res_num] == 49: # converts N-METHYL METHIONINE (MME) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MME -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'MME -> MET'))

                elif mod_pos[res_num] == 50: # converts N-methyl-L-alanine (MAA) to ALA

                    atom_site.setValue('ALA', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('ALA', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MAA -> ALA') not in mod_row:
                        mod_row.append((res_num, residue, 'MAA -> ALA'))

                elif mod_pos[res_num] == 51: # converts 3-FLUOROTYROSINE (YOF) to TYR

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'YOF -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'YOF -> TYR'))

                elif mod_pos[res_num] == 52: # converts (4R)-4-fluoro-L-proline (FP9) to PRO

                    atom_site.setValue('PRO', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PRO', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'FP9 -> PRO') not in mod_row:
                        mod_row.append((res_num, residue, 'FP9 -> PRO'))

                elif mod_pos[res_num] == 53: # converts N-formyl-L-valine (FVA) to VAL

                    atom_site.setValue('VAL', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('VAL', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'FVA -> VAL') not in mod_row:
                        mod_row.append((res_num, residue, 'FVA -> VAL'))

                elif mod_pos[res_num] == 54: # converts N-methyl-D-leucine (MLU) to LEU

                    atom_site.setValue('LEU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LEU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MLU -> LEU') not in mod_row:
                        mod_row.append((res_num, residue, 'MLU -> LEU'))

                elif mod_pos[res_num] == 55: # converts (betaR)-3-chloro-beta-hydroxy-L-tyrosine (OMY) to TYR

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'OMY -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'OMY -> TYR'))

                elif mod_pos[res_num] == 56: # converts GAMMA-D-GLUTAMIC ACID (FGA) to GLU

                    atom_site.setValue('GLU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('GLU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'FGA -> GLU') not in mod_row:
                        mod_row.append((res_num, residue, 'FGA -> GLU'))

                elif mod_pos[res_num] == 57: # converts N-METHYLPHENYLALANINE (MEA) to PHE

                    atom_site.setValue('PHE', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('PHE', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MEA -> PHE') not in mod_row:
                        mod_row.append((res_num, residue, 'MEA -> PHE'))

                elif mod_pos[res_num] == 58: # converts S-(METHYLMERCURY)-L-CYSTEINE (CMH) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CMH -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CMH -> CYS'))

                elif mod_pos[res_num] == 59: # converts D-HISTIDINE (DHI) to HIS

                    atom_site.setValue('HIS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('HIS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'DHI -> HIS') not in mod_row:
                        mod_row.append((res_num, residue, 'DHI -> HIS'))

                elif mod_pos[res_num] == 60: # converts SELENOCYSTEINE (SEC) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SEC -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'SEC -> CYS'))

                elif mod_pos[res_num] == 61: # converts (betaR)-3-CHLORO-BETA-HYDROXY-D-TYROSINE (OMZ) to TYR

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'OMZ -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'OMZ -> TYR'))

                elif mod_pos[res_num] == 62: # converts S-ACETYL-CYSTEINE (SCY) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SCY -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'SCY -> CYS'))

                elif mod_pos[res_num] == 63: # converts S-OXYMETHIONINE (MHO) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MHO -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'MHO -> MET'))

                elif mod_pos[res_num] == 64: # converts D-METHIONINE (MED) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MED -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'MED -> MET'))

                elif mod_pos[res_num] == 65: # converts S-DIMETHYLARSINOYL-CYSTEINE (CAF) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'CAF -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'CAF -> CYS'))

                elif mod_pos[res_num] == 66: # converts META-NITRO-TYROSINE (NIY) to TYR

                    atom_site.setValue('TYR', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('TYR', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'NIY -> TYR') not in mod_row:
                        mod_row.append((res_num, residue, 'NIY -> TYR'))

                elif mod_pos[res_num] == 67: # converts O-ACETYLSERINE (OAS) to SER

                    atom_site.setValue('SER', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('SER', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'OAS -> SER') not in mod_row:
                        mod_row.append((res_num, residue, 'OAS -> SER'))

                elif mod_pos[res_num] == 68: # converts S-METHYL-THIO-CYSTEINE (SCH) to CYS

                    atom_site.setValue('CYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('CYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SCH -> CYS') not in mod_row:
                        mod_row.append((res_num, residue, 'SCH -> CYS'))

                elif mod_pos[res_num] == 69: # converts 2-methyl-L-norleucine (MK8) to LEU

                    atom_site.setValue('LEU', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LEU', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'MK8 -> LEU') not in mod_row:
                        mod_row.append((res_num, residue, 'MK8 -> LEU'))

                elif mod_pos[res_num] == 70: # converts METHIONINE SULFOXIDE (SME) to MET

                    atom_site.setValue('MET', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('MET', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'SME -> MET') not in mod_row:
                        mod_row.append((res_num, residue, 'SME -> MET'))

                elif mod_pos[res_num] == 71: # converts 5-HYDROXYLYSINE (LYZ) to LYS

                    atom_site.setValue('LYS', 'label_comp_id', atom_row)

                    try:
                        atom_site.setValue('LYS', 'auth_comp_id', atom_row)
                    except IndexError:
                        pass

                    if (res_num, residue, 'LYZ -> LYS') not in mod_row:
                        mod_row.append((res_num, residue, 'LYZ -> LYS'))

        ########################################################################
        ## Notify user about modification made to cif data                    ##
        ########################################################################

        if DEBUG_MODE:
            mod_model_num = len([ msg for msg in cha_row if msg[2] == 'model_num' ])
            mod_ins_code = len([ msg for msg in cha_row if msg[2] == 'ins_code' ])
            mod_group_PDB = len([ msg for msg in cha_row if msg[2] == 'group_PDB' ])
            
            if mod_model_num != 0:
                print ('! {p} {c}: modified atom_site.pdbx_PDB_model_num for {cr} residues to 1.'.format(
                        p = pdb_code,
                        c = chain,
                        cr = mod_model_num))

            if mod_ins_code != 0:
                print ('! {p} {c}: modified atom_site.pdbx_PDB_ins_code for {cr} residues to "?".'.format(
                        p = pdb_code,
                        c = chain,
                        cr = mod_ins_code))

            if mod_group_PDB != 0:
                print ('! {p} {c}: modified atom_site.group_PDB for {cr} residues to "ATOM".'.format(
                        p = pdb_code,
                        c = chain,
                        cr = mod_group_PDB))

            for residue in reversed(mod_row):
                print ('! {p} {c}: modified cif pos {cr} ({nr}).'.format(
                    p = pdb_code, 
                    c = chain,
                    cr = residue[0],
                    ca = residue[1],
                    nr = residue[2]))


            for residue in reversed(rem_row):
                print ('! {p} {c}: removed cif pos {cr} ({ca})'.format(
                    p = pdb_code, 
                    c = chain,
                    cr = residue[0],
                    ca = residue[1]))

        cif_edits[pdb_code] = block

    # write modified pir to disk
    for pdb_code in cif_edits:
        out = open(os.path.join(output_path, pdb_code + '.cif'), 'w')
        writer = PdbxWriter(out)
        writer.writeContainer(cif_edits[pdb_code])

    # Delete missing entries from the last template sequence to the first
    for row in reversed(del_row):
        template_grid.del_row(row)

    return template_grid

def remove_self_alignment(template_grid, query_name):
    """ Removes a self alignment from the final pir alignment to prevent clashes with MODELLER """
    
    to_delete = list()
    
    for row in range(template_grid.get_grid_height()):
        if template_grid._pdb_code[row] == query_name:
            to_delete.append(row)

    for row in reversed(to_delete):
        template_grid.del_row(row)

    return True

def write_to_file(line_list, fname):
    """ Writes the final pir file """
    
    with open(fname, 'w+') as fout:
        for line in line_list:
            fout.write(line + "\n")

def arg():
    import argparse
    description = """Creates a MODELLER alignment (*.pir) from a HHSearch results file (*.hhr)."""
    epilog= '2016 Harald Voehringer.'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help = 'results file from HHsearch with hit list and alignment', metavar = 'FILE')    
    parser.add_argument('cifs', help = 'path to the folder containing cif files', metavar = 'DIR')
    parser.add_argument('pir', help = 'output file (PIR-formatted multiple alignment)', metavar = 'FILE')
    parser.add_argument('output', help = 'path to the folder where modified cif files should be written to', metavar = 'DIR')

    parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose mode')
    parser.add_argument('-m', nargs = '+', help = 'pick hits with specified indices (e.g. -m 2 5)', metavar = 'INT')
    parser.add_argument('-e', type = float, help = 'maximum E-Value threshold (e.g. -e 0.001)', metavar = 'FLOAT')
    parser.add_argument('-r', type = float, help = 'residue ratio (filter alignments that have contribute at least residues according to the specified ratio).',
        default = 0, metavar = 'FLOAT')

    parser.add_argument('-c', help = 'convert non-canonical residues (default = True)', action = 'store_true', default = True)
    

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])
    
    global DEBUG_MODE
    
    if args.verbose:
        DEBUG_MODE = True

    query_name, query_chain = get_query_name(args.input)

    data = read_result(args.input)
    selected_templates = list()

    if args.m and not args.e:
        selection = map(lambda x: int(x), args.m)
        print ('Selected templates {st}.'.format(st = ', '.join(args.m)))
        
        for i in selection:
            tmp_info = str(data[i - 1].template_info.split('>')[1])
            print ('{i}: {t}'.format(
                i = i,
                t = tmp_info[0:80]))

            selected_templates.append(data[i - 1])

    elif args.e and not args.m:
        print ('Selected templates satisfying E-val <= {e}'.format(e = args.e))
        
        e_values = { float(j.evalue):i for i, j in enumerate(data) }
        selection = sorted([ val for key, val in e_values.items() if key <= args.e ])

        for i in selection:
            tmp_info = str(data[i - 1].template_info.split('>')[1])
            print ('{i}: {t}'.format(
                i = i + 1,
                t = tmp_info[0:80]))

            selected_templates.append(data[i - 1])

    elif args.m and args.e:
        print ('! Please do not use option -m and -e at the same time ! Exiting.')
        sys.exit()
    else:
        selected_templates = data
        
        print ('Creating pir file using all templates ({n})'.format(
            n = len(selected_templates)))

    query_grid = create_query_grid(selected_templates) # load query grid
    print ('query_grid')
    print(query_grid)
    gapless_query_grid = create_gapless_grid(query_grid) # remove gaps
    print ('gapless_query_grid')
    print(gapless_query_grid)
    processed_query_grid = process_query_grid(query_grid, gapless_query_grid) # insert gaps
    ##processed_query_grid = process_query_grid(query_grid, query_grid) # insert gaps
    print ('processed_query_grid')
    print (processed_query_grid)
    glob_seq = derive_global_seq(processed_query_grid, query_name, query_chain) # derive query sequence
    template_grid = create_template_grid(selected_templates) # create template grid
    print ('template_grid')
    print (template_grid)
    processed_template_grid = process_template_grid(query_grid, template_grid) # insert gaps to template sequnces
    print ('processed_query_grid')
    print (processed_query_grid)
    print ('hzhu processed_template_grid')
    print (processed_template_grid)
    final_grid = compare_with_cifs(processed_template_grid, args.cifs, args.output, args.c, args.r) # compare with atom section of cifs
    remove_self_alignment(final_grid, query_name) # remove self alignment if any
    write_to_file([glob_seq, final_grid.display()], args.pir)


if __name__ == "__main__":
    main()
