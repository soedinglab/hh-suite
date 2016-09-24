#!/usr/bin/env python

"""
Created on Mon Jun 15 21:49:32 2015

@author: Harald Voehringer
"""

import sys, os, glob, textwrap, itertools
from optparse import OptionParser
from collections import defaultdict
from os.path import splitext
from pdbx.reader.PdbxReader import PdbxReader
from multiprocessing import Pool

DEBUG_MODE = False
MIN_SEQ_LEN = None
SCOP_LIBRARY = False

THREE2ONE = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 
        'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 
        'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE': 'M',
        'HYP': 'P', 'MLY': 'K', 'SEP': 'S', 'TPO': 'T', 'CSO': 'C', 'PTR': 'Y', 'KCX': 'K',
        'CME': 'C', 'CSD': 'A', 'CAS': 'C', 'MLE': 'L', 'DAL': 'A', 'CGU': 'E', 'DLE': 'L',
        'FME': 'M', 'DVA': 'V', 'OCS': 'C', 'DPR': 'P', 'MVA': 'V', 'TYS': 'Y', 'M3L': 'K',
        'SMC': 'C', 'ALY': 'K', 'CSX': 'C', 'DCY': 'C', 'NLE': 'L', 'DGL': 'E', 'DSN': 'S',
        'CSS': 'C', 'DLY': 'K', 'MLZ': 'K', 'DPN': 'F', 'DAR': 'R', 'PHI': 'F', 'IAS': 'D',
        'DAS': 'D', 'HIC': 'H', 'MP8': 'P', 'DTH': 'T', 'DIL': 'I', 'MEN': 'N', 'DTY': 'Y',
        'CXM': 'M', 'DGN': 'G', 'DTR': 'W', 'SAC': 'S', 'DSG': 'N', 'MME': 'M', 'MAA': 'A',
        'YOF': 'Y', 'FP9': 'P', 'FVA': 'V', 'MLU': 'L', 'OMY': 'Y', 'FGA': 'E', 'MEA': 'F',
        'CMH': 'C', 'DHI': 'H', 'SEC': 'C', 'OMZ': 'Y', 'SCY': 'C', 'MHO': 'M', 'MED': 'M',
        'CAF': 'C', 'NIY': 'Y', 'OAS': 'S', 'SCH': 'C', 'MK8': 'L', 'SME': 'M', 'LYZ': 'K'
    }

CANONICAL_RESIDUES = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
      'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

class CIF2FASTA(object):
    
    def __init__(self, cif_path):
        self.cif_path = cif_path
        self.block = self.open_cif()
            
    def open_cif(self):
        """ Assumes a mmCif file and returns a data block used for subsequent procedures. """
        # The "usual" procedure to open a mmCIF with pdbX/mmCIF

        with open(self.cif_path) as cif_fh:
            data = []
            reader = PdbxReader(cif_fh)
            reader.read(data)
            if len(data) == 0:
                return None
            else:
                return data[0]

    def is_valid(self):
        return self.block != None
 
    def chain_to_seq(self):
        """Extracts the sequence of the cif from entity_poly.pdbx_seq_one_letter_code"""
        
        cif_chain_to_seq = dict()
        non_polypeptide_chains = list()         
        
        try:
            entity_poly = self.block.getObj('entity_poly')
        except AttributeError:

            if DEBUG_MODE > 0:
                print ('! {pdb} Could not extract entity_poly table.'.format(
                    pdb = self.pdb_entry()))

            return False

        try:
            total_rows = entity_poly.getRowCount()
        except AttributeError:
            print ('! {pdb} Could not extract rows from entity_poly.'.format(
                pdb = self.pdb_entry()))

            return False

        for row in range(0, total_rows):

            if entity_poly.getValue('type', row) == 'polypeptide(L)':
                seq = entity_poly.getValue('pdbx_seq_one_letter_code', row)                                                            
                parsed_seq = parse_seq(seq) # removes special amino acids and newlines
                
                try:
                    chains = entity_poly.getValue('pdbx_strand_id', row)
                except ValueError:
                    
                    if total_rows == 1:
                        print ('! {pdb} Only one polypeptide chain, but no chain identifiers, setting it to ".".'.format(
                            pdb = self.pdb_entry()))
                        
                        cif_chain_to_seq['.'] = parsed_seq
                        
                        return cif_chain_to_seq

                    print ('! {pdb} Could not extract pdbx_strand_id from entity_poly table (polypeptide).'.format(
                        pdb = self.pdb_entry()))

                    return False

                chain_list = chains.split(',')
            
                for chain in chain_list:     
                    cif_chain_to_seq[chain] = parsed_seq
            
            else:
                
                try:
                    chains = entity_poly.getValue('pdbx_strand_id', row)
                except ValueError:
                    print ('! {pdb} Could not extract pdbx_strand_id from entity_poly table (non-polypeptide).'.format(
                        pdb = self.pdb_entry()))
                    
                    return False

                non_polypeptide_chains.append(chains)

        chains = list(cif_chain_to_seq.keys())
        # remove chains that contain only unknown residues
        for chain in chains:
            # this is a very odd way to check whether a string contains only a single char
            tmp_set = set(cif_chain_to_seq[chain])
            if len(tmp_set) == 1 and 'X' in tmp_set:
                print ('! Removing {pdb}_{chain} (contains only unknown residues).'.format(
                    pdb = self.pdb_entry(),
                    chain = chain))
                
                del cif_chain_to_seq[chain]
                continue
            
            if len(cif_chain_to_seq[chain]) < MIN_SEQ_LEN:
                print ('! Removing {pdb}_{chain} (sequence length < {min_len}).'.format(
                    pdb = self.pdb_entry(),
                    chain = chain,
                    min_len = MIN_SEQ_LEN))
                
                del cif_chain_to_seq[chain]

        if len(cif_chain_to_seq) != 0:
            
            if DEBUG_MODE > 1:
                print ('- Extracted chains of {pdb} {chains}.'.format(
                    pdb = self.pdb_entry(),
                    chains = ' '.join( str(chain) + ' (' + str(len(cif_chain_to_seq[chain])) + ')' for chain in cif_chain_to_seq.keys())))
                
                if len(non_polypeptide_chains) != 0:                
                    print ('- Following chains were non polypeptide chains {chains} no polypeptide chains were found.'.format(
                        chains = ', '.join(non_polypeptide_chains)))
            
            return cif_chain_to_seq
        
        else:
            if DEBUG_MODE > 0:
                print ('! {pdb} No polypeptide chains were found.'.format(
                    pdb = self.pdb_entry()))
            
        return False

    def chain_ratios(self, chain_to_seq):
        """ Tries to extract Sequence from the atom section """

        # chain_to_seq = self.chain_to_seq()

        if chain_to_seq != False:
            chain_ratios = dict()
            
            # compute the lengths of sequences found in _entity_poly
            entity_length = { chain : float(len(seq)) for chain, seq in chain_to_seq.items() }
            entity_chains = entity_length.keys()

            # load the atomsite and set up dictionary to keep track of sequences
            atom_site = self.block.getObj('atom_site')
            atom_seq = defaultdict(str)

            current_residue = 0
            # Iterate through the atomsection of the cif file
            for atom_row in range(0, atom_site.getRowCount()):

                # NMR structures contain many confomers
                try:
                    model_num = int(atom_site.getValue('pdbx_PDB_model_num', atom_row))
                except ValueError:
                    model_num = 1

                if model_num > 1:
                    continue

                atom_chain = atom_site.getValue('label_asym_id', atom_row)

                # get the alternative chain identifier too
                try:
                    alt_chain = atom_site.getValue('auth_asym_id', atom_row)
                except ValueError:
                    alt_chain = None

                # handle cases where there are no chains but only one structure
                if atom_chain == '.' and entity_chains[0] == '.':
                    atom_chain = '.'
                
                # get the residue and the residue number 
                try:
                    res_num = int(atom_site.getValue("label_seq_id", atom_row))
                except ValueError:
                    continue

                if res_num != current_residue:
                    residue = atom_site.getValue('label_comp_id', atom_row)
                    
                    try:
                        residue = THREE2ONE[residue]
                    except KeyError:
                        residue = 'X'

                    # try to get the chain identifier from alt_chain first, if this does not work use label_asym_id
                    if alt_chain != None:
                    	atom_seq[alt_chain] += residue

                    # sometimes we find the right chain identifier not in the alt_chain
                    if not (atom_chain in atom_seq.keys()) and atom_chain != None:
                    	atom_seq[atom_chain] += residue
                        
                    current_residue = res_num

            for chain in entity_length.keys():
                if chain in atom_seq.keys():
                    chain_ratios[chain] = len(atom_seq[chain]) / entity_length[chain]
                else:
                    chain_ratios[chain] = 0

            return chain_ratios
        else:
            return False

    def pdb_entry(self):
        """Extracts the PDB entry information of a cif file."""
        try:
            entry = self.block.getObj('entry')
            entry_id = entry.getValue('id')
            return entry_id.replace('\n', ' ')
        
        except AttributeError:
            if DEBUG_MODE > 0:
                print ('! {pdb} Could not extract id from entry.'.format(
                    pdb = self.pdb_entry()))

    def protein_description(self):
        """Extracts the protein description annotated in struct.pdbx_descriptor of the cif file."""
        try:
            # Get struct table which contains the protein description
            struct = self.block.getObj('struct')
            # Get the pdbx description and make format it appropritaly
            protein_description = struct.getValue('pdbx_descriptor')
            protein_description = protein_description.replace('\n', ' ')
            protein_description = protein_description.replace(';', ' ') # to prevent parsing errors
            
            if len(protein_description.split(' ')) >= 5: 
                protein_description = ' '.join(protein_description.split(' ')[0:5]) # maximum of 5 words in header
            
            return protein_description.strip(',')
        
        except AttributeError:
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract pdbx_descriptor from struct table.'.format(
                    pdb = self.pdb_entry()))
            
            return False

    def compounds(self):
        """ Extracts all compounds annotated in the HETATM section of the atom
            struct table if the compound appears at least 10 times and is not water 
            (HOH)."""    
        
        atom_site = self.block.getObj('atom_site')
        compounds = {}
        
        for row in range(0, atom_site.getRowCount()):
            if atom_site.getValue('group_PDB', row) == 'HETATM':

                label_comp_id = atom_site.getValue('label_comp_id', row)
                
                if label_comp_id not in compounds.keys():
                    compounds[label_comp_id] = 1
                else:
                    compounds[label_comp_id] += 1
    
        filtered_compounds = set()

        for compound in compounds.keys():
            if compounds[compound] >= 10 and compound != 'HOH':
                filtered_compounds.add(compound)

        if len(filtered_compounds) == 0:
            return False
        else:
            return ', '.join(filtered_compounds).replace('\n', ' ')
    
    def resolution(self):
        """Extracts the resolution of the mmCIF."""

        try:
            refine = self.block.getObj('refine')
            resolution = refine.getValue('ls_d_res_high')
            
            try: 
                resolution = float(resolution)
            except ValueError:
                return False
            
            return resolution
        
        except AttributeError:
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract ls_d_res_high from refine table.'.format(
                    pdb = self.pdb_entry()))            

        try:
            reflns = self.block.getObj('reflns')
            # Extract the resolution of the crystal
            resolution = reflns.getValue('d_resolution_high')
            
            try:
                resolution = float(resolution)
            except ValueError:
                return False
            
            return resolution
        
        except AttributeError:
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract d_resolution_high from reflns table.'.format(
                    pdb = self.pdb_entry()))

        # This is true for some Electron Microscopy structures
        try:
            em_3d = self.block.getObj('em_3d_reconstruction')
            resolution = em_3d.getValue('resolution')
            
            try:
                resolution = float(resolution)
            except ValueError:
                return False
            
            return resolution
        except AttributeError:
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract resolution from em_3d_reconstruction table.'.format(
                    pdb = self.pdb_entry())) 

        return False

    def experimental_method(self):
        """Extracts the experimental method of the mmCIF."""

        try:
            reflns = self.block.getObj('exptl')
            method = reflns.getValue('method')
            return method.replace('\n', ' ')
        
        except AttributeError:
        
            if DEBUG_MODE > 1:
                print ('! Could not extract text from exptl table.'.format(
                    pdb = self.pdb_entry()))
            
        return False

    def keywords(self):
        """Extracts the keywords of the mmCIF."""
        try:
            reflns = self.block.getObj('struct_keywords')
            keywords = reflns.getValue('text')
            # perform some string modifications
            keywords = keywords.replace('\n', ' ')
            keywords = keywords.replace(';', ' ')

            if len(keywords.split(' ')) >= 5:
                keywords = ' '.join(keywords.split(' ')[0:5]) 

            return keywords.rstrip(',')
        
        except AttributeError:
        
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract text from struct_keywords table.'.format(
                    pdb = self.pdb_entry()))
            
        return False   

    def organism(self):
        """Extracts the organism of the mmCIF."""
        try:
            entity_src_nat = self.block.getObj('entity_src_nat')
            organsim_scientific = entity_src_nat.getValue('pdbx_organism_scientific')
            return organsim_scientific.replace('\n', ' ')
        
        except AttributeError:
        
            if DEBUG_MODE > 1:
                print ("! {pdb} Could not extract from pdbx_organism_scientific from entity_src_gen table (".format(
                    pdb = self.pdb_entry()))
            pass
            
        try:
            entity_src_gen = self.block.getObj("entity_src_gen")
            src_scientific = entity_src_gen.getValue("pdbx_gene_src_scientific_name")
            return src_scientific.replace("\n", " ")
        
        except AttributeError:
        
            if DEBUG_MODE > 1:
                print ('! {pdb} Could not extract from pdbx_gene_src_scientific_name from entity_src_gen table'.format(
                    pdb = self.pdb_entry()))
            
        return False
        
    def r_free(self):

        try:
            refine = self.block.getObj('refine')
            r_free = refine.getValue('ls_R_factor_R_free')

            try: 
                r_free = float(r_free)
            except ValueError:
                return False

        except AttributeError:
            if DEBUG_MODE > 2:
                print ('! Could not extract R free factor ({pdb})'.format(
                    pdb = self.pdb_entry()))

        except ValueError:
            if DEBUG_MODE > 2:
                print ('! R free factor is not annotated ({pdb})'.format(
                    pdb = self.pdb_entry()))

            
        return False

# Helper functions

def parse_seq(orginal_seq):
    """Parses the cif fasta sequence and replaces non-canonical residues with their canonical counterparts."""
    seq = orginal_seq
    
    while seq.find('(') != -1: 
        start_pos = seq.find('(')
        stop_pos = seq.find(')')
        residue = seq[start_pos + 1:stop_pos]
        print(residue)

        try:
            canonical = THREE2ONE[residue]
        except KeyError:
            canonical = 'X'

        if DEBUG_MODE > 1:
            print ('! Replaced non canonical residue {nc} with {c}'.format(
                nc = residue,
                c = canonical))
                
        if start_pos == 0:
            seq = canonical + seq[stop_pos+1:]
        elif stop_pos == len(seq):
            seq = seq[0:start_pos] + canonical 
        else:
            pre_seq = seq[0:start_pos]
            post_seq = seq[stop_pos+1:]
            seq = pre_seq + canonical + post_seq
    
    seq = seq.replace('\n', '')
    seq_array = []
    
    for c in seq:
        if c in CANONICAL_RESIDUES:
            seq_array.append(c)
        else:
            seq_array.append('X')

    seq = ''.join(seq_array)

    return seq
    
    
def get_paths(in_folder, out_folder):
    """Searches a directory and all its subdirectories for files ending 
    with a specific suffix."""

    paths = list()
    
    for root, dirs, files in os.walk(in_folder):
        for fname in files:
            if fname.endswith(".cif"):
                in_path = os.path.join(root, fname)
                out_file = in_path.split('/')[-1].split('.')[0] + ".fasta"
                out_path = os.path.join(out_folder, out_file)
                #absolute_path = os.path.join(os.getcwd(), fpath)
                paths.append((in_path, out_path))
    
    return paths
        
def construct_header(cif2fasta):
    """Constructs a fasta header."""
    
    protein_description = cif2fasta.protein_description()
    keywords = cif2fasta.keywords()
    compounds = cif2fasta.compounds()
    resolution = cif2fasta.resolution()
    method = cif2fasta.experimental_method()
    organism = cif2fasta.organism()
    r_free = cif2fasta.r_free()

    # construct the header with the data given
    header_information = ''

    if protein_description:
        protein_description = protein_description.replace(';', ' ') # to prevent parsing errors
        protein_description = ' '.join(protein_description.split(' ')[0:5]) # maximum of 5 words in header
        header_information += 'DSC: ' + protein_description + '; '
    else:
        header_information += 'DSC: N/A; '
 
    # if keywords:
    #     header_information += keywords + '; '

    if method:
        header_information += 'MET: ' + method + '; '
    else:
        header_information += 'MET: N/A; '
    
    if resolution:
        header_information += 'RES: ' + resolution + '; '
    else:
        header_information += 'RES: N/A; '

    if r_free:
        header_information += 'RFR: ' + r_free + '; ' 
    else:
        header_information += 'RFR: N/A; ' 

    if organism:
        header_information += 'ORG: ' + organism + '; ' 
    else:
        header_information += 'ORG: N/A; ' 

    if compounds:
        header_information += 'HET: ' + compounds + '; '
    else:
        header_information += 'HET: N/A; '

    if SCOP_LIBRARY:
        
        try: 
            scop_idx = SCOP_LIBRARY[cif2fasta.pdb_entry()]
            if len(scop_idx) != 0:
                domains = ', '.join(scop_idx)
            else:
                domains = 'N/A'
        except KeyError:
            domains = 'N/A'
        
        header_information += 'SCOP: ' + domains + '; '

    return header_information.strip()

def create_fasta_entry(cif2fasta):
    """Combines information given in pdb_entry and strand_seq and yields a clean
    fasta file.
    """
    
    pdb_entry = cif2fasta.pdb_entry()
    header = construct_header(cif2fasta)
    chain_to_seq = cif2fasta.chain_to_seq()
    chain_ratios = cif2fasta.chain_ratios(chain_to_seq)
    #import pdb; pdb.set_trace()
    
    fasta_entry = ""
    
    if chain_to_seq and pdb_entry:
        for chain in sorted(chain_to_seq.keys()):
            if len(chain_to_seq[chain]) != 0:
                
                fasta_entry += '>{pdb}_{chain} {header} CMP: {rat:.2f}\n{seq}\n'.format(
                    pdb = pdb_entry,
                    chain = chain,
                    header = header,
                    rat = chain_ratios[chain],
                    seq = '\n'.join(textwrap.wrap(chain_to_seq[chain], 80)))
                
                # fasta_entry += '>' + pdb_entry + '_' + chain + ' ' + header + '\n' + '\n'.join(textwrap.wrap(chain_to_seq[chain], 80)) + '\n'
        return fasta_entry
    else:
        return None


def create_fasta_entry2(cif2fasta):
    """" Creates a fasta entry."""

    # Get all the information we need
    pdb_entry = cif2fasta.pdb_entry()
    protein_description = cif2fasta.protein_description()
    keywords = cif2fasta.keywords()
    compounds = cif2fasta.compounds()
    resolution = cif2fasta.resolution()
    method = cif2fasta.experimental_method()
    organism = cif2fasta.organism()
    r_free = cif2fasta.r_free()

    # get chains and sequences
    chain_to_seq = cif2fasta.chain_to_seq()
    chain_ratios = cif2fasta.chain_ratios(chain_to_seq)

    is_NMR = False
    if 'NMR' in method:
        is_NMR = True

    fasta_entry = ''
    pdb_filter = '' # If needed create entries that are parseable by pdb filter   
    
    if chain_to_seq and pdb_entry:
        for chain, seq in chain_to_seq.items():

            # check if we can find SCOP domains for the current chain
            domains = False
            
            if SCOP_LIBRARY:
                try: 
                    scop_idx = SCOP_LIBRARY[pdb_entry + '_' + chain]
                    if len(scop_idx) != 0:
                        domains = ', '.join(scop_idx)
                    else:
                        domains = False
                except KeyError:
                    domains = False

            if is_NMR:
                fasta_entry += '>{p}_{c}{d}{k}{h}{r}{{{o}}}{s}\n{seq}\n'.format(
                    p = pdb_entry,
                    c = chain,
                    d = ' ' + protein_description + ';' if protein_description else '',
                    k = ' ' + keywords + ';' if keywords else '',
                    h = ' HET: ' + compounds + ';' if compounds else '',
                    r = ' NMR ',
                    o = organism if organism else 'N/A',
                    s = ' SCOP: ' + domains if domains else '',
                    seq = '\n'.join(textwrap.wrap(seq, 80)))
                pdb_filter += '{p}_{c}\t{r}\t{f}\t{comp:.3f}\t{m}\n'.format(
                    p = pdb_entry,
                    c = chain,
                    m = method,
                    r = 'N/A',
                    f = round(r_free, 3) if r_free else 'N/A',
                    comp = round(chain_ratios[chain], 3))
            else:
                fasta_entry += '>{p}_{c}{d}{k}{h}{r}{{{o}}}{s}\n{seq}\n'.format(
                    p = pdb_entry,
                    c = chain,
                    d = ' ' + protein_description + ';' if protein_description else '',
                    k = ' ' + keywords + ';' if keywords else '',
                    h = ' HET: ' + compounds + ';' if compounds else '',
                    r = ' ' + str(resolution) + 'A ' if resolution else '',
                    o = organism if organism else 'N/A',
                    s = ' SCOP: ' + domains if domains else '',
                    seq = '\n'.join(textwrap.wrap(seq, 80)))
            # <pdb_entry>_<chain>\t<resolution>\t<r_free>\t<completness>
                pdb_filter += '{p}_{c}\t{r}\t{f}\t{comp:.3f}\t{m}\n'.format(
                    p = pdb_entry,
                    c = chain,
                    m = method,
                    r = round(resolution, 3) if resolution else 'N/A',
                    f = round(r_free, 3) if r_free else 'N/A',
                    comp = round(chain_ratios[chain], 3))

    return (fasta_entry, pdb_filter)

def parse_scop(scop_file):
    scop = defaultdict(set)

    with open(scop_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            
            scop_id, pdb_code, chain, scop_num = line.split('\t')[0:4]
            chain = chain.split(':')[0]
            entry = pdb_code.upper() + '_' + chain
            
            scop[entry].add(scop_num)

    return scop

def wrapper_function(paths):
    
    in_file = paths[0]
    out_file = paths[1]
    
    cif2fasta = CIF2FASTA(in_file)
    if cif2fasta.is_valid():
        fasta_entry = create_fasta_entry2(cif2fasta)
    else:
        print("Warning: Could not read %".format(in_file), file=sys.stderr)
        fasta_entry = None

    return fasta_entry

def write_to_file(line_list, fname, pdb_filter):
    """
    Input: A list containing all lines that have to be saved to the file (line_list).
    A filename (str) is the name of the file that is created.
    Output: A file containing the lines in line_list.
    """
    if pdb_filter:
        fasta_file = open(fname, 'w')
        pdb_filter = open(pdb_filter, 'w')

        for line in line_list:
            if line != None:
                fasta_file.write(line[0])
                pdb_filter.write(line[1])

        fasta_file.close()
        pdb_filter.close()
    else:
        fasta_file = open(fname, 'w')

        for line in line_list:
            if line != None:
                fasta_file.write(line[0])

def opt():
    # Initiate a OptionParser Class
    usage = "usage: cif2fasta.py -i cif_folder -o *.fasta -c num_cores -v"
    description = "cif2fasta.py takes a folder that contains cif files as input and outputs their sequences into fasta file."
    parser = OptionParser(usage = usage, description = description)
    # Call add_options to the parser
    parser.add_option("-i", help = "input cif folder.", dest = "input_files", metavar = "DIR")
    parser.add_option("-o", help = "output fasta file.", dest = "output_files", metavar = "FILE")  
    parser.add_option("-p", help = "output PDB filter file (optional).", dest= "pdb_filter", default = False, metavar = "FILE")   
    parser.add_option("-s", help = "SCOP annotation.", dest = "scop", default = False, metavar = "FILE")   
    parser.add_option('-c', help = "number of cores (default = 1).", dest = "cores", type = int, default = 1, metavar = "INT")
    parser.add_option('-l', help = "Remove chains with a length < X (default = 30).", dest = "seq_len", type = int, default = 30, metavar = "INT")
    parser.add_option('-v', help = 'Verbose Mode (quiet = 0, full verbosity = 2).', dest = 'bool', default = 0, type = int, metavar = "INT")

    return parser
 

def main():
    parser = opt()
    (options, argv) = parser.parse_args()    
    
    global DEBUG_MODE
    
    if options.bool:
        DEBUG_MODE = options.bool

    global SCOP_LIBRARY
    
    if options.scop:
        SCOP_LIBRARY = parse_scop(options.scop)

    global MIN_SEQ_LEN

    MIN_SEQ_LEN = options.seq_len

    paths = get_paths(options.input_files, options.output_files)
    print ("Found: " + str(len(paths)) + " files.")

    if options.cores > 1:
        pool = Pool(options.cores)
        fastas = pool.map(wrapper_function, paths)
    else:
        fastas = map(wrapper_function, paths)

    write_to_file(fastas, options.output_files, options.pdb_filter)

if __name__ == "__main__":
    main()

