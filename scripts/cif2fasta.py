#!/usr/bin/env python

"""
Created on Mon Jun 15 21:49:32 2015

@author: Sagar
"""

import itertools
from sys import argv, exit
import sys
from optparse import OptionParser
from os.path import splitext
from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
from pdbx.reader.PdbxContainers import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from multiprocessing import Pool
import os

import glob
import textwrap

DEBUG_MODE = True


class CIF2FASTA(object):
    
    def __init__(self, cif_path):
        self.cif_path = cif_path
        self.block = self.open_cif()
            
    def open_cif(self):
        """ Assumes a mmCif file and returns a data block used for subsequent procedures. """
        # The "usual" procedure to open a mmCIF with pdbX/mmCIF
        try:
            with open(self.cif_path) as cif_fh:
                data = []
                reader = PdbxReader(cif_fh)
                reader.read(data)
                return data[0]
        except:
            print "! Unexpected error during open_cif:", sys.exc_info()[0], in_file
            return None
                    
 
    def chain_to_seq(self):
        """Extracts the sequence of the cif from entity_poly.pdbx_seq_one_letter_code"""
        
        cif_chain_to_seq = dict()
        non_polypeptide_chains = list()         
        
        try:
            entity_poly = self.block.getObj("entity_poly")
            for row in range(0, entity_poly.getRowCount()):
                if entity_poly.getValue("type", row) == "polypeptide(L)":
                    seq = entity_poly.getValue('pdbx_seq_one_letter_code', row)                                                            
                    parsed_seq = parse_seq(seq) # removes special amino acids and newlines

                    chains = entity_poly.getValue('pdbx_strand_id', row)
                    chain_list = chains.split(',')
            
                    for chain in chain_list:     
                        cif_chain_to_seq[chain] = parsed_seq
                else:
                    chains = entity_poly.getValue('pdbx_strand_id', row)
                    non_polypeptide_chains.append(chains)
        except AttributeError: 
            if DEBUG_MODE:
                print "- Could not extract sequences from entity_poly (" + str(self.cif_path.split("/")[-1]) + ")." 
            return None
        except:
            print "! Unexpected error during chain_to_seq (entity_poly):", sys.exc_info()[0], in_file
            return None 

        if len(cif_chain_to_seq) != 0:
            if DEBUG_MODE:
                print "- Extracted chains of " + str(self.cif_path.split("/")[-1]) + " " + " ".join( str(chain) + " (" + str(len(cif_chain_to_seq[chain])) + ")" for chain in cif_chain_to_seq.keys() ) + "."
                if len(non_polypeptide_chains) != 0:                
                    print "- Following chains were non polypeptide chains " + ", ".join(non_polypeptide_chains) + " no polypeptide chains were found."
            return cif_chain_to_seq
        else:
            if DEBUG_MODE:
                print "- No polypeptide chains were found in (" + str(self.cif_path.split("/")[-1]) + ")." 
            
            return None

    def pdb_entry(self):
        """Extracts the PDB entry information of a cif file."""
        try:
            entry = self.block.getObj("entry")
            entry_id = entry.getValue("id")
            return entry_id.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract id from entry (" + str(self.cif_path.split("/")[-1]) + ")."   
        except:
            print "! Unexpected error during pdb_entry (entry):", sys.exc_info()[0], in_file
            return None 


    def protein_description(self):
        """Extracts the protein description annotated in struct.pdbx_descriptor of the cif file."""
        try:
            # Get struct table which contains the protein description
            struct = self.block.getObj('struct')
            # Get the pdbx description 
            protein_description = struct.getValue("pdbx_descriptor")
            return protein_description.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract pdbx_descriptor from struct table (" + str(self.cif_path.split("/")[-1]) + ")."
            return None
        except:
            print "! Unexpected error during protein_description (struct):", sys.exc_info()[0], in_file
            return None 

    def compounds(self):
        """ Extracts all compounds annotated in the HETATM section of the atom
            struct table if the compound appears at least 10 times and is not water 
            (HOH)."""    
        
        atom_site = self.block.getObj('atom_site')

        compounds = {}
        for row in range(0, atom_site.getRowCount()):
            if atom_site.getValue("group_PDB", row) == "HETATM":

                label_comp_id = atom_site.getValue("label_comp_id", row)
                if label_comp_id not in compounds.keys():
                    compounds[label_comp_id] = 1
                else:
                    compounds[label_comp_id] += 1
    
        filtered_compounds = set()

        for compound in compounds.keys():
            if compounds[compound] >= 10 and compound != "HOH":
                filtered_compounds.add(compound)

        if len(filtered_compounds) == 0:
            return None
        else:
            return ", ".join(filtered_compounds).replace("\n", " ")
    
    def resolution(self):
        """Extracts the resolution of the mmCIF."""
        try:
            reflns = self.block.getObj('reflns')
            # Extract the resolution of the crystal
            resolution = reflns.getValue("d_resolution_high")
            return resolution.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract d_resolution_high from reflns table (" + str(self.cif_path.split("/")[-1]) + ")."
            return None
        except:
            print "! Unexpected error during resolution (reflns):", sys.exc_info()[0], in_file
            return None

    def experimental_method(self):
        """Extracts the experimental method of the mmCIF."""

        try:
            reflns = self.block.getObj('exptl')
            method = reflns.getValue("method")
            return method.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract text from exptl table (" + str(self.cif_path.split("/")[-1]) + ")."
            return None
        except:
            print "! Unexpected error during experimental_method (exptl):", sys.exc_info()[0], in_file
            return None    


    def keywords(self):
        """Extracts the keywords of the mmCIF."""
        try:
            reflns = self.block.getObj("struct_keywords")
            keywords = reflns.getValue("text")
            return keywords.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract text from struct_keywords table (" + str(self.cif_path.split("/")[-1]) + ")."
            return None
        except:
            print "! Unexpected error during keywords (struct_keywords):", sys.exc_info()[0], in_file
            return None    

    def organism(self):
        """Extracts the organism of the mmCIF."""
        try:
            entity_src_nat = self.block.getObj('entity_src_nat')
            organsim_scientific = entity_src_nat.getValue("pdbx_organism_scientific")
            return organsim_scientific.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract from pdbx_organism_scientific from entity_src_gen table (" + str(self.cif_path.split("/")[-1]) + ")."
            pass
        except:
            print "! Unexpected error during get_organsm (pdbx_organism_scientific):", sys.exc_info()[0], in_file
            return None
            
        try:
            entity_src_gen = self.block.getObj("entity_src_gen")
            src_scientific = entity_src_gen.getValue("pdbx_gene_src_scientific_name")
            return src_scientific.replace("\n", " ")
        except AttributeError:
            if DEBUG_MODE:
                print "- Could not extract from pdbx_gene_src_scientific_name from entity_src_gen table (" + str(self.cif_path.split("/")[-1]) + ")."
            return None
        except:
            print "! Unexpected error during organism (pdbx_gene_src_scientific_name):", sys.exc_info()[0], in_file
            return None
# Helper functions

def parse_seq(orginal_seq):
    seq = orginal_seq
    
    while seq.find('(') != -1: 
        start_pos = seq.find('(')
        stop_pos = seq.find(')')
                
        if start_pos == 0:
            seq = 'X' + seq[stop_pos+1:]
        elif stop_pos == len(seq):
            seq = seq[0:start_pos] + 'X' 
        else:
            pre_seq = seq[0:start_pos]
            post_seq = seq[stop_pos+1:]
            seq = pre_seq + 'X' + post_seq

    return seq.replace('\n', '')
    
    
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

    # construct the header with the data given
    header_information = ""

    if protein_description:
        header_information += protein_description + "; "
 
    if keywords:
        header_information += keywords + "; "

    if compounds:
        header_information += "HET: " + compounds + "; "
    
    if resolution:
        header_information += resolution +"A "
    
    elif method:
        header_information += method + " "

    if organism:
        header_information += "{" + organism + "} "

    return header_information.strip()

def create_fasta_entry(cif2fasta):
    """Combines information given in pdb_entry and strand_seq and yields a clean
    fasta file.
    """
    
    pdb_entry = cif2fasta.pdb_entry()
    header = construct_header(cif2fasta)
    chain_to_seq = cif2fasta.chain_to_seq()
    
    fasta_entry = ""
    
    if chain_to_seq and pdb_entry:
        try:
            for chain in chain_to_seq.keys():
                if len(chain_to_seq[chain]) != 0:
                    fasta_entry += ">" + pdb_entry + "_" + chain + " " + header  + "\n" + "\n".join(textwrap.wrap(chain_to_seq[chain], 80)) + "\n"
            return fasta_entry
        except:
            fasta_entry = None
            print "Unexpected error during create_fasta_entry:", sys.exc_info()[0], in_file
            return None
    else:
        return None

            
def wrapper_function(paths):
    
    in_file = paths[0]
    out_file = paths[1]
    
    try:
        cif2fasta = CIF2FASTA(in_file)
        fasta_entry = create_fasta_entry(cif2fasta)
    except:
        fasta_entry = None
        print "Unexpected error:", sys.exc_info()[0], in_file
        pass        

    if fasta_entry:
        with open(out_file, 'w') as ofile:
            ofile.write(fasta_entry)
    
    

def opt():
    # Initiate a OptionParser Class
    parser = OptionParser()
    # Call add_options to the parser
    parser.add_option("-i", dest="input_files",
                      help="Input mmCIF folder.", metavar="GLOB")
    parser.add_option("-o", dest="output_files", 
                      help="Output fasta file", metavar="FILE")   
    parser.add_option('-c', dest="cores", type=int, default=1,
                      help="How many cores should be used?")
    return parser
 

def main():
    parser = opt()
    # parse the the parser object and save the data into options and args
    # opts contains all values received via command line
    # argv contains a list of all commited argument 
    (options, argv) = parser.parse_args()    
    
    paths = get_paths(options.input_files, options.output_files)
    print "Found: " + str(len(paths)) + " files."

    pool = Pool(options.cores)
    pool.map(wrapper_function, paths)

if __name__ == "__main__":
    main()

