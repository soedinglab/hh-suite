#!/usr/bin/env python

"""
Created in Jun 2016

@author: Harald Voehringer
"""


from itertools import groupby
from collections import namedtuple, defaultdict
import textwrap

DEBUG = False
FASTA = namedtuple('FASTA', ['header', 'seq'])
PDBDATA = namedtuple('PDBDATA', ['entry', 'res', 'rfr', 'comp', 'met'])

def as_pairs(fasta):
    """ Reads in fasta headers as defined """

    for header, group in groupby(fasta, lambda x: x[0] == '>'):
        if header:
            line = next(group)
            identifier = line.split(' ')[0].split('>')[1]
            header_line = str(line.strip())

        else:
            sequence = ''.join(line.strip() for line in group).upper()
            data = FASTA(header_line, sequence)

            yield identifier, data
 
def read_fasta(fasta_file):
    """ Reads in fasta sequences with help of the as_pairs function."""
    fasta_dic = dict()
    duplicate = 0
    
    with open(fasta_file) as fasta:
        for title, sequence in as_pairs(fasta):
            if title not in fasta_dic:
                fasta_dic[title] = sequence
            else:
                duplicate +=1
    
    if duplicate != 0:
        print ('! Warning found duplicate {num} fasta.'.format(
            num = duplicate))
        
    return fasta_dic

def read_fasta_annotations(fname):
    """ Reads the information about resolution, R-free and completness which 
    be outputed by cif2fasta."""

    annotations = dict()

    with open(fname) as fh:
        for line in fh:
            if len(line) > 0 and line[0] == '#':
                continue

            identifier, res, r_free, comp, method = line.strip().split('\t')

            try:
                res = float(res)
            except ValueError:
                res = None

            try:
                r_free = float(r_free)
            except ValueError:
                r_free = None

            try:
                comp = float(comp)
            except ValueError:
                comp = None

            annotations[identifier] = PDBDATA(identifier, res, r_free, comp, method)

    return annotations

def read_cluster(cluster_file):
    """ Reads in the clusters (this is the output of MMseqs). """
    cluster = defaultdict(set)

    with open(cluster_file) as fh:
        for line in fh:
            exemplar, node = line.split()
            
            if node in cluster[exemplar]:
                raise Exception('! Please check the clustering procedure: {n} was found twice in cluster {e}.'.format(
                    n = node,
                    e = exemplar))
            else:
                cluster[exemplar].add(node)

    return cluster

def read_pdblist(in_file):

    pdb_list = set()
    
    with open(in_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            strip_line = line.strip()
            pdb_code = strip_line.split('_')[0]

            # this is a very very basic check
            if len(pdb_code) != 4:
                print ('! Warning: {line} seems to be an incorrect identifer. Skipping it.'.format(
                    line = strip_line))
                continue
            
            pdb = strip_line.upper()
            pdb_list.add(pdb)

    return pdb_list

def select_sequences(clusters, annotations):

    selected_sequences = set()

    for idx, cluster in enumerate(clusters):
        
        nodes = [ annotations[entry] for entry in clusters[cluster] ]
        
        if DEBUG:
            print ('Processing Cluster {c} ({i}): {m}'.format(
                c = cluster,
                i = idx, 
                m = ', '.join([x.entry for x in nodes])))

        # select the best entries of the cluster
        best_res = float('inf')
        best_entry_res = None

        best_rfr = float('inf')
        best_entry_rfr = None

        best_comp = -float('inf')
        best_entry_comp = None

        # iterate through each entry in nodes while selecting the representative sequence
        for node in nodes:
            
            if (node.res is not None) and (node.res < best_res):
                best_res = node.res
                best_entry_res = node.entry

            if (node.rfr is not None) and (node.rfr < best_rfr):
                best_rfr = node.rfr
                best_entry_rfr = node.entry

            if (node.comp is not None) and (node.comp > best_comp):
                best_comp = node.comp
                best_entry_comp = node.entry

        if best_entry_res is not None:
            selected_sequences.add(best_entry_res)
            
            if DEBUG:
                print (' - Selected {n} (best resolution = {r}).'.format(
                    n = best_entry_res,
                    r = best_res))

        if best_entry_rfr is not None:
            selected_sequences.add(best_entry_rfr)
            
            if DEBUG:
                print (' - Selected {n} (best R-free = {r}).'.format(
                    n = best_entry_rfr,
                    r = best_rfr))
        
        if best_entry_comp is not None:
            selected_sequences.add(best_entry_comp)
            
            if DEBUG:
                print (' - Selected {n} (best completness = {r}).'.format(
                    n = best_entry_comp,
                    r = best_comp))    

        if best_entry_res == None and best_entry_rfr == None and best_entry_comp == None:
            print ('! Warning: Did not find any representative entry for cluster {c}.'.format(
                c = cluster))
    
    return selected_sequences

def write_sequences(out_file, fasta_db, selected_sequences):
    """ Writes selected sequences to a fasta file."""
    with open(out_file, 'w') as fh:
        for seq in selected_sequences:
            
            fasta_entry = '{h}\n{seq}\n'.format(
                h = fasta_db[seq].header,
                seq = '\n'.join(textwrap.wrap(fasta_db[seq].seq, 80)))
            
            fh.write(fasta_entry)

def arg():
    import argparse
    description = """
    pdbfilter.py selects from sequence clusters (determined by MMSeqs) the sequences 
    which have the best resolution, R-free factor and/or completness and writes them to a fasta file. 
    """.replace('\n', '')
    
    epilog = '2016 Harald Voehringer'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('fasta', help = 'input fasta file (created by cif2fasta.py)', metavar = 'FILE')    
    parser.add_argument('cluster', help = 'sequence clusters (MMseqs)', metavar = 'FILE')
    parser.add_argument('annotations', help = 'annotations file (created by cif2fasta using the -p flag, contains information about resolution, R-free and completness of sequences).', metavar = 'FILE')
    parser.add_argument('out_file', help = 'output fasta file', metavar = 'FILE')
    parser.add_argument('-i', '--include', help = 'include PDB chains', metavar = 'FILE')
    parser.add_argument('-r', '--remove', help = 'exclude PDB chains', metavar = 'FILE')

    parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose mode')


    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    global DEBUG
    
    if args.verbose:
        DEBUG = True

    fasta = read_fasta(args.fasta)
    clu70 = read_cluster(args.cluster)
    annot = read_fasta_annotations(args.annotations)

    print ("Found {i} clusters.".format(
        i = len(clu70.keys())))

    # choose representative sequences from clusters
    selection = select_sequences(clu70, annot)

    # make sure that pdbs specified in the argument are included
    if args.include:
        to_include = read_pdblist(args.include)
        for pdb in to_include:
            if pdb in fasta.keys():
                if pdb not in selection:
                    if DEBUG:
                        print ('Adding {p}.'.format(
                            p = pdb))
                    selection.add(pdb)
            else:
                print ('! Warning: {p} was not found in input fasta.'.format(
                    p = pdb))

    # removes entries
    if args.remove:
        to_remove = read_pdblist(args.remove)
        for pdb in to_remove:
            if pdb in selection:
                if DEBUG:
                    print ('Removing {p}.'.format(
                        p = pdb))
                selection.remove(pdb)
        
    # write them to file
    write_sequences(args.out_file, fasta, selection)


if __name__ == "__main__":
    main()
