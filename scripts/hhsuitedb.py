#!/usr/bin/env python

"""
    hhbsuite.py
    Creates HH-suite database files from A3M and HHM files 
    Usage: Usage: python hhsuite_db.py -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]

    HHsuite version 3.0.0 (15-03-2015)

    Reference: 
    Remmert M., Biegert A., Hauser A., and Soding J.
    HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
    Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

    (C) Johannes Soeding, 2012

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de
"""


from optparse import OptionParser
from glob import glob
from subprocess import check_call
import ffindex
import a3m
import tempfile
import os
import sys
import shutil


def write_glob_to_file(glob_expr, filename):
  fh = open(filename, "w")
  for f in glob(glob_expr):
    fh.write(f+"\n")
  fh.close()


def write_set_to_file(set, filename):
  fh = open(filename, "w")
  for f in set:
    fh.write(f+"\n")
  fh.close()


def remove_files_from_index(files_file, database_index_path):
  if(os.path.exists(database_index_path) and os.path.getsize(database_index_path) > 0):
    check_call(" ".join(["ffindex_modify -us -f ", files_file, database_index_path]), shell=True)


def build_ffindex_database(files_file, data_path, index_path):
  check_call(" ".join(["ffindex_build","-s", data_path, index_path, "-f", files_file]), shell=True)


def optimize_database(database_data_path, database_index_path):
  print(os.path.getsize(database_data_path))
  if os.path.exists(database_data_path) and os.path.exists(database_index_path) and os.path.getsize(database_data_path) > 0 and os.path.getsize(database_index_path) > 0: 
    check_call(" ".join(["ffindex_build", "-as", database_data_path+".optimized", database_index_path+".optimized", "-d", database_data_path, "-i", database_index_path]), shell=True)
    shutil.move(database_data_path+".optimized", database_data_path)
    shutil.move(database_index_path+".optimized", database_index_path)


def get_large_a3ms(a3m_base_path):
  entries = ffindex.read_index(a3m_base_path+"_a3m.ffindex")
  data = ffindex.read_data(a3m_base_path+"_a3m.ffindex")
  
  large_alignments = set()
  for entry in entries:
    lines = ffindex.read_lines(entry, data)
    alignment = a3m.A3M_Container()
    try:
      alignment.read_a3m_from_lines(lines)
      
      if alignment.get_number_sequences() > 50:
        large_alignments.add(entry.name)
    except:
      sys.stderr.write("Warning: A3M "+entry.name+" is corrupted!\n")

  return large_alignments


def calculate_hhm(threads, a3m_base_path, hhm_base_path):
  tmp_dir = tempfile.mkdtemp()
  
  large_a3ms = get_large_a3ms(a3m_base_path)
  large_a3m_index = os.path.join(tmp_dir, "large.ffindex")

  a3m_index = read_ffindex(a3m_base_path+".ffindex")
  write_subset_index(a3m_index, large_a3ms, large_a3m_index)
  
  #TODO at the moment we are building all hmms...
  check_call(" ".join(["mpirun", "-np", threads, "ffindex_apply_mpi", a3m_base_path+".ffdata", large_a3m_index, "-d", hhm_base_path+".ffdata", "-i", hhm_base_path+".ffindex", "--", "hhmake", "-v", str(0), "-i", "stdin", "-o" ,"stdout"]), shell=True)


def calculate_cs219(threads, a3m_db_base, cs_db_base):
  hhlib_environment = os.environ['HHLIB']
  if(hhlib_environment):
    #TODO: check
    check_call(" ".join(["cstranslate", "-A ", os.path.join(hhlib_environment, "data/cs219.lib"), "-D", os.path.join(hhlib_environment, "data/context_data.lib"), "-x 0.3 -c 4 --ffindex", "-i", a3m_db_base, "-o", cs_db_base, "-I a3m -b"]), env=dict(os.environ, OMP_NUM_THREADS=threads), shell=True)
  else:
    sys.err("ERROR: HHLIB environment variable not set! See manual!\n")
    exit(1)
    #TODO throw error


def merge_databases(source_data_path, source_index_path, dest_data_path, dest_index_path):
  check_call(["ffindex_build", "-as", "-d", source_data_path, "-i", source_index_path, dest_data_path, dest_index_path])


def sort_database(data_path, index_path):
  check_call(["ffindex_build", "-as", data_path, index_path])


def read_ffindex(file):
  fh = open(file, "r")
  
  index = []
  for line in fh:
    index.append(line.rstrip().split())
  
  fh.close()
  
  return index


def write_subset_index(index, subset, output_file):
  fh = open(output_file, "w")
  for entry in index:
    if entry[0] in subset:
      fh.write("{name:.64}\t{offset}\t{length}\n".format(name=entry[0], offset=entry[1], length=entry[2]))
  fh.close()


def is_sorted(index):
  for i in range(len(index)-1):
    if(index[i][0] > index[i+1][0]):
      return False
  return True


def get_duplicates(index):
  duplicates = set()
  entries = set()
  
  for entry in index:
    if(entry[0] in entries):
      duplicates.add(entry[0])
      
    entries.add(entry[0])

  return duplicates


def get_missing(index, cmp_index):
  missing = set()
  
  index_set = set()
  for entry in index:
    index_set.add(entry[0])
    
  cmp_index_set = set()
  for entry in cmp_index:
    cmp_index_set.add(entry[0])
  
  for key in cmp_index_set:
    if key not in index_set:
      missing.add(key)
  
  return missing


def get_overhead(index, cmp_index):
  overhead = set()
  
  index_set = set()
  for entry in index:
    index_set.add(entry[0])
    
  cmp_index_set = set()
  for entry in cmp_index:
    cmp_index_set.add(entry[0])
  
  for key in index_set:
    if key not in cmp_index_set:
      overhead.add(key)
  
  return overhead


def handle_duplicates(suffix, calculate, threads, db_basename, force_mode):
  index = read_ffindex(db_basename+"_"+suffix+".ffindex")
  duplicates = get_duplicates(index)
  
  if(suffix == "a3m" and len(duplicates) > 0):
    sys.stderr.write("ERROR: "+db_basename+"_a3m.ffindex contains duplicates!\n")
    sys.stderr.write("ERROR: Your database is broken!\n")
    exit(1);
    
  if(len(duplicates) == 0):
    return
  
  for duplicate in duplicates:
    sys.stderr.write("WARNING: "+db_basename+"_"+suffix+".ffindex contains duplicate "+duplicate+"!\n")

  
  if force_mode:
    tmp_dir = tempfile.mkdtemp()
    
    try:
      sys.stderr.write("WARNING: remove duplicates and recalculate them from their corresponding a3m's!\n")
      a3m_index = read_ffindex(db_basename+"_a3m.ffindex")
      
      duplicates_a3m_base = os.path.join(tmp_dir, "duplicates_a3m")
      duplicates_a3m_index_file = os.path.join(tmp_dir, "duplicates_a3m.ffindex")
      duplicates_a3m_data_file = os.path.join(tmp_dir, "duplicates_a3m.ffdata")
      write_subset_index(a3m_index, duplicates, duplicates_a3m_index_file)
      os.symlink(db_basename+"_a3m.ffdata", duplicates_a3m_data_file)
      
      duplicates_new_base = os.path.join(tmp_dir, "duplicates_new")
      duplicates_new_index_file = os.path.join(tmp_dir, "duplicates_new.ffindex")
      duplicates_new_data_file = os.path.join(tmp_dir, "duplicates_new.ffdata")
    
      duplicates_index_file = os.path.join(tmp_dir, "duplicates.dat")
      write_set_to_file(duplicates, duplicates_index_file)
      remove_files_from_index(duplicates_index_file, db_basename+"_"+suffix+".ffindex")
      
      calculate(threads, duplicates_a3m_base, duplicates_new_base)
      sort_database(duplicates_new_data_file, duplicates_new_index_file)
      
      merge_databases(duplicates_new_data_file, duplicates_new_index_file, 
        db_basename+"_cs219.ffdata", db_basename+"_cs219.ffindex")
      sort_database(db_basename+"_cs219.ffdata", db_basename+"_cs219.ffindex")
      optimize_database(db_basename+"_cs219.ffdata", db_basename+"_cs219.ffindex")
    finally:
      shutil.rmtree(tmp_dir)
  else:
    sys.stderr.write("You may try to use the option --force to fix the database!\n")

    
def handle_unsorted(suffix, db_basename, force_mode):
  index = read_ffindex(db_basename+"_"+suffix+".ffindex")
  if(not is_sorted(index)):
    sys.stderr.write("Index "+db_basename+"_"+suffix+".ffindex is unsorted!\n")
    if force_mode:
      sys.stderr.write("Try to sort unsorted index "+db_basename+"_"+suffix+".ffindex!\n")
      sort_database(db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
    else:
      sys.stderr.write("You may try to use the option --force to fix the database!\n")

    
def handle_missing(suffix, calculate, db_basename, threads, force_mode):
  index = read_ffindex(db_basename+"_"+suffix+".ffindex")
  a3m_index = read_ffindex(db_basename+"_a3m.ffindex")
  
  missing = get_missing(index, a3m_index)
  
  if(len(missing) == 0):
    return
  
  for f in missing:
    sys.stderr.write("WARNING: Missing entry "+f+" in "+db_basename+"_"+suffix+".ff{data,index}!\n")
  
  if force_mode:
    sys.stderr.write("WARNING: Try to calculate missing entries!\n")
    tmp_dir = tempfile.mkdtemp()
    
    try:
      missing_base = os.path.join(tmp_dir, "missing_"+suffix)
      missing_index_file = os.path.join(tmp_dir, "missing_"+suffix+".ffindex")
      missing_data_file = os.path.join(tmp_dir, "missing_"+suffix+".ffdata")
    
      missing_a3m_base = os.path.join(tmp_dir, "missing_a3m")
      missing_a3m_index_file = os.path.join(tmp_dir, "missing_a3m.ffindex")
      missing_a3m_data_file = os.path.join(tmp_dir, "missing_a3m.ffdata")
      write_subset_index(a3m_index, missing, missing_a3m_index_file)
      os.symlink(os.path.abspath(db_basename+"_a3m.ffdata"), missing_a3m_data_file)
      
      calculate(threads, missing_a3m_base, missing_base)
    
      merge_databases(missing_data_file, missing_index_file, db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
      optimize_database(db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
    finally:
      shutil.rmtree(tmp_dir)
  else:
    sys.stderr.write("You may try to use the option --force to fix the database!\n")


def handle_overhead(suffix, db_basename, force_mode):
  index = read_ffindex(db_basename+"_"+suffix+".ffindex")
  a3m_index = read_ffindex(db_basename+"_a3m.ffindex")

  overhead = get_overhead(index, a3m_index)

  #delete overhead cs219 files
  if(len(overhead) == 0):
    return

  for f in overhead:
    sys.stderr.write("WARNING: Entry "+f+" from "+db_basename+"_"+suffix+".ff{data,index} has no corresponding entry in the a3m database!\n")

  if force_mode:
    sys.stderr.write("WARNING: Try to fix overhead entries!\n")
    tmp_dir = tempfile.mkdtemp()
    
    try:
      index_file = os.path.join(tmp_dir, "to_delete.dat")
      write_set_to_file(overhead, index_file)
      remove_files_from_index(index_file, db_basename+"_"+suffix+".ffindex")
      optimize_database(db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
    finally:
      shutil.rmtree(tmp_dir)
  else:
    sys.stderr.write("You may try to use the option --force to fix the database!\n")


def check_a3m_format(db_basename, force_mode):
  entries = ffindex.read_index(db_basename+"_a3m.ffindex")
  data = ffindex.read_data(db_basename+"_a3m.ffdata")
  
  corrupted_alignments = set()
  for entry in entries:
    lines = ffindex.read_lines(entry, data)
    alignment = a3m.A3M_Container()
    try:
      alignment.read_a3m_from_lines(lines)
    except:
      corrupted_alignments.add(entry.name)
      sys.stderr.write("Warning: A3M "+entry.name+" is corrupted!\n")
  
  if len(corrupted_alignments) == 0:
    return
  
  if force_mode:
    tmp_dir = tempfile.mkdtemp()
    
    try:
      sys.stderr.write("WARNING: remove corrupted a3m's!\n")
      
      corrupted_index_file = os.path.join(tmp_dir, "corrupted.dat")
      write_set_to_file(corrupted_alignments, corrupted_index_file)
      
      for suffix in ["a3m", "cs219", "hhm"]:
        remove_files_from_index(corrupted_index_file, db_basename+"_"+suffix+".ffindex")
        sort_database(db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
        optimize_database(db_basename+"_"+suffix+".ffdata", db_basename+"_"+suffix+".ffindex")
    finally:
      shutil.rmtree(tmp_dir)
  else:
    sys.stderr.write("You may try to use the option --force to fix the database!\n")
    

def check_database(output_basename, threads, force_mode):
  if not os.path.exists(output_basename+"_a3m.ffindex") or not os.path.exists(output_basename+"_a3m.ffdata"):
    sys.stderr.write("Error: No a3m database found!")
    exit(1)
    
  if not os.path.exists(output_basename+"_cs219.ffindex") or os.path.exists(output_basename+"_cs219.ffdata"):
    open(output_basename+"_cs219.ffindex", 'a').close()
    open(output_basename+"_cs219.ffdata", 'a').close()
    
  if not os.path.exists(output_basename+"_cs219.ffindex") or os.path.exists(output_basename+"_cs219.ffdata"):
    open(output_basename+"_hhm.ffindex", 'a').close()
    open(output_basename+"_hhm.ffdata", 'a').close()

  check_a3m_format(output_basename, force_mode)
  #TODO check hhm's...
  #TODO check cs219...
  
  handle_unsorted("a3m", output_basename, force_mode)
  handle_unsorted("hhm", output_basename, force_mode)
  handle_unsorted("cs219", output_basename, force_mode)
  
  handle_duplicates("a3m", calculate_hhm, threads, output_basename, force_mode)
  handle_duplicates("hhm", calculate_hhm, threads, output_basename, force_mode)
  handle_duplicates("cs219", calculate_cs219, threads, output_basename, force_mode)

  handle_missing("cs219", calculate_cs219, output_basename, threads, force_mode)
  
  handle_overhead("hhm", output_basename, force_mode)
  handle_overhead("cs219", output_basename, force_mode)


def add_new_files(globular_expression, suffix, output_basename):
  tmp_dir = tempfile.mkdtemp()

  files = set(glob(globular_expression))
  
  if(len(files) == 0):
    return

  try:
    files_index = os.path.join(tmp_dir, "files.dat")
    write_set_to_file(files, files_index)
    
    new_base = os.path.join(tmp_dir, "new")
    new_index_file = new_base + ".ffindex"
    new_data_file = new_base + ".ffdata"
    build_ffindex_database(files_index, new_data_file, new_index_file)
  
    output_index_file = output_basename+"_"+suffix+".ffindex"
    output_data_file = output_basename+"_"+suffix+".ffdata"
    remove_files_from_index(files_index, output_index_file)
    merge_databases(new_data_file, new_index_file, output_data_file, output_index_file)
    optimize_database(output_data_file, output_index_file)
    
    if(suffix == "a3m"):
      for other_suffix in ["hhm", "cs219"]:
        output_index_file = output_basename+"_"+other_suffix+".ffindex"
        output_data_file = output_basename+"_"+other_suffix+".ffdata"
        if os.path.exists(output_index_file) and os.path.exists(output_data_file) and os.path.getsize(output_data_file) > 0:
          remove_files_from_index(files_index, output_index_file)
          optimize_database(output_data_file, output_index_file)
      
  finally:
    shutil.rmtree(tmp_dir)


def opt():
  parser = OptionParser()
  parser.add_option("--ia3m", dest="a3m_files_glob",
    help="Glob for a3m files", metavar="<GLOB>")
  parser.add_option("--ics219", dest="cs_files_glob",
    help="Glob for cs219 files", metavar="<GLOB>")
  parser.add_option("--ihhm", dest="hhm_files_glob",
    help="Glob for hhm files", metavar="<GLOB>")
  parser.add_option("-o", metavar="FILE", dest="output_basename", 
    help="Output hhsuite database basename")
  parser.add_option("--cpu", metavar="INT", dest="nr_cores", 
    help="Number of threads that shall used")
  parser.add_option("--force", 
    action="store_true", dest="force_mode", default=False,
    help="Try to fix problems with the database")


  return parser
  

def check_options(options, parser):
  if not options.nr_cores:
    sys.stderr.write("Please use --cpu!\n")
    parser.print_help()
    sys.exit(1)
    
  if not options.output_basename:
    sys.stderr.write("Please use -o!\n")
    parser.print_help()
    sys.exit(1)


def main():
  parser = opt()
  (options, args) = parser.parse_args()
  check_options(options, parser)
  
  #Important to do a3m's first... deleting out of date hhm's and cs219's
  if(options.a3m_files_glob):
    add_new_files(options.a3m_files_glob, "a3m", options.output_basename)
    
  if(options.hhm_files_glob):
    add_new_files(options.hhm_files_glob, "hhm", options.output_basename)

  if(options.cs_files_glob):
    add_new_files(options.cs_files_glob, "cs219", options.output_basename)

  check_database(options.output_basename, options.nr_cores, options.force_mode)
  

if __name__ == "__main__":
  main()
  
