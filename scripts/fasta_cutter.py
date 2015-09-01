#!/usr/bin/env python

from optparse import OptionParser


def opt():
  parser = OptionParser()
  parser.add_option("-i", dest="input_file",
                    help="Input fasta file", metavar="FILE")
  parser.add_option("-o", dest="output_file", 
                    help="Output fasta file", metavar="FILE")
  parser.add_option("--max_length", dest="max_length", type=int,
                    help="Max length of fasta sequence", metavar="INT")
  parser.add_option("--overlap", dest="overlap", type=int,
                    help="Length of overlap of cutted sequences", metavar="INT")

  parser.set_default(max_length, 14999)
  parser.set_default(overlap, 100)
  
  return parser


def cut_sequence(header, sequence, max_length, overlap):
  subsequences = []
  if len(sequence) == 0:
    return subsequences
  
  if len(sequence) <= max_length:
    subsequences.append((header, sequence))
    return subsequences
  
  total_length = len(sequence)
  number = 0
  offset = 0
  while offset < total_length:
    subsequence = sequence[max(0, offset - overlap) : min(total_length, offset - overlap + max_length)]
    offset = offset - overlap + max_length

    number += 1
    new_header = adjust_header(header, number)

    subsequences.append((new_header, subsequence))
    
  return subsequences


def adjust_header(header, number):
  header_tokens = header.split()
  name = header_tokens[0]
  name = adjust_name(name, number)
  header_tokens[0] = name
  return " ".join(header_tokens)

  
def adjust_name(name, number):
  short_name = name
  
  if name.find("|") != -1:
    name_tokens = name.split("|")
    short_name = name_tokens[1]
    short_name = short_name+"_"+str(number)
    name_tokens[1] = short_name
    return "|".join(name_tokens)
  else:
    return short_name+"_"+str(number)

  
def write_parts(out, parts):
  for part in parts:
    out.write(part[0])
    out.write("\n")
    out.write(part[1])
    out.write("\n")


def cut_sequences(input_file, max_length, overlap, output_file):
  out = open(output_file, "w")
  
  with open(input_file, "r") as fh:
    sequence = ""
    header = ""
    
    for line in fh:
      line = line.strip()
      if line[0] == ">":
        subsequences = cut_sequence(header, sequence, max_length, overlap)
        write_parts(out, subsequences)
        
        sequence = ""
        header = line
      else:
        sequence += line
        
    subsequences = cut_sequence(header, sequence, max_length, overlap)
    write_parts(out, subsequences)
    
  out.close()


def main():
  parser = opt()
  (options, argv) = parser.parse_args()
  
  cut_sequences(options.input_file, options.max_length, options.overlap, options.output_file)    
    

if __name__ == "__main__":
  main()
