#!/usr/bin/env python

import getopt
import sys
from collections import namedtuple

hhr_alignment = namedtuple('hhr_alignment', ['query_id', 'query_length', 'query_neff', 'template_id', 'template_info', 'template_neff', 'query_ali', 'template_ali', 'start', 'end', 'probability', 'evalue', 'score', 'aligned_cols', 'identity', 'similarity', 'sum_probs'])

def get_sequence_name(header):
  name = header.replace(">", "").split()[0]
  return name

def read_a3m(lines):
  header = "";
  seq = "";
  
  l = []
  
  for line in lines:
    if(len(line) > 0 and line[0] == '#'):
      continue
    if(len(line) > 0 and line[0] == '>'):
      if(len(header) != 0):
        l.append((header, seq))
      header = line
      seq = ""
    else:
      seq += line.rstrip()

  if(len(header) != 0):
    l.append((header, seq))

  return l

def parse_result(lines):
  results = []
	
  query_id = ""
  query_length = 0
  query_neff = 0
  query_seq = ""
  template_id = ""
  template_seq = ""
  template_info = ""
  query_start = 10**6
  query_end = 0
  template_start = 10**6
  template_end = 0
  probability = 0
  evalue = 0
  score = 0
  identity = 0
  similarity = 0
  template_neff = 0
  sum_probs = 0
  aligned_cols = 0
  
  skipped_ali_tags = ["ss_dssp", "ss_pred", "Consensus"]
	
  isAlignmentSection = False
	
  for line in lines:
    if(line.startswith("Query")):
      query_id = line.split()[1]
    elif(line.startswith("Match_columns")):
      query_length = int(line.split()[1])
    elif(line.startswith("Neff")):
      query_neff = float(line.split()[1])
    elif(isAlignmentSection and (line.startswith("No") or line.startswith("Done!"))):
      if query_start != 10**6:
      	result = hhr_alignment(query_id, query_length, query_neff, template_id, template_info, template_neff, query_seq, template_seq, (query_start, template_start), (query_end, template_end), probability, evalue, score, aligned_cols, identity, similarity, sum_probs)
      	results.append(result)
      template_id = ""
      template_info = ""
      query_seq = ""
      template_seq = ""
      
      query_start = 10**6
      query_end = 0
      template_start = 10**6
      template_end = 0
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
      isAlignmentSection = True
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
      	e = sys.exc_info()[0]
      	print (tokens[2])
      	print("Error while converting token 2 (query).")

      query_start = min(query_start, token_2)

      try:
      	token_4 = tokens[4].replace("(", "").replace(")", "")
      	token_4 = int(token_4)
      except:
      	e = sys.exc_info()[0]
      	print (tokens[4])
      	print ("Error while converting token 4 (query).")

      query_end = max(query_end, token_4)

      query_seq += tokens[3]
    elif(line.startswith("T")):
      tokens = line.split()
      if(tokens[1] in skipped_ali_tags):
        continue
      template_seq += tokens[3]

      try:
      	token_2 = tokens[2].replace("(", "").replace(")", "")
      	token_2 = int(token_2)
      except:
      	e = sys.exc_info()[0]
      	print (tokens[2])
      	print ("Error while converting token 2 (template).")

      template_start = min(template_start, token_2)

      try:
      	token_4 = tokens[4].replace("(", "").replace(")", "")
      	token_4 = int(token_4)
      except:
      	e = sys.exc_info()[0]
      	print (tokens[4])
      	print ("Error while converting token 4 (template).")

      template_end = max(template_end, token_4)
      
  if(len(template_id) > 0 and query_start != 10**6):
    result = hhr_alignment(query_id, query_length, query_neff, template_id, template_info, template_neff, query_seq, template_seq, (query_start, template_start), (query_end, template_end), probability, evalue, score, aligned_cols, identity, similarity, sum_probs)
    results.append(result)
			
  return results


def read_result(inputFile):
	with open(inputFile) as fh:
		lines = fh.readlines()
		return parse_result(lines)

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:")
	except getopt.GetoptError as err:
		print(str(err)) # will print something like "option -a not recognized"
		sys.exit(2)

	for o, a in opts:
		if o in ("-i"):
			input_file = a

	counter = 0
	for result in read_result(input_file):
		print("Alignment "+str(counter)+"\t evalue: "+str(result.evalue)+"\t probability: "+str(result.probability))
		print(result.query_id+"\t"+str(result.start[0])+"\t"+result.query_ali+"\t"+str(result.end[0]))
		print(result.template_id+"\t"+str(result.start[1])+"\t"+result.template_ali+"\t"+str(result.end[1]))
		counter += 1

if __name__ == "__main__":
	main()

