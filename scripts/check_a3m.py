#!/usr/bin/env python

from a3m import A3M_Container
import fileinput
import sys


def check_a3m(filename):
  a3m = A3M_Container()
  
  if(filename.lower() == "stdin"):
    fh = fileinput.input()
  else:
    fh = open(filename, "r")
    
  try:
    a3m.read_a3m(fh)
  except A3MFormatError as e:
    print(e)
    exit(1)
  

def main():
  filename = sys.argv[1]
  check_a3m(filename)


if __name__ == "__main__":
  main()

