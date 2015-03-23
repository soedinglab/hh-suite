#!/usr/bin/env python

from a3m import A3M_Container
import fileinput
import sys


def main():
  filename = sys.argv[1]

  a3m = A3M_Container()
  
  if(filename.lower() == "stdin"):
    fh = fileinput.input()
  else:
    fh = open(filename, "r")
    
  try:
    a3m.read_a3m(fh)
    size = a3m.get_number_sequences()
    print(size)
  except A3MFormatError as e:
    sys.stderr.write()
    exit(1)


if __name__ == "__main__":
  main()

