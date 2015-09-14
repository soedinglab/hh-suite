#!/usr/bin/env python

from a3m import A3M_Container
from a3m import A3MFormatError
import sys


def main():
  filename = sys.argv[1]

  a3m = A3M_Container()
  
  if(filename.lower() == "stdin"):
    fh = sys.stdin
  else:
    fh = open(filename, "r")

  threshold = int(sys.argv[2])
    
  try:
    a3m.read_a3m(fh)
    size = a3m.get_number_sequences()
    if size >= threshold:
      exit(0)
    else:
      exit(1)
  except A3MFormatError as e:
    sys.stderr.write(e)
    exit(1)


if __name__ == "__main__":
  main()

