#!/usr/bin/env python
import sys
sys.path.append('.')
import re

if __name__ == '__main__':

 fp = sys.stdin
 for line in fp:
   match = re.findall('[\d\_.e-]+',line);
#   print match[0] + "\t" + match[1] + ' ' + match[2]
   one = match[0]
   two = match[1]
   three = match[2]


#   print line.rstrip('\n')
