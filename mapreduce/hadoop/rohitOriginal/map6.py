#!/usr/bin/env python
import sys
sys.path.append('.')
import re

if __name__ == '__main__':

 fp = sys.stdin
 for line in fp:
   match = re.findall('[\d\_.e-]+',line)
   try:
     print match[0] + '\t' + match[1] + ' ' + match[2]
   except: continue

#   print line.rstrip('\n')
