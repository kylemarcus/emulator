import re
import sys

if __name__ == '__main__':

  fp = sys.stdin
  line = fp.readline()

  my_str = re.split('\,|\[|\]|\n',line)
#  print my_str
  flag = 0
  base = ''

  for pos in my_str:

#    if pos == '\n':
#	sys.exit(1)

    if pos == '':
	continue

    reg = re.compile('^[a-z].*')
    if reg.match(pos):
        if flag != 0:
	  base = pos
	  flag = 0
	else:
	  flag = 0
	  if base != '':
	    print base
	  base = pos

    elif pos.find('-') >= 0:
	my_range = re.search('(\d+)-(\d+)',pos)
	start = int(my_range.group(1))
	stop = int(my_range.group(2))
	for num in range(start,stop+1):
	    print base+'%02d'%num
	flag = 1

    else:
	print base + pos
	flag = 1

  if flag == 0 and base != '':
	print base
