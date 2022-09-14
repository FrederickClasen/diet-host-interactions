import sys

fin = open(sys.argv[1],'r')

for line in fin:
	line = line.strip()
	line = line.replace('=>','').strip()
	line = line.replace(' ','').strip()
	print(' => '+line)
