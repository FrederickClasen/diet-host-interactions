import sys

fin = open(sys.argv[1],'r')

for line in fin:
	line = line.strip()
	line = line.split()
	counter = 0
	for item in line:
		try:
			item = float(item)
		except ValueError:
			if counter > 0:
				print(item)
		counter += 1
