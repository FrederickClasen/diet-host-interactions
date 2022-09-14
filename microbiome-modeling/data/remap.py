import sys

fin = open(sys.argv[1],'r')
map_dict = {}
for line in fin:
	line = line.strip()
	line = line.split()
	map_dict[line[1]] = line[0]
fin.close()
print(map_dict)
