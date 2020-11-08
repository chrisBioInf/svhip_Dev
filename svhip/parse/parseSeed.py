import sys

filename = sys.argv[1]
arr = []

with open(filename, 'r') as seed:
	for line in seed.readlines():
		block = [None, None]
		block[0] = line[:25]
		block[1] = line[26:]
		arr.append(block)

with open(filename + '.fasta', 'w') as outf:
	print('H')
	for block in arr:
		outf.write('>' + block[0].replace(' ','') + '\n')
		outf.write(block[1].replace(' ',''))

print('Finished.')
