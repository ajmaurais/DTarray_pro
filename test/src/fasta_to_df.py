
import sys
import argparse
import re

def main():
	parser = argparse.ArgumentParser(description = 'Extract Uniprot IDs from a fasta file.')
	parser.add_argument('fasta_path')
	
	args = parser.parse_args()
	
	uniprot_id_re = '^[>](sp|tr)\\|([A-Z0-9-]+)\\|'
	inF = open(args.fasta_path, 'r')
	lines = inF.readlines()
	ids = list()
	sequences = list()
	nLines = len(lines)
	i = 0
	while i < nLines:
		match = re.search(uniprot_id_re, lines[i])
		if(match):
			ids.append(match.group(2))
			sequences.append(lines[i+1].strip())
			i += 2
			continue
		else:
			i += 1
	
	outF = open('sequences.tsv', 'w')
	outF.write('id\tseq\n')
	for i, s in zip(ids, sequences):
		outF.write('{}\t{}\n'.format(i, s))
	
if __name__ == "__main__":
	main()
