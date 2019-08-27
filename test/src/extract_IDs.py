
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
	for x in lines:
		match = re.search(uniprot_id_re, x)
		if(match):
			ids.append(match.group(2))
	
	outF = open('ids.tsv', 'w')
	outF.write('id\n')
	for l in ids:
		outF.write('{}\n'.format(l))
	
if __name__ == "__main__":
	main()

