#!/usr/bin/env python3

import argparse
import sys

def parse_arguments():

	parser = argparse.ArgumentParser(description='Compair different TE classification methods and create a consensus.')

	parser.add_argument('-i', '--input_file', help='Path to input file.')
	parser.add_argument('-o', '--output_file', help='Path to output file.')

	return parser.parse_args()

def main():
	args = parse_arguments()

	with open(args.input_file, "r") as f , open(args.output_file, "w") as o:
	
	# Loop over the input lines.
		for line in f:
			line=line.rstrip()
			if line.startswith("TE_ID,EarlGrey,DeepTE,TEsorter"):
				continue
			columns=line.split(",")
			
			if len(columns) !=4:
				sys.exit()
			
			if columns[1].lower() == columns[2].lower() and columns[1].lower() == columns[3].lower():
				print(f'{columns[0]}\t{columns[1]}', file=o)
			
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Gypsy" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}', file=o)

			else:
				print(line)

if __name__ == '__main__':
    main()
