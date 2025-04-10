#!/usr/bin/env python3

import argparse
import sys

def parse_arguments():

	parser = argparse.ArgumentParser(description='Compair different TE classification methods and create a consensus.')

	parser.add_argument('-i', '--input_file', help='Path to input file.')
	parser.add_argument('-d', '--deepte_domain', help='Path to input file of DeepTE domain file.')
	parser.add_argument('-o', '--output_file', help='Path to output file.')

	return parser.parse_args()

def main():
	args = parse_arguments()
	deeptedomains = {}

	with open(args.deepte_domain, "r") as d:
		for line in d:
			line=line.rstrip()
			columns=line.split()
			deeptedomains[columns[0]] = columns[1]
			continue

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
			
			#LTRs
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Gypsy" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassI LTR Gypsy" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter', file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Gypsy" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}', file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown', file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassII DNA hAT MITE" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown',file=o)
			
			#TIR
			elif columns[1] == "DNA/PIF-Harbinger" and columns[2] == "ClassII DNA hAT nMITE" and columns[3].startswith("TIR/PIF_Harbinger"):
				print(f'{columns[0]}\t{columns[3]}',file=o)

			#The CACTA
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassII DNA CACTA") and columns[3].startswith("TIR/EnSpm_CACTA/"):#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\t{columns[3]}', file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and 'CACTA' not in columns[2] and columns[2].startswith("ClassII DNA ") and columns[3].startswith("TIR/EnSpm_CACTA/"):#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\t{columns[3]}', file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassII DNA CACTA") and columns[3] == 'Unknown':#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\tTIR/EnSpm_CACTA/unknown', file=o)
			elif columns[1].startswith("DNA/hAT-") and columns[2].startswith("ClassII DNA hAT ") and (columns[3].startswith("TIR/hAT/") or columns[3] == "Unknown"):
				if columns[3] == "Unknown":
					col3sub = 'unknown'
				else:
					col3sub = columns[3].split("/")[2]
				
				col1sub = columns[1].split("-")[1]
				col2sub= columns[2].split(" ")[3]
				
				if(col1sub == col2sub and col2sub == col3sub):
					newclassif='TIR/hAT/'+col1sub
				elif(col1sub == col2sub and col3sub == 'unknown'):
					newclassif='TIR/hAT/'+col1sub
				elif(col1sub != col3sub and col1sub != col2sub):
					newclassif='TIR/hAT/unknown'
				elif(col1sub != col3sub and col3sub == col2sub):
					newclassif='TIR/hAT/'+col2sub
				elif(col1sub == col2sub and col1sub != col3sub):
					newclassif='TIR/hAT/'+col1sub
				else:
					newclassif='TIR/hAT/unknown'
				print(f'{columns[0]}\t{newclassif}', file=o)
			elif columns[1] == 'Unknown' and columns[3] == 'Unknown' and (columns[2].startswith("ClassII DNA") and columns[2].endswith("MITE")):
				if 'hAT' in columns[2]:
					print(f'{columns[0]}\tTIR/hAT/unknown *MITE', file=o) #TODO revisar se é correcto, MITe vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'TcMar' in columns[2]:
					print(f'{columns[0]}\tTIR/Tc1/Mariner *MITE', file=o) #TODO revisar se é correcto, MITe vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'Mutator' in columns[2]:
					print(f'{columns[0]}\tTIR/MuDR/Mutator *MITE', file=o) #TODO revisar se é correcto, MITe vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'Harbinger' in columns[2]:
					print(f'{columns[0]}\tTIR/PIF/Harbinger *MITE', file=o) #TODO revisar se é correcto, MITe vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183							
			elif columns[1] == 'Unknown' and columns[3] == 'Unknown' and columns[2].startswith("ClassI LTR") and columns[0] not in deeptedomains.keys():
					print(f'{columns[0]}\tUnknown', file=o)
			else:
				print(line)

if __name__ == '__main__':
    main()
