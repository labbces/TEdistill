#!/usr/bin/env python3

import argparse
import sys

def parse_arguments():

	parser = argparse.ArgumentParser(description='Compare different TE classification methods (Earl Grey, DeepTE and TEsorter) and create a consensus on the classification.')


	parser.add_argument('-i', '--input_file', help='Path to input file containing TE classifications table.', required=True)
	parser.add_argument('-d', '--deepte_domain', help='Path to input file of DeepTE domain classifications file.', required=True)
	parser.add_argument('-o', '--output_file', help='Path to output file with the consensus classification.', required=True)
	parser.add_argument('-l', '--delimiter', default=',', help='Delimiter used in the input file. Default is ",".')

	return parser.parse_args()

def main():
	args = parse_arguments()
	deeptedomains = {}

	# Read the DeepTE domain classifications into a dictionary (first column as key(TE_ID), second column as value (Domain)).
	with open(args.deepte_domain, "r") as d:
		for line in d:
			line=line.rstrip()
			columns=line.split()
			id=columns[0].split("#")[0]
			deeptedomains[id] = columns[1]
			continue

	with open(args.input_file, "r") as f , open(args.output_file, "w") as o:
	
	# Loop over the input lines.
		for line in f:
			line=line.rstrip()
			if line.startswith("TE_ID,EarlGrey,DeepTE,TEsorter"):
				continue
			columns=line.split(args.delimiter)
			
			if len(columns) !=4:
				sys.exit()

			#TODO: Add a check to see if the TE_ID is in the DeepTE domain classifications dictionary.
			# Check if the TE_ID is in the DeepTE domain classifications dictionary.
			if columns[0] in deeptedomains:
				domain = deeptedomains[columns[0]]
			else:
				domain = "Unknown"

			if columns[1].lower() == columns[2].lower() and columns[1].lower() == columns[3].lower():
				print(f'{columns[0]}\t{columns[1]}\tAll agree', file=o)
			
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
			elif columns[1] == "Unknown" and columns[2] == "ClassI LTR" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Gypsy" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter', file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown\tEarl Grey and DeepTE', file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassII DNA hAT MITE" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown\tEarl Grey and DeepTE',file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassI" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Copia" and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown\tEarl Grey and DeepTE',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2].startswith("ClassII ") and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LTR/Copia" and columns[2].startswith("ClassII ") and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassI LTR") and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassI LTR") and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tDeepTE and TESorter',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI nLTR LINE L1" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI nLTR LINE L1" and columns[3].startswith("LTR/Gypsy"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI nLTR DIRS" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "unknown" and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassII ") and columns[3].startswith("LTR/Copia"):
				print(f'{columns[0]}\t{columns[3]}\tOnly TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassI LTR") and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassI" and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)
			elif columns[1] == "LTR/Copia" and columns[2] == "ClassI LTR Gypsy" and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\tClassI/LTR/unknown\tEarl Grey and DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tOnly TESorter',file=o)
			elif columns[1] == "LTR/Gypsy" and columns[2] == "ClassI LTR Copia" and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\tClassI/LTR/unknown\tEarl Grey and DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tOnly TESorter',file=o)
			elif columns[1].startswith("LTR") and columns[2].startswith("ClassII DNA") and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\tUnknown\tAll disagree',file=o)
			elif columns[1].startswith("Unknown") and columns[2].startswith("ClassI LTR") and columns[3].startswith("TIR/"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tAll disagree',file=o)
			elif columns[1].startswith("Unknown") and columns[2].startswith("ClassII DNA") and columns[3].startswith("LTR/"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tAll disagree',file=o)

			#PLE
			elif columns[1] == "Unknown" and columns[2] == "ClassI nLTR PLE" and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tDomain found by DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)

			#LINE and SINE
			elif columns[1] == "LINE/L1" and columns[2].startswith("ClassI ") and columns[3].startswith("LINE/unknown"):
				print(f'{columns[0]}\t{columns[1]}/unknown\tEarl Grey and TESorter',file=o)
			elif columns[1] == "LINE/RTE-BovB" and columns[2].startswith("ClassI ") and columns[3].startswith("LINE"):
				print(f'{columns[0]}\t{columns[1]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassI LTR") and columns[3].startswith("LINE/unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tClassI/unknown\tDeepTE and TESorter, but domain not found by DeepTE',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassI nLTR") and columns[3].startswith("Unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)
			elif columns[1].startswith("LINE") and columns[2].startswith("ClassI LTR") and columns[3].startswith("Unknown"): 
				if columns[0] in deeptedomains.keys():
					domains = deeptedomains[columns[0]].split(',')
					unique_domains = set(d.strip() for d in domains if d.strip())
					if len(unique_domains) >= 2:
						print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE', file=o)
					else:
						domain_str = ', '.join(unique_domains) if unique_domains else 'None'
						print(f'{columns[0]}\t{columns[1]}\tOnly Earl Grey', file=o)
				else:
					print(f'{columns[0]}\tClassI/unknown\tEarl Grey and DeepTE, but domain not found by DeepTE', file=o)
			elif columns[1].startswith("LINE") and columns[2].startswith("ClassII DNA") and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\tUnknown\tAll disagree',file=o)
			
			#Helitron
			elif columns[1] == "RC/Helitron" and columns[2].startswith("ClassII DNA") and columns[3].startswith("Helitron"):
				print(f'{columns[0]}\tClassII/{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassIII Helitron") and columns[3].startswith("LTR/Copia"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\t{columns[3]}\tOnly TESorter',file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassIII Helitron" and columns[3].startswith("Unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\tClassII/Helitron\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)
			elif columns[1] == "RC/Helitron" and columns[2].startswith("ClassI LTR") and columns[3].startswith("Unknown"):
					print(f'{columns[0]}\tUnknown\tAll disagree',file=o)

			#TIR
			elif columns[1] == "DNA/PIF-Harbinger" and columns[2] == "ClassII DNA hAT nMITE" and columns[3].startswith("TIR/PIF_Harbinger"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/PIF-Harbinger" and columns[2].startswith("ClassII DNA ") and columns[3].startswith("TIR/PIF_Harbinger"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/MULE-MuDR" and columns[2].startswith("ClassII DNA") and columns[3].startswith("TIR/MuDR_Mutator"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/PIF-Harbinger" and columns[2].startswith("ClassI LTR") and columns[3].startswith("TIR/PIF_Harbinger"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/hAT-Ac" and columns[2] == "ClassII DNA TcMar MITE" and columns[3].startswith("TIR/hAT"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/hAT-Ac" and columns[2].startswith("ClassI LTR") and columns[3].startswith("TIR/hAT"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/hAT-Tip100" and columns[2].startswith("ClassI LTR") and columns[3].startswith("TIR/hAT"):
				print(f'{columns[0]}\t{columns[1]}\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/hAT-Tag1" and columns[2].startswith("ClassII DNA") and columns[3].startswith("TIR/hAT"):
				print(f'{columns[0]}\t{columns[3]}-Tag1\tEarl Grey and TESorter',file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassII DNA") and columns[3].startswith("Unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tEarl Grey and DeepTE',file=o)
				else:
					print(f'{columns[0]}\tClassII/unknown\tEarl Grey and DeepTE, but domain not found by DeepTE',file=o)
			
			#MITEs
			elif columns[1] == 'Unknown' and columns[3] == 'Unknown' and (columns[2].startswith("ClassII DNA") and columns[2].endswith("MITE")):
				if 'hAT' in columns[2]:
					print(f'{columns[0]}\tTIR/hAT/unknown\tOnly DeepTE *MITE', file=o) #TODO revisar se é correcto, MITE vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'TcMar' in columns[2]:
					print(f'{columns[0]}\tTIR/Tc1/Mariner\tOnly DeepTE *MITE', file=o) #TODO revisar se é correcto, MITE vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'Mutator' in columns[2]:
					print(f'{columns[0]}\tTIR/MuDR/Mutator\tOnly DeepTE *MITE', file=o) #TODO revisar se é correcto, MITE vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
				elif 'Harbinger' in columns[2]:
					print(f'{columns[0]}\tTIR/PIF/Harbinger\tOnly DeepTE *MITE', file=o) #TODO revisar se é correcto, MITE vs nMITE https://academic.oup.com/bioinformatics/article/36/15/4269/5838183
			elif columns[1] == "Unknown" and columns[2] == "ClassII MITE" and columns[3] == "Unknown":
				print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
			#"Stowaway MITEs are derived from and mobilised by elements from the mariner superfamily. (https://pmc.ncbi.nlm.nih.gov/articles/PMC6881990/)"
			elif columns[1] == "DNA/TcMar-Stowaway" and columns[2] == "ClassII DNA hAT nMITE" and columns[3].startswith("TIR/Tc1_Mariner"):
				print(f'{columns[0]}\t{columns[1]}\tEarl Grey and TESorter',file=o)
			#"DNA transposons (TIR) and LTR can  be a template for satDNA formation. (https://link.springer.com/article/10.1007/s10577-015-9483-7)"
			elif columns[1] == "Satellite" and 'MITE' in columns[2] and columns[2].startswith("ClassII DNA") and columns[3].startswith("LTR"):
				print(f'{columns[0]}\t{columns[1]}\tOnly Earl Grey', file=o)
			elif columns[1] == "Satellite" and 'MITE' in columns[2] and columns[2].startswith("ClassII DNA") and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}\tOnly Earl Grey', file=o)
			elif columns[1] == "Satellite" and columns[2].startswith("ClassI ") and columns[3].startswith("Unknown"):
				print(f'{columns[0]}\t{columns[1]}\tOnly Earl Grey', file=o)
			elif columns[1] == "Unknown" and columns[2] == "ClassII nMITE" and columns[3].startswith("Unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\tClassII/unknown\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tUnknown\tEarl Grey and TESorter, domain not found by DeepTE',file=o)
			elif columns[1] == "Unknown" and columns[2].startswith("ClassII DNA") and columns[3].startswith("Unknown"):
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE',file=o)
				else:
					print(f'{columns[0]}\tClassII/unknown\tOnly DeepTE, but domain not found by DeepTE',file=o)

			elif columns[1] == "DNA/CMC-EnSpm" and columns[2] == "ClassII DNA hAT MITE" and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE *MITE',file=o)
				else:
					print(f'{columns[0]}\tClassII/unknown\tOnly DeepTE, but domain not found by DeepTE',file=o)
			#"To date, no MITE families have been reported to be related to the Maverick/Polinton, Transib, CACTA, or Merlin superfamilies.
			#However, some MITE families show certain characteristic features that may identify them in these superfamilies."
			#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3290980/
			
			#CACTA
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassII DNA CACTA") and columns[3].startswith("TIR/EnSpm_CACTA/"):#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and 'CACTA' not in columns[2] and columns[2].startswith("ClassII DNA") and columns[3].startswith("TIR/EnSpm_CACTA/"):#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\t{columns[3]}\tAll agree', file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassII DNA CACTA") and columns[3] == 'Unknown':#EnSpm_CACTA is the same as CMC-EnSpm following https://www.jstage.jst.go.jp/article/ggs/94/6/94_18-00024/_html/-char/en
				print(f'{columns[0]}\tTIR/EnSpm_CACTA/unknown\tEarl Grey and DeepTE', file=o)
			elif columns[1] == "DNA/CMC-EnSpm" and columns[2].startswith("ClassI LTR") and columns[3].startswith("TIR/EnSpm_CACTA/"):
				print(f'{columns[0]}\t{columns[3]}\tEarl Grey and TESorter', file=o)
			elif columns[1].startswith("DNA/") and columns[2].startswith("ClassI LTR") and columns[3] == "Unknown":
				if columns[0] in deeptedomains.keys():
					print(f'{columns[0]}\t{columns[2]}\tOnly DeepTE', file=o)
				else:
					print(f'{columns[0]}\t{columns[1]}\tOnly Earl Grey', file=o)
			elif columns[1].startswith("DNA/hAT-") and columns[2].startswith("ClassII DNA hAT") and (columns[3].startswith("TIR/hAT/") or columns[3] == "Unknown"):
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
				print(f'{columns[0]}\t{newclassif}\tAll agree', file=o)
							
			# elif columns[1] == 'Unknown' and columns[3] == 'Unknown' and columns[2].startswith("ClassI LTR") and columns[0] not in deeptedomains.keys():
					# print(f'{columns[0]}\tUnknown\tNone', file=o)
			#elif columns[1] == 'Unknown' and columns[3] == 'Unknown':
			#		print(f'{columns[0]}\t{columns[2]}\tDMRP', file=o)
			else:
				print(line)

if __name__ == '__main__':
	main()
