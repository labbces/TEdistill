import argparse

def parse_arguments():
	parser = argparse.ArgumentParser(description="Padronize files classification nomenclature.")
	parser.add_argument('--Onthology', type=str, required=True, help='Sequence onthology file with all nomenclature from selected tools')
	parser.add_argument('--EarlGrey', type=str, required=True, help='Earl Grey output file path')
	parser.add_argument('--DeepTE', type=str, required=True, help='DeepTE output file path')
	parser.add_argument('--TEsorter', type=str, required=True, help='TEsorter output file path')
	parser.add_argument('--TEtrimmer', type=str, required=True, help='TEtrimmer output file path')
	parser.add_argument('--RepClassifier', type=str, required=True, help='RepeatClassifier output file path')
	parser.add_argument('--MITEtracker', type=str, required=True, help='MITEtracker output file path')
	parser.add_argument('--output', type=str, required=True, help='Output file path')

	args = parser.parse_args()
	return args

def normalize_classification():

	# Define dictionary for nomenclature normalization using the seq_ontology.txt
	classification = {
			ClassI: {eg(), dte(), tes(), rc(), mt()},
			ClassII: {eg(), dte(), tes(), tet(), rc(), mt()}
			}

def main():
	args = parse_arguments()
	normalize_classification()

	# Further processing and file handling would go here
