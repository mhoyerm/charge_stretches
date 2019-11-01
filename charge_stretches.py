import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Search for stretches with certain charges in proteins', usage='[options]')
parser.add_argument('-F', '--fasta', type=str, required=True, help='fasta archive containing ORFs list')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (csv format)')
parser.add_argument('-L', '--length', type=int, default=30, help='length of the stretch')
parser.add_argument('-C', '--charge', type=int, required=True, help='Insert wanted charge value')
parser.add_argument('-P', '--finalpos', type=int, default=999999, help='Position of the last aminoacid of each protein to be searched for')
parser.add_argument('-I', '--initialpos', type=int, default=0, help='Position of the last aminoacid of each protein to be searched for')
args = parser.parse_args()

def percent(index, genome_length): # Iterface updating the percentage of the ORFs read
	percentage = 100*index/genome_length
	interface = int(percentage)
	print str(interface) + '% complete\r',
import sys


def get_charge(nt1, nt2, nt3):
	if nt1 == 'G' and nt2 == 'A':
		return -1

	if nt1 == 'C' and nt2 == 'G':
		return 1

	if nt1 == 'A':
		if nt2 == 'A' or nt2 == 'G':
			if nt3 == 'A' or nt3 == 'G':
				return 1

	return 0

################################################################################################

def evaluate(prote, fout, prot_name):
	netCharge = ''
	last = -2
	#for each character (does not include the last field size characteres)

	end_pos = len(prote)/3 - args.length + 1
	if end_pos > args.finalpos:
		end_pos = args.finalpos

	if args.initialpos + args.length > end_pos:
		return

	for i in range(args.initialpos, end_pos):
		charge = 0
		#read the next field size caracteres
		for j in range(args.length): #compute the netCharge
			nt1 = prote[(i+j)*3]
			nt2 = prote[(i+j)*3 + 1]
			nt3 = prote[(i+j)*3 + 2]
			charge += get_charge(nt1, nt2, nt3)

		relative = float(i)/float(len(prote))

		if charge >= args.charge and args.charge > 0:
			if i > last + 1:
				fout.write("\n" + str(prot_name) + ";" + str(i*3+1) + ";" + str(relative) + ';' + str(charge))
			last = i
		if charge <= args.charge and args.charge < 0:
			if i > last + 1:
				fout.write("\n" + str(prot_name) + ";" + str(i*3+1) + ";" + str(relative) + ';' + str(charge))
			last = i
	return



###################################################################################
def main():
	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()

	out_file = open(args.output, 'w') # Output file
	
	actual_path = os.getcwd() # Get file path
	
	header_basic = 'Search for stretches with certain charges in proteins\n'
	header_fasta = 'fasta file:;' + actual_path + '/' + args.fasta + '\n' # Prepare fasta path
	header_length = 'length of window:;' + str(args.length) + '\n'
	header_difference = 'charge wanted:;' + str(args.charge) + '\n'
	header_gene_stretch = 'position of aminoacids in the protein:;' + str(args.initialpos) + ' - ' + str(args.finalpos) + '\n\n\n\n\n'
	header_data = 'name;position (in nucleotides);relative position;charge'
	
	file_header = header_basic + header_fasta + header_length + header_difference + header_gene_stretch + header_data # Build header without difference value
	out_file.write(file_header) # Write header

	proteome_list = in_file.replace('>UniRef100_','>')
	proteome_list = str.split(in_file, '>') # Split genes by '>'

	for i in range(len(proteome_list)): # Repeat for each gene
		percent(i, (len(proteome_list)))

		protein_code = ''

		proteome_list[i] = proteome_list[i].replace('\r', '')
		protein_by_line = str.split(proteome_list[i], '\n') # Split each line
		header = protein_by_line[0] # The first line is the header
		name = str.split(header,' ')[0]
		for j in range(1, len(protein_by_line)):
			protein_code += protein_by_line[j] # Every line (except the first) is part of the code

		protein_aminoacids = list(protein_code) # List of all aminoacids in a gene
		evaluate(protein_aminoacids, out_file, name)
	print '\r100% complete!'

main()