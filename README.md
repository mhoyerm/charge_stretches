Charge stretches

Searches for stretches with certain charges in proteins.

Usage:

python charge_stretches.py -F proteome.fasta -O output.csv -L 30 -C 10 -P 999999 -I 0

-F [required]: Fasta archive containing ORFs list

-O [required]: Output file in csv format

-L [optional]: Lenght of the stretch, default is 30

-C [required]: Insert wanted charge value

-P [optional]: Position of the last amino acid of each protein to be searched for

-I [optional]: Position of the last amino acid of each protein to be searched for
