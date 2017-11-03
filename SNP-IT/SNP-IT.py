#programme to report the lineage/subspecies of a TB sample alligned to NC000962
#useage: python fasta_typer9.py [guuid] [name of outfile]
#or for lots of samples: cat [list of samples] | parallel -j[no. of threads] python fasta_typer9.py {} {}.out
#Output1: is to file containing absolute and % hits for all subspecies
#Output2: only the top call, is to standard out - redirect to file if you want eg 1>calls.log
import sys
from Bio import SeqIO
import gzip, csv, operator, collections
from sys import argv
script, query, out = argv

csv_file=open(out, 'wb')
try:
	query_file = gzip.open(query)
except IOError:
	print query, '\t', '0', 'No_file'
	sys.exit()
#library is list of subspecies to compare, change this address to point to folder where library and library files are contained
library = open('./library')
out_dic={}
for record in library:
	query_dic = {}
	line = record.strip()
	with open(line) as pos_file:
		positions = []		
		reader=csv.reader(pos_file, delimiter= '\t')            
                for x in reader:
                        positions.append(x[0])

		
		query_file = gzip.open(query)
		record=SeqIO.read(query_file, "fasta")
		
		for x in positions:
			pos = int(x)
			python_pos=(pos - 1)
			nuc = record[python_pos]
			query_dic[pos]=nuc
		
	ref_dic={}
	with open(line) as ref_file:
		reader=csv.reader(ref_file, delimiter = '\t')
		for record in reader:
			pos = record[0]
			call = record[1]
			intpos = int(pos)
			ref_dic[intpos] = call
	shared_items = set(query_dic.items()) & set(ref_dic.items())
	csv_writer = csv.writer(csv_file, delimiter='\t')
	
	shared = float(len(shared_items))
	ref = float(len(ref_dic))
	pc = ((shared / ref) * 100) 

			

	out_dic[line]=(pc)
		
	csv_writer.writerow([query, line.rstrip(), shared, pc])

call_dic = collections.defaultdict(list)
made = sorted(out_dic.values(), reverse=True)

for k, v in out_dic.items():
	call_dic[v].append(k)

for key, value in out_dic.items():
	if made[0]==value:
        	key_0 = key
        if made[1]==value:
        	key_1 = key

#if two samples match the 6% threshold, we want to select the sublineage eg. xtype and not the parent lineage eg. lineage 4
if made[0] > 10 and made[1] > 10:
	for key, value in out_dic.items():
        	if made[0]==value:
                	key_0 = key
       		if made[1]==value:
                	key_1 = key

	if 'lineage4' in key_0 and 'lineage4' in key_1:
		print query, '\t', made[0], '\t', key_0
	elif 'lineage4' not in key_0 and 'lineage4' not in key_1:		
		print query, '\t', made[0], '\t', key_0
	elif 'lineage4' not in  key_0:
					
		print query, '\t', made[0], '\t', key_0
	else:
		
		print query, '\t', made[1], '\t', key_1

elif made[0] > 10 and made[1] <10:
	#if only one sublineage meets the threshold then this is our call
	for key, value in out_dic.items():
		if made[0]==value:
			print query, '\t', made[0], '\t', key
else:
	print query, '\t',  made[0],  '\t', 'none'
#if no sublineages meet the threshold, then we have no call


