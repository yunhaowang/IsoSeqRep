#!/usr/bin/env python
import sys,time,re,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_repeat_bed = parse_repeat_bed(args.input_bed)
	output_gpd = args.output
	p = Pool(processes=args.cpu)
	csize = 100
	results = p.imap(func=stat_repeat,iterable=generate_tx(args.input_gpd,args.min_len,args.min_rep_pct,args.min_iso_pct),chunksize=csize)
	for res in results:
		if not res: continue
		output_gpd.write(res+"\n")
	output_gpd.close()
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_repeat_bed(input_bed):
	global dic_repeat_bed
	dic_repeat_bed = {}
	for line in input_bed: # sort -k1,1 -k6,6 -k2,2n -k3,3n
		chr,start,end,name,family,strand = [str(i) for i in line.strip().split("\t")[:6]]
		chr_strand = chr+"&"+strand
		if chr_strand not in dic_repeat_bed.keys():
			dic_repeat_bed[chr_strand] = []
			dic_repeat_bed[chr_strand].append("\t".join([start,end,name,family]))
		else:
			dic_repeat_bed[chr_strand].append("\t".join([start,end,name,family]))
	input_bed.close()
	return dic_repeat_bed

def calculate_iso_len(iso_exon_number,iso_exon_start,iso_exon_end):
	iso_len = 0
	exon_start_list = [int(i) for i in iso_exon_start.split(",")[:-1]]
	exon_end_list = [int(i) for i in iso_exon_end.split(",")[:-1]]
	for i in range(0,int(iso_exon_number)):
		iso_len += (exon_end_list[i]-exon_start_list[i])
	return iso_len

def calculate_overlap(rep_start,rep_end,iso_exon_number,iso_exon_start,iso_exon_end):
	overlap_len = 0
	exon_start_list = [int(i) for i in iso_exon_start.split(",")[:-1]]
	exon_end_list = [int(i) for i in iso_exon_end.split(",")[:-1]]
	for i in range(0,int(iso_exon_number)):
		if exon_end_list[i] <= int(rep_start):
			break
		elif exon_start_list[i] >= int(rep_end):
			continue
		else:
			overlap_len += (min([int(rep_end),exon_end_list[i]])-max([int(rep_start),exon_start_list[i]]))
	return overlap_len

def generate_tx(inf,min_overlap_len,min_rep_pct,min_iso_pct):
	z = 0
	for line in inf:
		z += 1
		yield (line,z,min_overlap_len,min_rep_pct,min_iso_pct)

def stat_repeat(inputs):
	(line,z,min_overlap_len,min_rep_pct,min_iso_pct) = inputs
	gene_id,isoform_id,chr,strand,tss,tts,full_length_LR_count,LR_count,exon_number,exon_start,exon_end = line.strip().split("\t")[:11]
	chr_strand = chr+"&"+strand
	if chr_strand not in dic_repeat_bed.keys():
		output_line = "\t".join([line.strip(),"NA","NA"])
	else:
		rep_set_list = []
		iso_len = calculate_iso_len(exon_number,exon_start,exon_end)
		for rep_info in dic_repeat_bed[chr_strand]:
			rep_start,rep_end,name,family = rep_info.split("\t")
			if int(tss) >= int(rep_end):
				continue
			elif int(tts) <= int(rep_start):
				break
			else:
				rep_len = int(rep_end) - int(rep_start)
				overlap_len = calculate_overlap(rep_start,rep_end,exon_number,exon_start,exon_end)
				if (overlap_len >= min_overlap_len):
					if (float(overlap_len)/rep_len) >= min_rep_pct and (float(overlap_len)/iso_len) >= min_iso_pct:
						rep_set_list.append("|".join([name,family,str(overlap_len),str(round(float(overlap_len)/rep_len,2)),str(round(float(overlap_len)/iso_len,2))]))
		
		if rep_set_list == []:
			output_line = "\t".join([line.strip(),"NA","NA"])
		else:
			output_line = "\t".join([line.strip(),",".join(rep_set_list),str(len(rep_set_list))])
	return output_line

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. For repeat element analysis, repeat set split by ','. For each set split by '|': Repeat name; Repeat family; Overlap length between repeat and isoform sequence; Overlap percentage over repeat element; Overlap percentage over isoform sequence
21. For repeat element analysis, number of repeat elements overlapping with this isoform'''

	parser = argparse.ArgumentParser(description="Function: generate final output file (gpd format)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input_gpd',type=argparse.FileType('r'),required=True,help="Input: constructed isoform gpd file generated by 'py_isoseqcon_generate_output.py'")
	parser.add_argument('-b','--input_bed',type=argparse.FileType('r'),required=True,help="Input: repeat element file (BED file, 'chr,start,end,repeat_name,repeat_family,strand'). Must be 'sort -k1,1 -k6,6 -k2,2n -k3,3n'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform with Repeat analysis (gpd file)")
	parser.add_argument('--min_len',type=int,default=100,help="Minimal overlap length between repeat element and isoform sequence (bp)")
	parser.add_argument('--min_rep_pct',type=float,default=0.0,help="Minimal overlap percentage over the repeat element")
	parser.add_argument('--min_iso_pct',type=float,default=0.0,help="Minimal overlap percentage over the isoform sequence")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of threads")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
