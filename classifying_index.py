import pysam 

def rc(s): return s.translate(s.maketrans('ATGC', 'TACG'))[::-1]

bar = ['TATAGCCTGATATAGCCT', 'CAGGACGTGACAGGACGT', rc('cgagtaatgacgagtaat'.upper()), rc('agcttcaggaagcttcag'.upper())] #501 507 701 707 index sequence
score_cutoff = 0.5

input_file = 'single.sorted.bam' 
output_file1 = '501_701.bam'
output_file2 = '507_707.bam'

score_dict = {0: [], 1:[], 2:[], 3:[]}

for i in range(len(bar[0]) - 5):
	for n in range(4):
		b = bar[n][i: i + 5]
		if b not in score_dict[n]:
			score_dict[n].append(b)


with pysam.AlignmentFile(input_file, 'rb') as f, pysam.AlignmentFile(output_file, 'wb', template=f) as fw_501, pysam.AlignmentFile(output_file2, 'wb', template=f) as fw_507:
	
	ref_len = int(f.lengths[0])
	cnt = {'all': 0, 'align': 0, 0: 0, 1: 0}

	fw = [fw_501, fw_507]
	#for i in fw:
	#	fw.write(f.header)

	for read in f.fetch():
		cnt['all'] += 1
		if read.is_unmapped == True:
			continue

		#if read.pos > 20 or read.qend < ref_len - 20:
		if read.pos > 20:
			continue
		
		end_p = int(read.pos)
		for cigar in read.cigartuples:
			if cigar[0] == 0:
				end_p += cigar[1]
			if cigar[0] == 2:
				end_p += cigar[1]
		
		if end_p < ref_len - 20:
			continue 

		cnt['align'] += 1
		seq = read.query_sequence
		if read.cigar[0][0] == 4 and read.cigar[0][1] >= 18:
			pos = int(read.cigar[0][1])
			score = {0: 0, 1: 0}
			for i in range(20 - 5):
				for x in range(2):
					if seq[pos - i - 5: pos - i] in score_dict[x]:
						score[x] += 1

			max_val = sorted(score.items(), key = lambda x: x[1], reverse=True)[0]
			if score[0] != score[1] and max_val[1] > len(bar[max_val[0]]) * 0.5:
				fw[max_val[0]].write(read)
				cnt[max_val[0]] += 1
				continue

		if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= 18:
			pos = int(read.cigar[-1][1])
			score = {0: 0, 1: 0}
			for i in range(20 - 5):
				for x in range(2):
					if seq[-pos: -pos + i] in score_dict[x + 2]:
						score[x] += 1

			max_val = sorted(score.items(), key = lambda x: x[1], reverse=True)[0]
			if score[0] != score[1] and max_val[1] > len(bar[max_val[0]]) * 0.5:
				fw[max_val[0]].write(read)
				cnt[max_val[0]] += 1
				continue
	
	print(cnt)

