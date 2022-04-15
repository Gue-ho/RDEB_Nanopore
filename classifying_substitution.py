import pysam, os

fn = '07_707.sorted.bam' #sorted_file
ntp = 2638 #position
nt1 = 'C' #reference_nt
nt2 = 'T' #substitution_nt

sub_fn = fn[:fn.find('.')]

with pysam.AlignmentFile(fn, 'rb') as bamf, pysam.AlignmentFile('{0}_{1}_{2}.bam'.format(sub_fn, ntp, nt1), 'wb', template=bamf) as fw1, pysam.AlignmentFile('{0}_{1}_{2}.bam'.format(sub_fn, ntp, nt2), 'wb', template=bamf) as fw2:

	ref_len = int(bamf.lengths[0])
	
	cnt = {'all': 0, 'nt1': 0, 'nt2': 0, 'other': 0}
	nt_list_frame = ['-' for i in range(ref_len)]

	for read in bamf.fetch():

		cnt['all'] += 1
		
		seq = read.query_sequence
		nt_list = list(nt_list_frame)
		read_p = 0
		ref_p = int(read.pos) 
		s = ''

		for cigar in read.cigartuples:
			if cigar[0] == 0:
				for i in range(cigar[1]):
					nt_list[ref_p + i] = seq[read_p + i]
				ref_p += cigar[1]
				read_p += cigar[1]
			elif cigar[0] == 1:
				read_p += cigar[1]
			elif cigar[0] == 2:
				ref_p += cigar[1]
			elif cigar[0] == 4:
				read_p += cigar[1]

		if nt_list[ntp - 1] == nt1:
			fw1.write(read)
			cnt['nt1'] += 1
		elif nt_list[ntp - 1] == nt2:
			fw2.write(read)
			cnt['nt2'] += 1
		else:
			cnt['other'] += 1

	print(cnt)

for i in [nt1, nt2]:
	os.system('samtools sort {0}_{1}_{2}.bam > {0}_{1}_{2}.sorted.bam'.format(sub_fn, ntp, i))
	os.system('samtools index {0}_{1}_{2}.sorted.bam'.format(sub_fn, ntp, i))
