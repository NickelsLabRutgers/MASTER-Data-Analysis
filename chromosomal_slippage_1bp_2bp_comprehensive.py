#!/usr/bin/env python

import sys
import numpy
asciiseq = " !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
def seq_match(seqa, seqb):
	length = min(len(seqa), len(seqb))
	for i in xrange(-1, -length - 1, -1):
		if seqa[i] != seqb[i]:
			return(False)
	return True


def genSlippageStats(tss_fn, ptn_fn, fastq_fn_list):
	ptn_list = []
	with open(ptn_fn, 'r') as ptn_file:
		for line in ptn_file:
			ptn = line.strip()
			if len(ptn) != 0:
				ptn_list.append(ptn)

	qscore_cutoff = asciiseq[25]
	tss_dict = {}
	tss_file = open(tss_fn, 'r')
	for line in tss_file:
		lsl = line.strip().split()
		tss = lsl[2].upper()
		if tss in tss_dict:
			tss_dict[tss] += 1
		else:
			tss_dict[tss] = 1

	seq_mat_dict = {}
	for tss in tss_dict:
		for i in xrange(12):
			tss_sl = tss[i:]
			if tss_sl not in seq_mat_dict:
				seq_mat_dict[tss_sl] = [tss]
			else:
				if tss not in seq_mat_dict[tss_sl]:
					if tss[5:12] != seq_mat_dict[tss_sl][0][5:12]:
						del seq_mat_dict[tss_sl]

	n7_dict = {}
	for tss in tss_dict:
		n7 = tss[5:12]
		for i in xrange(2,6):
			if n7[i] == n7[i+1]:
				ln7 = n7 + '-H' + str(i+4)
				if ln7 not in n7_dict:
					n7_dict[ln7] = [0,0,[0]*12, 0]
					#number TSS, match, slp, all
				n7_dict[ln7][0] += tss_dict[tss]

			if n7[i] != n7[i+1]:
				ln7 = n7 + '-T' + str(i+4)
				if ln7 not in n7_dict:
					n7_dict[ln7] = [0,0,[0]*12, 0]
				n7_dict[ln7][0] += tss_dict[tss]

			if n7[i-1] != n7[i]:
				ln7 = n7 + '-I' + str(i+4)
				if ln7 not in n7_dict:
					n7_dict[ln7] = [0,0,[0]*12, 0]
				n7_dict[ln7][0] += tss_dict[tss]
	
	for ifn in fastq_fn_list:
		ifile = open(ifn)
		line1 = ifile.readline()
		seq = ifile.readline()
		line3 = ifile.readline()
		qscore = ifile.readline()
		
		while line1 != '' and seq != '' and line3 != '' and qscore != '':
			for q in qscore[:38]:
				if q < qscore_cutoff:
					line1 = ifile.readline()
					seq = ifile.readline()
					line3 = ifile.readline()
					qscore = ifile.readline()
					break
				
			seq_tss_list = None
			mat_found = False
			for i in xrange(12):
				for j in xrange(i, 12):
					if seq[i:38-j+i] in seq_mat_dict:
						seq_tss_list = (seq_mat_dict[seq[i:38-j+i]][0], (i,j))
						mat_found = True
						break
				if mat_found:
					break

			if seq_tss_list:
				tss = seq_tss_list[0]
				n7 = tss[5:12]
				seq_start = seq_tss_list[1][0]
				tss_start = seq_tss_list[1][1]
				rna_seq = seq[:12 - tss_start + seq_start]
				if len(rna_seq) > 12:
					print "RNA seq len error %s" % rna_seq
					sys.exit(2)

				match_boolean = seq_match(rna_seq, tss[:12])

				for i in xrange(2,6):
					if n7[i] == n7[i+1]:
						ln7 = n7 + '-H' + str(i+4)
							#n7_dict[ln7] = [0,0,[0]*12, 0]
							#number TSS, match, slp, all
						n7_dict[ln7][3] += 1
						if match_boolean and len(rna_seq) == (7-i):
							n7_dict[ln7][1] += 1
						else:
							if (len(rna_seq) > (7-i)) and (rna_seq[:-(7-i)] == n7[i] * (len(rna_seq) - (7-i))) and seq_match(rna_seq[-(7-i):], n7):
								n7_dict[ln7][2][len(rna_seq) - 1] += 1

						if n7[i-1] != n7[i]:
							ln7 = n7 + '-I' + str(i+4)
							n7_dict[ln7][3] += 1
							if match_boolean and len(rna_seq) == (8-i):
								n7_dict[ln7][1] += 1
							else:
								if (len(rna_seq) > (8-i)) and (rna_seq[:-(7-i)] == n7[i-1] + (n7[i] * (len(rna_seq) - (7-i) - 1))) and seq_match(rna_seq[-(7-i):], n7):
									n7_dict[ln7][2][len(rna_seq) - 1] += 1
				
					if n7[i] != n7[i+1]:
						ln7 = n7 + '-T' + str(i+4)
						n7_dict[ln7][3] += 1
						if match_boolean and len(rna_seq) == (7-i):
							n7_dict[ln7][1] += 1
						elif len(rna_seq) > (7-i) and (len(rna_seq) - (7-i)) % 2 == 0:
							if (rna_seq[:-(7-i)] == n7[i:i+2] * ((len(rna_seq) - (7-i)) / 2)) and seq_match(rna_seq[-(7-i):], n7):
								n7_dict[ln7][2][len(rna_seq) - 1] += 1

				

			line1 = ifile.readline()
			seq = ifile.readline()
			line3 = ifile.readline()
			qscore = ifile.readline()
		ifile.close()
	r_n7_dict = {}
	for ln7 in n7_dict:
		count_list = n7_dict[ln7]
		if count_list[1] + sum(count_list[2]) > 0:
			r_n7_dict[ln7] = [0, 0, [0]*12, 0, 0, 0]
			m_tnc = count_list[1]
			s_tnc = sum(count_list[2])
			all_tnc = count_list[3]
			r_n7_dict[ln7][0] = count_list[0]
			r_n7_dict[ln7][1] = m_tnc
			r_n7_dict[ln7][2] = count_list[2]
			r_n7_dict[ln7][3] = all_tnc
			r_n7_dict[ln7][4] = 100 * (float(s_tnc) / (s_tnc + m_tnc))
			r_n7_dict[ln7][5] = 100 * (float(s_tnc) / all_tnc)

	out_n7_dict = {}
	for ln7 in r_n7_dict:
		for ptn in ptn_list:
			start_pos = int(ptn[-1]) - 4
			ptn_match = True
			for i in xrange(len(ptn)):
				if ptn[i] != 'N':
					if ptn[i] == 'Z':
						if ln7[i] == ptn[start_pos]:
							ptn_match = False
							break
					else:
						if ln7[i] != ptn[i]:
							ptn_match = False
							break
			if ptn_match:
				count_list = r_n7_dict[ln7]
				if ln7[:7] not in out_n7_dict:
					out_n7_dict[ln7[:7]] = [count_list[0], count_list[3], 0, [0] * 12, [], [], set()]
					#num TSS, all, mat, slp list, sr1, sr2, slp types
					#0        1    2    3         4    5    6
				elif ln7 in out_n7_dict[ln7[:7]][6]:
					continue

				if ln7[-2:] != 'I8':
					out_n7_dict[ln7[:7]][2] += count_list[1]
					
				out_n7_dict[ln7[:7]][3][0] += count_list[2][0]
				out_n7_dict[ln7[:7]][3][1] += count_list[2][1]
				out_n7_dict[ln7[:7]][3][2] += count_list[2][2]
				out_n7_dict[ln7[:7]][3][3] += count_list[2][3]
				out_n7_dict[ln7[:7]][3][4] += count_list[2][4]
				out_n7_dict[ln7[:7]][3][5] += count_list[2][5]
				out_n7_dict[ln7[:7]][3][6] += count_list[2][6]
				out_n7_dict[ln7[:7]][3][7] += count_list[2][7]
				out_n7_dict[ln7[:7]][3][8] += count_list[2][8]
				out_n7_dict[ln7[:7]][3][9] += count_list[2][9]
				out_n7_dict[ln7[:7]][3][10] += count_list[2][10]
				out_n7_dict[ln7[:7]][3][11] += count_list[2][11]
				out_n7_dict[ln7[:7]][4].append(count_list[4])
				out_n7_dict[ln7[:7]][5].append(count_list[5])
				out_n7_dict[ln7[:7]][6].add(ln7)

#	for n7tuple in itertools.product('ACGT', repeat = 7):
#		in7 = ''.join(n7tuple)
#		if in7 not in out_n7_dict:
#			out_n7_dict[in7] = []

	summary_list = []
	for on7 in out_n7_dict:
		o_count_list = out_n7_dict[on7]
		if o_count_list == []:
			summary_list.append(on7 + '\t' + '\t'.join(['0'] * 14 + ['NaN'] * 4) + '\t')
			continue

		a_tnc = o_count_list[1]
		m_tnc = o_count_list[2]
		s_tnc = o_count_list[3]

		if sum(s_tnc) + m_tnc > 0:
			slp3 = float(sum(s_tnc)) / (sum(s_tnc) + m_tnc) * 100
		else:
			slp3 = 'NaN'

		if a_tnc > 0:
			slp4 = float(sum(s_tnc)) / a_tnc * 100
		else:
			slp4 = 'NaN'
		summary_list.append('\t'.join(map(str, [on7, o_count_list[0], a_tnc, m_tnc, 
												'\t'.join(map(str, s_tnc)), 
												numpy.mean(o_count_list[4]), 
												numpy.mean(o_count_list[5]), 
												slp3, slp4, ','.join(o_count_list[6])])))

	stats = '\n'.join(sorted(summary_list, key = lambda i: i[:7]))
	return stats

if __name__ == "__main__":
	argv = sys.argv
	usage = '''Usage:
%s <input TSS file b8a29> <pattern list file> <input FASTQ files>

TSS file is tab delimited file listing transcription start position in reference genome, strand (+ or -), and 5'-8bp-TSS-29bp-3' sequence.

Pattern list file lists slippage patterns, 1 pattern per line. Pattern example, NNZAAZN-H7:
N denotes A/C/G/T
H7 denotes slippage type homopolymeric tract start at position 7.
Z denotes nucleotides different from the base at slippage sequence start, C/G/T in this case.

Output to stdout.\
''' % argv[0]
	if len(argv) < 4:
		print usage
		sys.exit(2)
	tss_fn = argv[1]
	ptn_fn = argv[2]
	fastq_fn_list = argv[3:]
	stats = genSlippageStats(tss_fn, ptn_fn, fastq_fn_list)
	print stats