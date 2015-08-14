#!/usr/bin/env python

import sys
import numpy

#FASTQ quality score sequence. 2nd character ! is score 1.
asciiseq = " !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

#Check if 2 strings match from right to left.
def seq_match(seqa, seqb):
	length = min(len(seqa), len(seqb))
	for i in xrange(-1, -length - 1, -1):
		if seqa[i] != seqb[i]:
			return(False)
	return True

#Function to generate slippage summary for a slippage pattern.
def gen_ptn_summary(ptn, n7_dict):
	#Count of N7s included in the pattern.
	n7_c = 0
	#Count of TSS sites included in the pattern.
	tss_c = 0
	#Matched count.
	m_tnc = 0
	#Slippage count.
	s_tnc = [0] * 12
	#All count.
	all_tnc = 0
	#List slippage rate 1 of all N7s included in the pattern.
	slp1_list = []
	#List slippage rate 2 of all N7s included in the pattern.
	slp2_list = []

	#For each of the N7, check if it is included in the pattern.
	for n7 in n7_dict:
		ptn_match = True
		start_pos = int(n7[-1]) - 4
		for i in xrange(len(ptn)):
			if ptn[i] != 'N':
				if ptn[i] == 'Z':
					if n7[i] == ptn[start_pos]:
						ptn_match = False
						break
				else:
					if n7[i] != ptn[i]:
						ptn_match = False
						break
		#If the N7 is included in the pattern, sum its counts with pattern counts.
		if ptn_match:
			count_list = n7_dict[n7]
			n7_c += 1
			tss_c += count_list[0]
			m_tnc += count_list[1]
			s_tnc[0] += count_list[2][0]
			s_tnc[1] += count_list[2][1]
			s_tnc[2] += count_list[2][2]
			s_tnc[3] += count_list[2][3]
			s_tnc[4] += count_list[2][4]
			s_tnc[5] += count_list[2][5]
			s_tnc[6] += count_list[2][6]
			s_tnc[7] += count_list[2][7]
			s_tnc[8] += count_list[2][8]
			s_tnc[9] += count_list[2][9]
			s_tnc[10] += count_list[2][10]
			s_tnc[11] += count_list[2][11]
			all_tnc += count_list[3]
			slp1_list.append(count_list[4])
			slp2_list.append(count_list[5])
	#Calculate pattern slippage 1, 2, and 3.
	slp1 = numpy.mean(slp1_list)
	slp2 = numpy.mean(slp2_list)
	if sum(s_tnc) + m_tnc > 0:
		slp3 = float(sum(s_tnc)) / (sum(s_tnc) + m_tnc) * 100
	else:
		slp3 = 'NaN'
	#Return summary line for the pattern.
	summary_str = '\t'.join(map(str, [ptn, n7_c, tss_c, all_tnc, m_tnc, '\t'.join(map(str, s_tnc)), slp1, slp2, slp3]))
	return summary_str


def genSlippageStats(tss_fn, ptn_fn, fastq_fn_list):
	#Read pattern list file into a list.
	ptn_list = []
	with open(ptn_fn, 'r') as ptn_file:
		for line in ptn_file:
			ptn = line.strip()
			if len(ptn) != 0:
				ptn_list.append(ptn)

	#Set quality score cutoff.
	qscore_cutoff = asciiseq[25]

	#Read TSS file into a dictionary with schma:
	#{5'-8bp-TSS-29bp-3' : number of 5'-8bp-TSS-29bp-3' in TSS file}
	tss_dict = {}
	tss_file = open(tss_fn, 'r')
	for line in tss_file:
		lsl = line.strip().split()
		tss = lsl[2].upper()
		if tss in tss_dict:
			tss_dict[tss] += 1
		else:
			tss_dict[tss] = 1

	#Create sequence matching dictionary for testing RNA read substring matching to TSS sequence.
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

	stats = '\n'.join(map(lambda p: gen_ptn_summary(p, r_n7_dict), ptn_list)) + '\n'
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