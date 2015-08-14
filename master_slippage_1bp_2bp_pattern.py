#!/usr/bin/env python

import sys, numpy

BEG_SEQ = 'GCTTCGGCTCGTATAAT'
def seq_match(seqa, seqb):
	length = min(len(seqa), len(seqb))
	for i in xrange(-1, -length - 1, -1):
		if seqa[i] != seqb[i]:
			return(False)
	return True

def gen_ptn_summary(ptn, n7_dict):
	n7_c = 0
	m_tnc = 0
	s_tnc = [0] * 12
	all_tnc = 0
	slp1_list = []
	slp2_list = []
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
		if ptn_match:
			count_list = n7_dict[n7]
			n7_c += 1
			m_tnc += count_list[0]
			s_tnc[0] += count_list[1][0]
			s_tnc[1] += count_list[1][1]
			s_tnc[2] += count_list[1][2]
			s_tnc[3] += count_list[1][3]
			s_tnc[4] += count_list[1][4]
			s_tnc[5] += count_list[1][5]
			s_tnc[6] += count_list[1][6]
			s_tnc[7] += count_list[1][7]
			s_tnc[8] += count_list[1][8]
			s_tnc[9] += count_list[1][9]
			s_tnc[10] += count_list[1][10]
			s_tnc[11] += count_list[1][11]
			all_tnc += count_list[2]
			slp1_list.append(count_list[3])
			slp2_list.append(count_list[4])
	slp1 = numpy.mean(slp1_list)
	slp2 = numpy.mean(slp2_list)
	if sum(s_tnc) + m_tnc > 0:
		slp3 = float(sum(s_tnc)) / (sum(s_tnc) + m_tnc) * 100
	else:
		slp3 = 'NaN'
	summary_str = '\t'.join(map(str, [ptn, n7_c, all_tnc, m_tnc, '\t'.join(map(str, s_tnc)), slp1, slp2, slp3]))
	return summary_str

def gen_MASTER_slippage_stats(ptn_fn, ifn_list):
	ptn_list = []
	with open(ptn_fn, 'r') as ptn_file:
		for line in ptn_file:
			ptn = line.strip()
			if len(ptn) != 0:
				ptn_list.append(ptn)

	n7_dict = {}
	for ifn in ifn_list:
		ifile = open(ifn)
		for line in ifile:
			lsl = line.strip().split()
			n7 = lsl[0]
			rna_seq = lsl[1]
			tnc = int(lsl[5])
			match_boolean = seq_match(rna_seq, BEG_SEQ + n7)

			for i in xrange(2,6):
				if n7[i] == n7[i+1]:
					ln7 = n7 + '-H' + str(i+4)
					if ln7 not in n7_dict:
						n7_dict[ln7] = [0,[0]*12, 0]
					n7_dict[ln7][2] += tnc
					if match_boolean and len(rna_seq) == (7-i):
						n7_dict[ln7][0] += tnc
					else:
						if (len(rna_seq) > (7-i)) and (rna_seq[:-(7-i)] == n7[i] * (len(rna_seq) - (7-i))) and seq_match(rna_seq[-(7-i):], n7):
							n7_dict[ln7][1][len(rna_seq) - 1] += tnc
					if n7[i-1] != n7[i]:
						ln7 = n7 + '-I' + str(i+4)
						if ln7 not in n7_dict:
							n7_dict[ln7] = [0,[0]*12, 0]
						n7_dict[ln7][2] += tnc
						if match_boolean and len(rna_seq) == (8-i):
							n7_dict[ln7][0] += tnc
						else:
							if (len(rna_seq) > (8-i)) and (rna_seq[:-(7-i)] == n7[i-1] + (n7[i] * (len(rna_seq) - (7-i) - 1))) and seq_match(rna_seq[-(7-i):], n7):
								n7_dict[ln7][1][len(rna_seq) - 1] += tnc

				if n7[i] != n7[i+1]:
					ln7 = n7 + '-T' + str(i+4)
					if ln7 not in n7_dict:
						n7_dict[ln7] = [0,[0]*12, 0]
					n7_dict[ln7][2] += tnc
					if match_boolean and len(rna_seq) == (7-i):
						n7_dict[ln7][0] += tnc
					elif len(rna_seq) > (7-i) and (len(rna_seq) - (7-i)) % 2 == 0:
						if (rna_seq[:-(7-i)] == n7[i:i+2] * ((len(rna_seq) - (7-i)) / 2)) and seq_match(rna_seq[-(7-i):], n7):
							n7_dict[ln7][1][len(rna_seq) - 1] += tnc

		ifile.close()

	r_n7_dict = {}
	for ln7 in n7_dict:
		count_list = n7_dict[ln7]
		if count_list[0] + sum(count_list[1]) > 0:
			r_n7_dict[ln7] = [0, 0, 0, 0, 0]
			m_tnc = count_list[0]
			s_tnc = sum(count_list[1])
			all_tnc = count_list[2]
			r_n7_dict[ln7][0] = m_tnc
			r_n7_dict[ln7][1] = count_list[1]
			r_n7_dict[ln7][2] = all_tnc
			r_n7_dict[ln7][3] = 100 * (float(s_tnc) / (s_tnc + m_tnc))
			r_n7_dict[ln7][4] = 100 * (float(s_tnc) / all_tnc)

	stats = '\n'.join(map(lambda p: gen_ptn_summary(p, r_n7_dict), ptn_list)) + '\n'
	return stats




if __name__ == "__main__":
	argv = sys.argv
	if len(argv) < 3:
		print "Usage:\n%s <Pattern table> <input all sequences PhoneBook tables>\nOutput to stdout." %argv[0]
		sys.exit(2)

	ptn_fn = argv[1]
	ifn_list = argv[2:]

	stats = gen_MASTER_slippage_stats(ptn_fn, ifn_list)
	print stats