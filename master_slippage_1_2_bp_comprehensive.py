#!/usr/bin/env python

import sys, numpy, itertools

BEG_SEQ = 'GCTTCGGCTCGTATAAT'
def seq_match(seqa, seqb):
	length = min(len(seqa), len(seqb))
	for i in xrange(-1, -length - 1, -1):
		if seqa[i] != seqb[i]:
			return(False)
	return True

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
			r_n7_dict[ln7] = [0, [0]*12, 0, 0, 0]
			m_tnc = count_list[0]
			s_tnc = sum(count_list[1])
			all_tnc = count_list[2]
			r_n7_dict[ln7][0] = m_tnc
			r_n7_dict[ln7][1] = count_list[1]
			r_n7_dict[ln7][2] = all_tnc
			r_n7_dict[ln7][3] = 100 * (float(s_tnc) / (s_tnc + m_tnc))
			r_n7_dict[ln7][4] = 100 * (float(s_tnc) / all_tnc)


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
					out_n7_dict[ln7[:7]] = [count_list[2], 0, [0] * 12, [], [], set()]
				elif ln7 in out_n7_dict[ln7[:7]][5]:
					continue
				#all, mat, slp, slpR1, slpR2, slpR3, ptn
				if ln7[-2:] != 'I8':
					out_n7_dict[ln7[:7]][1] += count_list[0]
				out_n7_dict[ln7[:7]][2][0] += count_list[1][0]
				out_n7_dict[ln7[:7]][2][1] += count_list[1][1]
				out_n7_dict[ln7[:7]][2][2] += count_list[1][2]
				out_n7_dict[ln7[:7]][2][3] += count_list[1][3]
				out_n7_dict[ln7[:7]][2][4] += count_list[1][4]
				out_n7_dict[ln7[:7]][2][5] += count_list[1][5]
				out_n7_dict[ln7[:7]][2][6] += count_list[1][6]
				out_n7_dict[ln7[:7]][2][7] += count_list[1][7]
				out_n7_dict[ln7[:7]][2][8] += count_list[1][8]
				out_n7_dict[ln7[:7]][2][9] += count_list[1][9]
				out_n7_dict[ln7[:7]][2][10] += count_list[1][10]
				out_n7_dict[ln7[:7]][2][11] += count_list[1][11]
				out_n7_dict[ln7[:7]][3].append(count_list[3])
				out_n7_dict[ln7[:7]][4].append(count_list[4])
				out_n7_dict[ln7[:7]][5].add(ln7)
	
	for n7tuple in itertools.product('ACGT', repeat = 7):
		in7 = ''.join(n7tuple)
		if in7 not in out_n7_dict:
			out_n7_dict[in7] = []

	summary_list = []
	for on7 in out_n7_dict:
		o_count_list = out_n7_dict[on7]
		if o_count_list == []:
			#summary_list.append(on7 + '\t' + '\t'.join(['0'] * 14 + ['NaN'] * 4) + '\t')
			continue

		m_tnc = o_count_list[1]
		s_tnc = o_count_list[2]
		a_tnc = o_count_list[0]

		if sum(s_tnc) + m_tnc > 0:
			slp3 = float(sum(s_tnc)) / (sum(s_tnc) + m_tnc) * 100
		else:
			slp3 = 'NaN'

		if a_tnc > 0:
			slp4 = float(sum(s_tnc)) / a_tnc * 100
		else:
			slp4 = 'NaN'

		summary_list.append('\t'.join(map(str, [on7, o_count_list[0], m_tnc, '\t'.join(map(str, s_tnc)), numpy.mean(o_count_list[3]), numpy.mean(o_count_list[4]), slp3, slp4, ','.join(o_count_list[5])])))

	stats = '\n'.join(sorted(summary_list, key = lambda i: i[:7]))
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