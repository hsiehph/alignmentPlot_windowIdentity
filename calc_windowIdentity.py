import argparse
from Bio import AlignIO, SeqIO

def hamming_distance(seq1, seq2):
	return sum([base1 != base2 for base1, base2 in zip(seq1, seq2) if base1 != "-" and base2 != "-"])

def window(seq, width=500, step=100):
	seqlen = len(seq)
	for i in range(0, seqlen, step):
		if i + width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield (i, j)
		if j == seqlen:
			break

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_alignment")
	parser.add_argument("--window", nargs="?", const=1, default=1000)
	parser.add_argument("--step", nargs="?", const=1, default=100)
	parser.add_argument("--removeGap", nargs="?", const=1, default=True)
	args = parser.parse_args()

	# Load the sequence alignment.
	alignment = AlignIO.read(args.input_alignment, "fasta")

	for win in window(range(alignment.get_alignment_length()), width=int(args.window), step=int(args.step)):
		start, end = win
		if set(alignment[0][start : end]) != {'-'} and set(alignment[1][start : end]) != {'-'}:
			frac_identity = (int(args.window) - hamming_distance(alignment[0][start : end], alignment[1][start : end])) / float(args.window)
			print('\t'.join([str(x) for x in [start, end, frac_identity, "w%s_s%s" % (args.window, args.step)]]))
		
	

