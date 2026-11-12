import sys
import argparse
from collections import defaultdict

def reverse_complement(dna: str) -> str:
	"""
	Compute reverse complement of DNA string.
	"""
	complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
	return dna.translate(complement)[::-1]


class HaplotypeIndex:
	"""
	Store syncmer sequence for fast access.
	"""
	def __init__(self, syncmer_file, haplotype_len=-1, reverse=False):
		self._syncmers = []
		self._syncmer_to_pos = defaultdict(lambda: [None, None, None])
		self._reversed = reverse

		syncmer_idx = 0
		with open(syncmer_file, 'r') as infile:
			for line in infile:
				if line.startswith("syncmer_ID"):
					continue
				fields = line.strip().split()

				if reverse:
					# reverse orientation of the syncmer_id
					syncmer_id = fields[0][1:] if fields[0].startswith("-") else "-" + fields[0]
					self._syncmers.append(syncmer_id)
					assert not syncmer_id in self._syncmer_to_pos
					# a syncmer spanning [p, p+k) on the forward strand spans [L-p-k, L-p)
					# on the reverse-complement strand, so its start there is L-p-k, not L-p
					kmer_size = len(fields[3])
					if int(fields[1]) + kmer_size > haplotype_len:
						raise RuntimeError(
							f"HaplotypeIndex: syncmer {fields[0]} at position {fields[1]} with "
							f"length {kmer_size} extends beyond the given haplotype_len {haplotype_len}."
						)
					pos = haplotype_len - int(fields[1]) - kmer_size
					# reverse complement of syncmer sequence
					sequence = reverse_complement(fields[3])
					self._syncmer_to_pos[syncmer_id] = [syncmer_idx, pos, sequence]
				else:
					self._syncmers.append(fields[0])
					assert not fields[0] in self._syncmer_to_pos
					self._syncmer_to_pos[fields[0]] = [syncmer_idx, int(fields[1]), fields[3]]
				syncmer_idx += 1

		# if reverse is true, need to adjust syncmer_idx
		if reverse:
			self._syncmers = self._syncmers[::-1]
			for syncmer in self._syncmer_to_pos:
				self._syncmer_to_pos[syncmer][0] = self.get_length() - self._syncmer_to_pos[syncmer][0] - 1

	def is_reversed(self):
		"""
		Check if this haplotype has been reversed.
		"""
		return self._reversed


	def get_pos(self, syncmer_id):
		"""
		Get position of syncmer on haplotype.
		"""
		if not syncmer_id in self._syncmer_to_pos:
			return None
		else:
			return self._syncmer_to_pos[syncmer_id][1]


	def get_idx(self, syncmer_id):
		"""
		Get index of the syncmer in path.
		"""
		if not syncmer_id in self._syncmer_to_pos:
			return None
		else:
			return self._syncmer_to_pos[syncmer_id][0]


	def get_seq(self, syncmer_id):
		"""
		Get nucleotide sequence of syncmer.
		"""
		if not syncmer_id in self._syncmer_to_pos:
			return None
		else:
			return self._syncmer_to_pos[syncmer_id][2]


	def is_next_syncmer(self, syncmer_id):
		"""
		Check if there is a next syncmer.
		"""
		if not syncmer_id in self._syncmer_to_pos:
			return False
		else:
			syncmer_idx = self.get_idx(syncmer_id)
			largest_idx = len(self._syncmers) - 1
			if syncmer_idx < largest_idx:
				return True
			else:
				return False


	def get_syncmer(self, idx):
		"""
		Get the syncmer at index idx.
		"""
		if idx < len(self._syncmers):
			return self._syncmers[idx]
		else:
			raise RuntimeError("get_syncmer: idx is too large.")


	def get_next_syncmer(self, syncmer_id):
		"""
		Get the next syncmer on the path.
		"""
		if not syncmer_id in self._syncmer_to_pos:
			return None
		else:
			syncmer_idx = self.get_idx(syncmer_id)
			largest_idx = len(self._syncmers) - 1
			if syncmer_idx < largest_idx:
				return self._syncmers[syncmer_idx + 1]
			else:
				return None


	def get_length(self):
		"""
		Get length of syncmer path
		"""
		return len(self._syncmers)


	def contains_syncmer(self, syncmer_id):
		"""
		Check if the haplotype contains a given syncmer.
		"""
		return syncmer_id in self._syncmer_to_pos


def find_mem(haplotype_index1, haplotype_index2, syncmer_id):
	"""
	Find a MEM starting at syncmer_id
	"""
	mem = []
	mem_sequence = ""

	if not haplotype_index1.contains_syncmer(syncmer_id) or not haplotype_index2.contains_syncmer(syncmer_id):
		return mem, mem_sequence, -1
	else:
		mem.append(syncmer_id)
		assert haplotype_index1.get_seq(syncmer_id) == haplotype_index2.get_seq(syncmer_id)
		mem_sequence = haplotype_index1.get_seq(syncmer_id)

	# determine kmer size
	seq_len1 = len(haplotype_index1.get_seq(syncmer_id))
	seq_len2 = len(haplotype_index2.get_seq(syncmer_id))
	if seq_len1 != seq_len2:
		raise RuntimeError("find_mem: kmer sizes between both haplotypes do not match.")
		sys.exit(1)
	
	kmer_size = seq_len1

	cur_id = syncmer_id
	while haplotype_index1.is_next_syncmer(cur_id) and haplotype_index2.is_next_syncmer(cur_id):
		sync_h1 = haplotype_index1.get_next_syncmer(cur_id)
		sync_h2 = haplotype_index2.get_next_syncmer(cur_id)
		if (sync_h1 == sync_h2):

			# get syncmer sequences. Should be identical
			# in both haplotypes.
			seq1 = haplotype_index1.get_seq(sync_h1)
			seq2 = haplotype_index2.get_seq(sync_h1)
			assert seq1 == seq2
			
			# determine sequence overlap between preceding syncmers. Make
			# sure overlap sizes are consistent between both haplotypes.
			overlap1 = (haplotype_index1.get_pos(mem[-1]) + kmer_size) - haplotype_index1.get_pos(sync_h1)
			overlap2 = (haplotype_index2.get_pos(mem[-1]) + kmer_size) - haplotype_index2.get_pos(sync_h1)

			if (overlap1 != overlap2):
				# same syncmer IDs but different overlap size and thus different
				# sequences (can happen in repetitive regions), end of MEM
				break
			if overlap1 <= 0:
				# the syncmers are further than the kmer_size apart, end of MEM
				break

			# combine the sequences
			mem_sequence += seq1[overlap1:]

			# keep mem syncmer id
			mem.append(sync_h1)
			cur_id = sync_h1
		else:
			break
	return mem, mem_sequence, kmer_size


def generate_paf(syncmers_h1, syncmers_h2, outfile_name, chromosome1, chromosome2, lenH1, lenH2):
	"""
	Given two haplotypes, find MEMs and output them in PAF format.
	"""
	with open(outfile_name, 'w') as outfile:
		if (syncmers_h1.get_length() > 0) and (syncmers_h2.get_length() > 0):

			cur_syncmer = syncmers_h1.get_syncmer(0)

			while True:
				# look for MEM starting at cur_syncmer in both haplotypes
				mem, mem_sequence, kmer_size = find_mem(syncmers_h1, syncmers_h2, cur_syncmer)
				if len(mem) > 0:
					# write PAF line
					query_name = chromosome1
					query_len = lenH1
					if syncmers_h1.is_reversed():
						# get_pos() is in the reverse-complement frame (L - p - k); convert
						# back to the original forward-strand position for the PAF output
						query_start = lenH1 - syncmers_h1.get_pos(mem[-1]) - kmer_size
					else:
						query_start = syncmers_h1.get_pos(mem[0])
					query_end = query_start + len(mem_sequence)
					strand = "-" if syncmers_h1.is_reversed() else "+"
					target_name = chromosome2
					target_len = lenH2
					if syncmers_h2.is_reversed():
						target_start = lenH2 - syncmers_h2.get_pos(mem[-1]) - kmer_size
					else:
						target_start = syncmers_h2.get_pos(mem[0])
					target_end = target_start + len(mem_sequence)
					residue_matches = str(len(mem_sequence))
					alignment_len = residue_matches
					mapping_qual = "255"
					cigar = "cg:Z:" + str(len(mem_sequence)) + "M"
					outfile.write("\t".join([query_name, str(query_len), str(query_start), str(query_end), strand, target_name, str(target_len), str(target_start), str(target_end), residue_matches, alignment_len, mapping_qual, cigar]) + "\n")
					# update current syncmer to the one after the MEM (if exists)
					if syncmers_h1.is_next_syncmer(mem[-1]):
						cur_syncmer = syncmers_h1.get_next_syncmer(mem[-1])
					else:
						# reached end of haplotype
						break
				else:
					# update current syncmer to next in haplotype
					if syncmers_h1.is_next_syncmer(cur_syncmer):
						cur_syncmer = syncmers_h1.get_next_syncmer(cur_syncmer)
					else:
						# reached end of haplotype
						break



if __name__== "__main__":
	parser = argparse.ArgumentParser(prog='generate-syncmer-paf.py', description=__doc__)
	parser.add_argument('-syncmersH1', metavar='SYNCMERS_H1', required=True, help='List of interval syncmers of haplotype 1 output by query command.') 
	parser.add_argument('-syncmersH2', metavar='SYNCMERS_H2', required=True, help='List of interval syncmers of haplotype 2 output by query command.')
	parser.add_argument('-chromosomeH1', metavar="CHROMOSOME_H1", required=True, help='Which chromosome in sample 1 the interval is located on.')
	parser.add_argument('-chromosomeH2', metavar="CHROMOSOME_H2", required=True, help='Which chromosome in sample 2 the interval is located on.')
	parser.add_argument('-lenH1', metavar="LEN_H1", type=int, required=True, help='Needs the total length of the respective chromosome of haplotype 1.')
	parser.add_argument('-lenH2', metavar="LEN_H2", type=int, required=True, help='Needs the total length of the respective chromosome of haplotype 2.')
	parser.add_argument("--reverse", action='store_true', default=False, help='Reverse first haplotype.')
	parser.add_argument('-outfile', metavar="OUTFILE", required=True, help="Name of the output file.")
	args = parser.parse_args()

	# parse the syncmer files
	if args.reverse:
		syncmers_h1 = HaplotypeIndex(args.syncmersH1, args.lenH1, True)
	else:
		syncmers_h1 = HaplotypeIndex(args.syncmersH1)
	syncmers_h2 = HaplotypeIndex(args.syncmersH2)

	# find MEMs and output PAF
	generate_paf(syncmers_h1, syncmers_h2, args.outfile, args.chromosomeH1, args.chromosomeH2, args.lenH1, args.lenH2)
