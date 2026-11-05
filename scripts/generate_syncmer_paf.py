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

		syncmer_idx = 0
		for line in open(syncmer_file, 'r'):
			if line.startswith("syncmer_ID"):
				continue
			fields = line.strip().split()

			if reverse:
				# reverse orientation of the syncmer_id
				syncmer_id = fields[0][1:] if fields[0].startswith("-") else "-" + fields[0]
				self._syncmers.append(syncmer_id)
				assert not syncmer_id in self._syncmer_to_pos
				# total length minus position
				pos = haplotype_len - int(fields[1])
				# reverse complement of syncmer sequence
				sequence = reverse_complement(fields[3])
				assert not syncmer_id in self._syncmer_to_pos
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
		return mem, mem_sequence
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
			assert overlap1 == overlap2
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
	return mem, mem_sequence


if __name__== "__main__":
	parser = argparse.ArgumentParser(prog='generate-syncmer-paf.py', description=__doc__)
	parser.add_argument('-syncmersH1', metavar='SYNCMERS_H1', required=True, help='List of interval syncmers of haplotype 1 output by query command.') 
	parser.add_argument('-syncmersH2', metavar='SYNCMERS_H2', required=True, help='List of interval syncmers of haplotype 2 output by query command.')
	parser.add_argument('-chromosome', metavar="CHROMOSOME", required=True, help='Which chromosome the interval is located on.')
	parser.add_argument('-sample', metavar="SAMPLE", required=False, default="", help='Name of the sample (to be specified in output)')
	parser.add_argument('-reverse', metavar="LEN_H1", type=int, required=False, help='Reverse first haplotype. Needs the total length of the respective chromosome of haplotype 1.')
	args = parser.parse_args()

	# parse the syncmer files
	if args.reverse:
		syncmers_h1 = HaplotypeIndex(args.syncmersH1, args.reverse, True)
	else:
		syncmers_h1 = HaplotypeIndex(args.syncmersH1)
	syncmers_h2 = HaplotypeIndex(args.syncmersH2)

	if (syncmers_h1.get_length() > 0) and (syncmers_h2.get_length() > 0):

		cur_syncmer = syncmers_h1.get_syncmer(0)
		cur_idx = 0

		while syncmers_h1.is_next_syncmer(cur_syncmer):
			# look for MEM starting at cur_syncmer in both haplotypes
			mem, mem_sequence = find_mem(syncmers_h1, syncmers_h2, cur_syncmer)
			if len(mem) > 0:
				# write PAF line
				query_name = args.sample + "#hap1#" + args.chromosome
				query_len = len(mem_sequence)
				query_start = syncmers_h1.get_pos(mem[0])
				query_end = query_start + query_len
				strand = "+"
				target_name = args.sample + "#hap2" + args.chromosome
				target_len = len(mem_sequence)
				target_start = syncmers_h2.get_pos(mem[0])
				target_end = target_start + target_len
				residue_matches = len(mem_sequence)
				alignment_len = residue_matches
				mapping_qual = "255"
				cigar = "cg:Z:" + len(mem_sequence) + "M"
				print("\t".join([query_name, query_len, query_start, query_end, strand, target_name, target_len, target_start, target_end, residue_matches, alignment_len, mapping_qual, cigar]))




