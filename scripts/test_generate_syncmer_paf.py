import unittest
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
DATA_FOLDER = HERE / "data"

from generate_syncmer_paf import reverse_complement 
from generate_syncmer_paf import HaplotypeIndex
from generate_syncmer_paf import find_mem
from generate_syncmer_paf import generate_paf

class TestReverseComplement(unittest.TestCase):

	def test_reverse_complement(self):
		seq1 = "ATGCGCTTGAGCTTCTCTCTTTTGaaaagtcTGA"
		seq1_rev = "TCAgacttttCAAAAGAGAGAAGCTCAAGCGCAT"
		self.assertTrue(reverse_complement(seq1) == seq1_rev)

		seq2 = "aatagatcgctcgctcgctcgatatcgcgattGGTAGATCgt"
		seq2_rev = "acGATCTACCaatcgcgatatcgagcgagcgagcgatctatt"
		self.assertTrue(reverse_complement(seq2) == seq2_rev)

		seq3 = "GGGGTTTTTTGATAGAGTTcccccNNNNnncgctnn"
		seq3_rev = "nnagcgnnNNNNgggggAACTCTATCAAAAAACCCC"
		self.assertTrue(reverse_complement(seq3) == seq3_rev)


class TestHaplotypeIndex(unittest.TestCase):

	def setUp(self):
		self.syncmer_file = DATA_FOLDER / "test_syncmers_hap1.tsv"
		self.hap_len = 34985800

	def test_haplotype_index_forward(self):
		"""
		Test HaplotypeIndex on example file.
		"""
		haplotype_index = HaplotypeIndex(self.syncmer_file)
		self.assertTrue(haplotype_index.get_length() == 30)

		expected_ids = []
		expected_positions = []
		expected_seqs = []
		for line in open(self.syncmer_file, 'r'):
			if line.startswith('syncmer_ID'):
				continue
			fields = line.strip().split()
			expected_ids.append(fields[0])
			expected_positions.append(fields[1])
			expected_seqs.append(fields[3])

		index = 0
		while index < haplotype_index.get_length():
			syncmer_id = haplotype_index.get_syncmer(index)
			self.assertTrue(syncmer_id == expected_ids[index])
			self.assertTrue(haplotype_index.get_pos(syncmer_id) == int(expected_positions[index]))
			self.assertTrue(haplotype_index.get_idx(syncmer_id) == index)
			self.assertTrue(haplotype_index.get_seq(syncmer_id) == expected_seqs[index])
			self.assertTrue(haplotype_index.contains_syncmer(syncmer_id))

			if index < haplotype_index.get_length() - 1:
				self.assertTrue(haplotype_index.is_next_syncmer(syncmer_id))
				self.assertTrue(haplotype_index.get_next_syncmer(syncmer_id) == expected_ids[index + 1])
			index += 1

	def test_haplotype_index_reverse(self):
		"""
		Test HaplotypeIndex on example file with reverse option
		"""
		haplotype_index = HaplotypeIndex(self.syncmer_file, self.hap_len, True)
		self.assertTrue(haplotype_index.get_length() == 30)

		expected_ids = []
		expected_positions = []
		expected_seqs = []
		for line in open(self.syncmer_file, 'r'):
			if line.startswith('syncmer_ID'):
				continue
			fields = line.strip().split()
			expected_ids.append(fields[0][1:] if fields[0].startswith('-') else "-" + fields[0] )
			expected_positions.append(str(self.hap_len - int(fields[1])))
			expected_seqs.append(reverse_complement(fields[3]))
		expected_ids = expected_ids[::-1]
		expected_positions = expected_positions[::-1]
		expected_seqs = expected_seqs[::-1]

		index = 0
		while index < haplotype_index.get_length():
			syncmer_id = haplotype_index.get_syncmer(index)
			self.assertTrue(syncmer_id == expected_ids[index])
			self.assertTrue(haplotype_index.get_pos(syncmer_id) == int(expected_positions[index]))
			self.assertTrue(haplotype_index.get_idx(syncmer_id) == index)
			self.assertTrue(haplotype_index.get_seq(syncmer_id) == expected_seqs[index])
			self.assertTrue(haplotype_index.contains_syncmer(syncmer_id))

			if index < haplotype_index.get_length() - 1:
				self.assertTrue(haplotype_index.is_next_syncmer(syncmer_id))
				self.assertTrue(haplotype_index.get_next_syncmer(syncmer_id) == expected_ids[index + 1])
			index += 1

	def test_haplotype_index_exception(self):
		"""
		Test come corner cases.
		"""
		haplotype_index = HaplotypeIndex(self.syncmer_file)
		self.assertTrue(haplotype_index.get_pos("nonexistent") is None)
		self.assertTrue(haplotype_index.get_idx("nonexistent") is None)
		self.assertTrue(haplotype_index.get_seq("nonexistent") is None)
		self.assertTrue(not haplotype_index.is_next_syncmer("nonexistent"))
		self.assertTrue(haplotype_index.is_next_syncmer("61105397"))
		self.assertTrue(not haplotype_index.is_next_syncmer("118792640"))
		self.assertTrue(haplotype_index.get_next_syncmer("118792640") is None)
		self.assertRaises(RuntimeError, haplotype_index.get_syncmer, 30)
		self.assertFalse(haplotype_index.contains_syncmer("nonexistent"))


class TestFindMem(unittest.TestCase):

	def setUp(self):
		self.syncmer_file1 = DATA_FOLDER / "test_syncmers_hap1.tsv"
		self.hap_len1 = 34985800
		self.syncmer_file2 = DATA_FOLDER / "test_syncmers_hap2.tsv"
		self.hap_len2 = 34976900
		self.expected_mem = DATA_FOLDER / "expected-mem.txt"

	def test_find_mem1(self):

		# syncmer is only present in first HT, so no result should be empty
		haplotype_index1 = HaplotypeIndex(self.syncmer_file1)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file2)
		mem1, mem_sequence1 = find_mem(haplotype_index1, haplotype_index2, "61105364")
		self.assertTrue(mem1 == [])
		self.assertTrue(mem_sequence1 == "")

		# there should be a mem returned
		expected_mem2 = [
			"61105370",
			"-61105371",
			"-61105372",
			"61105373",
			"61105374",
			"61105375",
			"-61105376",
			"61105377",
			"-61105378",
			"61105379",
			"61105380",
			"61105381",
			"61105382",
			"61105383",
			"61105384",
			"61105385",
			"61105386",
			"-61105387"
		]
		expected_mem_sequence2 = ""
		for line in open(self.expected_mem, 'r'):
			expected_mem_sequence2 += line.strip()
		mem2, mem_sequence2 = find_mem(haplotype_index1, haplotype_index2, "61105370")

		self.assertTrue(mem2 == expected_mem2)
		self.assertTrue(mem_sequence2 == expected_mem_sequence2)

		mem3, mem_sequence3 = find_mem(haplotype_index1, haplotype_index2, "-61105393")
		expected_mem3 = [
			"-61105393",
			"-61105394",
			"-61105395",
			"61105396",
			"61105397",
			"118792640"
		]
		expected_mem_sequence3 = "gatgtgggtaaggtgatcctcaaaaccggaaaaatggagctgagtccaaaaaaatctttgaatttcttctgtaagaataaaatgtttgttatatcctggattcatcttagtgaaaggcttggggcaatggatgataaagaatgatgctgtggcaggtcttctgtcgaacctccccctcttccccgctgccatgggttttcaaccatggatagtactacaccctgtgggaggcaattggaaatttgggagggaggagtttggggaggtcacaaca"
		self.assertTrue(mem3 == expected_mem3)
		self.assertTrue(mem_sequence3 == expected_mem_sequence3)


	def test_find_mem_reverse(self):

		haplotype_index1 = HaplotypeIndex(self.syncmer_file1, self.hap_len1, True)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file2, self.hap_len2, True)
		mem1, mem_sequence1 = find_mem(haplotype_index1, haplotype_index2, "-61105364")
		self.assertTrue(mem1 == [])
		self.assertTrue(mem_sequence1 == "")

		# there should be a mem returned
		expected_mem2 = [
			"-61105370",
			"61105371",
			"61105372",
			"-61105373",
			"-61105374",
			"-61105375",
			"61105376",
			"-61105377",
			"61105378",
			"-61105379",
			"-61105380",
			"-61105381",
			"-61105382",
			"-61105383",
			"-61105384",
			"-61105385",
			"-61105386",
			"61105387"
		]
		expected_mem_sequence2 = ""
		for line in open(self.expected_mem, 'r'):
			expected_mem_sequence2 += line.strip()
		mem2, mem_sequence2 = find_mem(haplotype_index1, haplotype_index2, "61105387")

		self.assertTrue(mem2 == expected_mem2[::-1])
		self.assertTrue(mem_sequence2 == reverse_complement(expected_mem_sequence2))

		mem3, mem_sequence3 = find_mem(haplotype_index1, haplotype_index2, "-118792640")
		expected_mem3 = [
			"61105393",
			"61105394",
			"61105395",
			"-61105396",
			"-61105397",
			"-118792640"
		]
		expected_mem_sequence3 = "gatgtgggtaaggtgatcctcaaaaccggaaaaatggagctgagtccaaaaaaatctttgaatttcttctgtaagaataaaatgtttgttatatcctggattcatcttagtgaaaggcttggggcaatggatgataaagaatgatgctgtggcaggtcttctgtcgaacctccccctcttccccgctgccatgggttttcaaccatggatagtactacaccctgtgggaggcaattggaaatttgggagggaggagtttggggaggtcacaaca"
		self.assertTrue(mem3 == expected_mem3[::-1])
		self.assertTrue(mem_sequence3 == reverse_complement(expected_mem_sequence3))


class TestGeneratePAF(unittest.TestCase):

	def setUp(self):
		self.syncmer_file1 = DATA_FOLDER / "test_syncmers_hap1.tsv"
		self.hap_len1 = 34985800
		self.syncmer_file2 = DATA_FOLDER / "test_syncmers_hap2.tsv"
		self.hap_len2 = 34976900
		self.expected_mem = DATA_FOLDER / "expected-mem.txt"

	def test_generate_paf(self):
		haplotype_index1 = HaplotypeIndex(self.syncmer_file1)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file2)

		with tempfile.TemporaryDirectory() as tmp_dir:
			output_file = Path(tmp_dir) / "output.paf"
			generate_paf(haplotype_index1, haplotype_index2, output_file, "sample1", "chr1", self.hap_len1, self.hap_len2)

			# check if PAF file looks as expected
			expected_lines = []
			for line in open(DATA_FOLDER / "expected-paf.paf", 'r'):
				expected_lines.append(line.strip())
			computed_lines = []
			for line in open(output_file, 'r'):
				computed_lines.append(line.strip())
			self.assertTrue(expected_lines == computed_lines)


if __name__ == "__main__":
	unittest.main()
