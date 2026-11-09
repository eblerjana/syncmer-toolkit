import sys
import subprocess
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
		# must be >= (position + len(sequence)) of the last syncmer in the file (34985804)
		self.hap_len = 34985804

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
			# a syncmer spanning [p, p+k) on the forward strand spans [L-p-k, L-p) on the
			# reverse-complement strand, so its start there is L-p-k, not L-p
			expected_positions.append(str(self.hap_len - int(fields[1]) - len(fields[3])))
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

	def test_haplotype_index_reverse_length_too_small(self):
		"""
		A haplotype_len that does not cover the full extent (position + len(sequence))
		of the last syncmer in the file must raise, instead of silently producing a
		negative position.
		"""
		self.assertRaises(RuntimeError, HaplotypeIndex, self.syncmer_file, self.hap_len - 1, True)
		# self.hap_len is exactly the minimum valid value, so this must NOT raise
		HaplotypeIndex(self.syncmer_file, self.hap_len, True)


class TestFindMem(unittest.TestCase):

	def setUp(self):
		self.syncmer_file1 = DATA_FOLDER / "test_syncmers_hap1.tsv"
		# must be >= (position + len(sequence)) of the last syncmer in the file (34985804)
		self.hap_len1 = 34985804
		self.syncmer_file2 = DATA_FOLDER / "test_syncmers_hap2.tsv"
		self.hap_len2 = 34976900
		self.expected_mem = DATA_FOLDER / "expected-mem.txt"

		self.syncmer_file_repeat1 =  DATA_FOLDER / "test_syncmers_hap1_repeat.tsv"
		self.syncmer_file_repeat2 =  DATA_FOLDER / "test_syncmers_hap2_repeat.tsv"


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

	def test_find_mem_sequence_mismatch(self):
		"""
		If the two haplotypes disagree about the sequence of the (shared-ID) starting
		syncmer - which would also mean its length differs, since unequal-length
		strings can never be equal - find_mem must fail loudly rather than silently
		returning a wrong MEM. Note: the explicit "kmer sizes ... do not match"
		RuntimeError further down in find_mem is actually unreachable for this via
		this entry point under normal (assertion-enabled) execution, because the
		"assert ...get_seq(syncmer_id) == ...get_seq(syncmer_id)" check right above it
		always fires first for the starting syncmer - any length mismatch is also a
		content mismatch. What actually gets raised today is an AssertionError.
		"""
		with tempfile.TemporaryDirectory() as tmp_dir:
			f1 = Path(tmp_dir) / "h1.tsv"
			f2 = Path(tmp_dir) / "h2.tsv"
			with open(f1, 'w') as f:
				f.write("syncmer_ID\tposition\tcount\tsequence\n")
				f.write("999\t10\t1\tAAAAAAAAAA\n")
			with open(f2, 'w') as f:
				f.write("syncmer_ID\tposition\tcount\tsequence\n")
				f.write("999\t10\t1\tAAAAA\n")

			haplotype_index1 = HaplotypeIndex(f1)
			haplotype_index2 = HaplotypeIndex(f2)
			self.assertRaises(AssertionError, find_mem, haplotype_index1, haplotype_index2, "999")

	def test_find_mem_repeat(self):
		"""
        In repetitive areas, it can happen that two haplotypes traverse the same syncmeres,
        but the overlaps between neighboring syncmers differs. This test checks if
        the algorithm correctly breaks MEMs at such positions.
		"""

		haplotype_index1 = HaplotypeIndex(self.syncmer_file_repeat1)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file_repeat2)
		mem, mem_sequence = find_mem(haplotype_index1, haplotype_index2, "61117411")

		# only a single syncmer in MEM, as overlaps to next syncmer differs in both files
		expected_mem1 = ["61117411"]
		expected_mem_sequence1 = "aaatataaatatttataaaaatatttatataaatatttataaaaatatttatataaatattta"
		self.assertTrue(mem == expected_mem1)
		self.assertTrue(mem_sequence == expected_mem_sequence1)

		mem2, mem_sequence2 = find_mem(haplotype_index1, haplotype_index2, "61117414")
		expected_mem2 = ["61117414", "61117415"]
		expected_mem_sequence2 = "ttataaaaatatttatataaatatttataaaaatatttatataaatatttataaaaaatataaatattagccatgtttggtggctcacgcctgtaa"
		self.assertTrue(mem2 == expected_mem2)
		self.assertTrue(mem_sequence2 == expected_mem_sequence2)



class TestGeneratePAF(unittest.TestCase):

	def setUp(self):
		self.syncmer_file1 = DATA_FOLDER / "test_syncmers_hap1.tsv"
		# must be >= (position + len(sequence)) of the last syncmer in the file (34985804)
		self.hap_len1 = 34985804
		self.syncmer_file2 = DATA_FOLDER / "test_syncmers_hap2.tsv"
		self.hap_len2 = 34976900
		self.expected_mem = DATA_FOLDER / "expected-mem.txt"

		self.syncmer_file_repeat1 =  DATA_FOLDER / "test_syncmers_hap1_repeat.tsv"
		self.syncmer_file_repeat2 =  DATA_FOLDER / "test_syncmers_hap2_repeat.tsv"
		self.hap_len_repeat1 = 36270000
		self.hap_len_repeat2 = 36270000


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

	def test_generate_paf_reverse_strand(self):
		"""
		When haplotype 1 is reversed, the strand column must be "-", and (since matching
		requires the syncmer IDs, including sign, to agree) haplotype 2 must be reversed
		too for any MEM to be found at all.
		"""
		haplotype_index1 = HaplotypeIndex(self.syncmer_file1, self.hap_len1, True)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file2, self.hap_len2, True)

		with tempfile.TemporaryDirectory() as tmp_dir:
			output_file = Path(tmp_dir) / "output.paf"
			generate_paf(haplotype_index1, haplotype_index2, output_file, "sample1", "chr1", self.hap_len1, self.hap_len2)

			expected_lines = [
				"sample1#hap1#chr1\t34985804\t0\t274\t-\tsample1#hap2#chr1\t34976900\t38\t312\t274\t274\t255\tcg:Z:274M",
				"sample1#hap1#chr1\t34985804\t385\t904\t-\tsample1#hap2#chr1\t34976900\t423\t942\t519\t519\t255\tcg:Z:519M",
			]
			computed_lines = []
			for line in open(output_file, 'r'):
				computed_lines.append(line.strip())
			self.assertTrue(computed_lines == expected_lines)
			for line in computed_lines:
				self.assertTrue(line.split("\t")[4] == "-")

	def test_generate_paf_empty_haplotype(self):
		"""
		generate_paf must still create the output file (empty) when one of the
		haplotypes has no syncmers at all, instead of not creating it.
		"""
		haplotype_index1 = HaplotypeIndex(self.syncmer_file1)

		with tempfile.TemporaryDirectory() as tmp_dir:
			empty_file = Path(tmp_dir) / "empty.tsv"
			with open(empty_file, 'w') as f:
				f.write("syncmer_ID\tposition\tcount\tsequence\n")
			haplotype_index_empty = HaplotypeIndex(empty_file)
			self.assertTrue(haplotype_index_empty.get_length() == 0)

			output_file = Path(tmp_dir) / "output.paf"
			generate_paf(haplotype_index1, haplotype_index_empty, output_file, "sample1", "chr1", self.hap_len1, 0)

			self.assertTrue(output_file.exists())
			self.assertTrue(open(output_file, 'r').read() == "")

	def test_generate_paf_no_shared_syncmers(self):
		"""
		Two non-empty haplotypes that do not share any syncmer at all must still
		produce an (empty) output file, with generate_paf's main loop advancing
		all the way through haplotype 1 without ever finding a MEM.
		"""
		haplotype_index1 = HaplotypeIndex(self.syncmer_file1)

		with tempfile.TemporaryDirectory() as tmp_dir:
			disjoint_file = Path(tmp_dir) / "disjoint.tsv"
			with open(disjoint_file, 'w') as f:
				f.write("syncmer_ID\tposition\tcount\tsequence\n")
				f.write("999999991\t5\t1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
				f.write("999999992\t80\t1\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
			haplotype_index2 = HaplotypeIndex(disjoint_file)
			self.assertTrue(haplotype_index2.get_length() == 2)

			output_file = Path(tmp_dir) / "output.paf"
			generate_paf(haplotype_index1, haplotype_index2, output_file, "sample1", "chr1", self.hap_len1, 200)

			self.assertTrue(output_file.exists())
			self.assertTrue(open(output_file, 'r').read() == "")

	def test_generate_paf_repeat(self):
		"""
		The haplotypes traverse the same syncmers, but overlaps between
		syncmers is different and therefore underlying sequences differ.
		Check if MEM is broken into pieces in such cases.
		"""
		haplotype_index1 = HaplotypeIndex(self.syncmer_file_repeat1)
		haplotype_index2 = HaplotypeIndex(self.syncmer_file_repeat2)

		with tempfile.TemporaryDirectory() as tmp_dir:
			repeat_output = Path(tmp_dir) / "repeat.paf"
			generate_paf(haplotype_index1, haplotype_index2, repeat_output, "sample1", "chr9", self.hap_len_repeat1, self.hap_len_repeat2)

			self.assertTrue(repeat_output.exists())
			computed_lines = []
			for line in open(repeat_output, 'r'):
				computed_lines.append(line.strip())
			self.assertTrue(len(computed_lines) == 2)

			first_mem = "sample1#hap1#chr9\t36270000\t35262222\t35262285\t+\tsample1#hap2#chr9\t36270000\t35271097\t35271160\t63\t63\t255\tcg:Z:63M"
			second_mem = "sample1#hap1#chr9\t36270000\t35262282\t35262378\t+\tsample1#hap2#chr9\t36270000\t35271133\t35271229\t96\t96\t255\tcg:Z:96M"

			self.assertTrue(first_mem == computed_lines[0])
			self.assertTrue(second_mem == computed_lines[1])


class TestCLI(unittest.TestCase):
	"""
	Exercises the actual command-line entry point (__main__/argparse), rather than
	importing functions directly - this is the only way to catch wiring bugs in the
	argument parser itself (e.g. a typo'd parser method that only fails when the
	script is actually run, not when its functions are imported).
	"""

	def setUp(self):
		self.script = HERE / "generate_syncmer_paf.py"
		self.syncmer_file1 = DATA_FOLDER / "test_syncmers_hap1.tsv"
		self.syncmer_file2 = DATA_FOLDER / "test_syncmers_hap2.tsv"

	def test_cli_forward(self):
		with tempfile.TemporaryDirectory() as tmp_dir:
			output_file = Path(tmp_dir) / "output.paf"
			result = subprocess.run(
				[sys.executable, str(self.script),
				 "-syncmersH1", str(self.syncmer_file1),
				 "-syncmersH2", str(self.syncmer_file2),
				 "-chromosome", "chr1",
				 "-lenH1", "34985804",
				 "-lenH2", "34976900",
				 "-outfile", str(output_file),
				 "-sample", "sample1"],
				capture_output=True, text=True,
			)
			self.assertTrue(result.returncode == 0, msg=result.stderr)

			expected_lines = []
			for line in open(DATA_FOLDER / "expected-paf.paf", 'r'):
				expected_lines.append(line.strip())
			computed_lines = []
			for line in open(output_file, 'r'):
				computed_lines.append(line.strip())
			self.assertTrue(computed_lines == expected_lines)

	def test_cli_reverse_flag_does_not_crash(self):
		with tempfile.TemporaryDirectory() as tmp_dir:
			output_file = Path(tmp_dir) / "output.paf"
			result = subprocess.run(
				[sys.executable, str(self.script),
				 "-syncmersH1", str(self.syncmer_file1),
				 "-syncmersH2", str(self.syncmer_file2),
				 "-chromosome", "chr1",
				 "-lenH1", "34985804",
				 "-lenH2", "34976900",
				 "-outfile", str(output_file),
				 "--reverse"],
				capture_output=True, text=True,
			)
			self.assertTrue(result.returncode == 0, msg=result.stderr)
			self.assertTrue(output_file.exists())


if __name__ == "__main__":
	unittest.main()
