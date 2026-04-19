# tests/test_main.py
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from main import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    parse_blast_results,
    build_parser,
)



class TestDNASequence(unittest.TestCase):

    def test_check_alphabet_valid(self):
        """valid alphabet check"""
        seq = DNASequence("ATGC")
        self.assertTrue(seq.check_alphabet())

    def test_check_alphabet_invalid(self):
        """invalid alphabet check"""
        seq = DNASequence("AUGC")
        self.assertFalse(seq.check_alphabet())

    def test_transcribe(self):
        """transribe check"""
        dna = DNASequence("ATGCTT")
        rna = dna.transcribe()
        self.assertIsInstance(rna, RNASequence)
        self.assertEqual(str(rna), "AUGCUU")

    def test_reverse_complement(self):
        """reverse complement check"""
        dna = DNASequence("ATGC")
        self.assertEqual(str(dna.reverse_complement()), "GCAT")


class TestAminoAcidSequence(unittest.TestCase):

    def test_check_alphabet_valid(self):
        """valid aminoacid check"""
        aa = AminoAcidSequence("ACDEFG")
        self.assertTrue(aa.check_alphabet())

    def test_check_alphabet_invalid_raises_false(self):
        """invalid aminoacid check"""
        aa = AminoAcidSequence("ACDB")
        self.assertFalse(aa.check_alphabet())

    def test_hydrophobic_ratio(self):
        """hydrophobic ratio check"""
        # A, V, L — гидрофобные (3 из 4)
        aa = AminoAcidSequence("AVLK")
        self.assertAlmostEqual(aa.hydrophobic_ratio(), 0.75)



class TestParseBlastResults(unittest.TestCase):

    def _make_blast_file(self, content: str):
        """blast I/O check"""
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
        tmp.write(content)
        tmp.close()
        return tmp.name

if __name__ == "__main__":
    unittest.main(verbosity=2)