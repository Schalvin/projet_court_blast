import unittest
import blast


class TestBlast(unittest.TestCase):
    sequence = 'AGTKKLPMST'
    words = {'AGT': [0], 'GTK': [1], 'TKL': [2],
             'KLP': [3], 'LPM': [4], 'PMS': [5], 'MST': [6]}
    sequence2 = "KAGTKKLPMSTL"

    def test_seq_to_words(self):
        result = blast.seq_to_words(sequence)
        self.assertEqual(result, words)

    def test_double_hit_finder(self):
        result = blast.double_hit_finder(words, sequence2)
        self.assertEqual(result, words)
