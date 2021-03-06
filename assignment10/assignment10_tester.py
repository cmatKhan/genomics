import unittest
import call_orfs


class MyTestCase(unittest.TestCase):
    def test_complement(self):
        input_sequence =   'TCAATCAATCTAATCATGCATGCATGATTGATGACTAA'
        input_complement = 'AGTTAGTTAGATTAGTACGTACGTACTAACTACTGATT'
        function_complement = call_orfs.complement(input_sequence)
        self.assertEqual(input_complement, function_complement)

    def test_reverse(self):
        input_sequence = 'TCAATCAATCTAATCATGCATGCATGATTGATGACTAA'
        input_reverse = 'AATCAGTAGTTAGTACGTACGTACTAATCTAACTAACT'
        function_reverse = call_orfs.reverse(input_sequence)
        self.assertEqual(input_reverse, function_reverse)

    def test_reverseComplement(self):
        input_sequence = 'ATGC'
        input_reverse_complement = 'GCAT'
        function_reverse_complement = call_orfs.reverseComplement(input_sequence)
        self.assertEqual(input_reverse_complement, function_reverse_complement)

    def test_translateSequenceToAminoAcid(self):
        input_sequence = 'ATGATAGGGCACACCTATTAA'
        input_translated = 'MIGHTY'
        function_translated = call_orfs.translateSequenceToAminoAcid(input_sequence)
        self.assertEqual(input_translated, function_translated)

    def test_findLongestOpenReadingFrame(self):
        input_sequence = 'TCAATCAATCTAATCATGCATGCATGATTGATGACTAA'
        input_reverse_complement = call_orfs.reverseComplement(input_sequence)

        correct_return = ['ATGCATGCATGA']
        function_rf1 = call_orfs.findLongestOpenReadingFrame(input_sequence)
        self.assertEqual(correct_return, function_rf1)


        input_forward_rf1 = ['ATGCATGCATGA']
        function_rf1 = call_orfs.findLongestOpenReadingFrame(input_sequence)
        self.assertEqual(input_forward_rf1, function_rf1)

        input_forward_rf2 = ['ATGCATGATTGA']
        function_rf2 = call_orfs.findLongestOpenReadingFrame(input_sequence[1:])
        self.assertEqual(input_forward_rf2, function_rf2)

        input_forward_rf3 = ['ATGATTGATGACTAA']
        function_rf3 = call_orfs.findLongestOpenReadingFrame(input_sequence[2:])
        self.assertEqual(input_forward_rf3, function_rf3)

        input_reverse_rf1 = ['']
        function_reverse_rf1 = call_orfs.findLongestOpenReadingFrame(input_reverse_complement)
        self.assertEqual(input_reverse_rf1, function_reverse_rf1)

        input_reverse_rf2 = ['ATGCATGCATGA']
        function_reverse_rf2 = call_orfs.findLongestOpenReadingFrame(input_reverse_complement[1:])
        self.assertEqual(input_reverse_rf2, function_reverse_rf2)

        input_reverse_rf3 = 'ATGAAATAGATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGA'
        input_reverse_rf3_correct = ['ATGAAATAG', 'ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGA']
        function_reverse_rf3 = call_orfs.findLongestOpenReadingFrame(input_reverse_rf3)
        self.assertEqual(input_reverse_rf3_correct, function_reverse_rf3)


if __name__ == '__main__':
    unittest.main()
