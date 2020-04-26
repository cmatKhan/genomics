import unittest
import neutral_rate


class MyTestCase(unittest.TestCase):
    def test_parseAlignmentFile(self):
        known_skud = '--TTATTTTTTTGCAGCCAACTTCACTACAGTTTAACCTACCTTATAGAAAGAGGTGAAAAGAATAGTTCATTAAATGAACATATATTTCTGAAAAGTTTACAGGTAAAATCACGTCCTTCGCTTAACTGCAAAAAAAAAAAGAGAACTCTTGCAT-GTACTTTCCCTTTTCATCCATCA'
        known_smik = '------TTTTTAACATCTATCCTTGATGA----TAGCATGTTAT-TACGAAATAATACAATGGAAAGTCCAATATCCCAGTGCTGAGTTCTGAAA-GTTTAAAAACATA-TGTCAACAACAAAACATCTTTCATATGTGGAAAAGCA-TATTGCTT----CTTTAACGTTTTACCTG--A'
        known_scer = '--TACTTTTTTGATGGCTACAATCATAAT---TTGACATACCTCATATGAA--TATGATGTGG-TAGTAAATTAAATAGGC--TGAATCCTGATA-GTTTGAGCCTGGGTTCTGACTATCCGAGAAGCGCACAGCAACATAACACCACCTTCTTTT-AGTCTTTTAGCGAAAAAGTAGGA'
        known_sbay = 'TTTTCTTTTTTAAGCGCCGGTCTAATTATCGGTTATCCTGCTACGCCAAAAGTGATATTATAG--AGTCTAATAAAGCA-------------AGATGTTCATATACAGC----AATCCCTCATTTGACT---------GGATAAATAGTCTTGCATCGTGCTTTGTTGTCGCACTTAGCA'
        known_alignment = '      *****     *     *          *  * *          **    *         ***  * **                  * * ***                            *            *  *  *   *    *    ****        *      *'

        known_alignment_dict = alignment_dict = {'Skud': known_skud, 'Smik' : known_smik, 'Scer': known_scer, 'Sbay': known_sbay, 'alignment': known_alignment}

        test_alignment_file_path = '/home/chase/code/cmatkhan/genomics/assignment12/data/test_align.aln'
        function_alignment_dict = neutral_rate.parseAlignmentFile(test_alignment_file_path)

        self.assertDictEqual(function_alignment_dict, known_alignment_dict)

    def test_calculateNeutralWobbleRate(self):
        pass


if __name__ == '__main__':
    unittest.main()
