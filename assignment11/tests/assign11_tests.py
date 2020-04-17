import unittest
import filter_variants
import count_barcodes

class MyTestCase(unittest.TestCase):
    count_barcodes_test_data_path = '/home/chase/code/cmatkhan/genomics/assignment11/data/test_data/count_barcodes_test_data.txt'
    filtered_variant_barcode_path = '/home/chase/code/cmatkhan/genomics/assignment11/filtered_variant_to_barcode.tsv'
    fastq_test_data_path = '/home/chase/code/cmatkhan/genomics/assignment11/data/test_data/fastq_test_data.fq'

    def test_createLineGenerator(self):
        line_generator = filter_variants.createLineGenerator(self.count_barcodes_test_data_path)
        func_header = next(line_generator)
        known_header = 'Variant_ID\tREF_barcode\tALT_barcode\n'
        func_line = next(line_generator)
        known_line = 'rs1057902\tCCTCGCTGG:TTCTACACA:CAGGATGCT:TACTCAACC:GCAATTCAG:TAGTGCCTT:GGTATGAGG:ATTCACGAG:AGGTTGGTC\tACTGTCCAT:TAACACAGT:ACTGCGAAC:GAGAGGCAC:CCGTGTCTA:CTCAGCCGT\n'
        self.assertEqual(func_header, known_header)
        self.assertEqual(func_line, known_line)

    def test_createVariantBarcodeDict(self):
        line_generator = filter_variants.createLineGenerator(self.count_barcodes_test_data_path)
        func_dict = filter_variants.createVariantBarcodeDict(line_generator)

        known_dict = {'rs1057902':   ['CCTCGCTGG','TTCTACACA','CAGGATGCT','TACTCAACC','GCAATTCAG', 'TAGTGCCTT','GGTATGAGG',
                                      'ATTCACGAG','AGGTTGGTC','ACTGTCCAT', 'TAACACAGT','ACTGCGAAC', 'GAGAGGCAC','CCGTGTCTA', 'CTCAGCCGT'],
                      'rs116386160': ['CACGATGAC','TGTAGTTGG','AAGCATGAA','TCAAGGCAC','TGTAGATCG','CTCCTGATG','TAGGTCTAG',
                                       'TCTTGCCTC','TGACTATAA','GAATATCAC','CAGAGTAGG','CACTTCACG','CACCTTAAC','TGATAGATC',
                                       'GATGACGGA','ATCCGATAG','TAAGGTGAA','CATTAAGTA', 'TGGACAATT','TCGCTGTTG', 'AACGCAACC']}
        for key, value in func_dict.items():
            for list in value:
                for item in list:
                    self.assertIn(item, known_dict[key])
        self.assertEqual(len(func_dict['rs1057902']), 2)
        self.assertEqual(len(func_dict['rs116386160']), 2)

    def test_countBarcodes(self):
        correct_dict = {'CGTGCGAAT': 1}
        barcode_list = count_barcodes.createBarcodeList(self.filtered_variant_barcode_path)
        func_dict = count_barcodes.countBarcodes(barcode_list, self.fastq_test_data_path)
        self.assertDictEqual(correct_dict, func_dict)


if __name__ == '__main__':
    unittest.main()
