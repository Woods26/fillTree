#!/usr/bin/env python
# coding=utf-8
"""Testing suite for fudger.py"""

import unittest
from unittest import TestCase
import hashlib


class IndependentUnitTests(TestCase):
    """Unit tests that are not dependent on any other units"""

    def test_parse_fasta_filenames_from_directory(self):
        from fudger import parse_fasta_filenames_from_directory
        fasta_path = "input/orthoExon_fasta"
        calculated = parse_fasta_filenames_from_directory(fasta_path)
        expected = {'orth10136_1834-2121': 'input/orthoExon_fasta/orth10136_1834-2121.full.fasta',
                    'orth10223_1638-1932': 'input/orthoExon_fasta/orth10223_1638-1932.full.fasta',
                    'orth10262_1395-1632': 'input/orthoExon_fasta/orth10262_1395-1632.full.fasta',
                    'orth10262_574-1105': 'input/orthoExon_fasta/orth10262_574-1105.full.fasta',
                    'orth10297_2517-2710': 'input/orthoExon_fasta/orth10297_2517-2710.full.fasta',
                    'orth10315_1281-1727': 'input/orthoExon_fasta/orth10315_1281-1727.full.fasta',
                    'orth10339_505-906': 'input/orthoExon_fasta/orth10339_505-906.full.fasta',
                    'orth10362_263-601': 'input/orthoExon_fasta/orth10362_263-601.full.fasta',
                    'orth10375_1062-1241': 'input/orthoExon_fasta/orth10375_1062-1241.full.fasta',
                    'orth10375_1291-1484': 'input/orthoExon_fasta/orth10375_1291-1484.full.fasta',
                    'orth10394_1043-1328': 'input/orthoExon_fasta/orth10394_1043-1328.full.fasta',
                    'orth10395_1616-1808': 'input/orthoExon_fasta/orth10395_1616-1808.full.fasta',
                    'orth10410_311-493': 'input/orthoExon_fasta/orth10410_311-493.full.fasta',
                    'orth10410_543-761': 'input/orthoExon_fasta/orth10410_543-761.full.fasta',
                    'orth10425_453-875': 'input/orthoExon_fasta/orth10425_453-875.full.fasta',
                    'orth2374_1109-1760': 'input/orthoExon_fasta/orth2374_1109-1760.full.fasta',
                    'orth2652_1294-1523': 'input/orthoExon_fasta/orth2652_1294-1523.full.fasta',
                    'orth2652_1938-2132': 'input/orthoExon_fasta/orth2652_1938-2132.full.fasta',
                    'orth2782_175-1267': 'input/orthoExon_fasta/orth2782_175-1267.full.fasta'}
        self.assertEqual(calculated, expected)

    def test_parse_fasta_file(self):
        from fudger import parse_fasta_file
        file_name = 'input/orthoExon_fasta/orth10136_1834-2121.full.fasta'
        calculated = str(parse_fasta_file(file_name))
        expected = ("{'Bcur': SeqRecord(seq=Seq('ggcaaaagatcttgccggatgtgctgcagcatatcggcgtacggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bcur', name='Bcur', description='Bcur cds274_2', dbxrefs=[]), "
                    "'Bdor': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bdor', name='Bdor', description='Bdor cds10257_3', dbxrefs=[]), "
                    "'Bole': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcaacatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bole', name='Bole', description='Bole cds4643_2', dbxrefs=[]), "
                    "'Ccap': SeqRecord(seq=Seq('ggcaaaagatcttgccagacgtactacaacatatcggtgtgcgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Ccap', name='Ccap', description='Ccap cds6222_3', dbxrefs=[]), "
                    "'Afra': SeqRecord(seq=Seq('ggcaaaagatattgcctgatgtgctgcagcatatcggggtacgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Afra', name='Afra', "
                    "description='Afra TR39090|c2_g2_i1|m.26898', dbxrefs=[]), "
                    "'Asus': SeqRecord(seq=Seq('ggcaaaagatattgcctgatgtgctgcagcatatcggggtacgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Asus', name='Asus', "
                    "description='Asus TR39090|C2_G2_I1|M.26898_R0', dbxrefs=[]), "
                    "'Bcor': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bcor', name='Bcor', "
                    "description='Bcor TR9663|c1_g1_i1|m.7449', dbxrefs=[]), "
                    "'Bzon': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcaacatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bzon', name='Bzon', "
                    "description='Bzon TR38546|c0_g8_i1|m.23458', dbxrefs=[])}")
        self.assertEqual(calculated, expected)

    def test_create_unique_dir(self):
        from fudger import create_unique_dir
        import os
        path = "test_dir"
        os.mkdir(path)
        calculated = create_unique_dir(path)
        expected = path
        self.assertEqual(calculated, expected)

        with open(os.path.join(path, "test.txt"), "wt") as f:
            f.write("some text")
        calculated_01 = create_unique_dir(path)
        expected_01 = "test_dir_01"
        self.assertEqual(calculated_01, expected_01)

        with open(os.path.join(calculated_01, "test.txt"), "wt") as f:
            f.write("some text")
        with self.assertRaises(Exception):
            create_unique_dir(path, limit=1)

        os.remove(os.path.join(calculated_01, "test.txt"))
        os.rmdir(calculated_01)
        os.remove(os.path.join(path, "test.txt"))
        os.rmdir(path)


class DependentUnitTests(TestCase):
    """Unit tests that depend on other tested units"""

    def test_closest_available_relatives(self):
        from fudger import closest_available_relatives
        from fudger import parse_fasta_file
        from Bio import Phylo
        tree = Phylo.read("input/tapir_ref_959genes.tre", 'newick')
        terminals_dict = {leaf.name: leaf for leaf in tree.get_terminals()}
        fasta = parse_fasta_file('input/orthoExon_fasta/orth10136_1834-2121.full.fasta')
        avail_sp_set = set(fasta.keys())
        sp = 'Bmin'
        calculated = closest_available_relatives(sp, avail_sp_set, terminals_dict, tree)
        expected = [('Bole', 8.895795999999999), ('Bcur', 14.586389999999998),
                    ('Bdor', 17.182138000000002), ('Bzon', 17.185076000000002),
                    ('Bcor', 17.186880000000002), ('Ccap', 36.786696), ('Afra', 113.23144299999998),
                    ('Asus', 113.23283799999999)]

        self.assertEqual(calculated, expected)

    def test_write_fasta_file_with_substitutetions(self):
        from fudger import write_fasta_file_with_substitutetions
        from fudger import parse_fasta_file
        from Bio import Phylo
        import os

        tree = Phylo.read("input/tapir_ref_959genes.tre", 'newick')

        sp_list = [cl.name for cl in tree.get_terminals()]
        terminals_dict = {leaf.name: leaf for leaf in tree.get_terminals()}
        fasta = parse_fasta_file('input/orthoExon_fasta/orth10136_1834-2121.full.fasta')

        ortho = "orth10136_1834-2121"
        output_path = ""

        write_fasta_file_with_substitutetions(ortho, fasta, sp_list, output_path, terminals_dict, tree)

        calculated = str(parse_fasta_file(ortho + ".padded.fasta"))
        expected = ("{'Aobl': SeqRecord(seq=Seq('ggcaaaagatattgcctgatgtgctgcagcatatcggggtacgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Aobl', name='Aobl', "
                    "description='Aobl subbed-missing-with_Afra distance_14.227426000000001', dbxrefs=[]), "
                    "'Afra': SeqRecord(seq=Seq('ggcaaaagatattgcctgatgtgctgcagcatatcggggtacgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Afra', name='Afra', "
                    "description='Afra TR39090|c2_g2_i1|m.26898', dbxrefs=[]), "
                    "'Asus': SeqRecord(seq=Seq('ggcaaaagatattgcctgatgtgctgcagcatatcggggtacgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Asus', name='Asus', "
                    "description='Asus TR39090|C2_G2_I1|M.26898_R0', dbxrefs=[]), "
                    "'Bcor': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bcor', name='Bcor', "
                    "description='Bcor TR9663|c1_g1_i1|m.7449', dbxrefs=[]), "
                    "'Bzon': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcaacatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bzon', name='Bzon', "
                    "description='Bzon TR38546|c0_g8_i1|m.23458', dbxrefs=[]), "
                    "'Bdor': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bdor', name='Bdor', "
                    "description='Bdor cds10257_3', dbxrefs=[]), "
                    "'Blat': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Blat', name='Blat', "
                    "description='Blat subbed-missing-with_Bdor distance_0.048722999999999995', dbxrefs=[]), "
                    "'Bjar': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bjar', name='Bjar', "
                    "description='Bjar subbed-missing-with_Bdor distance_5.712655', dbxrefs=[]), "
                    "'Btry': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcagcatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Btry', name='Btry', "
                    "description='Btry subbed-missing-with_Bdor distance_5.711662', dbxrefs=[]), "
                    "'Bole': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcaacatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bole', name='Bole', "
                    "description='Bole cds4643_2', dbxrefs=[]), "
                    "'Bmin': SeqRecord(seq=Seq('ggcaaaagatcttgccagatgtgctgcaacatatcggcgtgcggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bmin', name='Bmin', "
                    "description='Bmin subbed-missing-with_Bole distance_8.895795999999999', dbxrefs=[]), "
                    "'Bcur': SeqRecord(seq=Seq('ggcaaaagatcttgccggatgtgctgcagcatatcggcgtacggccagccaaca...cag', "
                    "SingleLetterAlphabet()), id='Bcur', name='Bcur', "
                    "description='Bcur cds274_2', dbxrefs=[]), "
                    "'Ccap': SeqRecord(seq=Seq('ggcaaaagatcttgccagacgtactacaacatatcggtgtgcgaccagccaata...cag', "
                    "SingleLetterAlphabet()), id='Ccap', name='Ccap', "
                    "description='Ccap cds6222_3', dbxrefs=[])}")

        os.remove(ortho + ".padded.fasta")

        self.assertEqual(calculated, expected)


class RegressionTesting(TestCase):
    """Test of main. Checks that output matches historical expected values."""

    def test_main(self):
        from fudger import main
        import shutil
        import os

        fasta_path = "input/orthoExon_fasta"
        output_path = "main_test_dir/"
        tree_fn = "input/tapir_ref_959genes.tre"
        seed = 1159072149485953569
        sort_tree = True
        main(fasta_path, output_path, tree_fn, seed, sort_tree)

        files = [f for f in os.listdir(output_path) if os.path.isfile(os.path.join(output_path, f))]

        calculated = {file: filehash(os.path.join(output_path, file)) for file in files}
        expected = {
            'orth10136_1834-2121.padded.fasta': 'b370310637c146fda4d28af92b4df6391e0103aa5be440586b86c89249f8d1df',
            'orth10223_1638-1932.padded.fasta': '3e5e0f7bd029002a1d7259132f459ee0a5dd997b3363e9d20e91095cf9ae3dbf',
            'orth10262_1395-1632.padded.fasta': '37cd77e3c195ff7965860f9c947dd016437ee737ace76f82d8b534cb7a0184be',
            'orth10262_574-1105.padded.fasta': '223bc7e2bf9f2c6f1c3c2e4f1b0ad158d9c5e0b5f0fc12b2a33f1fd2bf2cf494',
            'orth10297_2517-2710.padded.fasta': 'ad421e61e94aeef2234ca113d013128503e5d502766b95b107672dfff183b10f',
            'orth10315_1281-1727.padded.fasta': '3e0045fe32c04019d669b5f0647cd0ce18014345213353756f0bfa3ee76557d2',
            'orth10339_505-906.padded.fasta': '64baf2c9b4300cd305dd0ae3ee178068a510b8a9b7b1dc401c36fc4da8592af6',
            'orth10362_263-601.padded.fasta': 'ef39ef36a82bd67fc918c28179a43cdbeaaf03cd3b5802e0efd3880b43aceda6',
            'orth10375_1062-1241.padded.fasta': 'a6774610a0b29fc9adb9c2ced7e2420e9f47150245ad5cc41b1d442e7c532246',
            'orth10375_1291-1484.padded.fasta': '7463f2be28d85d60805eac51e4186c869d6133325ee828723cd642c56ff4e09b',
            'orth10394_1043-1328.padded.fasta': 'a651bf00be10c92e43b8875029310789b27ecbe2fda9c45f2660a708fcb1b814',
            'orth10395_1616-1808.padded.fasta': '01f7637e05723c26d4960b4007c046a723d58543716f892a9a83c508ae8f7747',
            'orth10410_311-493.padded.fasta': '9f1d6d2d3a625bf46aac1fcee552086a495e990d8eea417ea20ca8dc7f910519',
            'orth10410_543-761.padded.fasta': '89fd46d5585a3bca0368c161325b1f9232a66b09ec45e9907c5364b5893b9226',
            'orth10425_453-875.padded.fasta': '4c933ebe49c8f6de4ce7e4e82c79f0fb8dc550a4736fdeb038b4558d4aeface4',
            'orth2374_1109-1760.padded.fasta': '71174408d061765d9489bff687a500f991684e9dc7380fb506b8f286d73a3873',
            'orth2652_1294-1523.padded.fasta': '0cda52976dc287fa55e4544bb4b44e5613f8be03427cb0876c16a690942307a8',
            'orth2652_1938-2132.padded.fasta': '18e7927e5e1728200090657c09217221583e294f03b2fa9ad1f6914ea34dc874',
            'orth2782_175-1267.padded.fasta': '1d57f77e383e87cfd1eec1aaf2dbfed8a0416c6788a3bcae47b6d2c8dd3c8b2b'}

        self.assertEqual(calculated, expected)

        shutil.rmtree(output_path)


def filehash(filepath, blocksize=4096):
    """ Return the hash hexdigest for the file `filepath', processing the file
    by chunk of `blocksize'.

    function found in tutorial at:
    https://thomassileo.name/blog/2013/12/12/tracking-changes-in-directories-with-python/

    :type filepath: str
    :param filepath: Path to file

    :type blocksize: int
    :param blocksize: Size of the chunk when processing the file

    """
    sha = hashlib.sha256()
    with open(filepath, 'rb') as fp:
        while True:
            data = fp.read(blocksize)
            if data:
                sha.update(data)
            else:
                break
    return sha.hexdigest()


if __name__ == '__main__':
    unittest.main()
