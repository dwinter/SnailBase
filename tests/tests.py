import unittest

import SnailBase as sb
from Bio import SeqIO

class TestSnailBase(unittest.TestCase)

    def setup(self):
    	d = sb.dataset_from_seqs(SeqIO.parse(open("spiders.fasta"), "fasta"))
    
    def test_read(self):
    """ test on reading the spider data """
    self.assertEqual(len(d), 40)
    self.assertEqual(d.species[7:11], 
                     ['katipo', 'hasseltii', 'hasseltii', 'atritus'])
    self.assertEqual(len(d.get_sequence("ITS")), 39)
    self.assertEqual(len(d.get_sequence("COI")), 40)
    
    def test_select(self):
    """ do the selection methods work """
    self.assertEqual(len(sb.select(d, "gene", "COI")), 40)
    self.assertEqual(len(sb.select(d, "species", "katipo")), 20)
    self.assertEqual(len(sb.select(d, "ngenes", 1)), 1)
    self.assertEqual(len(sb.seletct(d, "id", "La13")),1)
    
    def test_write(self):
    """ do all the writers work ? """
    self.assertEqual(sb.write_species(d, "BEST"), best_block)
    self.assertEqual(sb.write_species(d, "BEAST"), beast_block)
    self.assertEqual(sb.write_species(d, "gsi"), gsi_table)
    # other formats will need files


