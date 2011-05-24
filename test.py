import unittest, sys, filecmp, os
import SnailBase as sb

#setup
os.chdir("tests/")
d = sb.IO.read("spiders.fasta", "fasta")

d2 = sb.IO.read("spiders.fasta", "fasta")
d2.add_seqs("spidersCOI.fasta", "fasta")
      
files_to_delete = ["spiders.imap","spidersBEST.nex", "spiders.arp","spiders_its.fasta"]
      
class TestSnailBase(unittest.TestCase):

  def test_read(self):
    """ Can SnailBase read data and represent in correctly """
    self.assertEqual(len(d), 39)
    self.assertEqual(d.get_species()[::11], 
                    ['atritus', 'atritus', 'katipo', 'hasseltii'])
    self.assertEqual(len(d.get_sequences("ITS")), 39)

  def test_add(self):
    """ Does adding a sequence add specimens/sequences to dataset """
    self.assertEqual(len(d2), 40)
    self.assertEqual(len(d2.get_sequences("COI")), 39)

  def test_select(self):
    """ do the selection methods work """
    self.assertEqual(len(sb.select(d2, "species", "katipo")), 20)
    self.assertEqual(len(sb.select(d2, "ngenes", 2)), 38)
    self.assertEqual(len(sb.select(d2, "id", "La13")),1)


class TestWriters(unittest.TestCase):
  def tearDown(self):
    """ delete files only neede for testing (must be a pretteir way...) """
    for f in files_to_delete:
      if os.path.isfile(f):
        os.remove(f)
  
  def test_alignment(self):
    """ Can we write an alignment for one gene """
    sb.IO.write_alignment(d2, "ITS", "spiders_its.fasta", "fasta")
    self.assertTrue(filecmp.cmp("spiders_its.fasta", "test_ITS.fasta"))
  
  def test_arp(self):
    """ Can we write valid arlequin files """
    sb.IO.write_alignment(d2, "COI", "spiders.arp", "arp", sample="species")
    self.assertTrue(filecmp.cmp("spiders.arp", "test_arp.arp"))
  
  def test_multi(self):
    """ Can we write formats for multi-species stuff """
    sb.IO.write_multispecies(d2, "spidersBEST", "BEST")
    self.assertTrue(filecmp.cmp("spidersBEST.nex", "test_BEST.nex"))
    sb.IO.write_multispecies(d2, "spiders", "GSI")
    self.assertTrue(filecmp.cmp("spiders.imap", "test_GSI.imap"))

if __name__ == "__main__":
    unittest.main()   
