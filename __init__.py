import random
import IO
from Bio import AlignIO, SeqIO, Alphabet
from Bio.Nexus import Nexus

"""
This module is used to handle DNA sequences assoiciated with individual
samples in phylogeography and multi-locus phylogeny projects. The main aim is
to make it easy to manage and modify data, to subselect sequences based on the
species or the site they belong to and to create the specfic datafiles 
required for analysis.

This docstring provides an overview of the module, you should read teh docs of e
ach class and each function to get the full detials

SnailBase doesn't use persitant objects, instead it accepts a flat fasta DNA
sequence file with the description field carrying information on species, gene 
and site: 

|>Hsap1 human COI Dunedin
|ATAGCAGTAATGCT
|>Hsap2 human COI Mosgeil
|ATAFHCATAGGTAC

Data is read into memory with the function ``SnailBase.IO.read()``

>>> import SnailBase as sb
>>> spider_data = sb.IO.read("tests/spiders.fasta", "fasta")
>>> spider_data
< Dataset with 40 specimens >

The ``Dataset`` object acts like a list, with each element of the list being a 
``Specimen`` object.

>>> len(spider_data)
40
>>> spider_data[1]
>>> <Specimen 'lh5' with 2 sequences>

Each specimen contains the specimen id, species and site information as well 
as a dictionary mapping sequences to the gene name. 

>>> spider1 = spider_data[1]
>>> dir(spider1)[4:]
['add_seq', 'id', 'sequences', 'site', 'species']
>>> spider1.site
'unknown'
>>> spider1.sequences.keys()
["ITS", "COI"]

Because the ``Dataset`` object inherits from a Python list, you can use the the 
information in each ``Specimen`` to select from the list. Say you wanted to find
those specimens with less that two genes associated

>>> [s for s in spider_data if len(s.sequences.keys() < 2)]
<Dataset with 1 specimens>

Instead or writting a list expression every time you want to subselect a 
dataset you can use the SnailBase function ``select()``. Think of the  
``select()`` function as some pre-built list comprehensions:

>>> halseti_data = sb.select(spider_data, "species", "halseti")
>>> #won't work with the test data, whcih doesn't have sites
>>> region1 = sb.select(spider_data, "site", ["site1", "site2"])

Once you've cut the dataset down, you'll want to write it out. You can single
gene alignments for any gene in your dataset to any of the formats Biopython 
can write 

>>> spiders = d 
>>> sb.IO.write_alignment(d, "COI", "spiderCOI.nex", "nexus")

There are some more sophisticated writers too. If you want to do some 
phylogeography or population genetics, you can use ``write_alignment()`` to 
generate an Arlequin (.arp) file with either species or sites used to split
specimens into samples

>>>sb.IO.write_alignment(d, "COI", "spiders.arp", "arp", sample="species")
wrote records in 3 samples"

For species tree estimation you can make files for BEST and BEAST and you can 
write an imap file for GSI analysis using ``write_multispecies()``
>> sb.IO.write_multispecies(d, "spiders", "BEAST")
wrote spidersCOI.nex
wrote spidersITS.nex
>> sb.IO.write_multispecies(d, "spiders", "BEST")
Add the following to the MrBayes block in spiders.nex
begin MyBayes;
taxset atritus = 1 2 9 10 11 12 13 14 15 16 17 18
taxset katipo = 3 4 5 6 7 8 19 20 21 22 23 24 25 26 27 28 29 30 31 32
taxset hasseltii = 33 34 35 36 37 38 39
>>> sb.IO.write_multispecies(d, "spiders", "GSI")
wrote GSI imap file for 39 taxa

"""



class Region:
    """ Class that represents a region, a set of sites
    
    work in progress, for now this is just a list"""
    def __init__(self, L=[]):
        list.__init__(self, L)
    

class Specimen:
    """ A specimen in a dataset with sequences and site information """
    def __init__(self, seq):
        self.sequences = {}
        self.id, self.species, gene, self.site = seq.description.split()
        self.sequences[gene] = seq
    
    def __repr__(self):
        return "<Specimen '%s' with %i sequences>" % \
               (self.id, len(self.sequences.keys()))
    
    def add_seq(self, record):
        """ add a SeqRecord object to a Specimen """
        record.seq.alphabet = Alphabet.Gapped(Alphabet.IUPAC.ambiguous_dna, "-")
        id, sp, gene, site = record.description.split()
        #TODO
        #these should have better error messages
        assert sp == self.species, "conflict about species for" + id
        assert site == self.site, "conflict about site for" + id
        self.sequences[gene] = record

class Dataset(list):
    """ A whole dataset with multiple species and tricks to subselect """
    def __init__(self, L=[]):
        list.__init__(self, L)
    
    def __repr__(self):
        return "< Dataset with %i specimens >" % len(self)

    def _new_seq_ids(self, gene, new_ids):
        """ Returns sequences with changed id (used by various writers) """
        L = []
        for i, seq in enumerate(self.get_sequences(gene)):
            yield SeqRecord(new_ids[i], seq.seq)
        
    def add_seqs(self, fname, format):
        """ Add a records from a sequence file to Dataset"""
        for record in SeqIO.parse(fname, format):
            if record.id in [spec.id for spec in self]:
                #specimen already exists, add new sequence
                spec_index = [spec.id for spec in self].index(record.id)
                self[spec_index].add_seq(record)
            else:
                #start a new specimen with the record
                record.seq.alphabet = Alphabet.IUPAC.ambiguous_dna
                self.append(Specimen(record))
    
    def get_species(self):
        """ return a list of species names for specimens """
        return [s.species for s in self]
    
    def get_ids(self):
        """ return a list of specimen ids """
        return [s.id for s in self]
    
    def get_sites(self):
        """ return a list of sites for specimens in dataset """
        return [s.site for s in self]
    
    def get_sequences(self, seq_name):
        """ get sequences corresonding to a sequence name """
        return [s.sequences[seq_name] for s in self \
                if seq_name in s.sequences]

    def get_genes(self):
        """ get a list of names of loci in this dataset """
        genes = []
        for s in self:
         for g in s.sequences.keys():
          if g not in genes:
            genes.append(g)
        return genes

    def change_species(self, from_species, to_species):
        """ changen a give species name to something else"""
        for i, sp in enumerate(self.species()):
            if sp == from_species:
                self[i].species = to_species
    
    def randomize_species(self):
        """Randomly assign existing species names to specimens 
        
        You might want to do this for permutation tests, returns a new dataset
        """
        d = Dataset(self[:])
        sp = d.get_species()
        random.shuffle(sp)
        for name, specimen in zip(sp,d):
            specimen.species = name
        return d
            
    def sample_by_taxon(self, taxon_tuple):
        """ Randomly select samples from different taxa (returns new dataset)
        
        Dataset.random_tax([("venosa", 2), ("globosa", "2")] returns Dataset
        with 2 venosa and 2 globosa specimens
        """
        out = []
        for species, n in taxon_tuple:
            sp_pool = [s for s in self if s.species == species]
            #can't append because sample() returns a list
            out = out + random.sample(sp_pool, n)
        return Dataset(out)



def select(dataset, attr, values, match_all=True):
    """A tool to subselect datasets based on arrtibutes of specimens """ 
 
    if attr == "ngenes":
        return [d for d in dataset if d.sequences.keys > values]
    else:      
        #type checking!!! (but I want to be able to pass string or list here)
        if isinstance(values, basestring):
            values = [values]
        if match_all:
            return Dataset([d for d in dataset if getattr(d, attr) in values])
        else:
            L = []
            for v in values:
                for d in dataset:
                    if specimen.attr_map[attr] == v and d not in L:
                        L.append(d)
            return Dataset(L)
