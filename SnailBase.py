import random
from collections import defaultdict

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

SnailBase doesn't use persitant objects, instead it accepts a flat "fasta" DNA
sequence file with the description field carrying information on species, gene 
and site: 

|>Hsap1 human COI Dunedin
|ATAGCAGTAATGCT
|>Hsap2 human COI Mosgeil
|ATAFHCATAGGTAC

Data is read into memory with the function dataset_from_seq() 

>>> import SnailBase as sb
>>> spider_data = sb.dataset_from_seqs("tests/spiders.fasta", "fasta")
>>> spider_data
< Dataset with 40 specimens >

The Dataset object acts like a list, with each element of the list being a 
Specimen object.

>>> len(spider_data)
40
>>> spider_data[1]
>>> <Specimen xxx with 2 sequences>
>>> spider1 = spider_data[1]

Each specimen contains the specimen id, species and site information as well 
as a dictionary mapping sequences to the gene name. 

>>> spider1 = spider_data[1]
>>> spider1.site


Because the Dataset object
inherits from a Python list, you can use the the information in each Specimen 
to select from the list



>>> [s for s in spider_data if len(s.sequences.keys() < 2)]
<Dataset with 1 specimens>

Instead or writting a list expression every time you want to subselect a 
dataset you can use the SnailBase function select. Think of hte  select() 
function as some pre-built list comprehensions:

>>> halseti_data = sb.select(spider_data, "species", "halseti")
>>> region1 = sb.select(spider_data, "site", ["site1, "site2"])

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
        
    def add_seqs(self, seqs):
        """ Add a records from a sequence file to Dataset"""

        for record in seqs:
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
      genes = []
      for d in self:
        for g in d.sequences.keys():
          if g not in genes:
            genes.append(g)
      return genes

    def change_species(self, from_species, to_species):
        """ changen a give species name to something else"""
        for i, sp in enumerate(self.species()):
            if sp == from_species:
                self[i].species = to_species
    
    def randomize_species(self):
        """Randomly assign existing species names to specimens """
        sp = random.shuffle(self.get_species())
        for name in sp:
            self.species = name
            
    def sample_by_taxon(self, taxon_tuple):
        """ Randomly select samples from different taxa
        
        Dataset.random_tax([("venosa", 2), ("globosa", "2")] returns Dataset
        with 2 venosa and 2 globosa specimens
        """
        out = []
        for species, n in taxon_tuple:
            sp_pool = self.select_species(species)
            #can't append because sample() returns a list
            out = out + random.sample(sp_pool, n)
        return Dataset(out)


##_multi_map = {"bpp":Writers.bpp, "BEST": Writers.BEST, "BEAST": Writers.BEAST, 
##              "gsi": Writers.gsi, "nexus": Writers.nexus}
##

def select(dataset, attr, values, match_all=True):
    """A tool to subselect datasets based on arrtibutes of specimens """ 
 
    if attr == "ngenes":
        return [d for d in dataset if d.sequences.keys > values]
    else:      
        attr_map = {"species": d.species, "id": d.id, "site": d.site, 
                    "gene": d.sequences.keys(), }
        if type(values) == string:
            values = [values]
        if match_all:
            return Dataset([d for d in dataset if attr_map[attr] in values])
        else:
            L = []
            for v in values:
                for d in dataset:
                    if specimen.attr_map[attr] == v and d not in L:
                        L.append(d)
            return Dataset(L)


def write_alignment(dataset, gene, filename, format, **kwarks):
    """ 
    Writes a single alignment 
    
    For most formats, this wraps Bio.SeqIO.write(). The user provides a gene
    name, filename and file format. If format = "arp" for Arlequin then
    the key-word argument 'sample' can be used to split the sequences into 
    samples based on the 'sites' or 'species' attributes of the specimens
    """
    
    if  format == "arp":
        if sample == "sites":
            Writers.Arlequin(dataset, gene).write(
                    filename, sample_map = dataset.get_sites())
        elif sample == "species":
            Writers.Arlequin(dataset, gene).write( 
                             filename, dataset.get_species())
        else: #write all the sequences to one big sample
            Writers.Arlequin(dataset.get_sequences(gene), 
                             filename, ["OneBigOne"] * len(dataset))
    else:
        SeqIO.write(dataset.get_sequences(gene), handle, format)
                
def write_multilocus(dataset, filename, format):
    """ Writes a multilocus alignment for species tree estimation """
    #writer_class = _sptree_map[format]
    #writer_cass(filename).write_files(datasets)
    raise NotImplementedError
    
########
def write_beast(dataset, file_stem):
  """ writes a nexus file for each gene, adding species name to id 
        
  adds *_speciesID to each sequence, so BEAUTi can 'guess' the species
  in the 'traits' tab """
  genes = dataset.get_genes()
  for g in genes:
    seqs = [(d.species, d.sequences[g]) for d in dataset if g in d.sequences]
    for (sp,s) in seqs:
     s.id = "%s_%s" % (s.id,sp)
    fname = "%s_%s.nex" % (file_stem, g)
    SeqIO.write([seq for sp,seq in seqs], open(fname, "w"), "nexus")
        

def write_best(self, file_handle=None):
        """ write a MrBayes block for BEST species tree estimation """
        d = defaultdict(list)
        for sp, i in zip(self.species(),
                         [str(i) for i in xrange(1,len(self)+1)]):
          d[sp].append(i)
        contents = ["begin MyBayes;"]
        for species, OTUs in d.items():
            contents.append("taxset %s = % s" % (species, " ".join(OTUs)))
        if file_handle:
            file_handle.write("\n".join(contents))
        else:
            for line in contents:
                print line
    
    
def write_gsi(self, file_handle):
    """ Write a mapping file for genealogical sorting w/ R """    
    for id, sp in [(t.id, t.species) for t in self]:
        file_handle.write( '"%s" "%s"\n' % (id, sp))

def dataset_from_seq(sequences):
    """ take a SeqIO object (or list of sequences) and make a dataset """
    d = Dataset()
    d.add_seqs(sequences)
    return d

def nexify(sequences):
    """ set up a nexus file for one gene """
    n = Nexus.Nexus("""#NEXUS\nbegin data; dimensions ntax=0 nchar=0;
                   format datatype=DNA; end;""")
    n.alphabet = sequences[1].seq.alphabet
    for record in sequences:
        n.add_sequence(record.id, record.seq.tostring())
    return n



spiders = SeqIO.parse("tests/spiders.fasta", "fasta")
d = dataset_from_seq(spiders)
print d.get_ids()
