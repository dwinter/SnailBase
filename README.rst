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

| >Hsap1 human COI Dunedin
| ATAGCAGTAATGCT
| >Hsap2 human COI Mosgeil
| ATAFHCATAGGTAC

Data is read into memory with the function ``SnailBase.IO.read()``

``
>>> import SnailBase as sb
>>> spider_data = sb.IO.read("tests/spiders.fasta", "fasta")
>>> spider_data
< Dataset with 40 specimens >
``

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
>>> sb.IO.write_multispecies(d, "spiders", "BEAST")
wrote spidersCOI.nex
wrote spidersITS.nex
>>> sb.IO.write_multispecies(d, "spiders", "BEST")
Add the following to the MrBayes block in spiders.nex
begin MyBayes;
taxset atritus = 1 2 9 10 11 12 13 14 15 16 17 18
taxset katipo = 3 4 5 6 7 8 19 20 21 22 23 24 25 26 27 28 29 30 31 32
taxset hasseltii = 33 34 35 36 37 38 39
>>> sb.IO.write_multispecies(d, "spiders", "GSI")
wrote GSI imap file for 39 taxa

