import SnailBase as sb
from collections import defaultdict
from Bio import SeqIO
from Bio.Nexus import Nexus


def read(sequences, format):
  """ take a SeqIO object (or list of sequences) and make a dataset """
  d = sb.Dataset()
  d.add_seqs(sequences, format)
  return d
    

def write_alignment(dataset, gene, filename, format, **kwarks):
  """ 
  Writes an alignment for one gene
  
  For most formats, this wraps Bio.SeqIO.write(). The user provides a gene
  name, filename and file format. If format = "arp" for Arlequin then
  the key-word argument 'sample' can be used to split the sequences into 
  samples based on the 'sites' or 'species' attributes of the specimens
  """
  seqs = dataset.get_sequences(gene)  
  if format == "arp":
    if sample == "sites":
      sb.Arlequin.write(seqs, filename, dataset.get_sites())
    elif sample == "species":
      sb.Arlequin.write(seqs, filename, dataset.get_species())
    else: #write all the sequences to one big sample
      sb.Arlequin.write(seqs, filename,["OneBigOne"] * len(dataset)))
  else:
    SeqIO.write(dataset.get_sequences(gene), filename, format)

def write_dataset(dataset, filename, format):
  """
  Write out the whole data set
  
  If
  """
  if format in ["snailbase", "sb", "fasta"]:
    seqs = (s.sequences.values() for s in dataset)
    SeqIO.write(seqs, open(filename, "w"), "fasta")
  elif format == "nexus":
    nexi = []
    for g in dataset.get_genes():
      nexi.append( (g, _nexify( dataset.get_sequences(g)))) 
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=filename)
  else:
    raise ValueError, "Can't handle format %s for whole dataset" % format

def write_multispecies(dataset, filestem, format):
  """
  """
  multi_map = {"BEAST":_write_BEAST, "BEST":_write_BEST, "GSI":_write_GSI}
  if format in multi_map:
    multi_map[format](dataset, filestem)
  else: 
    raise ValueError, "Can't handle format '%s'" % format
    

def _write_GSI(dataset, filestem):
    """ Write a mapping file for genealogical sorting w/ R """    
    handle = open(filestem, 'w')
    counter = 0
    for id, sp in [(t.id, t.species) for t in dataset]:
        handle.write( '"%s" "%s"\n' % (id, sp))
        counter += 1
    print "wrote GSI imap file for %s taxa" % counter 


def _write_BEAST(dataset, filestem):
  """ Writes a nexus file for each gene, adding species name to id 
        
  adds *_speciesID to each sequence, so BEAUTi can 'guess' the species
  in the 'traits' tab 
  """
  genes = dataset.get_genes()
  for g in genes:
    seqs = [(d.species, d.sequences[g]) for d in dataset if g in d.sequences]
    for (sp,s) in seqs:
     s.id = "%s_%s" % (s.id,sp)
    fname = "%s_%s.nex" % (filestem, g)
    SeqIO.write([seq for sp,seq in seqs], open(fname, "w"), "nexus")
    print "wrote %" % fname
        

def _write_BEST(dataset, filestem):
  """ write a MrBayes block for BEST species tree estimation """
  fname = filestem + ".nex"
  #write a nexus file with partitions for each gene
  nexi = []
  for g in dataset.get_genes():
    nexi.append( (g, _nexify( dataset.get_sequences(g)))) 
  combined = Nexus.combine(nexi)
  combined.write_nexus_data(filename=fname)
  #then build a MrBayes blog for BEST
  d = defaultdict(list)
  for sp, i in zip(dataset.get_species(),
                   [str(i) for i in xrange(1,len(dataset)+1)]):
    d[sp].append(i)
  contents = ["begin MyBayes;"]
  for species, OTUs in d.items():
      contents.append("taxset %s = % s" % (species, " ".join(OTUs)))
  print "Add the following to the MrBayes block in %s" % fname
  for line in contents:
      print line

def _nexify(sequences):
    """ set up a nexus file for one gene (used for concatenation) """
    n = Nexus.Nexus("""#NEXUS\nbegin data; dimensions ntax=0 nchar=0;
                   format datatype=DNA; end;""")
    n.alphabet = sequences[1].seq.alphabet
    for record in sequences:
        n.add_sequence(record.id, record.seq.tostring())
    return n
