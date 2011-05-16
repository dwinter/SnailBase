def _one_to_many_dict(keys, values):
    """ creates a dictionary in which keys point to a list of values """
    d = dict()
    for key, value in zip(keys, values):
		d.setdefault(k, []).append(value)
	return d
	


class AlignmentWriter(GenericWriter):
	
	def __init__(self, dataset, gene, format):

	
	
class ArpWriter(GenricWriter):
	"""Writes an alignment to Arlequin file format """
	
	def __init__(self, sequences):
		GenericWriter.__init__(self, sequences)
		
		self.profile_text = """
		[Profile]

			Title="Automatically generated with biopython"
			NbSamples= %i
		    GenotypicData=0
		    DataType=DNA
		    LocusSeparator=NONE
		    MissingData="-"

		[Data]

		    [[Samples]]
		"""
		
		self.sample_text = """
		    SampleName="%s"
		    SampleSize= %i
		    SampleData={\n\n
		"""
		
    def write(self, filename, sample_map):
    """ writes sequences to arlequin file """

        if len(records) == 0:
            raise ValueError("Need at at least one sequence to write")

        samples = one_to_many_dict(sample_map, self.sequences)    
        hanle = open(filename, "w")
        handle.write(self.profile_text % len(samples.keys()))
        samples = _one_to_many_dict(sample_map, records)
        for sample, records in samples.items():
            handle.write(self.sample_text % (sample, len(records)))
            for rec in records:
                handle.write("\t\t%s 1 %s\n" % (rec.name, rec.seq))
            handle.write("}\n")


class BPPWriter:
	""" """
	def __init__(self, dataset, filestem):
		seqfile = open(filestem + '.phy', 'a')
		for gene in dataset.genes():
			seqfile.write(" %s %s\n" % (len(dataset), len(dataset[1].sequences[gene])))
			for d in dataset:
				seqfile.write("a^%s  %s\n" % (d.id[0:10], str(d.sequences[gene].seq)))
		seqfile.close()
		
		mapfile = open(filestem + '.imap', 'a')
		for d in dataset:
			mapfile.write("%s %s\n" % (d.id[0:10], d.species))
		mapfile.close()
		
		print "species counts for control file:"
		species = dataset.get_species()
		for sp in set(species):
			print sp, species.count(sp) 
	
	

class BEAST:
	"""Write files needed for analysis with BEAST
	
	I tried to make this write a tab delimted traits file for the species 
	but this doesn't seem to work. Instead we add a the species name to the
	identifiersof a nexus file after a '|' character. You can then use BEAUti
	to 'guess' the species adn set up a *BEAST run"""
	
	def _annotated(seqs):
		for i, seq in enumerate(dataset.get_sequences(gene)):
			seq.id = "%S|%S" % (seq.id, dataset[i].species)
			yield seq
	
	def write_file(dataset):
		for gene in dataset.sequences.keys():
			


class Gsi(GenericWriter):
	"""write a species map for gsi
	
	The Genetic Sorting Indec (gsi) is calculated from gene trees, to calculate
	the statistic using R you need a mapping file that relates the tip-labels of 
	the trees to species. That file is written by this writer
	"""
	
	def write_file():
     for id, sp in [(dataset.id, dataset.species) for dataset in self]:
             self.handle.write( '"%s" "%s"\n' % (id, sp))