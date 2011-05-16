from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment 
from Interfaces import AlignmentIterator, SequentialAlignmentWriter 


profile_text = """
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

sample_text = """
    SampleName="%s"
    SampleSize= %i
    SampleData={\n\n
"""



class ArlequinIO()
    """ Write an alignment in Arlequin format (work in progress...)
        
        Works like AlignIO in Biopython with the adition of "sample map" to map
        records to a sample in an arp file.
        
        >>> sequences = SeqIO.parse("primates.fasta", "fasta")
        >>> sp = ["H", "H", "G", "G", "O", "O"]
        >>> ArlequinIO.write(sequences, sp, open("primates.arp", "arp")
    """

	def _one_to_many_dict(keys, values):
        """ creates a dictionary in which keys point to a list of values """
        d = dict()
        for k,v in zip(keys, values):
            if k not in d.keys():
                d[k] = [v]
            else:
                d[k].append(v)
        return d
    
    def write(self, records, sample_map, handle):
    """ writes sequences to arlequin file """

        if len(records) == 0:
            raise ValueError("Need at at least one sequence to write")
        samples = one_to_many_dict(sample_map, records)            
        handle.write(profile_text % len(samples.keys()))
        samples = self._one_to_many_dict(sample_map, records)
        for sample, records in samples.items():
            handle.write(sample_text % (sample, len(records)))
            for rec in records:
                handle.write("\t\t%s 1 %s\n" % (rec.name, rec.seq))
            handle.write("}\n")
       
 
	
############################################################
##This works in Biopython, needs refactoring for SnailBase##
############################################################
	
class ArlequinIterator(AlignmentIterator):
    """Iterate over alignment of DNA sequences in Arlequin format
    """
    def _clean(self, line):
        """Remove comments from a line
        """
        return line.split("#")[0]

    def _get_attr(self, line):
        """get attributes from the header sections of file
        """
        value = line.split("=")[1].replace('"', "").replace("'", "")
        return value.strip()
        
    def _is_profile(self):
        line = self.handle.readline()
        while "[Profile]" not in self._clean(line):
            if not line:
                return 
            line = self.handle.readline()
        return True
        
    def _get_profile(self):
        line = self.handle.readline()
        while "[Data]" not in self._clean(line):
            line = self._clean(self.handle.readline())
            if "DataType" in line:
                dtype = self._get_attr(line)
                if dtype != 'DNA':
                    raise ValueError("At the moment the arlequin parser only "
                    "supports haplotypic datafiles, this file has %s data" %
                    dtype)
            if "MissingData" in line:
                self.alphabet = Gapped(IUPAC.ambiguous_dna,
                                       self._get_attr(line))
            if "NbSamples" in line:
                try:
                    self.nsamples = int(self._get_attr(line))
                except ValueError:
                   print "Can't coerce NbSamples from profile to an integer"

    def next(self):
        #is there a new profile to read?
        if not self._is_profile():
            raise StopIteration
        self._get_profile()
        alignment = Alignment(self.alphabet)
        for sample in xrange(0, self.nsamples):
            line = self._clean(self.handle.readline())
            while "SampleName" not in line:
                line = self._clean(self.handle.readline())
            sample_name = self._get_attr(line)
            line = self._clean(self.handle.readline())
            nseqs = int(self._get_attr(line))
            line = self.handle.readline()
            records = []
            for sequence in xrange(0, nseqs):
                line = self._clean(self.handle.readline())
                records.append(line.split())
            for id, freq, seq in records:
                alignment.add_sequence(id, seq)
                record = alignment.get_all_seqs()[-1]
                record.annotations["sample"] = sample_name
                record.annotations["frequency"] = int(freq)
        return alignment