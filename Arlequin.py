from collections import defaultdict
from Bio.Alphabet import IUPAC, Gapped
from Bio.AlignIO.Interfaces import AlignmentIterator


profile_text = """
[Profile]

	  Title="Automatically generated with Biopython"
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

def write(records, filename, sample_map):
  """ writes sequences to arlequin file """
  handle = open(filename, "w")
  if len(records) == 0:
      raise ValueError("Need at at least one sequence to write")
  samples = defaultdict(list)
  for r, s in zip(records,sample_map):
    samples[s].append(r)
  handle.write(profile_text % len(samples.keys()))
  for sample, records in samples.items():
      handle.write(sample_text % (sample, len(records)))
      for rec in records:
          handle.write("\t\t%s 1 %s\n" % (rec.name, rec.seq))
      handle.write("}\n")
  print "wrote records in %s samples" % samples.keys()
 
	
#####################################################################
##This works in Biopython, probably needs refactoring for SnailBase##
#####################################################################
	
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
        #return because there is only one alignment per file        
        return alignment
