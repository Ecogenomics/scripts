#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__prog_name__ = 'trimSeqs.py'
__prog_desc__ = 'trim ends of sequences in a multiple sequence alignment'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse

class Trim(object):
    def __init__(self):
      pass

    def readFasta(self, fasta_file):
      """Read sequences from fasta file.

      Parameters
      ----------
      fasta_file : str
          Name of fasta file to read.

      Returns
      -------
      dict : dict[seq_id] -> seq
          Sequences indexed by sequence id.
      """

      if not os.path.exists(fasta_file):
          raise InputFileException

      if os.stat(fasta_file).st_size == 0:
          return

      open_file = open
      if fasta_file.endswith('.gz'):
          open_file = gzip.open

      seqs = {}
      for line in open_file(fasta_file):
          # skip blank lines
          if not line.strip():
              continue

          if line[0] == '>':
              seq_id = line[1:].split(None, 1)[0]
              seqs[seq_id] = []
          else:
              seqs[seq_id].append(line[0:-1])

      for seq_id, seq in seqs.iteritems():
          seqs[seq_id] = ''.join(seq)

      return seqs

    def run(self, inputMSA, outputMSA, bRemoveIdentical, minPerTaxa, minBP):
      # read seqs
      seqs = self.readFasta(inputMSA)

      # filter identical seqs
      identicalSeqs = set()
      if bRemoveIdentical:
        print 'Filtering identical sequences.'

        seqIds = seqs.keys()
        for i in xrange(0, len(seqIds)):
          seqIdI = seqIds[i]

          if seqIdI in identicalSeqs:
            continue

          for j in xrange(i+1, len(seqIds)):
            seqIdJ = seqIds[j]
            if seqs[seqIdI] == seqs[seqIdJ]:
              print '  Seq %s and %s are identical.' % (seqIdI, seqIdJ)
              identicalSeqs.add(seqIdJ)

        print '  Identified %d of %d sequeces as identical.' % (len(identicalSeqs), len(seqs))

      # trim start and end columns to consensus alignment
      firstChar = []
      lastChar = []
      for seqId, seq in seqs.iteritems():
        if seqId in identicalSeqs:
          continue

        for i, ch in enumerate(seq):
          if ch != '.' and ch != '-':
            firstChar.append(i)
            break

        for i in xrange(len(seq)-1, -1, -1):
          if seq[i] != '.' and seq[i] != '-':
            lastChar.append(i)
            break

      firstChar.sort()
      lastChar.sort(reverse=True)

      trimIndex = int((len(firstChar) * minPerTaxa) + 1)

      start = firstChar[trimIndex]
      end = lastChar[trimIndex]

      print ''
      print 'Trimming seqs from %d to %d leaving a %dbp length alignment.' % (start, end, end-start+1)

      shortSeqFile = outputMSA + '.short'
      fout = open(outputMSA, 'w')
      foutShort = open(shortSeqFile, 'w')
      numFilteredSeq = 0
      for seqId, seq in seqs.iteritems():
        if seqId in identicalSeqs:
          continue

        validBP = 0
        for i in xrange(start, min(len(seq), end+1)):
          ch = seq[i]
          if ch != '.' and ch != '-':
            validBP += 1

        if validBP >= minBP:
          fout.write('>' + seqId + '\n')
          fout.write(seq[start:end+1] + '\n')
        else:
          print 'Filtering seq %s with %d of %d (%.1f%%) aligned bases.' % (seqId, validBP, (end-start+1), validBP*100.0/(end-start+1))
          numFilteredSeq += 1
          foutShort.write('>' + seqId + '\n')
          foutShort.write(seq[start:end+1] + '\n')
      fout.close()
      foutShort.close()

      print 'Filtered %d of %d sequences due to length.' % (numFilteredSeq, len(seqs) - len(identicalSeqs))
      print 'Short sequence written to: %s' % shortSeqFile

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_msa', help='input msa')
    parser.add_argument('output_msa', help='output msa')
    parser.add_argument('-r', '--remove_identical', help='remove identical sequences from file', action='store_true')
    parser.add_argument('-p', '--min_per_taxa', help='minimum percentage of taxa required to keep leading and trailing columns', type=float, default=0.7)
    parser.add_argument('-b', '--min_bp', help='minimum base pairs required to keep trimmed seqs', type=int, default=1000)

    args = parser.parse_args()

    try:
        trim = Trim()
        trim.run(args.input_msa, args.output_msa, args.remove_identical, args.min_per_taxa, args.min_bp)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
