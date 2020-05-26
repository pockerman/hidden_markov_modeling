"""
Utilities for working with bam files
"""
import math
import pysam
import logging
import numpy as np
import re
import array


from helpers import Window
from helpers import INFO, WARNING, DEBUG
from exceptions import Error


def extract_windows(chromosome, ref_filename, bam_filename, **args):

    """
    Extract the windows that couple the seq_file and ref_files
    for the given chromosome
    :param chromosome: chromosome name (str)
    :param ref_file: The reference file
    :param test_file: The sequence file
    :return:
    """

    windowcapacity = args["windowsize"]
    start_idx = args["start_idx"]
    end_idx = args["end_idx"]

    # the windows list
    windows = []

    with pysam.FastaFile(ref_filename) as fastafile:

        print("{0} Reference file: {1}".format(INFO, fastafile.filename))

        with pysam.AlignmentFile(bam_filename, "rb") as sam_file:
            print("{0} Sam file: {1} ".format(INFO, sam_file.filename))

            wcounter = 0
            while start_idx < end_idx:
              #print("{0} Window Start/End {1}/{2}.".format(INFO, start_idx, start_idx+windowcapacity))
              sam_output = window_sam_file(chromosome=chromosome,
                                       sam_file=sam_file,
                                       fastafile=fastafile,
                                       start=start_idx, end=start_idx+windowcapacity,
                                       **args)
              windows.append(Window(idx=wcounter,
                                    capacity=windowcapacity,
                                    samdata=sam_output))

              start_idx += windowcapacity

              #print("{0} Created window: {1}".format(INFO, wcounter))

              wcounter += 1
              #if wcounter == 1:
              #  break
    return windows


def window_sam_file(chromosome, sam_file, fastafile,
                    start, end, **kwargs):

  # store the refseq
  refseq = '' #array.array('u')

  # store for the sample
  samseq = array.array('u')
  pos = array.array('L') #[] #reference pos
  nseq = array.array('L') #[] #all reads
  nalign = array.array('L') #[] #quality limited reads
  head = 0 #whether read starts in window
  head1 = 0 #start of read 1 in pair
  head2 = 0 #start of read 2 in pair
  read1 = 0 #equivalent nalign just for read 1s in pairs
  read2 = 0 #as above for read 2s
  errorAlert = False

  for pcol in sam_file.pileup(chr=chromosome,
                         start=start, end=end,
                         truncate=kwargs["sam_read_config"]["truncate"],
                         ignore_orphans=kwargs["sam_read_config"]["ignore_orphans"],
                         fastafile=fastafile,
                         max_depth=kwargs["sam_read_config"]["max_depth"]):

    # if there is a quality threshold then use it
    if kwargs["sam_read_config"]["quality_threshold"] is not None:
      pcol.set_min_base_quality(kwargs["sam_read_config"]["quality_threshold"])

    #fill in start when there are no reads present
    while pcol.reference_pos != start:
            samseq.append('_')
            refseq += fastafile.fetch(chromosome,start,start+1)
            pos.append(start)
            nseq.append(0)
            nalign.append(0)
            start += 1

    #collect metrics for pileup column
    pos.append(start)

    #all reads over a pileup column
    nseq.append(pcol.nsegments)

    #above but quality limited
    nalign.append(pcol.get_num_aligned())
    for p in pcol.pileups:
            head += p.is_head #start of read present
            read1 += p.alignment.is_read1 #first of mate pair (from nalign)
            read2 += p.alignment.is_read2 #second of mate pair (from nalign)
            if p.is_head == 1 and p.alignment.is_read1 == 1: #above but is start of read
                head1 += 1
            if p.is_head == 1 and p.alignment.is_read2 == 1:
                head2 += 1

    #get bases from reads at pileup column (quality limited)
    try:
            if pcol.get_num_aligned() == 0:
                samseq.append('_')
                refseq += fastafile.fetch(chromosome, pcol.reference_pos,
                                              pcol.reference_pos+1)
            else:
                x = pcol.get_query_sequences(add_indels=kwargs["sam_read_config"]["add_indels"])
                x = [a.upper() for a in x]
                samseq.append(set(x)) #store as unique set
                refseq += fastafile.fetch(chromosome,pcol.reference_pos,pcol.reference_pos+1)
    except Exception as e: # may fail if large number of reads
            try:
                x = pcol.get_query_sequences(add_indels=kwargs["sam_read_config"]["add_indels"])
                x = [a.upper() for a in x]
                samseq.append(set(x)) #store as unique set
                refseq += fastafile.fetch(chromosome, pcol.reference_pos,
                                              pcol.reference_pos+1)
            except Exception as e:
                errorAlert = True
    start += 1

  #fill in end if no reads at end of window
  while start < end:
            samseq.append('_')
            refseq += fastafile.fetch(chromosome,start,start+1)
            pos.append(start)
            nseq.append(0)
            nalign.append(0)
            start+=1

  #Metrics for read depth
  #allmean = nseq
  rdmean = np.mean(nseq)
  rdstd  = np.std(nseq)
  rdmedian = np.median(nseq)
  rdsum = np.sum(nseq)

  #qseq = nalign
  qmean = np.mean(nalign)
  qstd = np.std(nalign)
  qmedian = np.median(nalign)
  qsum = np.sum(nalign)

    #GC content
  gcr = (refseq.count('G') + refseq.count('g') + refseq.count('C') + refseq.count('c'))/len(refseq)
  gapAlert = True if 'N' in refseq or 'n' in refseq else False

  altseq = array.array('u')
  for bs in samseq:
        talt = []
        for i in range(len(altseq)):
            for el in bs:
                alt = altseq[i] + el
                talt.append(alt)
        altseq = talt

  gcmax = 0
  gcmin = 1
  minLen = pos[-1] - pos[0] +1
  for alt in altseq:
        minLen = len(re.sub('[\+\-_Nn*]','',alt)) if len(re.sub('[\+\-_Nn*]','',alt)) < minLen else minLen

        try:
          gct = len(re.sub('[\+\-_Nn*AaTt]','',alt))/len(re.sub('[\+\-_Nn*]','',alt))
          gcmax = gct if gct > gcmax else gcmax
          gcmin = gct if gct < gcmin else gcmin
        except:
          pass

  output={'gcmax': gcmax,
            'gcmin': gcmin,
            'gcr': gcr,
            'gapAlert': gapAlert,
            #'allmean': nseq,
            'rdmean': rdmean,
            'rdstd': rdstd,
            'rdmedian': rdmedian,
            'rdsum':rdsum,
            #'qseq': nalign,
            'qmean':qmean,
            'qstd':qstd,
            'qmedian':qmedian,
            'qsum':qsum,
            'errorAlert':errorAlert,
            'head':head,
            'start':pos[0],
            'end': pos[-1],
            'minLen':minLen
            }

  return output

