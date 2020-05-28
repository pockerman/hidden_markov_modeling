import pysam
import numpy as np
from scipy import stats
import re
import time
import sys

def total_memory(windows):

  total = 0
  for window in windows:
    total += sys.getsizeof(window.keys())
    total += sys.getsizeof(window.values())
    total += sys.getsizeof(window)

  return total

def windowAna(chr,start,end,qual,fas,sam):
    refseq = '' #store for the refseq
    samseq = [] #stote for the sample
    pos = [] #reference pos
    nseq = [] #all reads
    nalign = [] #quality limited reads
    head = 0 #whether read starts in window
    head1 = 0 #start of read 1 in pair
    head2 = 0 #start of read 2 in pair
    read1 = 0 #equivalent nalign just for read 1s in pairs
    read2 = 0 #as above for read 2s
    errorAlert = False

    try:

      time_start = time.perf_counter()
      for pcol in sam.pileup(chr,start,end,
                             truncate=True,ignore_orphans=False,
                             max_depth=1000):
          pcol.set_min_base_quality(qual)

          #fill in start when there are no reads present
          while pcol.reference_pos != start:
              samseq.append('_')
              refseq += fas.fetch(chr,start,start+1)
              pos.append(start)
              nseq.append(0)
              nalign.append(0)
              start+=1

          #collect metrics for pileup column
          pos.append(start)
          nseq.append(pcol.nsegments) #all reads over a pileup column
          nalign.append(pcol.get_num_aligned()) #above but quality limited
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
                  refseq+= fas.fetch(chr,pcol.reference_pos,pcol.reference_pos+1)
              else:
                  x = pcol.get_query_sequences(add_indels=True)
                  x = [a.upper() for a in x]
                  samseq.append(set(x)) #store as unique set
                  refseq+= fas.fetch(chr,pcol.reference_pos,pcol.reference_pos+1)
          except Exception as e: # may fail if large number of reads
              try:
                  x = pcol.get_query_sequences(add_indels=True)
                  x = [a.upper() for a in x]
                  samseq.append(set(x)) #store as unique set
                  refseq+= fas.fetch(chr,pcol.reference_pos,pcol.reference_pos+1)
              except Exception as e:
                  errorAlert = True
          start+=1


      time_end = time.perf_counter()
      print("Time for loop over pcol {0}".format(time_end - time_start))
      sys.stdout.flush()

      time_start = time.perf_counter()

      #fill in end if no reads at end of window
      while start < end:
              samseq.append('_')
              refseq += fas.fetch(chr,start,start+1)
              pos.append(start)
              nseq.append(0)
              nalign.append(0)
              start+=1

      #Metrics for read depth
      allmean = np.mean(nseq) #metrics, may not need all these
      qmean = np.mean(nalign)
      allmedian = np.median(nseq)
      qmedian = np.median(nalign)
      allsum = np.sum(nseq)
      qsum = np.sum(nalign)

      #GC content
      gcr = (refseq.count('G') + refseq.count('g') + refseq.count('C') + refseq.count('c'))/len(refseq)
      gapAlert = True if 'N' in refseq or 'n' in refseq else False


      gcmax = 0
      gcmaxlen = 0
      gcmin = 0
      gcminlen = 0
      for bs in samseq:
          minelgc = None
          lenminelgc = None
          maxelgc = None
          lenmaxelgc = None
          for el in bs:
              el = el.split('-')
              ellen = len(re.sub('[\+\-_Nn*]','',el[0]))
              elgc = len(re.sub('[\+\-_Nn*AaTt]','',el[0]))
              if minelgc == None or elgc < minelgc:
                  minelgc = elgc
                  lenminelgc = ellen
              if maxelgc == None or elgc > maxelgc:
                  maxelgc = elgc
                  lenmaxelgc = ellen
          gcmax += maxelgc
          gcmaxlen += lenmaxelgc
          gcmin += minelgc
          gcminlen += lenminelgc

      gcmax = gcmax/gcmaxlen if gcmaxlen > 0 else None
      gcmin = gcmin/gcminlen if gcminlen > 0 else None

      time_end = time.perf_counter()
      print("Time for remaining function {0}".format(time_end - time_start))
      sys.stdout.flush()

      output={'gcmax':gcmax,
              'gcmin':gcmin,
              'gcr':gcr,
              'gapAlert':gapAlert,
              'allmean':allmean,
              'qmean':qmean,
              'allsum':allsum,
              'qsum':qsum,
              'errorAlert':errorAlert,
              'head':head,
              'start':pos[0],
              'end':pos[-1],
              'minLen':gcminlen
              }
      return output
    except MemoryError as e:
      total = sys.getsizeof(refseq)
      total += sys.getsizeof(samseq)
      total += sys.getsizeof(pos)
      total += sys.getsizeof(nseq)
      total += sys.getsizeof(nalign)

      print("Memory used by function {0} GB".format(total*1e-9))
      sys.stdout.flush()
      raise

if __name__ == '__main__':
  fas = pysam.FastaFile("/scratch/spectre/a/ag568/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
  sam = pysam.AlignmentFile("/scratch/spectre/a/ag568/m585_verysensitive_trim_sorted.bam", "rb")


  c = 'chr1'
  start = 1000000
  end = 1500000
  qual = 20

  counter = 0
  windows = []
  time_start = time.perf_counter()
  windows = []
  try:
    while start < end:
        met = windowAna(c,start,start+100,qual,fas,sam)
        windows.append(met)
        print("Created window: ", counter)
        print("Window pos: {0}/{1} ".format(met['start'], met['end']))

        counter += 1
        start = start + 100
    time_end = time.perf_counter()
    print("Time to create  {0} windows is {1} secs".format(counter, time_end - time_start))
    sys.stdout.flush()
  except MemoryError as e:
    time_end = time.perf_counter()
    print("Time to create  {0} windows is {1} secs (Exception case)".format(counter, time_end - time_start))
    sys.stdout.flush()
    total = total_memory(windows)
    print("MemoryError exception detected. "
          "Windows memory used is: {0} GB".format(total*1e-9))
    sys.stdout.flush()
    raise


