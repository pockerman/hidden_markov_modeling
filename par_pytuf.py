import pysam
import numpy as np
from scipy import stats
import re
import time
import sys
from multiprocessing import Queue
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import Manager

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
                             fastafile=fas, max_depth=1000):
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

      altseq = ['']
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
              'minLen':minLen
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


def process_worker(p, start, end, windows_dict):

    fas = pysam.FastaFile("/scratch/spectre/a/ag568/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
    sam = pysam.AlignmentFile("/scratch/spectre/a/ag568/m585_verysensitive_trim_sorted.bam", "rb")

    c = 'chr1'
    start = 1000000
    end = 2000000
    qual = 20


    windows = []
    try:

      time_start = time.perf_counter()
      while start < end:
         met = windowAna(c,start,start+100,qual,fas,sam)
         windows.append(met)
         start = start + 100
      time_end = time.perf_counter()
      windows_dict[p] = windows
      print("Process {0} finished. "
            "Total time for {1} windows {2} secs".format(p, len(windows), time_end - time_start))
      sys.stdout.flush()
    except MemoryError as e:

      total = total_memory(windows)
      print("MemoryError exception detected in process {0}. "
           "Windows memory used is: {1} GB".format(p, total*1e-9))
      sys.stdout.flush()
      return
    except Exception as e:
       print("An exception detected in process {0}. "
            "Exception is: {1}".format(p, str(e)))
       sys.stdout.flush()
       return

if __name__ == '__main__':


  n_procs = 5
  start = 1000000
  end = 2000000
  qual = 20

  manager = Manager()
  windows_dict = manager.dict()

  for p in range(n_procs):
    windows_dict[p] = []

  # list of processes we use
  procs = []

  load = (end - start )// n_procs
  chunks = []

  start_p = start
  end_p = start_p + load
  for p in range(n_procs):

    chunks.append((start_p, end_p ))
    start_p = end_p
    end_p += load


  print("Used chunks: {0}".format(chunks))
  sys.stdout.flush()


  time_start = time.perf_counter()

  for p in range(n_procs):
     procs.append(Process(target=process_worker, args=(p,
                                                       chunks[p][0], chunks[p][1],
                                                        windows_dict)))
     procs[p].start()

  # wait here and join the processes
  for p in range(n_procs):
    procs[p].join()


  counter = 0
  for p in windows_dict:
    print("Process {0} created"
          " {1} windows ".format(p, len(windows_dict[p])))
    sys.stdout.flush()
    counter += len(windows_dict[p])

  time_end = time.perf_counter()
  print("Time to create  {0} windows"
        " with {1} processes is {2} secs".format(counter, n_procs, time_end - time_start))
  sys.stdout.flush()



