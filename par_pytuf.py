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


      #time_end = time.perf_counter()
      #print("Time for loop over pcol {0}".format(time_end - time_start))
      #sys.stdout.flush()

      #time_start = time.perf_counter()

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
      #print("Time for remaining function {0}".format(time_end - time_start))
      #sys.stdout.flush()

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

      #print("Memory used by function {0} GB".format(total*1e-9))
      #sys.stdout.flush()
      raise


def process_worker(p, start, end, windows_dict, errors_dict, msg_dict):

    fas = pysam.FastaFile("/scratch/spectre/a/ag568/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
    sam = pysam.AlignmentFile("/scratch/spectre/a/ag568/m585_verysensitive_trim_sorted.bam", "rb")

    c = 'chr1'
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
      msg_dict[p] = ("Process {0} finished. "
      "Total time for {1} windows {2} secs").format(p, len(windows), time_end - time_start)

    except MemoryError as e:

      total = total_memory(windows)
      #print("MemoryError exception detected in process {0}. "
      #     "Windows memory used is: {1} GB".format(p, total*1e-9))

      errors_dict[p] = ("MemoryError exception "
      "detected in process {0}. Windows constructed {1}").format(p, len(windows))

      return
    except Exception as e:
        errors_dict[p] = ("An exception detected "
        "in process {0}. Exception is: {1}").format(p, str(e))

        return

if __name__ == '__main__':


  n_procs = 5
  start = 1000000
  end = 2000000
  qual = 20

  manager = Manager()
  windows_dict = manager.dict()
  errors_dict = manager.dict()
  msg_dict = manager.dict()

  for p in range(n_procs):
    windows_dict[p] = []
    errors_dict[p] = "No error"
    msg_dict[p] = "No msg"

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

  for p in range(n_procs-1):
     procs.append(Process(target=process_worker, args=(p,
                                                       chunks[p][0], chunks[p][1],
                                                        windows_dict,
                                                        errors_dict,
                                                        msg_dict)))
     procs[p].start()

  print("Created: {0} processes".format(n_procs))
  sys.stdout.flush()

  p = n_procs -1
  print("Master process is: {0} ".format(p))
  sys.stdout.flush()

  print("Master process doing its share")
  sys.stdout.flush()

  p = n_procs -1
  process_worker(n_procs-1,
                 chunks[p][0], chunks[p][1],
                 windows_dict,
                 errors_dict,
                 msg_dict)

  print("Process {0} msg: {1}".format(p, msg_dict[p]))
  sys.stdout.flush()
  print("Process {0} errs: {1}".format(p, errors_dict[p]))
  sys.stdout.flush()

  # wait here and join the processes
  for p in range(n_procs-1):
    procs[p].join()
    print("Process {0} msg: {1}".format(p, msg_dict[p]))
    sys.stdout.flush()
    print("Process {0} errs: {1}".format(p, errors_dict[p]))
    sys.stdout.flush()


  counter = 0
  for p in windows_dict:

    if errors_dict[p] == "No error":

      print("Process {0} created "
            "{1} windows ".format(p, len(windows_dict[p])))
      sys.stdout.flush()
      counter += len(windows_dict[p])
    else:
      print(errors_dict[p])

  time_end = time.perf_counter()
  print("Time to create  {0} windows"
        " with {1} processes is {2} secs".format(counter, n_procs, time_end - time_start))
  sys.stdout.flush()


