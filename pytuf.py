import pysam
import numpy as np
from scipy import stats
import re
import time

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

    for pcol in sam.pileup(chr,start,end,
                           truncate=True,ignore_orphans=False,
                           fastafile=fas, max_depth=1000):
        pcol.set_min_base_quality(qual)

        print("Start first while loop")
        time_start = time.perf_counter()
    #fill in start when there are no reads present
        while pcol.reference_pos != start:
            samseq.append('_')
            refseq += fas.fetch(chr,start,start+1)
            pos.append(start)
            nseq.append(0)
            nalign.append(0)
            start+=1

        time_end = time.perf_counter()
        print("Done first while loop. Execution time"
          " {0} secs".format(time_end - time_start))

    #collect metrics for pileup column
        pos.append(start)
        nseq.append(pcol.nsegments) #all reads over a pileup column
        nalign.append(pcol.get_num_aligned()) #above but quality limited

        time_start = time.perf_counter()
        print("Start first for loop")
        for p in pcol.pileups:
            head += p.is_head #start of read present
            read1 += p.alignment.is_read1 #first of mate pair (from nalign)
            read2 += p.alignment.is_read2 #second of mate pair (from nalign)
            if p.is_head == 1 and p.alignment.is_read1 == 1: #above but is start of read
                head1 += 1
            if p.is_head == 1 and p.alignment.is_read2 == 1:
                head2 += 1

        time_end = time.perf_counter()
        print("Done first for loop. Execution time"
          " {0} secs".format(time_end - time_start))

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

    #fill in end if no reads at end of window
    time_start = time.perf_counter()
    print("Start second while loop")
    while start < end:
            samseq.append('_')
            refseq += fas.fetch(chr,start,start+1)
            pos.append(start)
            nseq.append(0)
            nalign.append(0)
            start+=1
    time_end = time.perf_counter()
    print("Done second while loop. Execution time"
          " {0} secs".format(time_end - time_start))

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

    time_start = time.perf_counter()
    print("Start second for loop")
    altseq = ['']
    for bs in samseq:
        talt = []
        for i in range(len(altseq)):
            for el in bs:
                alt = altseq[i] + el
                talt.append(alt)
        altseq = talt
    time_end = time.perf_counter()
    print("Done second loop. Execution time"
          " {0} secs".format(time_end - time_start))

    gcmax = 0
    gcmin = 1

    time_start = time.perf_counter()
    print("Start third loop")
    for alt in altseq:
        t1 = re.sub('[\+\-_Nn*]','',alt)
        gct = len(re.sub('[\+\-_Nn*AaTt]','',alt))/len(re.sub('[\+\-_Nn*]','',alt))
        gcmax = gct if gct > gcmax else gcmax
        gcmin = gct if gct < gcmin else gcmin

    time_end = time.perf_counter()
    print("Done third loop. Execution time"
          " {0} secs".format(time_end - time_start))

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
            'end':pos[-1]
            }
    return output

if __name__ == '__main__':
  fas = pysam.FastaFile("/scratch/spectre/a/ag568/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
  sam = pysam.AlignmentFile("/scratch/spectre/a/ag568/m585_verysensitive_trim_sorted.bam", "rb")


  c = 'chr1'
  start = 1000000
  end = 2000000
  qual = 20
  wcounter = 0

  time_start = time.perf_counter()


  while start < end:
      print("==============================")
      met = windowAna(c,start,start+100,qual,fas,sam)
      print("Created window: {0}".format(wcounter))
      start = start + 100

      if wcounter == 10:
        break

  time_end = time.perf_counter()
  print("Total time for {0} is {1} secs".format(wcounter, time_end - time_start))