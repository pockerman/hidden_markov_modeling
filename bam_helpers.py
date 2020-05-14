"""
Utilities for working with bam files
"""
import math
import pysam
import logging
from collections import Counter
from helpers import Window, Observation, DUMMY_ID
from helpers import add_window_observation
from helpers import INFO, WARNING, DEBUG
from exceptions import FullWindowException
from exceptions import Error

DUMMY_BASE = "*"

def extract_windows(chromosome, ref_filename, test_filename, **args):

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
    quality_threshold = args.get("quality_threshold", None)

    with pysam.FastaFile(ref_filename) as ref_file:

        print("{0} Reference file: {1}".format(INFO, ref_file.filename))

        with pysam.AlignmentFile(test_filename, "rb") as test_file:
            print("{0} Alternative file: {1} ".format(INFO, test_file.filename))

            try:

                ref_list = ref_file.fetch(chromosome, 0, )

                logging.info("Working with test file: {0}".format(test_filename))

                bam_out, errors, adjusted =\
                  bam_strip(chromosome=chromosome,
                            file=test_file,
                            start=start_idx, stop=end_idx,
                            quality_threshold=quality_threshold,
                            fastadata=ref_list)

                print("{0} Number of errors: {1}".format(INFO, errors))
                print("{0} Number of adjusted: {1}".format(INFO, adjusted))
                print("{0} Bam length: {1}".format(INFO, len(bam_out)))

                if "debug" in args:
                  debug = args["debug"]
                  if "log_bam_for_debug" in debug and\
                    debug["log_bam_for_debug"]:
                      logging.info("Bam sequence: {0}".format(bam_out))

                print("{0} Extracting common bases".format(INFO))

                # get the common bases
                # TODO: explain what are we trying to do here
                bam_out=common_bases(bamdata=bam_out, fastadata=ref_list)

                if "debug" in args:
                  debug = args["debug"]
                  if "log_bam_for_debug" in debug and\
                    debug["log_bam_for_debug"]:
                      logging.info("Bam sequence (after common base): {0}".format(bam_out))

                # extract the windows
                windows = create_windows(bamlist=bam_out,
                                          indel_dict=None,
                                          fastdata=ref_list,
                                          windowcapacity=windowcapacity,
                                          start=start_idx,
                                          end=end_idx,
                                          fill_missing_window_data=args.get("fill_missing_window_data", False),
                                          fill_missing_window_data_factor=args.get("fill_missing_window_data_factor",0))

                return windows
            except Exception as e:
                print(str(e))
                raise


def bam_strip(chromosome, file, start, stop, **kwargs):

    """
    get the contents of the file starting at start position and ending
    at stop position for the given chromosome
    """
    errors = 0
    adjusted = 0
    bam_out = []

    quality_threshold = kwargs.get("quality_threshold", None)

    if quality_threshold is not None:
      print("{0} Using quality threshold {1}".format(INFO, quality_threshold))
    else:
      print("{0} Not using quality threshold".format(INFO))


    # move column-wise
    for pileupcolumn in file.pileup(chromosome, start, stop,
                                    truncate=True, ignore_orphans=False):


            # if there is a quality threshold then use it
            if quality_threshold is not None:
              pileupcolumn.set_min_base_quality(quality_threshold)


            bam_out, adjusted_tmp, errors_tmp = \
                get_query_sequences(pileupcolumn=pileupcolumn,
                                    bam_out=bam_out,
                                    use_indels=True,
                                    do_not_use_indels_on_error=True,
                                    quality_threshold=quality_threshold,
                                    fastadata=kwargs.get("fastadata", None))

            adjusted += adjusted_tmp
            errors += errors_tmp

    return bam_out, errors, adjusted


def create_windows(bamlist, indel_dict, fastdata,
                   windowcapacity, start, end, **kwargs):

    """
    Arrange the given bamlist into windows of size windowcapacity.
    Note that a window is not necessary that will have the
    windowcapacity items. windowcapacity simply indicates the
    maximum number of observations that should be
    added in a window
    :param bamlist:
    :param indel_dict: insertions/deletions directory
    :param fastdata: The reference sequence
    :return: a list of Window instances
    """

    if not bamlist:
        raise Error("No test sequence is provided")

    if not fastdata:
        raise Error("No reference sequence is provided")

    print("{0} Estimated number"
          " of windows: {1} ".format(INFO, len(bamlist)//windowcapacity))

    # the returned list of windows
    windows = []

    idstart = 0
    window = Window(idx=idstart, capacity=windowcapacity)
    previous_observation = None

    for idx, item in enumerate(bamlist):

        for base in item[2]:
          base.upper()



        # create an observation
        observation = Observation(position=int(item[0]),
                                  read_depth=item[1],
                                  base= item[2])

        if observation.position == 1005863:
          import pdb
          pdb.set_trace()

        if previous_observation is not None:

            # is this sequential?
            if int(observation.position) == int(previous_observation.position) + 1:

                # yes it is...nice add it to the window
                # and update the observation
                window = \
                  add_window_observation(window=window,
                                         windows=windows,
                                         observation=observation,
                                         windowcapacity=windowcapacity)

                previous_observation = observation
            else:

                import pdb
                pdb.set_trace()
                logging.info("For observation {0}"
                            " there is a gap. Next "
                            "observation is at {1}".format(previous_observation.position,
                                                           observation.position))


                gap_size = int(observation.position) - int(previous_observation.position)


                # there is a gap we cannot simply
                # add the observation as this may have
                # to be added to the next window depending
                # on the size of the gap.

                # fill in the missing info from the
                # reference file positions are adjusted
                # in _get_missing_gap_info to start +1, end +1
                # to access the referece section
                window_gaps = \
                  _get_missing_gap_info(start=int(previous_observation.position),
                                        end=int(observation.position),
                                        fastdata=fastdata)

                if len(window_gaps ) != gap_size:
                  raise Error("Invalid window_gaps. "
                              "Size {0} not equal to {1}".format(len(window_gaps ), gap_size))

                # after getting the missing info we try to add it
                # to the window. we may have accumulated so much info that
                # we exceed the window capacity. For example
                # a window already has windowcapacity - 2 items
                # and len(window_gaps) is 10. In this case we need
                # to create a new window

                for win_gap_item in window_gaps:
                    dummy_observation = Observation(position=win_gap_item[0],
                                                    read_depth=win_gap_item[1],
                                                    base=win_gap_item[2])

                    window = add_window_observation(window=window, windows=windows,
                                                    observation=dummy_observation,
                                                    windowcapacity=windowcapacity)

                # add also the current observation
                # that led us here
                window = add_window_observation(window=window, windows=windows,
                                                observation=observation,
                                                windowcapacity=windowcapacity)

                previous_observation = observation
        else:

            # that's the first observation
            window = add_window_observation(window=window, windows=windows,
                                            observation=observation,
                                            windowcapacity=windowcapacity)
            previous_observation = observation

    # catch also the last window. The last
    # window may not be using all its capacity
    # as this depends on the partitioning. Optionally
    # we fill in the missing data if that was requested
    if len(window) != window.capacity:

      print("{0} Window {1} size {2} is"
            " not equal capacity {3} ".format(WARNING,
                                              window.idx,
                                              len(window),
                                              window.capacity))

      # fill in missing data if this is requested
      if kwargs.get("fill_missing_window_data", False):

        miss_factor =kwargs["fill_missing_window_data_factor"]
        print("{0} Window size {1} is not \
              equal capacity {2} ".format(WARNING, miss_factor))

        while window.has_capacity():
          window.add(observation=Observation(position=DUMMY_ID,
                                  read_depth=miss_factor,
                                  base= [DUMMY_BASE]))

      windows.append(window)

    return windows


def get_query_sequences(pileupcolumn, bam_out,
                        use_indels=True,
                        do_not_use_indels_on_error=True,
                        **kwargs):

    """
    Given a pysam.PileupColumn instance, it updates the bam_out
    list with the relevant entries corresponding to the given column.
    That is it appends to bam_out the following sublist
    [pileupcolumn.reference_pos,
     pileupcolumn.n,
     [list of bases mapping to the reference_pos]]

    If an error occurs, the sublist of bases mapping to the
    reference position is substituted
    with the string "ERROR".
    Client code then should decide how to fill this entry or disregard it altogether

    :param pileupcolumn: the pysam.PileupColumn
    :param bam_out: list of lists. Each sublist has the
                    following structure [pileupcolumn.reference_pos,
                                         pileupcolumn.n,
                                         [list of bases mapping to the reference_pos]]

    :param use_indels: boolean whether indels should be used with the query
    :param do_not_use_indels_on_error: boolean if the query that uses indels fails
                                        and this flag is True allow for a fall over by
                                        performing the same query with use_indels = False

    :return: updated list, number of adjustments, number of errors
    """
    temp = []
    adjusted = 0
    errors = 0

    # add the reference position
    temp.append(pileupcolumn.reference_pos)

    """
    if int(pileupcolumn.reference_pos) == 1005897:

      import pdb
      pdb.set_trace()
      print("{0} Position {1} has"
            " been reached ".format(DEBUG, pileupcolumn.reference_pos))
    """
    # if the count is zero then we consult the reference
    # at this position if we have a reference file
    if pileupcolumn.nsegments == 0 and "fastadata" in kwargs:

      # when we consult the reference no RD
      temp.append(0)

      # plus 1 since bam is zero-bsed and FASTA is 1 based
      temp.append([kwargs["fastadata"][pileupcolumn.reference_pos + 1]])
    elif pileupcolumn.nsegments == 0:
      # we cannot do anything log the error and return
      logging.error("At position: {0} read \
                    count is zero and cannot use ref file.".format(pileupcolumn.reference_pos))
      return bam_out, 0, 1
    else:


      # get the read qualities for the column
      quality = pileupcolumn.get_query_qualities()
      quality_threshold = kwargs.get("quality_theshold", None)
      try:
          # append bases only if the quality of the read
          # satisfies our threshold

          # if
          query_sequences = \
            pileupcolumn.get_query_sequences(add_indels=use_indels)

          if len(query_sequences) != len(quality):
            logging.error("len(query_sequences) not equal\
                          to len(quality). pysam error?")
          else:

            if quality_threshold:
              filtered_bases = [ query_sequences[i] for i, q in enumerate(quality)
                                if q >= quality_threshold]


              # From the pysam documentation:
              # pileupcolumn.nsegments
              # ignores the base quality filter
              # Assume that the
              # number of segments is always larger
              # and subtract the
              nsegments = pileupcolumn.nsegments #- len(filtered_bases)
              temp.append(nsegments)
              temp.append(filtered_bases)
            else:
              temp.append(pileupcolumn.nsegments)
              temp.append(query_sequences)
            bam_out.append(temp)
      except Exception as e:

          # try a fall out extra step only if it makes sense
          if do_not_use_indels_on_error and use_indels:
              try:

                  query_sequences = pileupcolumn.get_query_sequences(add_indels=False)

                  if len(query_sequences) != len(quality):
                    logging.error("len(query_sequences) not equal to \
                                  len(quality). pysam error?")
                  else:

                    if quality_threshold:
                       filtered_bases = [ query_sequences[i]
                                         for i, q in enumerate(quality)
                                         if q >= quality_threshold]

                       nsegments = pileupcolumn.nsegments #- len(filtered_bases)
                       temp.append(nsegments)
                       temp.append(filtered_bases)
                    else:
                      temp.append(pileupcolumn.nsegments)
                      temp.append(query_sequences)

                    # flag to show indels not assessed.
                    # TODO: Do we really need that?
                    temp.extend("*")
                    bam_out.append(temp)

                    # once in here we always adjust
                    adjusted += 1
              except:
                  errors += 1
                  logging.error("Error at position {0}"
                                " at query.".format(pileupcolumn.reference_pos))
                  return bam_out, adjusted, errors
          else:
              errors += 1
              logging.error("Error at position {0}"
                            " at query.".format(pileupcolumn.reference_pos))
              return bam_out, adjusted, errors


    return bam_out, adjusted, errors


def create_indels_dictionary(chromosome, samfile, start, stop):
    """
    Insertion/deletions for the given position
    :param chromosome: the region to use
    :param samfile: BAM  file instance
    :param start:
    :param stop:
    :return:
    """

    indels = {}

    for i, pileupcolumn in enumerate(samfile.pileup(chromosome, start, stop,
                                                    truncate=True, ignore_orphans=False)):
        indel_count = 0
        for pileupread in pileupcolumn.pileups:

            # start counting indels when first
            # isertion/deletion  is encountered
            if (pileupread.indel != 0):
                indel_count += 1

            if indel_count > 0:
                indels.update({str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))})

    return indels


def common_bases(bamdata, fastadata):

    """
    Fill in the information in bamdata by using the
    reference sequence given in fastadata. The bamdata
    is a list containing lists of the form
    # [reference_pos, count_number, [list of bases.  usually a single base]]
    :param bamdata: The bamdata to fill in
    :param fastadata: The reference data
    :return:
    """

    for x in bamdata:

        try:
            # the bamdata should be [ref_pos, count, [bases]]
            # TODO:
            if x[1] == 1 and len(x) < 3:
                # appending the corresponding base for
                # this position from hg38 ref fasta
                #x.append((fastadata[x[0] + 1]))
                logging.error(" Bam item does not have the correct format")
            elif x[1] != 0  and x[2] == '':
                logging.warning("An entry of the form {0} found. Turnd RD to zero".format(str(x)))

                #  we miss a base so set the RD to zero
                x[1] = 0

                # and consult the refernce
                # (plus 1 since fasta is zero based  )
                x[2] = [fastadata[x[0] + 1]]
            else:
                # provides the element which reaches the
                # highest occurrence first.
                common_count = Counter(x[2]).most_common(1)

                try:
                    # when the most common mapped base to the position is an
                    # indel then all elements of the string are appended
                    # to the list (see screenshot on windows staff account).
                  if len(common_count) != 0:

                    # we will use the most common
                    del (x[2:])
                    indel = common_count[0][0]
                    x.extend([indel[0]])

                    """
                    if len(common_count[0][0]) > 1:
                      del (x[2:])
                      indel = common_count[0][0]
                      x.extend([indel[0]])
                    else:

                      del (x[2:])
                      # extend adds the new elements into the list,
                      # not into the list as a separate list.
                      x.extend([common_count[0][0]])
                    """
                  elif x[1] != 0 and x[2] == []:

                      logging.warning(" Found a delete marking it")
                      logging.warning(" x looked at is {0}".format(x))
                      # this is a delete mark is as such
                      x[2] = ["-"]

                  else:
                      logging.warning(" No common bases found don't know what to do")
                      logging.warning(" x looked at is {0}".format(x))

                except Exception as e:
                   print("Common count is: {0}".format(common_count))
                   print("x is: {0}" .format(x))
                   raise
        except Exception as e:
            raise Error("An error occurred whilst extracting\
                        common_bases {0}".format(str(e)))

    return bamdata


def _get_insertions_and_deletions_from_indel(indel, insertions, deletions):

    """
    For every base in the given indel check whether it contains
    "+" or "-". If yes append this  to the appropriate array
    :param indel:
    :param insertions: list representing insertions
    :param deletions: list representing deletions
    :return:
    """

    for base in indel:
        if "+" in base:
            insertions.append(base)
        if "-" in base:
            # for each element that contains one or more indels, the value (list of bases)
            # is iterated, deletions are represented in the read with a "+" in pysam, whilst
            # insertions are represented with a "-". Any indels that map to the position
            # are appended to deletions e.g. [T+2N] and insertions lists e.g. [A-3GTT]
            deletions.append(base)


def _get_missing_gap_info(start, end, fastdata):

    # plus 1 as this fasta is one based
    fastdata_range = fastdata[start + 1: end + 1 ]

    # the position of the skipped element.
    skipped_pos = start
    window_gaps = []

    for base in fastdata_range:

        skipped_temp = []
        skipped_temp.append(skipped_pos)

        # 0 to represent the read depth for this position.
        skipped_temp.append(0)
        skipped_temp.append(base)
        window_gaps.append(skipped_temp)

        # increase by one the reference position
        skipped_pos += 1

    return window_gaps
