"""
Utilities for working with bam files
"""
import math
import statistics
from collections import Counter, OrderedDict
from helpers import Window, Observation
from exceptions import FullWindowException
from exceptions import Error


def extract_windows(chromosome, ref_file, test_file, **args):
    """
    Extract the windows that couple the seq_file and ref_files
    for the given chromosome
    :param chromosome: chromosome name (str)
    :param ref_file: The reference file
    :param test_file: The sequence file
    :return:
    """

    windowsize = args["windowsize"]
    start_idx = args["start_idx"]
    end_idx = args["end_idx"]

    try:

        ref_list = ref_file.fetch(chromosome, 0, )

        for item in ref_list:
            print(item)

        bam_out, errors, adjusted = bam_strip(chromosome=chromosome,
                                              file=test_file,
                                              start=start_idx, stop=end_idx)

        print("\t Number of erros: ", errors)
        print("\t Number of adjusted: ", adjusted)
        print("\t bam output: ", len(bam_out))

        # find insertions and deletions
        indel_dict = get_indels(chromosome=chromosome, samfile=test_file,
                                start=start_idx, stop=end_idx)

        # get the common bases
        # TODO: explain what are we trying to do here
        common_bases(bamdata=bam_out, fastadata=ref_list)

        # extract the windows
        windows = create_windows(bamlist=bam_out,
                                  indel_dict=indel_dict,
                                  fastdata=ref_list,
                                  windowsize=windowsize,
                                  start=start_idx,
                                  end=end_idx)

        return windows
    except Exception as e:
        print(str(e))
        raise


def bam_strip(chromosome, file, start, stop):

    """
    get the contents of the file starting at start position and ending
    at stop position for the given chromosome
    """
    errors = 0
    adjusted = 0
    bam_out = []

    counter = 0
    # move column-wise
    for pileupcolumn in file.pileup(chromosome, start, stop,
                                        truncate=True, ignore_orphans=False):

            bam_out, adjusted_tmp, errors_tmp = \
                get_query_sequences(pileupcolumn=pileupcolumn, bam_out=bam_out)
            adjusted += adjusted_tmp
            errors += errors_tmp

    return bam_out, errors, adjusted


def get_indels(chromosome, samfile, start, stop):

    """
    find  insertions/deletions  and put in dictionary
    """
    indels = {}

    for i, pileupcolumn in enumerate(samfile.pileup(chromosome, start, stop,
                                                    truncate=True, ignore_orphans=False)):
        indel_count = 0
        for pileupread in pileupcolumn.pileups:

            # start counting indels when
            # first indel is encountered

            if (pileupread.indel != 0):
                indel_count += 1

            # TODO: Do we really need all the if-else stuff here?
            if (indel_count == pileupcolumn.n) and pileupcolumn.n > 1:  # homozygous indels.
                indel_1 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_1)

            elif indel_count >= 0.5 * pileupcolumn.n and pileupcolumn.n > 1:  # heterozygous indels.
                indel_2 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_2)

            elif (indel_count > 0):  # spontaneous indels.
                indel_3 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_3)

    return indels

def accumulate_windows(bamlist, windowsize):
    pass


def add_window_observation(window, windows, observation, windowsize):

    if window.has_capacity():
        window.add(observation=observation)
    else:
        windows.append(window)
        window = Window(capacity=windowsize)
        window.add(observation=observation)
    return window

def create_windows(bamlist, indel_dict, fastdata, windowsize, start, end):

    """
    Arrange the given bamlist into windows of size windowsize
    :param bamlist:
    :param indel_dict:
    :param fastdata:
    :return:
    """

    if not bamlist:
        raise Error("No test sequence is provided")

    if not fastdata:
        raise Error("No reference sequence is provided")

    # the returned list of windows
    tmp_windows = []
    windows = []

    window = Window(capacity=windowsize)
    previous_observation = None

    for idx, item in enumerate(bamlist):

        # create an observation
        observation = Observation(position=item[0],
                                  read_depth=item[1],
                                  base= item[2])

        if previous_observation is not None:

            # is this sequential?
            if int(observation.position) == int(previous_observation.position) + 1:

                # yes it is...nice add it to the window
                # and update the observation
                window = add_window_observation(window=window, windows=windows,
                                                observation=observation, windowsize=windowsize)

                previous_observation = observation
            else:
                # there is a gap we cannot simply
                # add the observation as this may have
                # to be added to the next window depending
                # on the size of the gap.

                # fill in the missing info from the
                # reference file
                window_gaps = _get_missing_gap_info(start=int(previous_observation.position)-1,
                                                    end=int(observation.position)-1, fastdata=fastdata)

                # after getting the missing info we try to add it
                # to the window we have accumulated so much info that
                # we may exceed the window capacity. For example
                # a window already has windowcapacity - 2 items
                # and len(window_gaps) is 10. In this case we need
                # to create a new window

                for win_gap_item in window_gaps:
                    dummy_observation = Observation(position=win_gap_item[0],
                                                    read_depth=win_gap_item[1],
                                                    base= win_gap_item[2])

                    window = add_window_observation(window=window, windows=windows,
                                                    observation=dummy_observation, windowsize=windowsize)

                # add also the current observation that led us
                # heer
                window = add_window_observation(window=window, windows=windows,
                                                observation=observation, windowsize=windowsize)

                previous_observation = observation
        else:

            # that's the first observation
            window = add_window_observation(window=window, windows=windows,
                                            observation=observation, windowsize=windowsize)
            previous_observation = observation

    return windows


def get_query_sequences(pileupcolumn, bam_out, use_indels=True, do_not_use_indels_on_error=True):
    """
    Given a pysam.PileupColumn update the bam_out list with the relevant entries
    corresponding to the given column. That is append to bam_out the following
    sublist [pileupcolumn.reference_pos, pileupcolumn.n, [list of bases mapping to the reference_pos]]
    If an error occurs, the sublist of bases mapping to the reference position is substituted
    with the string "ERROR". Client code then should decide how to fill this entry or disregard it
    altogether

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

    temp.append(pileupcolumn.reference_pos)
    # number of reads mapping to this column
    temp.append(pileupcolumn.n)

    try:
        # append bases
        temp.append(pileupcolumn.get_query_sequences(add_indels=use_indels))
        bam_out.append(temp)
    except:

        # try a fall out extra step only if it makes sense
        if do_not_use_indels_on_error and use_indels:
            try:

                # appending bases
                temp.append(pileupcolumn.get_query_sequences(add_indels=False))

                # flag to show indels not assessed.
                # TODO: Do we really need that?
                temp.extend("*")
                bam_out.append(temp)

                # once in here we always adjust
                adjusted += 1
            except:
                errors += 1
                temp.append(["ERROR"])
                bam_out.append(temp)
        else:
            errors += 1
            temp.append(["ERROR"])
            bam_out.append(temp)

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
    reference sequence given in fastadata
    :param bamdata: The bamdata to fill in
    :param fasta: The reference data
    :return:
    """

    for i, x in enumerate(bamdata):

        # x is expected to be of the form
        # [reference_pos, count_number, [list of bases.  usually a single base]]
        try:
            if x[1] == 1 and len(x) < 3:
                # appending the corresponding base for
                # this position from hg38 ref fasta
                x.append((fastadata[x[0] + 1]))

            # change to searching fasta and appending
            # corresponding base in fasta file instead.
            elif x[2] == '':
                del (x[2])

                # plus one because of the 0 indexing of fasta
                # list compared to the 1 indexing of the bam positions.
                x.insert(2, (fastadata[x[0] + 1]))

            else:
                # provides the element which reaches the
                # highest occurrence first.
                common_count = Counter(x[2]).most_common(1)

                # delete from the original bam entry what is present
                # from the entry position 2 and backwards

                del (x[2:])

                # when the most common mapped base to the position is an
                # indel then all elements of the string are appended
                # to the list (see screenshot on windows staff account).
                if len(common_count[0][0]) > 1:

                    indel = common_count[0][0]
                    x.extend(indel[0])
                else:

                    # extend adds the new elements into the list,
                    # not into the list as a separate list.
                    x.extend(common_count[0][0])

        except Exception as e:
            raise Error("An error occurred whilst extracting common_bases")


def get_winsize(data):

    adj_winsize = 0
    for x in data:

        # if the base for the position = "N"
        if x[2] != "N" or x[2] != "n":
            adj_winsize += 1
        else:
            return
    return adj_winsize


# use the fasta sequence to alter the GC counts:
# if gap or no reads/bases then use the bases in the fasta
# file for the skipped positions to calculate the gc content of window.







def sum_rd(data, expected_size, net_indel):

    """
    Adjust the given data with respect to net_indel.
    Only adjust for indels not for reference mapped gaps
    :param data:
    :param expected_size: The expected size of the window
    :param net_indel: total sum of indels
    :return:
    """
    sums = []
    winsize = 0

    for x in data:

        # if the base for the position = "N"
        if x[2] != "N" or x[2] != "n":

            # new method which excludes any
            # positions marked N from the calculation,
            # allowing the GC average (here)
            # and sum RD for a window (sum_rd function)
            # to be adjusted.
            sums.append(x[1])

            # we do not adjust winsize for any skipped positions
            # that mapped to a base in the reference fasta as
            # these may be bases that are supressed
            # by TUF so scaling up the rd of the window would make them seem regular.
            # similarly its important to scale up windows
            # with Ns as these are unmapped positions, that
            # will lower the rd of windows and make finding TUF too difficult.
            winsize += 1

    if net_indel is None or net_indel == 0.0:
        return sum(sums)
    else:
        # we use winsize to adjust for the indels
        # but not compensating for the gap_size
        # always divide by winsize as we do not want
        # to compensate for reference mapped gaps,
        # these could represent tuf regions/cores.

        adjusted_rd = (round(sum(sums) / (winsize + net_indel)) * expected_size)
        return adjusted_rd


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


def _uniquefy_insertions_and_deletions(insertions, deletions):

    insertions_set = set(insertions)
    deletions_set = set(deletions)

    # grab the third character in the string, which is he number of inserted/deleted bases.
    unique_insertions = [int(x[2]) for x in insertions_set]
    unique_deletions = [int(x[2]) for x in deletions_set]
    return unique_insertions, unique_deletions


def _get_insertion_deletion_difference(unique_insertions, unique_deletions):
        """
        Returns the absolute difference
        between the two given lists
        :param unique_insertions:
        :param unique_deletions:
        :return:
        """
        sum_unique_insert = sum(unique_insertions)
        sum_unique_delete = sum(unique_deletions)
        net_indel = math.fabs(sum_unique_delete - sum_unique_insert)
        return net_indel


def _get_missing_gap_info(start, end, fastdata):

    fastdata_range = fastdata[start:end]

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

    return window_gaps


def _get_skip_insert_positions(cur_pos, prev_pos, prev_end,
                               fastadata, win_gaps, windowsize):

    # the position of the skipped element.
    skipped_pos = (int(prev_pos) + 1)

    # last two digits of prev_position
    # before the gap, + 1 (first missing position)
    insert_pos = (prev_end + 1)

    # used to insert skipped_temp into the right position in cur_win later on.
    # iterating the data_list to pull bases corresponding to the skipped positions

    for base in fastadata[(int(prev_pos) + 1):(int(cur_pos) + 1)]:

        # in the bam file and filling the skipped
        # positions in cur_win with [pos, rd, base].
        if int(cur_pos) > ((int(prev_pos) + windowsize) - prev_end):

            if len(win_gaps) != 0:
                win_gaps = []

        # single wins retain win_gaps until the
        # end of their cur_win iteration
        # whilst win_gaps is emptied everytime a
        # new cross_win gap has been dealt with.

        skipped_temp = []
        skipped_temp.append(skipped_pos)

        # 0 to represent the read depth for this position.
        skipped_temp.append(0)
        skipped_temp.append(base)
        skipped_temp.append(insert_pos)
        win_gaps.append(skipped_temp)

        # -1 so that cur_pos does not end up in cur_win twice.
        if skipped_pos == (int(cur_pos) - 1):
            break

        skipped_pos += 1
        insert_pos += 1

    return win_gaps, skipped_pos, insert_pos







