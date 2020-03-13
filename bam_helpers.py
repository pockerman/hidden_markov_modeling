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

        bam_out, errors, adjusted = bam_strip(chromosome=chromosome,
                                              file=test_file,
                                              start=start_idx, stop=end_idx)

        print("\t Number of erros: ", errors)
        print("\t Number of adjusted: ", adjusted)
        print("\t bam output: ", len(bam_out))

        indel_dict = get_indels(chromosome=chromosome, samfile=test_file,
                                start=start_idx, stop=end_idx)

        # get the common bases
        #common_bases(bamdata=bam_out, fastadata=None)

        windows  = create_windows(bamlist=bam_out,
                                  indel_dict=indel_dict,
                                  fastdata=None,
                                  windowsize=windowsize,
                                  start=start_idx,
                                  end=end_idx)

        return  windows
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
    find all indels and put in dictionary
    """
    indels = {}

    for i, pileupcolumn in enumerate(samfile.pileup(chromosome, start, stop,
                                                    truncate=True, ignore_orphans=False)):
        indel_count = 0
        for pileupread in pileupcolumn.pileups:

            # start counting indels whenfirst indel is encountered
            if (pileupread.indel != 0):
                indel_count += 1

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

def create_windows(bamlist, indel_dict, fastdata, windowsize, start, end):

    """
    Arrange the given bamlist into windows of size windowsize
    :param bamlist:
    :return:
    """

    # the returned list of windows
    windows = []

    # list to hold all skipped temps (gaps) that are to be added to this window,
    # cannot append the gaps directly to cur_win whilst cur_win is being iterated
    # otherwise an infinite loop is created where the inserted
    # gap_data is constantly appended/iterated.
    win_gaps = []

    # setting what should be the start position for every window,
    # this is used for checking if gaps occur at the start cur_pos and prev_pos
    # will always be 1 apart and so if cur_pos does not equal ....01, the gap from
    # 0 to cur_pos can be calculated.
    start_counter = start

    #cur_win = Window(capacity=windowsize)
    #old_window = None

    cur_win = []
    chr_dict = OrderedDict()

    # binomial distribution
    r_len = 1000000
    iter_r_len = 1000000
    rvalues = find_r(bamlist, r_len)
    counter=0
    r = rvalues[0]

    for idx, item in enumerate(bamlist):

        # these two lists contain the lengths (as integers) of all the UNIQUE indels
        # that are discovered during the iteration of individual cur_win states.
        unique_insertions = []
        unique_deletions = []

        # create an observation to be added
        # to the current window
        #observation = Observation(position=item[0],
        #                          read_depth=item[1],
        #                          base=item[2])

        # this will throw when attempting to add to a window
        # that is already filled
        #try:
        #    cur_win.add(observation=observation)
        #except FullWindowException as e:
        #    old_window = cur_win
        #    cur_win = Window(capacity=windowsize)
        #    cur_win.add(observation=observation)

        cur_win.append(item)

        # if the position is divisible by the window size
        if item[0] % windowsize == 0:

            cw_start = 0
            cw_end = windowsize

            # try to find the gaps for the current window
            gap_counter = 0
            window_gaps = []
            cont = False
            for i, observation in enumerate(cur_win):

                # insertions and deletions are temporary
                # lists that exist are reset for each element in cur_win.
                insertions = []
                deletions = []

                if not cont:
                    # this code block solves the problem of when the gap is
                    # appended to cur_win, the index is shifted and so
                    if len(window_gaps) != 0:

                        gap_counter += 1
                        if gap_counter != len(window_gaps) + 1:
                            cont = True
                            continue
                    else:
                        cont = True

                # checking to see if the position of
                # the first element in cur_win
                # matches the start position of
                # cur_win (start_counter). If that's not
                # the case we have a gap at the beginning
                if i == 0 and cur_win[i][0] != int(start_counter):
                    raise Error("Gap detected at the start of the window...")
                else:

                    # for all elements except the first in cur_win,
                    # the position of the current element is
                    # compared to the position of the last element to
                    # identify skipped position in cur_win.
                    cur_pos = observation.position

                    # that won't work if we just created a new window
                    # we need to request it from the old window if
                    # there are no as many entries in the current window
                    #if len(cur_win) == 1:
                    #    prev_pos = old_window[windowsize-1].position
                    #else:
                    prev_pos = cur_win[i - 1].position

                    # changing prev_pos to a string to slice the last two characters of
                    # the position and assigning them to a variable (prev_end)
                    # which is used in calculations later on.
                    prev_pos = str(prev_pos)
                    prev_end = prev_pos[-2:]
                    prev_end = int(prev_end)

                    # if the current position is in the indel
                    # dictionary
                    if str(cur_pos) in indel_dict:

                        # get the value (a list of the bases/indels)
                        # of the key (the position)
                        indel = indel_dict.get(str(cur_pos))

                        # get the insertions and deletions for the indel
                        _get_insertions_and_deletions_from_indel(indel=indel,
                                                                 insertions=insertions,
                                                                 deletions=deletions)

                        # because insertions/deletions lists exist for only
                        # one element at a time it is likely the same
                        # indel will occur several times in the base lists as it may occur
                        # in several reads that map to this element's position.
                        # In order to not count the same indel several times the temp lists are
                        # changed to sets
                        unique_insertions, unique_deletions = \
                            _uniquefy_insertions_and_deletions(insertions=insertions, deletions=deletions)

                    if int(cur_pos) != (int(prev_pos) + 1):

                        # we have skipped positions or gaps

                        # calculate skip gap positions
                        _get_skip_insert_positions(cur_pos=cur_pos, prev_pos=prev_pos,
                                                   prev_end=prev_end,
                                                   data_list=fastdata, win_gaps=win_gaps,
                                                   windowsize=windowsize)

                        # once data for the gap has been pulled
                        # from the fastadata file and appended to win_gaps above
                        # the position of the current element (cur_pos)
                        # and the last element (prev_pos) are checked
                        # to see if the gap spans multiple windows.
                        if int(cur_pos) > ((int(prev_pos) + windowsize) - prev_end):

                            # get the absolute difference between the
                            # accumulted insertions and deletions for the
                            # working window
                            net_indel = _get_window_insertion_deletion_difference(unique_insertions=unique_insertions,
                                                                                  unique_deletions=unique_deletions)

                            # inserting the skipped gap data into the
                            # cur_win when the gap spans across more than one window.
                            # TODO: fix this
                            for gap_data in win_gaps:
                                cur_win.insert((gap_data[3] - 1), gap_data[:3])

                            # generating data for windows containing
                            # gaps that span out of the window.
                            window_range = cur_win[cw_start: cw_end]
                            win_sum = sum_rd(window_range, windowsize, net_indel)

                            # the first win the gap occurs in is represented
                            # by cur_win[:-1], which would end at prev_pos
                            win_gc = gc_count(cur_win[cw_start:(cw_end)])
                            window = cat_win(cur_win[cw_start:(cw_end)], r=r, gc=win_gc, sum=win_sum)

                            # the skipped positions are appended to cur_win which
                            # allows the window to be generated
                            # for the window the gap occurs in. if the gap does not
                            # skip an entire window (checked for below),
                            # e.g. if the distance from the end of the window the gap occurs in
                            # and the end of the gap is >= 100, then the itteration resumes at
                            # the position that is the end of the gap + 1 and unless anoter cross_win gap occurs
                            # it will continue onto the if cur_pos == cur_win[-1][0] code block.
                            position = {str(window[0]): {"rd": window[1], "pt": window[2],
                                                         "nb": window[3], "win_length": window[4],
                                                         "gc_count": window[5]}}

                            chr_dict.update(position)

                            # each time a window is produced, cw_start/end
                            # needs to += windowsize, until the iteration
                            cw_start += windowsize

                            # of cur_win is over, at which point cw_start/end
                            # need to be reset to index the new cur_win.
                            cw_end += windowsize

                            # checking for skipped windows between the gap start/end
                            # and appending data for skipped windows.
                            # end of the window the gap occurs in.
                            gap_start_win = int(prev_pos) + (windowsize - int(prev_pos[-2:]))

                            # the difference between the end of the gap, and the end of the
                            gap_end = ((cur_pos - gap_start_win) + cw_start)

                            # if the no. of skipped bases between the end of
                            # the first window in the gap and the next recorded position is > 100
                            if (int(cur_pos) - int(gap_start_win)) >= windowsize:

                                skip_count = 0
                                for skipped_win in cur_win[cw_start: gap_end: windowsize]:


                                    # adding 'empty' data for skipped windows to the dictionary.
                                    # Is it correct to put 0 for the pt/nb?
                                    # is it equivalent to rd in that respect.
                                    window = [skipped_win[0], 0, 0, 0, 100, '-']

                                    position = {str(window[0]): {"rd": window[1], "pt": window[2],
                                                                 "nb": window[3], "win_length": (window[4], "**"),
                                                                 "gc_count": window[5]}}

                                    chr_dict.update(position)
                                    skip_count += 1

                                cw_start += (windowsize * skip_count)

                                # when exiting the for loop, cw_start/end need to be adjusted for the number of windows skipped,
                                # so that any remaining windows in cur_win can be calculated correctly.
                                cw_end += (windowsize * skip_count)

                            start_counter += 100

                            if i < (iter_r_len + 1):
                                r = rvalues[counter]
                            else:
                                iter_r_len = (iter_r_len + r_len)
                                if counter < (len(rvalues) - 1):
                                    counter += 1

                    elif int(cur_pos) == cur_win[-1][0]:

                        net_indel = _get_window_insertion_deletion_difference(unique_insertions=unique_insertions,
                                                                              unique_deletions=unique_deletions)

                        # insert_pos = gap_data[3] # not 0 indexing on the range function.
                        for gap_data in win_gaps:
                            cur_win.insert(int(gap_data[3] - 1), gap_data[:3])
                            win_gaps = []

                        win_sum = sum_rd(cur_win[cw_start : cw_end], windowsize, net_indel)

                        # the first win the gap occurs in is represented
                        # by cur_win[:-1], which would end at prev_pos.
                        win_gc = gc_count(cur_win[cw_start:cw_end])
                        window = cat_win(cur_win[cw_start:cw_end], r=r, gc=win_gc, sum=win_sum)

                        # if no gaps were present in range of
                        # cur_win (since win_gaps is reset each time a new window is produced.)
                        if not win_gaps:
                            position = {str(window[0]): {"rd": window[1], "pt": window[2],
                                                         "nb": window[3], "win_length": window[4],
                                                         "gc_count": window[5]}}
                        else:
                            # if gaps were present in this range of cur_win,
                            # add a marker '*' to the winlength to show
                            # that positons were skipped and sequence was pulled from fasta.
                            position = {str(window[0]): {"rd": window[1], "pt": window[2],
                                                         "nb": window[3], "win_length": (window[4], "*"),
                                                         "gc_count": window[5]}}
                        chr_dict.update(position)
                        cw_start += windowsize
                        cw_end += windowsize
                        start_counter += windowsize


                        cur_win = []

                        if i < (iter_r_len + 1):
                            r = rvalues[counter]
                        else:
                            iter_r_len = (iter_r_len + r_len)
                            if counter < (len(rvalues) - 1):
                                counter += 1
                        break

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

    indels = {}

    for i, pileupcolumn in enumerate(samfile.pileup(chromosome, start, stop,
                                                    truncate=True, ignore_orphans=False)):
        indel_count = 0
        for pileupread in pileupcolumn.pileups:

            # start counting indels when first indel is encountered
            if (pileupread.indel != 0):
                indel_count += 1

            # homozygous indels.
            if (indel_count == pileupcolumn.n) and pileupcolumn.n > 1:
                indel_1 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_1)

            # heterozygous indels.
            elif indel_count >= 0.5 * pileupcolumn.n and pileupcolumn.n > 1:
                indel_2 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_2)

            # spontaneous indels.
            elif (indel_count > 0):
                indel_3 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                indels.update(indel_3)

    return indels


def common_bases(bamdata, fastadata):

    """
    Fill in the information in bam_data by using the
    reference sequence described in fasta_data
    :param data:
    :param fasta:
    :return:
    """

    from collections import Counter

    for i, x in enumerate(bamdata):

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
                # indexing list to get tuple, indexing tuple to get only first element (base).

        except Exception as e:
            print(e, x)


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

def gc_count(data):

    gc_count = 0
    at_count = 0

    for win_element in data:  # CALCULATING GC CONTENT FOR NORMAL

        if win_element[2] == "C" or win_element[2] == "c"\
                or win_element[2] == "G" or win_element[2] == "g":
            gc_count += 1
        elif win_element[2] != "N" or win_element[2] != "n":
            at_count += 1
            # winsize = (winsize - 1) # new method which excludes any
            # positions marked N from the calculation,
            # allowing the GC average (here) and sum RD for a window (sum_rd function) to be adjusted.

    if gc_count == 0:
        return 0
    elif at_count == 0:
        return 1

    return (gc_count / (gc_count + at_count))


def cat_win(data, r, gc=None, sum=None):

    """
    concatenate data from bam_list/cur_win into windows,
    data = cur_win. cw_start/end = cur_win start/end.
    If gaps contained to one win, cw_start/end need no counter.
    when cur_win (data) is input to function, the range of
    the window is specified is specified,
    but the function will work without specifying a range.
    """
    window = []

    # the end position of the window.
    window.insert(0, data[-1][0])

    if sum != None:
        window.insert(1, (sum))

    window.insert(2, poisson_transformation(window[1], decimals=2))
    window.insert(3, nb_transformation(window[1], r, decimals=2))

    # using the winsize function that removes any Ns from winsize.
    window.insert(4, get_winsize(data))

    if gc != None:
        window.insert(5, round(gc, 2))

    return window

def find_r(data, r_len):

    nbtemp = []
    r_values = []
    r_num = 1

    for i, x in enumerate(data):
        nbtemp.append(x[1])
        # using i, not x[0] because there are
        # too many gaps, e.g. pos 1Mb = i 900kb
        if i == (r_len * r_num):
            r_values.append((statistics.mean(nbtemp)**2)/(statistics.stdev(nbtemp)**2 - statistics.mean(nbtemp))) # cannot use numpy, error with package taking too long to fix.
            nbtemp = []
            r_num += 1
    return r_values

def nb_transformation(x, r, decimals):
    # x = window
    return round(float(2*(r**float(0.5))*math.log(((float(x)+0.25)/((100*r)-0.5))**0.5+(1+((float(x)+0.25)/((100*float(r))-0.5)))**0.5)), decimals)


def poisson_transformation(x, decimals):
    return round(float(2*(((float(x)+0.25)/100)**0.5)), decimals)

def partition_data(bamfile, fastafile):
    pass


def fetch_data(chromosome, bamfile, start, stop):

    counter =0
    for read in bamfile.fetch(chromosome, start, stop):

        if counter < 11:
            print(read)
            counter += 1


def sum_rd(data, expected_size, net_indel):  # we only adjust for indels, not for reference mapped gaps.

    """
    Adjust the given data with respect to net_indel
    :param data:
    :param expected_size:
    :param net_indel:
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
            # that mapped to a base in the refenece fasta as
            # these may be bases that are supressed
            # by tuf so scaling up the rd of the window would make them seem regular.
            # similarly its important to scale up windows
            # with Ns as these are unmapped positions, that
            # will lower the rd of windows and make finding TUF too dificult.
            winsize += 1

    if net_indel == None:
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


def _get_window_insertion_deletion_difference(unique_insertions, unique_deletions):
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


def _get_skip_insert_positions(cur_pos, prev_pos, prev_end,
                               data_list, win_gaps, windowsize):

    # the position of the skipped element.
    skipped_pos = (int(prev_pos) + 1)

    # last two digits of prev_position
    # before the gap, + 1 (first missing position)
    insert_pos = (prev_end + 1)

    # used to insert skipped_temp into the right position in cur_win later on.
    # iterating the data_list to pull bases corresponding to the skipped positions

    for base in data_list[(int(prev_pos) + 1):(int(cur_pos) + 1)]:

        # in the bam file and filling the skipped positions in cur_win with [pos, rd, base].
        if int(cur_pos) > ((int(prev_pos) + windowsize) - prev_end):

            if len(win_gaps) != 0:
                win_gaps = []

        # this code is executed regardless of
        # whether a gap spans several windows or not,
        # allows the indexing code above to work as
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

    return skipped_pos, insert_pos







