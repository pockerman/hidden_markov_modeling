"""
Utilities for working with bam files
"""
import math
from helpers import Window, Observation
from exceptions import FullWindowException



def extract_windows(chromosome, ref_file, test_file, **args):
    """
    Extract the windows that couple the seq_file and ref_files
    for the given chromosome
    :param chromosome: chromosome name (str)
    :param ref_file: The reference file
    :param test_file: The sequence file
    :return:
    """

    windowsize = args["window_size"]
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

        windows  = create_windows(bamlist=bam_out, indel_dict=indel_dict,
                                  windowsize=windowsize)

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

            adjusted_tmp, errors_tmp = get_query_sequences(pileupcolumn=pileupcolumn, bam_out=bam_out)
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

def create_windows(bamlist, indel_dict, windowsize):

    """
    Arrange the given bamlist into windows of size windowsize
    :param bamlist:
    :return:
    """

    # the returned list of windows
    windows = []

    # list to hold all skipped temps (gaps) that are to be added to this window,
    # cannot append the gaps directly to cur_win whilst cur_win is being iterated
    # otherwise an infinite loop is created where the inserted gap_data is constantly appended/iterated.
    win_gaps = []

    cur_win = Window(capacity=windowsize)
    old_window = None
    for idx, item in enumerate(bamlist):

        # these two lists contain the lengths (as integers) of all the UNIQUE indels
        # that are discovered during the itteration of indiviual cur_win states.
        unique_insertions = []
        unique_deletions = []

        # create an observation to be added
        # to the current window
        observation = Observation()
        observation.position = item[0]
        observation.read_depth = item[1]
        observation.base = item[2][0]

        # this will throw when attempting to add to a window
        # that is already filled
        try:
            cur_win.add(observation=observation)
        except FullWindowException as e:
            old_window = cur_win
            cur_win = Window(capacity=windowsize)
            cur_win.add(observation=observation)

        # if the position is divisible
        # by the window size
        if item[0] % windowsize == 0:

            cur_window_start = 0
            cur_window_end = windowsize

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

                if i != 0:

                    # for all elements except the first in cur_win,
                    # the position of the current element is
                    # compared to the position of the last element to
                    # identify skipped position in cur_win.
                    cur_pos = observation.position

                    # that won't work if we just created a new window
                    # we need to request it from the old window if
                    # there are no as many entries in the current window
                    prev_pos = cur_win[i - 1].position

                    # changing prev_pos to a string to slice the last two characters of
                    # the position and assigning them to a variable (prev_end)
                    # which is used in calculations later on.
                    prev_pos = str(prev_pos)
                    prev_end = prev_pos[-2:]
                    prev_end = int(prev_end)

                    # checking if the current element is
                    # present in the indel dictionary
                    if str(cur_pos) in indel_dict:

                        # get the value (a list of the bases/indels)
                        # of the key (the position)
                        indel = indel_dict.get(str(cur_pos))

                        # iterating through the value (base, indel list)
                        # for this element in cur_win.
                        _get_insert_delete_from_indel(indel=indel,
                                                      insertions=insertions,
                                                      deletions=deletions)

                        # because insertions/deletions lists exist for only one element at a time it is likely the same
                        # indel will occur several times in the base lists as it may occur
                        # in several reads that map to this element's position.
                        # In order to not count the same indel several times the temp lists are
                        # changed to sets which cannot contain element duplicates.
                        unique_insertions, unique_deletions = \
                            _uniquefy_insertions_and_deletions(insertions=insertions, deletions=deletions)

                    # execute code block when skipped positions/gaps occurs
                    if int(cur_pos) != (int(prev_pos) + 1):
                        _get_skip_insert_positions(cur_pos=cur_pos, prev_pos=prev_pos,
                                                   prev_end=prev_end,
                                                   data_list=None, win_gaps=win_gaps, windowsize=windowsize)

                        # once data for the gap has been pulled from the fasta file and appended to win_gaps above
                        # the position of the current element (cur_pos) and the last element (prev_pos) are checked
                        # to see if the gap spans multiple windows.
                        if int(cur_pos) > ((int(prev_pos) + windowsize) - prev_end):
                            net_indel = _get_window_indel_sum(unique_insertions=unique_insertions,
                                                          unique_deletions=unique_deletions)

                            # inserting the skipped gap data into the
                            # cur_win when the gap spans across more than one window.
                            for gap_data in win_gaps:
                                cur_win.insert((gap_data[3] - 1), gap_data[:3])

                            # generating data for windows containg
                            # gaps that span out of the window.
                            win_sum = sum_rd(cur_win[cw_start:(cw_end)], 100, net_indel)

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
                                                             "nb": window[3], "win_length": window[4], "gc_count": window[5]}}

                            chr1_dictionary.update(position)

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

                            for skipped_win in cur_win[cw_start:gap_end: windowsize]:
                                skip_count = 0

                                # adding 'empty' data for skipped windows to the dictionary.
                                # Is it correct to put 0 for the pt/nb?
                                # is it equivalent to rd in that respect.
                                window = [skipped_win[0], 0, 0, 0, 100, '-']

                                position = {str(window[0]): {"rd": window[1], "pt": window[2],
                                                             "nb": window[3], "win_length": (window[4], "**"),
                                                             "gc_count": window[5]}}

                                chr1_dictionary.update(position)
                                skip_count += 1

                                cw_start += (windowsize * skip_count)

                                # when exiting the for loop, cw_start/end need to be adjusted for the number of windows skipped,
                                # so that any remaining windows in cur_win can be calculated correctly.
                                cw_end += (windowsize * skip_count)

                            start_counter += 100

                            if i < (iter_r_len + 1):
                                r = r_values[counter]
                            else:
                                iter_r_len = (iter_r_len + r_len)
                                if counter < (len(r_values) - 1):
                                    counter += 1
                    elif int(cur_pos) == cur_win[-1][0]:

                        net_indel = _get_window_indel_sum(unique_insertions=unique_insertions,
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
                        chr1_dictionary.update(position)
                        cw_start += windowsize
                        cw_end += windowsize
                        start_counter += windowsize


                        cur_win = []

                        if i < (iter_r_len + 1):
                            r = r_values[counter]
                        else:
                            iter_r_len = (iter_r_len + r_len)
                            if counter < (len(r_values) - 1):
                                counter += 1
                        break


def get_query_sequences(pileupcolumn, bam_out):

    temp = []
    temp.append(pileupcolumn.reference_pos)
    adjusted = 0
    errors = 0

    # number of reads mapping to this column
    temp.append(pileupcolumn.n)

    try:
        # append bases
        temp.append(pileupcolumn.get_query_sequences(add_indels=True))
        bam_out.append(temp)
    except:
        try:

            # appending bases
            temp.append(pileupcolumn.get_query_sequences(add_indels=False))

            # flag to show indels not assessed.
            temp.extend("*")
            bam_out.append(temp)

            if len(temp) == 4:
                adjusted += 1
        except:
            errors += 1

    return adjusted, errors


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


def common_bases(data, fasta):

    from collections import Counter

    for i, x in enumerate(data):

        try:
            if x[1] == 1 and len(x) < 3:
                # appending the corresponding base for
                # this position from hg38 ref fasta
                x.append((fasta[x[0] + 1]))
                a = 1

            # change to searching fasta and appending
            # corresponding base in fasta file instead.
            elif x[2] == '':
                del (x[2])

                # plus one because of the 0 indexing of fasta
                # list compared to the 1 indexing of the bam positions.
                x.insert(2, (fasta[x[0] + 1]))
                a = 2
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
                a = 3
        except Exception as e:
            print(e, x)


def partition_data(bamfile, fastafile):
    pass


def fetch_data(chromosome, bamfile, start, stop):

    counter =0
    for read in bamfile.fetch(chromosome, start, stop):

        if counter < 11:
            print(read)
            counter += 1


def _get_insert_delete_from_indel(indel, insertions, deletions):

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


def _get_window_indel_sum(unique_insertions, unique_deletions):

        sum_unique_insert = sum(unique_insertions)
        sum_unique_delete = sum(unique_deletions)
        net_indel = math.fabs(sum_unique_delete - sum_unique_insert)
        return net_indel

def _get_skip_insert_positions(cur_pos, prev_pos, prev_end,
                              data_list, win_gaps, windowsize):

    # the position of the skipped element.
    skipped_pos = (int(prev_pos) + 1)

    # last two digits of prev_position before the gap, + 1 (first missing position)
    insert_pos = (prev_end + 1)

    # used to insert skipped_temp into the right position in cur_win later on.
    # iterating the data_list to pull bases corresponding to the skipped positions

    for base in data_list[(int(prev_pos) + 1):(int(cur_pos) + 1)]:

        # in the bam file and filling the skipped positions in cur_win with [pos, rd, base].
        if int(cur_pos) > ((int(prev_pos) + windowsize) - prev_end):

            if len(win_gaps) != 0:
                win_gaps = []



        # this code is executed regardless of whether a gap spans several windows or not, allows
        # the indexing code above to work as single wins retain win_gaps until the end of their cur_win itteration
        # whilst win_gaps is emptied everytime a new cross_win gap has been dealt with.
        skipped_temp = []
        skipped_temp.append(skipped_pos)
        skipped_temp.append(0)  # 0 to represent the read depth for this position.
        skipped_temp.append(base)
        skipped_temp.append(insert_pos)
        win_gaps.append(skipped_temp)

        if skipped_pos == (int(cur_pos) - 1):  # -1 so that cur_pos does not end up in cur_win twice.
            break

        skipped_pos += 1
        insert_pos += 1

    return skipped_pos, insert_pos







