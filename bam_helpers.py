"""
Utilities for working with bam files
"""


def extract_windows(chromosome, ref_file, test_file, **args):
    """
    Extract the windows that couple the seq_file and ref_files
    for the given chromosome
    :param chromosome: chromosome name (str)
    :param ref_file: The reference file
    :param test_file: The sequence file
    :return:
    """

    """
    if "start_test" not in args.keys():
        raise ValueError("Start position for test file has not been specified")
    else:
        start_test = int(args['start_test'])

    if "end_test" not in args.keys():
        raise ValueError("End position for test file has not been specified")
    else:
        stop_test = int(args['stop_test'])
    """

    #window_size = args["window_size"]
    #window_start = 0
    #window_end = window_size

    try:

        fetch_data(chromosome=chromosome, bamfile=test_file, start=0, stop=12000)

        """
        bam_out, errors, adjusted = bam_strip(chromosome=chromosome,
                                              file=test_file,
                                              start=0, stop=20000)

        print("\t Number of erros: ", errors)
        print("\t Number of adjusted: ", adjusted)
        print("\t bam output: ", len(bam_out))
        """
    except Exception as e:
        print(str(e))
        raise
    # extract common bases
    #common_bases(bam_out, ref_file)


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

            if counter < 11:
                print(bam_out)
                counter += 1

            adjusted += adjusted_tmp
            errors += errors_tmp

    return bam_out, errors, adjusted

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

    for read in bamfile.fetch(chromosome, start, stop):
        print(read)




