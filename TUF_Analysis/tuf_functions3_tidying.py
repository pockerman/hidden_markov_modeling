import statistics, math

def bam_strip(Chr, file, fasta, start, stop): # file must be SAM/BAM format #op = operation

    errors = 0
    adjusted = 0
    bam_out = [] # Chr, file, op must be type str # Chr must be name of chrom in bam file
    window = []
    for pileupcolumn in file.pileup(Chr ,start, stop, truncate=True,
                                       ignore_orphans=False):
        temp = []
        temp.append(pileupcolumn.reference_pos)
        temp.append(pileupcolumn.n)
        try:
            temp.append(pileupcolumn.get_query_sequences(add_indels=True)) # appending bases
        except:
            try:
                temp.append(pileupcolumn.get_query_sequences(add_indels=False)) # appending bases
                temp.append("*") # flag to show indels not assessed.
                #print("adjusted base ", "\n", temp)
                if len(temp) == 4:
                    adjusted += 1
            except:
                errors += 1
                #(print(pileupcolumn.reference_pos))
        window.append(temp)

        if (pileupcolumn.reference_pos % 100) == 0:
            common_bases(window, fasta)
            bam_out.append(window)
            window = []

    print("errors = ", errors, "adjusted = ", adjusted)
    return bam_out

# TESTED ////////////////////////////////////////////////////////////////////////////////////////////

def common_bases(data, fasta):
    
    from collections import Counter

    for i, x in enumerate(data):

        try:
            # if x[0] == 122524712:
            #     print(x)
            #     exit()
            if len(x) < 3:
                x.append((fasta[x[0] + 1])) # appending the corresponding base for this position from hg38 ref fasta
                a = 1
                #print("filling the list")
            elif x[2] == '': # change to searching fasta and appending corresponding base in fasta file instead.
                del(x[2])
                x.insert(2, (fasta[x[0] + 1])) # plus one because of the 0 indexing of fasta list compared to the 1 indexing of the bam positions. 
                #print("adding fasta bases", fasta[x[0]])
                a = 2
            else:
                common_count = Counter(x[2]).most_common(1) # provides the element which reaches the highest occurance first. 
                del(x[2])
                if len(common_count[0][0]) > 1: # when the most common mapped base to the position is an indel then all elements of the string are appended to the list (see screenshot on windows staff account).
                    indel = common_count[0][0]
                    x.insert(2, indel[0])
                else:
                    x.insert(2, common_count[0][0]) # extend adds the new elements into the list, not into the list as a separate list.
                a = 3                            # indexing list to get tuple, indexing tuple to get only first element (base).
        except Exception as e:
            print(e, x)
            #print(a)
    return

# ///////////////////////////////////////////////////////////////////////////////////////////////////

def indels(Chr, samfile, fastafile, start, stop): # find all indels and put in dictionary

    indels = {}

    try:
    
        for i, pileupcolumn in enumerate(samfile.pileup(Chr ,start, stop, 
                                        truncate=True, ignore_orphans=False)):
            indel_count = 0
            for pileupread in pileupcolumn.pileups:
                if (pileupread.indel != 0): #start counting indels whenfirst indel is encountered
                    indel_count +=1

                if (indel_count == pileupcolumn.n) and pileupcolumn.n > 1: # homozygous indels.
                    indel_1 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                    indels.update(indel_1)

                elif indel_count >= 0.5* pileupcolumn.n and pileupcolumn.n > 1: # heterozygous indels.
                    indel_2 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                    indels.update(indel_2)
                    
                elif (indel_count > 0): # spontaneous indels.
                    indel_3 = {str(pileupcolumn.pos):
                               (pileupcolumn.get_query_sequences(add_indels=True))}
                    indels.update(indel_3)
    
    except Exception as e:
        print(e, "\n", i)
    return indels

# ///////////////////////////////////////////////////////////////////////////////////////////////////

def gap_fill(bam_list, fasta_list):

    prev_window = None

    for window in bam_list:
        #print(window) # window before gap_fill
        index = 0
        for gap_data in window:
            if index != 0: # for all elements except the first in the window, the position of the current element is
                           # compared to the position of the last element to identify skipped positions.
                cur_pos = window[index][0]
                prev_pos = window[index-1][0]
                if cur_pos == window[-1][0]: # remvoes list index out of range error when the end of the itteration is reached.
                    #print(window) # window after gap_fill
                    prev_window = window[-1][0]
                    break
                if int(cur_pos) != (int(prev_pos) + 1):
                    #print("BLOCK 1", "prev_pos = ", prev_pos, "cur_pos = ", cur_pos)
                    position = int(prev_pos) + 1
                    for fasta_base in fasta_list[(int(prev_pos) + 1):(int(cur_pos))]:
                        skipped_temp = []
                        skipped_temp.append(int(position))
                        skipped_temp.append(0)
                        skipped_temp.append(fasta_base)
                        skipped_temp.append("*")
                        window.insert(index, skipped_temp)
                        position += 1
                        index += 1
                else:
                    index += 1 # code block tested
# ----------------------------------------------------------------------------------------------------             
            elif index == 0 and prev_window != None and window[index][0] > (prev_window + 101):
                 # if there are skipped window between the current window and the last window.
                sw_pos = (prev_window + 1) # skipped win pos
                print("BLOCK 3")
                #print(window)
                for skipped_win in fasta_list[(prev_window + 1):window[index][0]]:
                    skipped_temp = []
                    skipped_temp.append(int(sw_pos))
                    skipped_temp.append(0)
                    skipped_temp.append(fasta_base)
                    skipped_temp.append("**")
                    window.insert(index, skipped_temp)
                    sw_pos += 1
                    index += 1
                #print(window)
                                # code block tested
# -----------------------------------------------------------------------------------------------------
            #win_start = (cur_pos - 1), which is the same as str(window[index][0]) - 1))[:-2]
            elif index == 0 and str(window[index][0])[-2:] != "01": # for all elements except the first in the window, the position of the current element is
                                                                    # compared to the position of the last element to identify skipped positions.
                cur_pos = window[index][0]
                cur_end = str(cur_pos)[-2:]
                #print("BLOCK 2", "intended start = ", (cur_pos - int(cur_end)), "cur_pos = ", cur_pos)
                #print(int(cur_end))
                position = ((int(cur_pos) - int(cur_end)) + 1)
                for fasta_base in fasta_list[((cur_pos - int(cur_end)) + 1):cur_pos]:
                    skipped_temp = []
                    skipped_temp.append(int(position))
                    skipped_temp.append(0)
                    skipped_temp.append(fasta_base)
                    skipped_temp.append("*")
                    window.insert(index, skipped_temp)
                    position += 1
                    index += 1 # code block tested
            else:
                index += 1

# TESTED ////////////////////////////////////////////////////////////////////////////////////////////

def net_indels(window, fasta_list, indel_dict, sums=False):

    deletions = []
    unique_deletions = []
    insertions = []
    unique_insertions = []

    for pos_data in window:
        position = pos_data[0]
        if str(position) in indel_dict:  # checking if the current element is present in the indel dictionary.
            indel = indel_dict.get(str(position)) # get the value (a list of the bases/indels) of the key (the position)
            for base in indel: # itterating through the value (base, indel list) for this element in cur_win.
                if "+" in base:
                    insertions.append(base) 
                if "-" in base:
                    deletions.append(base) # for each element that contains one or more indels, the value (list of bases) 
                                           # is itterated, deletions are represented in the read with a "+" in pysam, whilst
                                           # insertions are represented with a "-". Any indels that map to the position
                                           # are appended to deletions e.g. [T+2N] and insertions lists e.g. [A-3GTT]
            insertions_set = set(insertions) 
            deletions_set = set(deletions) # because insertions/deletions lists exist for only one element at a time it is likely the same
                                           # indel will occur several times in the base lists as it may occur in several reads that map to this element's position.
                                           # In order to not count the same indel several times the temp lists are changed to sets which cannot contain element duplicates.
            if sums == True:
                for x in list(insertions_set):
                    unique_insertions.append(int(x[2])) # grab the third character in the string, which is he number of inserted/deleted bases.
                for y in list(deletions_set):
                    unique_deletions.append(int(y[2]))
            
            else:
                for x in list(insertions_set):
                    unique_insertions.append((x)) # grab the whole string, which is the indel
                for y in list(deletions_set):
                    unique_deletions.append((y))

            insertions = []
            deletions = []

    if sums == True:
    # finding the sum of indels for the window 
        if unique_insertions != None or unique_deletions != None:
            # print(unique_insertions, unique_deletions)
            # print(sum(unique_insertions))
            if sum(unique_insertions) > sum(unique_deletions):
                net_indel = (sum(unique_insertions)-sum(unique_deletions)) # value for adjusted mean
            elif sum(unique_deletions) > sum(unique_insertions):
                net_indel = ((sum(unique_deletions) - (2*sum(unique_deletions)))+sum(unique_insertions)) # (sum(unique_deletions) - (2*sum(unique_deletions)), this makes the integer negative
            else:
                net_indel = None

        return net_indel

    elif sums == False:
    # producing a tuple with all the unique indels for each positions in the window.
        all_indels = ()
        all_indels.append(tuple(unique_insertions))
        all_indels.append(tuple(unique_deletions))
        return all_indels

# TESTED /////////////////////////////////////////////////////////////////////////////////////////////

def sum_rd(data, expected_size, net_indel): # we only adjust for indels, not for reference mapped gaps.

    sums = []
    winsize = 0

    for x in data:
        if x[2] != "N" or x[2] != "n": #if the base for the position = "N"
            sums.append(x[1]) # new method which excludes any positions marked N from the calculation, allowing the GC average (here) and sum RD for a window (sum_rd function) to be adjusted.
            winsize += 1   # we do not adjust winsize for any skipped positions that mapped to a base in the refenece fasta as these may be bases that are supressed by tuf so scaling up the rd of the window would make them seem regular.
                            # similarly its important to scale up windows with Ns as these are unmapped positions, that will lower the rd of windows and make finding TUF too dificult.
        
    if net_indel == None:
        return sum(sums)
    else:
        if net_indel != None: # we use winsize to adjust for the indels but not compensating for the gap_size
            adjusted_rd = (round(sum(sums)/(winsize + net_indel))*expected_size) # always divide by winsize as we do not want to compensate for reference mapped gaps, these could represent tuf regions/cores.
            return adjusted_rd

# TESTED ////////////////////////////////////////////////////////////////////////////////////////////

def gc_count(window): #gc_counts takes a list as the itterable/data
    
    gc_count = 0
    at_count = 0
    n_count = 0

    for pos_base in window: # CALCULATING GC CONTENT FOR NORMAL WINDOWS
    #    print(x, x[2])
        if pos_base[2] == "C" or pos_base[2] == "c" or pos_base[2] == "G" or pos_base[2] == "g":
            gc_count += 1
        #    print("gc_count =", gc_count)
        elif pos_base[2] != "N" or pos_base[2] != "n":
            at_count += 1
            # winsize = (winsize - 1) # new method which excludes any positions marked N from the calculation, allowing the GC average (here) and sum RD for a window (sum_rd function) to be adjusted.
    #print(window, "\n", "gc = ", gc_count, "at = ", at_count)
    if gc_count == 0:
        return 0
    elif at_count == 0:
        return 1
    else:
        return (gc_count/(gc_count+at_count))

# TESTED ////////////////////////////////////////////////////////////////////////////////////////////

def get_winsize(data):
    
    adj_winsize = 0
    for x in data:
        if x[2] != "N" or x[2] != "n": #if the base for the position = "N"
            adj_winsize += 1
        #print(x[2])
        else:
            print("N")
            return
    return adj_winsize

def nb_transformation(x, r, decimals): # x = window
    if r == 0:
        nb = 0
        return nb
    else:    
        try:
            nb = round(float(2*(r**float(0.5))*math.log(((float(x)+0.25)/((100*r)-0.5))**0.5+(1+((float(x)+0.25)/((100*float(r))-0.5)))**0.5)), decimals)
            return nb
        except:
            print("no nb")
        #    print("1", 2*(r**float(0.5)))
        #    print("2", math.log(((float(x)+0.25)/((100*r)-0.5))**0.5+(1+((float(x)+0.25)/((100*float(r))-0.5)))**0.5))

    #return round(float(2*(r**float(0.5))*math.log(((float(x)+0.25)/((100*r)-0.5))**0.5+(1+((float(x)+0.25)/((100*float(r))-0.5)))**0.5)), decimals)

def pt_transformation(x, decimals):
    return round(float(2*(((float(x)+0.25)/100)**0.5)), decimals)

def cat_win(window, r, gc=None, sums=None): 

    pt_of_sum = pt_transformation(sums, decimals=2)
    nb_of_sum = nb_transformation(sums, r, decimals=2)
    win_length = get_winsize(window)
    #print(window[-1][0])

    win_position = {(window[-1][0]):{"rd" : sums, "pt" : pt_of_sum,
    "nb" : nb_of_sum, "win_length" : win_length, "gc_count" : gc}}
    
    return win_position

# ///////////////////////////////////////////////////////////////////////////////////////////////////

def find_r(data, r_len):

    nbtemp = []
    r_values = []

    for index, true_win in enumerate(data):
        for position in true_win:
            nbtemp.append(position[1])
            if (position[0]%1000000) == 0: # using i, not x[0] because there are too many gaps, e.g. pos 1Mb = i 900kb
                print("r_value")
                print(nbtemp[:1000000:100000])
                try:
                    r_values.append((statistics.mean(nbtemp)**2)/(statistics.stdev(nbtemp)**2 - statistics.mean(nbtemp))) # cannot use numpy, error with package taking too long to fix.
                    # print("1", statistics.mean(nbtemp)**2)
                    # print("2", (statistics.stdev(nbtemp)**2 - statistics.mean(nbtemp)))
                    nbtemp = []

                except Exception as e:
                    print(e)
                    if sum(nbtemp) == 0:
                        r_values.append(0)
                        print("Exception solved") 
                    # print("1", statistics.mean(nbtemp)**2)
                    # print("2", (statistics.stdev(nbtemp)**2 - statistics.mean(nbtemp)))

    return r_values

# ///////////////////////////////////////////////////////////////////////////////////////////////////

def cat_json(output_filename, input_filenames):
    
    def mangle(s): # String, e.g. infile name.
        return s.strip()[1:-1]

    with open(output_filename, "w") as outfile:
        first = True
        for infile_name in input_filenames:
            with open(infile_name) as infile:
                if first:
                    outfile.write('{')
                    first = False
                else:
                    outfile.write(',')
                #print(infile[1:-1])
                outfile.write(mangle(infile.read()))
        outfile.write('}')