import statistics, math

def bam_strip(Chr, file, start, stop): # file must be SAM/BAM format #op = operation
    errors = 0
    adjusted = 0
    bam_out = [] # Chr, file, op must be type str # Chr must be name of chrom in bam file
    for pileupcolumn in file.pileup(Chr ,start, stop, truncate=True,
                                       ignore_orphans=False):
        temp = []
        temp.append(pileupcolumn.reference_pos)
        temp.append(pileupcolumn.n) # number of reads mapping to this column
        try:
            temp.append(pileupcolumn.get_query_sequences(add_indels=True)) # appending bases
            bam_out.append(temp)
        except:
            try:
                temp.append(pileupcolumn.get_query_sequences(add_indels=False)) # appending bases
                temp.extend("*") # flag to show indels not assessed.
                bam_out.append(temp)
                #print("adjusted base ", "\n", temp)
                if len(temp) == 4:
                    adjusted += 1
            except:
                errors += 1
                #(print(pileupcolumn.reference_pos))
    print("errors = ", errors, "adjusted = ", adjusted)
    return bam_out

def indels(Chr, samfile, start, stop): # find all indels and put in dictionary

    indels = {}
    
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
                
    return indels


def is_mapped(start, stop, fasta): # need to make sure all skpped windows are put into cur_wins.

    n_count = 0
    #loop_control = "on"
    for x in fasta[int(start):int(stop)]: # string.count() may be faster.
        if x == "N":
            n_count +=1
            if n_count >= 0.5*(int(start)-int(stop)): # arbritrarily chosen value, can be altered.
                #loop_control = "off"
                #print(f"referene map for gap {start}:{stop} = False")
                return False #print(f"referene map for gap {start}:{stop} = False")
    #print(f"referene map for gap {start}:{stop} = True")
    return True 

#use the fasta sequence to alter the GC counts: if gap or no reads/bases then use the bases in the fasta file for the skipped positions to calculate the gc content of window.

def gc_count(data): #gc_counts takes a list as the itterable/data
    
    gc_count = 0
    at_count = 0
    n_count = 0

    for win_element in data: # CALCULATING GC CONTENT FOR NORMAL WINDOWS
    #    print(x, x[2])
        if win_element[2] == "C" or win_element[2] == "c" or win_element[2] == "G" or win_element[2] == "g":
            gc_count += 1
        #    print("gc_count =", gc_count)
        elif win_element[2] != "N" or win_element[2] != "n":
            at_count += 1
            # winsize = (winsize - 1) # new method which excludes any positions marked N from the calculation, allowing the GC average (here) and sum RD for a window (sum_rd function) to be adjusted.

    if gc_count == 0:
        return 0
    elif at_count == 0:
        return 1
    else:
        return (gc_count/(gc_count+at_count))

def common_bases(data, fasta):
    
    from collections import Counter

    for i, x in enumerate(data):

        try:
            if x[1] == 1 and len(x) < 3:
                x.append((fasta[x[0] + 1])) # appending the corresponding base for this position from hg38 ref fasta
                a = 1
                print("filling the list")
            elif x[2] == '': # change to searching fasta and appending corresponding base in fasta file instead.
                del(x[2])
                x.insert(2, (fasta[x[0] + 1])) # plus one because of the 0 indexing of fasta list compared to the 1 indexing of the bam positions. 
                #print("adding fasta bases")
                a = 2
            else:
                common_count = Counter(x[2]).most_common(1) # provides the element which reaches the highest occurance first. 
                del(x[2:])
                if len(common_count[0][0]) > 1: # when the most common mapped base to the position is an indel then all elements of the string are appended to the list (see screenshot on windows staff account).
                    indel = common_count[0][0]
                    x.extend(indel[0])
                else:
                    x.extend(common_count[0][0]) # extend adds the new elements into the list, not into the list as a separate list.
                a = 3                            # indexing list to get tuple, indexing tuple to get only first element (base).
        except Exception as e:
            print(e, x)
            print(a)
    return         


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

def find_r(data, r_len):

    nbtemp = []
    r_values = []
    r_num = 1

    for i, x in enumerate(data):
        nbtemp.append(x[1])
        if i == (r_len * r_num): # using i, not x[0] because there are too many gaps, e.g. pos 1Mb = i 900kb
            r_values.append((statistics.mean(nbtemp)**2)/(statistics.stdev(nbtemp)**2 - statistics.mean(nbtemp))) # cannot use numpy, error with package taking too long to fix.
            nbtemp = []
            r_num += 1  
    return r_values


def nb_transformation(x, r, decimals): # x = window
    return round(float(2*(r**float(0.5))*math.log(((float(x)+0.25)/((100*r)-0.5))**0.5+(1+((float(x)+0.25)/((100*float(r))-0.5)))**0.5)), decimals)


def poisson_transformation(x, decimals):
    return round(float(2*(((float(x)+0.25)/100)**0.5)), decimals)

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


def cat_win(data, r, gc=None, sum=None): #concatenate data from bam_list/cur_win into windows, gc and sum = optional parameters.
    # data = cur_win. cw_start/end = cur_win start/end. If gaps contained to one win, cw_start/end need no counter. 
    # when cur_win (data) is input to function, the range of the window is specified is specified, but the fucntion will work without specifying a range.
    window = []

    window.insert(0, data[-1][0]) # the end position of the window.
    if sum != None:
        window.insert(1, (sum))
    window.insert(2, poisson_transformation(window[1], decimals=2))
    window.insert(3, nb_transformation(window[1], r, decimals=2))
    window.insert(4, get_winsize(data)) # using the winsize function that removes any Ns from winsize.
    if gc != None:
        window.insert(5, round(gc, 2))
    
    return window