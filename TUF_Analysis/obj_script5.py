import pysam, json, pprint, time #, numpy as np
from collections import Counter, OrderedDict
from TUF_functions2 import bam_strip, is_mapped, indels, gc_count, common_bases, sum_rd
from TUF_functions2 import nb_transformation, poisson_transformation, find_r, get_winsize, cat_win

# sequence file for individual
samfile = pysam.AlignmentFile('m605_verysensitive_trim_sorted.bam',"rb")

# reference chain
fastafile = pysam.FastaFile('GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
fasta_list = list(fastafile.fetch("chr1", 0, )) # pysam.FastaFile not itterable.

start = 100001 # 180039201 #known gap in area
end = 1400001 # 180939401
start_counter = start  # setting the for what should be the start position for every window, 
					   # this is used for checking if gaps occur at the start cur_pos and prev_pos 
					   # will always be 1 apart and so if cur_pos does not equal ....01, the gap from 
					   # 0 to cur_pos can be calculated.

bam_list = bam_strip("chr1", samfile, start, end) # itterating the bam file, pulling out the 
print("bam_list generated") 					  # positions, counts/coverage, and bases.

indel_dict = indels("chr1", samfile, start, end) # itterating the bam file and putting all 
print("indels found") # positions containg at least 1 read with an indel into a dictionary
					  # where the position is the key and the mapped bases to the position are the vals.
common_bases(bam_list, fasta_list) # finding the most common base that is mapped to each given position and
					   # appending that base to bam_list; e.g. [pos, rd, most common base]



chr1_dictionary = OrderedDict()
cur_win = [] # bam_list is itterated and each element is appended to cur_win until the position 
			 # in the element (bam_data[0]), e.g. [POS, rd, base] is divisible by 100. cur_win
			 # can then contain more than one window's worth of positions if cross-win gaps occur.
# window size
winsize = 100

# parameters for generating negative
# binaomial dist
r_len = 1000000 # always match
iter_r_len = 1000000

r_values = find_r(bam_list, r_len)
counter = 0
r = r_values[counter]

print("itterating bam_list...")

for ind, x in enumerate(bam_list):
	# if (ind%1000) == 0:
	# 	print(ind)
	cur_win.append(x)
	unique_insertions = []
	unique_deletions = [] # these two lists contain the lengths (as integers) of all the UNIQUE indels
						  # that are discovered during the itteration of indiviual cur_win states.
	if (x[0]%100) == 0:
		cw_start = 0 # Because all the skipped positions from the fasta reference are added to cur_win
					 # below, all windows will be len 100 until the win get_winsize function is used.
					 # therefore when more than one window is present in cur_win (explained above)
					 # the ranges used for each window can be calculate with the two counters: cw_start, cw_end.
		cw_end = 100 # these counters must be reset before the itteration of a new cur_win.
		win_gaps = [] # list to hold all skipped temps (gaps) that are to be added to this window, 
					  # cannot append the gaps directly to cur_win whilst cur_win is being itterated 
					  # otherwise an infinite loop is created where the inserted gap_data is constantly appended/itterated.
		gap_counter = 0
		go = "no"
		for i, bam_data in enumerate(cur_win):
			if go == "no":
				if len(win_gaps) != 0: # this code block solves the problem of when the gap is appended to cur_win, the index is shifted and so 
					# print("win_gaps is not empty, length = ", len(win_gaps))
					gap_counter += 1 # if the gap was from 103851 (index 50) - 103907 (index 51), the next itteration after win_gap appending and window generation 
					if gap_counter != len(win_gaps) + 1: # would be (index 52), which shows as 103853 at the moment (even thhough 103852 is present in the gap_win updated cur_win data.) 
						go = "yes"
						continue
				else:
					go = "yes"	
			# if cur_win[i][0] >= 103801:
			# 	print("ACTUAL index for next itteration", i)
			# 	print("back up here", cur_win[i][0])
				pass
			if not cur_win: # quick fix for a problem that shouldn't exist
				break # cur_win is being reset to [] at end of itteration (bottom of script)
			          # and is empty up here, however bam_data still holds a value even
				      # though the variable is not used elsewhere, so the loop does not break
				      # as it should when the end of the range is reached.
			insertions = [] 
			deletions = [] # insertions and deletions are temporary lists that exist are reset for
						   # each element in cur_win.
			if i == 0: # checking to see if the position of the first element in cur_win matches the start position of cur_win (start_counter).
				# print(cur_win[i][0], "start counter = ", start_counter)
				# if(cur_win[i][0]) == 103801:
				# 	print("bam_data[0] == 103801", "\n", cur_win)
				if cur_win[i][0] != int(start_counter): # if first pos in win does not equal the start of win (e.g. 105 instead of 100).
					# print("start counter = ", start_counter)
					print("gap at the start of a windw", bam_data)
					print(cur_win)
					exit()

					# this code never executes due to the exit() above
					pos_counter = int(start_counter) # pos_counter = the position of the skipped element, it is used later on when the elements
													 # in win_gaps need to be inserted into cur_win in their correct positions.
					for fasta_base in fasta_list[int(start_counter):bam_data[0]]: # itterating the fasta file to append [pos, rd, base] for any gaps that occur 
																				  # at the start of the window, in the range of the expected win_start (start_counter)
																				  # to the observed win_start (bam_data[0]).
						skipped_temp = [] 
						insert_pos = str(pos_counter)
						insert_pos = insert_pos[-2:]
						skipped_temp.append(pos_counter)
						skipped_temp.append(0) # 0 to represent the read depth for this position.
						skipped_temp.append(fasta_base)
						skipped_temp.append(int(insert_pos))
						win_gaps.append(skipped_temp)
						pos_counter += 1
			else: # for all elements except the first in cur_win, the position of the current element is
				  # compared to the position of the last element to identify skipped position in cur_win.
				cur_pos = cur_win[i][0]	
				prev_pos = cur_win[i-1][0]
				prev_pos = str(prev_pos) # changing prev_pos to a string to slice the last two characters of 
										 # the position and assigning them to a variable (prev_end)
										 # which is used in calculations later on.
				prev_end = prev_pos[-2:]
				prev_end = int(prev_end)

				if str(cur_pos) in indel_dict: # checking if the current element is present in the indel dictionary.
					indel = indel_dict.get(str(cur_pos)) # get the value (a list of the bases/indels) of the key (the position)
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
												   # indel will occur several times in the base lists as it may occur
												   # in several reads that map to this element's position.
												   # In order to not count the same indel several times the temp lists are
												   # changed to sets which cannot contain element duplicates.
					for x in list(insertions_set):
						unique_insertions.append(int(x[2])) # grab the third character in the string, which is he number of inserted/deleted bases.
					for y in list(deletions_set):
						unique_deletions.append(int(y[2]))

				if int(cur_pos) != (int(prev_pos) + 1): # execute code block when skipped positions/gaps occurs
					skipped_pos = (int(prev_pos) + 1) # the position of the skipped element.
					insert_pos = (prev_end + 1) # last two digits of prev_position before the gap, + 1 (first missing position)
												# used to insert skipped_temp into the right position in cur_win later on.
					for fasta_base in fasta_list[(int(prev_pos)+1):(int(cur_pos)+1)]: # itterating the fasta_list to pull bases corresponsing to the skipped positions
																					  # in the bam file and filling the skipped positions in cur_win with [pos, rd, base].
						if int(cur_pos) > ((int(prev_pos) + winsize)-prev_end):
							if len(win_gaps) != 0:
								win_gaps = [] 
						# this code is executed regardless of whether a gap spans several windows or not, allows 
						# the indexing code above to work as single wins retain win_gaps until the end of their cur_win itteration
						# whilst win_gaps is emptied everytime a new cross_win gap has been dealt with.
						skipped_temp = []
						skipped_temp.append(skipped_pos)
						skipped_temp.append(0) # 0 to represent the read depth for this position.
						skipped_temp.append(fasta_base)
						skipped_temp.append(insert_pos)
						win_gaps.append(skipped_temp)
						if skipped_pos == (int(cur_pos)-1): # -1 so that cur_pos does not end up in cur_win twice.
							break
						skipped_pos += 1
						insert_pos += 1
	
					if int(cur_pos) > ((int(prev_pos) + winsize)-prev_end): # once data for the gap has been pulled from the fasta file and appended to win_gaps above
																			# the position of the current element (cur_pos) and the last element (prev_pos) are checked 
																			# to see if the gap spans multiple windows.
						# finding the sum of indels for the window 
						if unique_insertions != None or unique_deletions != None:
							if sum(unique_insertions) > sum(unique_deletions):
								net_indel = (sum(unique_insertions)-sum(unique_deletions)) # value for adjusted mean
							elif sum(unique_deletions) > sum(unique_insertions):
								net_indel = (sum(unique_deletions)-sum(unique_insertions))
							else:
								net_indel = None
							unique_insertions = []
							unique_deletions = []

						# inserting the skipped gap data into the cur_win when the gap spans across more than one window.
						for gap_data in win_gaps:
							cur_win.insert((gap_data[3]-1), gap_data[:3]) # insert_pos = gap_data[3] # weird list indexing??? not 0 indexing on the range function.
						# print(cur_win[cw_start:cw_end])

						#generating data for windows containg gaps that span out of the window.
						win_sum = sum_rd(cur_win[cw_start:(cw_end)], 100, net_indel)
						win_gc = gc_count(cur_win[cw_start:(cw_end)]) # the first win the gap occurs in is represented by cur_win[:-1], which would end at prev_pos.
						window = cat_win(cur_win[cw_start:(cw_end)], r=r, gc=win_gc, sum=win_sum)
						position = {str(window[0]):{"rd" : window[1], "pt": window[2], 
						"nb": window[3], "win_length" : window[4], "gc_count": window[5]}} # the skipped positions are appended to cur_win which allows the window to be generated 
																			     		   # for the window the gap occurs in. if the gap does not skip an entire window (checked for below), e.g. if the distance from the end of the window the gap occurs in
																			     		   # and the end of the gap is >= 100, then the itteration resumes at the position that is the end of the gap + 1 and unless anoter cross_win gap occurs
																			     		   # it will continue onto the if cur_pos == cur_win[-1][0] code block.

						chr1_dictionary.update(position)
						cw_start += 100 # each time a window is produced, cw_start/end needs to +=100, until the itteration
										# of cur_win is over, at which point cw_start/end need to be reset to index the new cur_win.
						cw_end += 100

						# checking for skipped windows between the gap start/end and appending data for skipped windows.
						gap_start_win = int(prev_pos) + (100 - int(prev_pos[-2:])) # end of the window the gap occurs in.
						gap_end = ((cur_pos - gap_start_win) + cw_start) # the difference between the end of the gap, and the end of the 			   						   			 # window the gap occurs in, added onto the current position of the index of cur_win (cw_start).
						if (int(cur_pos) - int(gap_start_win)) >= 100: # if the no. of skipped bases between the end of
							                                           # the first window in the gap and the next recorded position is > 100
							#print("skipped window")
							#print(gap_start_win, cur_pos)
							for skipped_win in cur_win[cw_start:gap_end:100]:
								skip_count = 0
								window = [skipped_win[0], 0, 0, 0, 100, '-'] # adding 'empty' data for skipped windows to the dictionary. Is it correct to put 0 for the pt/nb? is it equivalent to rd in that respect.
								#print(window)
								position = {str(window[0]):{"rd" : window[1], "pt": window[2], 
								"nb": window[3], "win_length" : (window[4], "**"), "gc_count": window[5]}}
								chr1_dictionary.update(position)
								skip_count += 1
								#print("skipped windows ", positions)
							cw_start += (100*skip_count)
							cw_end += (100*skip_count) # when exiting the for loop, cw_start/end need to be adjusted for the number of windows skipped, 
													   # so that any remaining windows in cur_win can be calculated correctly.
							#print(position)
						start_counter += 100


						if i < (iter_r_len + 1):
							r = r_values[counter]
						else:
							iter_r_len = (iter_r_len + r_len)
							if counter < (len(r_values) - 1):
								counter += 1
						# print("index after the window was produced ", i)
						# print("INTENDED index for next itteration", (i + len(win_gaps)))
									
				# code exectues when the cur_pos = the end of cur_win, this executes when no gaps occur in cur_win,
				# when gaps are contained within a single window in cur_win, or after exiting the final cross-win gap in cur_win.
				elif int(cur_pos) == cur_win[-1][0]:

					if unique_insertions != None or unique_deletions != None:
						if sum(unique_insertions) > sum(unique_deletions):
							net_indel = (sum(unique_insertions)-sum(unique_deletions)) # value for adjusted mean
						elif sum(unique_deletions) > sum(unique_insertions):
							net_indel = ((sum(unique_deletions)-sum(unique_insertions)) * -1) # * -1 to convert to negative for deletions.
						else:
							net_indel = None
						unique_insertions = []
						unique_deletions = []
					
					for gap_data in win_gaps:
						cur_win.insert(int(gap_data[3]-1), gap_data[:3]) # insert_pos = gap_data[3] # not 0 indexing on the range function.
						win_gaps = []

					#print(cur_win[cw_start:cw_end])
					win_sum = sum_rd(cur_win[cw_start:cw_end], 100, net_indel)
					win_gc = gc_count(cur_win[cw_start:cw_end]) # the first win the gap occurs in is represented by cur_win[:-1], which would end at prev_pos.
					window = cat_win(cur_win[cw_start:cw_end], r=r, gc=win_gc, sum=win_sum)
					if not win_gaps: # if no gaps were present in range of cur_win (since win_gaps is reset each time a new window is produced.)
						position = {str(window[0]):{"rd" : window[1], "pt": window[2], 
						"nb": window[3], "win_length" : window[4], "gc_count": window[5]}}
					else: # if gaps were present in this range of cur_win, add a marker '*' to the winlength to show 
						  # that positons were skipped and sequence was pulled from fasta.
						position = {str(window[0]):{"rd" : window[1], "pt": window[2], 
						"nb": window[3], "win_length" : (window[4], "*"), "gc_count": window[5]}}
					chr1_dictionary.update(position)
					cw_start += 100
					cw_end += 100
					start_counter += 100
					#print("normal window, may contain gaps ", position)
					cur_win = []

					
					if i < (iter_r_len + 1):
						r = r_values[counter]
					else:
						iter_r_len = (iter_r_len + r_len)
						if counter < (len(r_values) - 1):
							counter += 1
					break


# with open("new_m605_248Mb.json", "w") as data_file:
#     json.dump(chr1_dictionary, data_file)

pprint.pprint(chr1_dictionary)
key_vals = list(chr1_dictionary.items()) # explicitly convert to a list, in case it's Python 3.x
print(key_vals)
print(key_vals[0:10])
#print(key_vals[0]) # get first inserted element 
#print(key_vals[-1]) # get last inserted element 