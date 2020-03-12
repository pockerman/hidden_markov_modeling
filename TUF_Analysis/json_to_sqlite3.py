# import pymysql, sqlite3, os, json # tutorial at: https://blog.softhints.com/python-read-validate-and-import-csv-json-file-to-mysql/#pymysql
# # reading json data
# file_path = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/m585_0-40Mb.json') # setting the file path of the json file we are working with
# json_data = open(file_path).read() # reading in json file.
# json_obj = json.loads(json_data) # decoding json file for use in python.
# #print(json_obj.keys())

# #connecting to sqlite and insert json data with pymysql
# con = sqlite3.connect('json_trial') # #con = pymysql.connect(host = 'localhost', user = 'ap718', passwd = '', db = 'json_trial')
# cur = con.cursor()

# cur.execute('''CREATE TABLE json_table(
# POSITION 	INT PRIMARY KEY	NOT NULL,
# READ_DEPTH 	INT		NOT NULL,
# POISSON_TRANSFORMATION 	INT		NOT NULL,
# NEGATIVE_BINOMIAL 	INT		NOT NULL
# )''')

# # parsing the json data into SQL inserts
# for k, v in json_obj.items():
# 	#data = (json_obj[k].get("data", None))
# 	for k2, v2 in json_obj[k].items():
# 		if k2 == "data":
# 			position = k
# 			rd = json_obj[k][k2].get("rd", None) 
# 			pt = (json_obj[k][k2].get("pt", None)) # none = val to be returned if the key is not found.
# 			nb = (json_obj[k][k2].get("nb", None))
# 			if (int(position)%100000) == 0: 
# 				print("Position = ", position, "\n", "read depth = ", rd, "Poisson Transformation = ", pt, 
# 					"Negative Binomial = ", nb)
# 			cur.execute("""INSERT INTO json_table (POSITION, READ_DEPTH, NEGATIVE_BINOMIAL, 
# 				POISSON_TRANSFORMATION) VALUES (?, ?, ?, ?)""", (position, rd, nb, pt))
# con.commit()

# excellent sqlite3 tutorial at: https://www.tutorialspoint.com/sqlite/sqlite_python.htm

#---------------------------------------------------------------------------------------------------

import pymysql, sqlite3, os, json
from tuf_functions3_tidying import cat_json

#set the wd all the json files are located in before concatenating the files.
#store file paths in an itterable object type.
#this saves on specifying longer file paths.

#file_path = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/new_m585_120Mb.json') # setting the file path of the json file we are working with

#concatenating json data
input_files = []
fp_1 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_40Mb.json')
input_files.append(fp_1)
fp_2 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_40-80Mb.json')
input_files.append(fp_2)
fp_3 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_80-120Mb.json')
input_files.append(fp_3)
fp_4 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_120-160Mb.json')
input_files.append(fp_4)
fp_5 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_160-200Mb.json')
input_files.append(fp_5)
fp_6 = os.path.abspath('/scratch/spectre/a/ap718/YoungAdam/fromColin/wf_m585_200Mb-end.json')
input_files.append(fp_6)

#chr1_m585 = open("chr1_m585.json", "w")
cat_json("chr1_m585.json", input_files)

# reading json data

with open("chr1_m585.json", "r") as json_data:
	json_obj = json.load(json_data)
#json_data = open(file_path).read() # reading in json file., g_cou
#json_obj = json.loads(json_data) # decoding json file for use in python.

#print(json_obj.keys())

#connecting to sqlite and insert json data with pymysql
con = sqlite3.connect('wf_TUF_chr1') # #con = pymysql.connect(host = 'localhost', user = 'ap718', passwd = '', db = 'json_trial')
cur = con.cursor()

cur.execute("""CREATE TABLE chr1_m585(
POSITION 	INT PRIMARY KEY NOT NULL,
READ_DEPTH 	INT 	NOT NULL,
POISSON_TRANSFORMATION 	INT 	NOT NULL,
NEGATIVE_BINOMIAL 	INT 	NOT NULL,
WIN_LENGTH 	INT 	NOT NULL,
GC_COUNT 	INT 	NOT NULL
)""")

exceptions = 0

# parsing the json data into SQL inserts
for k, v in json_obj.items():
	position = k
	rd = json_obj[k].get("rd", None) 
	pt = (json_obj[k].get("pt", None)) # none = val to be returned if the key is not found.
	nb = (json_obj[k].get("nb", "N/A"))
	win_length = (json_obj[k].get("win_length", None))
	gc_count = (json_obj[k].get("gc_count", None))
	#if (int(position)%100000) == 0: 
		#print("Position = ", position, "\n", "read depth = ", rd, "Poisson Transformation = ", pt, 
		#	"Negative Binomial = ", nb, "win_length = ", win_length, "gc_count = ", gc_count)
	try:
		cur.execute("""INSERT INTO chr1_m585 (POSITION, READ_DEPTH, POISSON_TRANSFORMATION, 
			NEGATIVE_BINOMIAL, WIN_LENGTH, GC_COUNT) VALUES (?, ?, ?, ?, ?, ?)""", (position, rd, pt,
				nb, win_length, gc_count))
	except Exception as e:
		print(e, "\n", position, rd, pt, nb, win_length, gc_count)
		exceptions += 1
print(exceptions)
con.commit()
