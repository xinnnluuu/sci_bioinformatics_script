###
#calculate mutation rate from REDItools mutation result file 
###

import re
import sys

infilename = sys.argv[1]
infile = open(infilename,'rt')
count_dict = {
    "AC":0,
    "AG":0,
    "AT":0,
    "CA":0,
    "CG":0,
    "CT":0,
    "GA":0,
    "GC":0,
    "GT":0,
    "TA":0,
    "TC":0,
    "TG":0
}
for line in infile:
    if not re.match("^chr",line):
        continue
    info = re.split("\t",line)
    if info[2] == "N":
        continue
    num_list = re.findall("[0-9]+",info[6])
    if info[2] == "A":
        count_dict["AC"] += int(num_list[1])
        count_dict["AG"] += int(num_list[2])
        count_dict["AT"] += int(num_list[3])
    elif info[2] == "C":
        count_dict["CA"] += int(num_list[0])
        count_dict["CG"] += int(num_list[2])
        count_dict["CT"] += int(num_list[3])
    elif info[2] == "G":
        count_dict["GA"] += int(num_list[0])
        count_dict["GC"] += int(num_list[1])
        count_dict["GT"] += int(num_list[3])
    elif info[2] == "T":
        count_dict["TA"] += int(num_list[0])
        count_dict["TC"] += int(num_list[1])
        count_dict["TG"] += int(num_list[2])
infile.close()
print(count_dict)
