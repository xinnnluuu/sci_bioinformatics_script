###
#    extract labeled reads containing T>C mutation
#    this script is wrote for AMUC-seq
###

import sys
import re

mismatch = sys.argv[1]
sam = sys.argv[2]
outfile = sys.argv[3]

reverse_dict = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C'
}

#get mismatch dict
mismatch_file = open(mismatch,'rt')
mismatch_dict = dict()
for mis_line in mismatch_file:
    mis_line_info = re.split("\t",mis_line)
    mismatch_dict[mis_line_info[0] + "\t" + mis_line_info[1]] = mis_line_info[7]
mismatch_file.close()

#extract read contain mismatch
extract_file = open(outfile,'wt')
sam_file = open(sam,'rt')
CIGAR_pattern = "[0-9]+[MIDN]+"              #not contain S,H 
MD_pattern = "[0-9]+[ATGC\^]"
for sam_line in sam_file:
    if re.match('^@',sam_line):                #skip head line
        extract_file.write(sam_line)
        continue
    if re.match('NM:i:0',sam_line):             #skip no-mismatch line
        continue
    
    linelist = re.split("\t",sam_line)              #get line info
    if re.match('KI|GL|MT',linelist[2]):
        continue
    flag_value = int(linelist[1])               #find strand type
    flag_value_bin = bin(flag_value)
    if len(flag_value_bin) >= 8 and (flag_value_bin[-7:-4] == "110" or flag_value_bin[-8:-4] == "1001"):
        strandswitch = 1
    elif len(flag_value_bin) >= 8 and (flag_value_bin[-7:-4] == "101" or flag_value_bin[-8:-4] == "1010"):
        strandswitch = 0
    else:
        continue
    if not re.match("chr",linelist[2]):
        chr_name = "chr" + linelist[2]
    else:
        chr_name = linelist[2]
    start_pos = int(linelist[3])
    cigar = linelist[5]
    align_seq = linelist[9]
    align_quality = linelist[10]
    
    CIGAR_list = re.findall(CIGAR_pattern,cigar)             #get CIGAR info
    if len(CIGAR_list) > 5:
        continue
    MD_list = []                #get MD info
    for tag in linelist[11:]:
        if re.match('^MD',tag):
            MD_list = re.findall(MD_pattern,tag)

    if not MD_list:             #get mismatch site
        continue
    mismatch_list = []
    mismatch_ref_pos_list = []
    mis_pos = 0

    for mis in MD_list:
        mis_pos += (int(mis[0:-1]) + 1)             #mis postion in read seq
        if mis[-1] == "^":              #deal with del base
            mis_pos += -1
            continue
        if mis[-1] == "N" or align_seq[mis_pos - 1] == "N":
            continue
        if strandswitch == 0:
            mismatch_type = mis[-1] + align_seq[mis_pos - 1]              #ref_site + read_site  same as mismatch_dict
        elif strandswitch == 1:
            mismatch_type = reverse_dict[(mis[-1])] + reverse_dict[align_seq[mis_pos - 1]]
        if mismatch_type != "TC":              #only get TC mismatch
                continue
        mismatch_list.append(mismatch_type)
        if len(CIGAR_list) == 1:                #get mismatch ref pos
            site_pos = start_pos + mis_pos - 1
        elif int(CIGAR_list[0][0:-1]) >= mis_pos - 1:
            site_pos = start_pos + mis_pos - 1
        elif CIGAR_list[1][-1] in {'N','D'}:
            if int(CIGAR_list[0][0:-1]) + int(CIGAR_list[2][0:-1]) >= mis_pos - 1:
                site_pos = start_pos + mis_pos - 1 + int(CIGAR_list[1][0:-1])
            elif CIGAR_list[3][-1] in  {'N','D'}:
                site_pos = start_pos + mis_pos - 1 + int(CIGAR_list[1][0:-1]) + int(CIGAR_list[3][0:-1])
            else:
                site_pos = start_pos + mis_pos - 1 + int(CIGAR_list[1][0:-1]) - int(CIGAR_list[3][0:-1])
        elif CIGAR_list[1][-1] == 'I':
            if int(CIGAR_list[0][0:-1]) + int(CIGAR_list[2][0:-1]) >= mis_pos - 1 :
                site_pos = start_pos + mis_pos - 1 - int(CIGAR_list[1][0:-1])
            elif CIGAR_list[3][-1] in  {'N','D'}:
                site_pos = start_pos + mis_pos - 1 - int(CIGAR_list[1][0:-1]) + int(CIGAR_list[3][0:-1])
            else:   
                site_pos = start_pos + mis_pos - 1 - int(CIGAR_list[1][0:-1]) - int(CIGAR_list[3][0:-1])
        mismatch_ref_pos_list.append(chr_name + "\t" + str(site_pos))              #ref pos format same as mismatch_dict
    outputswitch = 0                #affirm if output read
    i = 0

    while i < len(mismatch_list):
        if mismatch_ref_pos_list[i] in mismatch_dict.keys():
            if re.match(mismatch_list[i],mismatch_dict[mismatch_ref_pos_list[i]]):
                outputswitch = 1
        i += 1
    if outputswitch == 1:
        extract_file.write(sam_line)


sam_file.close()
extract_file.close()
