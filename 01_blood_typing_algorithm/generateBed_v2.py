#!/usr/bin/env python
# coding=utf-8
# @Function : Generate bed file of reference database
# @Date     : 2021-07-07
# @Author   : JIANG YUAN
# @Version  : v2


import os,re
import argparse


parser = argparse.ArgumentParser(description="This procedure is to generate bed file for blood antigen reference database ")
parser.add_argument("--db", help = "reference database",required = True )
parser.add_argument("-o", "--outputDir", help = "output directory", required = True)

args = parser.parse_args()
db = args.db
outputDir = args.outputDir

# Set parameters
blood = re.search(r"(\w+)(_Typer_database_hg38-hg19_v2_jy)",db).group(1)
hg19_position_field = 10
hg38_position_field = 11

print(blood)
print(db)

chromose={}
chromose["LU"]="chr19"
chromose["002_MNS"]="chr4"
chromose["YT"]="chr7"
chromose["XK"]="chrX"
chromose["VEL"]="chr1"
chromose["SID"]="chr17"
chromose["SC"]="chr1"
chromose["RHAG"]="chr6"
chromose["RAPH"]="chr11"
chromose["P1PK"]="chr22"

chromose["RHD"]="chr1"
chromose["RHCE"]="chr1"

chromose["OK"]="chr19"
chromose["MAM"]="chr19"
chromose["LW"]="chr19"
chromose["LAN"]="chr2"
chromose["KN"]="chr1"
chromose["KEL"]="chr7"
chromose["KANNO"]="chr20"
chromose["JMH"]="chr15"
chromose["JK"]="chr18"
chromose["I"]="chr6"
chromose["IN"]="chr11"
chromose["H"]="chr19"
chromose["GLOB"]="chr3"
chromose["GIL"]="chr9"
chromose["GE"]="chr2"
chromose["FY"]="chr1"
chromose["FROS"]="chr9"
chromose["DO"]="chr12"
chromose["DI"]="chr17"
chromose["CTL2"]="chr19"
chromose["CROM"]="chr1"
chromose["CO"]="chr7"
chromose["CHRG"]="chr6"
chromose["CD59"]="chr11"
chromose["AUG"]="chr6"
chromose["JR"]="chr4"
chromose["LE"]="chr19"

chromose["GPIIIa"]="chr17"
chromose["GPIbα"]="chr17"
chromose["GPIIb"]="chr17"
chromose["GPIa"]="chr17"
chromose["GPIbβ"]="chr17"
chromose["CD109"]="chr17"


def generateBed(db, outputDir, hgType):
    all_pos = {}

    if (hgType == "hg19"):
        hg_position_field = 11
    elif (hgType == "hg38"):
        hg_position_field = 10
    
    if (not(os.path.exists(outputDir + hgType + "_bed/"))):
        os.mkdir(outputDir + hgType + "_bed/")
    bedFile = open((outputDir + hgType + "_bed/" + blood + "_bed"),"w")
    db_NY =  open((db + "_NY" + hgType),"w")

    with open(db) as f:
        # step 1: deduplication
        lines = f.readlines()
        for line in lines[1:]:
            s = line.split("\t")
            pos = s[hg_position_field]
            for p in pos.split(";"):
                if ((p.strip() not in all_pos.keys()) & (p.strip() != "-")): #del null
                    all_pos[p.strip()] = len(all_pos)  
        # step 2: generate bed
        for p in all_pos.keys():
            chr = chromose[blood]
            if (re.findall(r"del", p) != []): # type 1: del
                ori_start = re.search(r"(g.)(\d+)(_)(\d+)(del)", p).group(2)
                ori_end = re.search(r"(g.)(\d+)(_)(\d+)(del)", p).group(4)
                if (p[-3:] == "del"):
                    ori_seq = "*"
                else:    
                    ori_seq = re.search(r"(g.)(\d+)(_)(\d+)(del)(\w+)", p).group(6)
                sub_start = "*"
                sub_end = "*"
                sub_seq = "-"
            elif (re.findall(r"_", p) != []): # type 2: cnv
                ori_start =  re.search(r"(g.)(\d+)(_)(\d+)(.*)(g.)(\d+)(_)(\d+)", p).group(2)
                ori_end = re.search(r"(g.)(\d+)(_)(\d+)(.*)(g.)(\d+)(_)(\d+)", p).group(4)
                sub_start = re.search(r"(g.)(\d+)(_)(\d+)(.*)(g.)(\d+)(_)(\d+)", p).group(7)
                sub_end = re.search(r"(g.)(\d+)(_)(\d+)(.*)(g.)(\d+)(_)(\d+)", p).group(9)
                if (p[-1].isalpha()):
                    ori_seq = re.search(r"(g.)(\d+)(_)(\d+)(\w+)(>)(g.)(\d+)(_)(\d+)(\w+)", p).group(5)
                    sub_seq = re.search(r"(g.)(\d+)(_)(\d+)(\w+)(>)(g.)(\d+)(_)(\d+)(\w+)", p).group(11)
                elif (p[-1].isdigit()):
                    ori_seq = "*"
                    sub_seq = "*"
                else:
                    print("wrong input {}, please check it".format(p))
            else: # type 3: snp
                ori_start = re.search(r"(g.)(\d+)(\w)(>)(\w)", p).group(2)
                ori_end = re.search(r"(g.)(\d+)(\w)(>)(\w)", p).group(2)
                ori_seq = re.search(r"(g.)(\d+)(\w)(>)(\w)", p).group(3)
                sub_start = "*"
                sub_end = "*"
                sub_seq = re.search(r"(g.)(\d+)(\w)(>)(\w)", p).group(5)
            bedFile.write(chr + "\t" + ori_start + "\t" + ori_end + "\t" + ori_seq +
                          "\t" + sub_start + "\t" + sub_end +"\t" + sub_seq +"\n")
        bedFile.close() 
        
    # step 3: generate NY information
        db_NY.write(lines[0].strip() + "\t" + hgType +"_NY_information" + "\n")
        for line in lines[1:]:
            s = line.split("\t")
            pos = s[hg_position_field]
            line_append = ["N"] * (len(all_pos))
            for p in pos.split(";"):
                if (p != "-"):
                    line_append[all_pos[p.strip()]] = "Y"
            line_append = "".join(line_append)
            db_NY.write(line.strip() + "\t" + line_append + "\n")
        db_NY.close()  



# Main
generateBed(db, outputDir, "hg38")
generateBed((db+"_NYhg38"), outputDir, "hg19")
os.remove(db + "_NYhg38")

