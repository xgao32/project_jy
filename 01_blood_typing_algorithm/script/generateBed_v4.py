#!/usr/bin/env python
# coding=utf-8
# @Function : Generate bed file of reference database
# @Date     : 2021-11-18
# @Author   : JIANG YUAN
# @Version  : v4


import os,re
import argparse


parser = argparse.ArgumentParser(description="This procedure is to generate bed file for blood antigen reference database ")
parser.add_argument("--db", help = "reference database",required = True )
parser.add_argument("-o", "--outputDir", help = "output directory", required = True)

args = parser.parse_args()
db = args.db
outputDir = args.outputDir
print(db)

chromose={}
chromose["ABO"]="chr9"
chromose["MNS"]="chr4"
chromose["P1PK"]="chr22"
chromose["RH"]="chr1"
chromose["LU"]="chr19"
chromose["KEL"]="chr7"
chromose["LE"]="chr19"
chromose["FY"]="chr1"
chromose["JK"]="chr18"
chromose["DI"]="chr17"
chromose["YT"]="chr7"
chromose["XG"]="chrX"
chromose["SC"]="chr1"
chromose["DO"]="chr12"
chromose["CO"]="chr7"
chromose["LW"]="chr19"
chromose["CHRG"]="chr6"
chromose["H"]="chr19"
chromose["XK"]="chrX"
chromose["GE"]="chr2"
chromose["CROM"]="chr1"
chromose["KN"]="chr1"
chromose["IN"]="chr11"
chromose["OK"]="chr19"
chromose["RAPH"]="chr11"
chromose["JMH"]="chr15"
chromose["I"]="chr6"
chromose["GLOB"]="chr3"
chromose["GIL"]="chr9"
chromose["RHAG"]="chr6"
chromose["FORS"]="chr9"
chromose["JR"]="chr4"
chromose["LAN"]="chr2"
chromose["VEL"]="chr1"
chromose["CD59"]="chr11"
chromose["AUG"]="chr6"
chromose["KANNO"]="chr20"
chromose["SID"]="chr17"
chromose["CTL2"]="chr19"
chromose["PEL"]="chr13"
chromose["MAM"]="chr19"
chromose["EMM"]="chr4"

chromose["GPIIIa"]="chr17"
chromose["GPIbα"]="chr17"
chromose["GPIIb"]="chr17"
chromose["GPIa"]="chr17"
chromose["GPIbβ"]="chr17"
chromose["CD109"]="chr17"


def generateBed(db, outputDir, hgType):
    all_pos = {}
    blood = re.search(r"(.*)(_Typer_database_hg38-hg19)(.*)",os.path.basename(db)).group(1) 
    if (hgType == "hg19"):
        hg_position_field = 11
    elif (hgType == "hg38"):
        hg_position_field = 10
    
    if (not(os.path.exists(outputDir + hgType + "_bed/"))):
        os.mkdir(outputDir + hgType + "_bed/")
    bedFile = open((outputDir + hgType + "_bed/" + blood + "_" + hgType + "_bed"),"w")
    db_NY =  open((outputDir + os.path.basename(db) + "_NY" + hgType),"w")

    with open(db) as f:
        # step 1: deduplication
        lines = f.readlines()
        for line in lines[1:]:
            s = line.split("\t")
            pos = s[hg_position_field]
            for p in pos.split(";"):
                if ((p.strip() not in all_pos.keys()) & (p.strip() != "-") & (p.strip() != "")): #del null
                    all_pos[p.strip()] = len(all_pos) 
        # step 2: generate bed
        for p in all_pos.keys():
            chr = chromose[blood[4:]]
            if (("del" in p) & ("ins" in p)): # type 1: delins
                record = re.search(r"g.(\d+)_(\d+)del(.*)ins(\w+)", p)
                ori_start = record.group(1)
                ori_end = record.group(2)
                if (record.group(3) == ""):
                    ori_seq = "*"
                else:
                    ori_seq = record.group(3)
                sub_start = "*"
                sub_end = "*"
                sub_seq = record.group(4)
            elif ("del" in p): # type 2: del
                record = re.search(r"g.(.*)del(.*)", p)
                if ("_" in record.group(1)):
                    ori_start = re.search(r"(\d+)_(\d+)", record.group(1)).group(1)
                    ori_end = re.search(r"(\d+)_(\d+)", record.group(1)).group(2) 
                else:
                    ori_start = ori_end = record.group(1)
                if (record.group(2) == ""):
                    ori_seq = "*"
                else:    
                    ori_seq = record.group(2)
                sub_start = "*"
                sub_end = "*"
                sub_seq = "-"
            elif ("ins" in p): # type 3: ins
                record = re.search(r"g.(\d+)ins(\w+)", p)
                ori_start = ori_end = record.group(1)
                ori_seq = "-"
                sub_start = sub_end = "*"
                sub_seq = record.group(2)
            elif ("dup" in p): # type 4: dup
                record = re.search(r"g.(.*)dup(.*)", p)
                if ("_" in record.group(1)):
                    ori_start = ori_end = re.search(r"(\d+)_(\d+)", record.group(1)).group(2)
                else:
                    ori_start = ori_end = record.group(1)
                ori_seq = "-"
                sub_start = sub_end = "*"
                sub_seq = record.group(2)
            elif (re.search(r"g.(\d+)_(\d+)(.*)>(\d+)_(\d+)(.*)", p)): # type 5: replace
                record = re.search(r"g.(\d+)_(\d+)(.*)>(\d+)_(\d+)(.*)", p)
                ori_start =  record.group(1)
                ori_end = record.group(2)
                sub_start = record.group(4)
                sub_end = record.group(5)
                if (record.group(3) == ""):
                    ori_seq = "*"
                    sub_seq = "*"
                else:
                    ori_seq = record.group(3)
                    sub_seq = record.group(6)
            elif(re.search(r"g.(\d+)(\w)>(\w)", p)): # type 6: snp
                record = re.search(r"g.(\d+)(\w)>(\w)", p)
                ori_start = ori_end = record.group(1)
                ori_seq = record.group(2)
                sub_start = sub_end = "*"
                sub_seq = record.group(3)
            elif(p == "-"):
                continue
            else:
                print("wrong type:please check it" + p)
                continue
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
                if ((p != "-") & (p != "")):
                    line_append[all_pos[p.strip()]] = "Y"
            line_append = "".join(line_append)
            db_NY.write(line.strip() + "\t" + line_append + "\n")
        db_NY.close()  



# Main
generateBed(db, outputDir, "hg38")
generateBed((outputDir + os.path.basename(db) + "_NYhg38"), outputDir, "hg19")
os.remove(outputDir + os.path.basename(db) + "_NYhg38")

