#!/usr/bin/env python
# coding=utf-8

import os, re
from interval import interval

exonfile = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg38_jy.txt"
input = "/root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43-7_cnv_result/"
files = os.listdir(input)
exon = {
    "GYPA":{},
    "GYPB":{},
    "GYPE":{}
}

with open(exonfile) as f:
    lines = f.readlines()[1:]
    for line in lines:
        s = line.split("\t")
        exon[s[1]]["exon" + s[10]] = [int(s[8]), int(s[9])]


for file in files:
    print(file)
    sample = re.search(r"(.*)(.deduped_cnv.cnr)", file).group(1)
    file = input + file
    with open(file) as f, open((sample + "_MNS_cnv.cnr"),"w") as fw:
        for line in f.readlines():
            s = line.strip().split("\t")
            if ((s[3] == "GYPA") or (s[3] == "GYPB") or (s[3] == "GYPE")):
                for n in exon[s[3]].keys():
                    start = int(s[1])
                    end = int(s[2])
                    if ((interval(exon[s[3]][n]) & interval[start, end]) != interval()):
                        fw.write("\t".join(s[0:4]) + "\t" + n + "\t" +"\t".join(s[4:]) + "\n")


