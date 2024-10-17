#!/usr/bin/env python
# coding=utf-8

import re, os

inputDir = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/vcf_bloodRecord/"
files = os.listdir(inputDir)


result = {}
sampleList = []
#name = ""
#file = "./{sample}_vcf_bloodRecord.hg38_multianno.txt".format(sample = sample)
for file in files:
    name = re.search(r"(.*)_vcf_bloodRecord.hg38_multianno.txt",file)
    if (not name):
        continue
    sample = name.group(1)
    sampleList.append(sample)
    with open(file) as f:
        for line in f.readlines():
            s = line.strip().split("\t")
            if (s[0] == "chr"):
                continue
            key = ",".join(s[0:5]) 
            keyanno = s[9] + "," + s[10] + "," + s[12]
            if key not in result.keys():
                result[key] = {"keyanno":keyanno}
            result[key][sample] = s[5]
sampleList.sort()
with open("./vcf_bloodRecord.hg38_multianno_all.txt","w") as fw:
    fw.write("sample" + "\t"*7)
    for sample in sampleList:
        fw.write("\t" + sample)
    fw.write("\n")
    for key in result.keys():
        mut = "\t".join(key.split(","))
        mutanno = "\t".join(result[key]["keyanno"].split(","))
        fw.write(mut + "\t" + mutanno )
        for sample in sampleList:
            if (sample in result[key].keys()):
                fw.write("\t" + result[key][sample])
            else:
                fw.write("\t" + "-")
        fw.write("\n")

