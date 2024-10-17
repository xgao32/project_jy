#!/usr/bin/env python
# coding=utf-8

import re, os

serofile = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/vcf_bloodRecord/temp/temp.txt"
serology = {}

with open(serofile) as fs:
    for line in fs.readlines()[1:]:
        s = line.split("\t")
        if (s[1] in serology.keys()):
            serology[s[1]].append(s[0])
        else:
            serology[s[1]] = []
            serology[s[1]].append(s[0])

result = {}
name = ""
for p in serology.keys():
    for sample in serology[p]:
        if (name == ""):
            name = name + sample
        else:
            name = name + "-" + sample
        file = "./{sample}_vcf_bloodRecord.hg38_multianno.txt".format(sample = sample)
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

with open(("./temp/"+ name + ".txt"),"w") as fw:
    fw.write("sample" + "\t"*7)
    for p in serology.keys():
        for sample in serology[p]:
            fw.write("\t" + "QYW_" + sample + "-" + p)
    fw.write("\n")

    for key in result.keys():
        mut = "\t".join(key.split(","))
        mutanno = "\t".join(result[key]["keyanno"].split(","))
        fw.write(mut + "\t" + mutanno )
        for p in serology.keys():
            for sample in serology[p]:
                if (sample in result[key].keys()):
                    fw.write("\t" + result[key][sample])
                else:
                    fw.write("\t" + "-")
        fw.write("\n")

