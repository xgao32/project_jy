#!/usr/bin/env python
# coding=utf-8
import os

path = "/oss/biobc/Result/05.GATK/snp-indel_vcf"
files = os.listdir(path)
outf = open("./MNS_dose_effect.txt","w")
outf.write("sample\t144120554\t\t144120555\t\t144120567\t\t143880509\n")
sites = ["144120554", "144120555", "144120567", "143880509"]
for file in files:
    print(file)
    record = {}
    outf.write(file + "\t")
    with open(path + "/" + file) as f:
        line = f.readline()
        while(line):
            s = line.strip().split("\t")
            if ((s[0] == "chr4") and (s[1] == "144120554")):
                record[s[1]] = [s[3] + ">" +s[4], s[9].split(":")[1]]
            elif((s[0] == "chr4") and (s[1] == "144120555")):
                record[s[1]] = [s[3] + ">" +s[4], s[9].split(":")[1]]
            elif((s[0] == "chr4") and (s[1] == "144120567")):
                record[s[1]] = [s[3] + ">" +s[4], s[9].split(":")[1]]
            elif ((s[0] == "chr4") and (s[1] == "143880509")):
                record[s[1]] = [s[3] + ">" +s[4], s[9].split(":")[1]]
            line = f.readline()
        for site in sites:
            if (site in record.keys()):
                outf.write("\t".join(record[site]) + "\t")
            else:
                outf.write("-\t-\t")
        outf.write("\n")


