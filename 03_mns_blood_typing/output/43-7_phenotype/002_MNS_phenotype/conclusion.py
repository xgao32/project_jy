#!/usr/bin/env python
# coding=utf-8
import os,re
import pandas as pd
input = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/43-7_phenotype/002_MNS_phenotype/"
phenofiles = os.listdir(input)

dict = {}

for file in phenofiles:
    if ("_MNS_phenotype_result" in file):
        sample = re.search(r"(.*)(_MNS_phenotype_result)", file).group(1)
        dict[sample] = []
        with open(input + file) as f:
            lines = f.readlines()
            NY_record = re.search(r"(NY_result:)(\w+)" , lines[0].strip()).group(2)
            phenotype = ""
            for line in lines[1:]:
                s = line.strip().split("\t")
                if (s[-1] == "perfect_match"):
                    if (phenotype == ""):
                        phenotype = re.search(r"(MNS:)(.*)", s[1]).group(2)
                    else:
                        phenotype = phenotype + ";" + re.search(r"(MNS:)(.*)", s[1]).group(2)
            dict[sample].append(phenotype)
            for n in NY_record:
                dict[sample].append(n)

dataframe = pd.DataFrame(dict)
dataframe.to_csv(r"./MNS_phenotype_result.csv", sep= ",")

