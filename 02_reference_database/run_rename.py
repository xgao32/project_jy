#!/usr/bin/env python
# coding=utf-8
import os
path = "/root/projects/projects_jy/WGS_blood/02_reference_database/db_withNY/"
files = os.listdir(path)
for file in files:
    if(file[-6:] != "NYhg19" ):
        continue
    filepath = path + file
    outfile = "/root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/db_ref/" + file[1:-18]
    command = "cp {input} {output}".format(input = filepath, output = outfile)
    os.system(command)
