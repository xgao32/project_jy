#!/usr/bin/env python
# coding=utf-8
# @Content : pre-process bam file, get variant calling file
# @Date    : 2021-09-30
# @Author  : JIANG YUAN
# @Version : v1
# @Require : Snakefile_multi_Capture_wholegene

import os, re, sys, time, datetime
import argparse

parser = argparse.ArgumentParser(description="data pre-processing",
                                 prog="data_preprocess.py",
                                 usage = "python data_preprocess.py -i <input bam file> -o <output directory> -g <genome>")
parser.add_argument("-g", "--genome", help = "human reference genome", required = True)
parser.add_argument("-i", "--inputBam", help = "input the bam file",required = True )
parser.add_argument("-o", "--output", help = "output directory", required = True)
args = parser.parse_args()
genome = args.genome
inputBam = args.inputBam
outputDir = os.path.abspath(args.output) + "/data_preprocess"
samples = re.search(r"(.*).bam", os.path.basename(inputBam)).group(1)
if (not os.path.exists(outputDir)):
    os.mkdir(outputDir)
if (not os.path.exists(outputDir + "/Outputs")):
    os.mkdir(outputDir + "/Outputs")

configfile = outputDir + "/config.yaml"
with open(configfile, "w") as fw:
    fw.write(
        "samples: " + samples + "\n" + \
        "input_file: " + inputBam + "\n" + \
        "outputDir: " + outputDir + "\n" + \
        "genome: " + genome + "\n"
        )
cmd = "cp Snakefile_multi_Capture_wholegen " + outputDir
os.system(cmd)
cmd = "snakemake -c4 -s " + outputDir + "/Snakefile_multi_Capture_wholegene"
os.system(cmd)