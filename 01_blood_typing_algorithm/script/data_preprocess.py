#!/usr/bin/env python
# coding=utf-8
# @Content : pre-process bam file, get variant calling file
# @Date    : 2021-09-30
# @Author  : JIANG YUAN
# @Version : v1
# @Require : Snakefile_multi_Capture_wholegene
"""
The Python script data_preprocess.py is designed to pre-process a BAM file and generate a variant calling file. It utilizes the Snakemake workflow management system to execute a series of steps defined in a separate Snakefile.

Here's a breakdown of the script:

Argument Parsing:

The script uses argparse to handle command-line arguments:
-g/--genome: Specifies the human reference genome to be used.
-i/--inputBam: Path to the input BAM file.
-o/--output: Path to the output directory.
Directory and File Setup:

It extracts the sample name from the input BAM file name.
Creates an output directory data_preprocess if it doesn't exist.
Creates a subdirectory Outputs within the output directory.
Generates a config.yaml file containing the sample name, input BAM file path, output directory, and genome information.
Snakemake Integration:

Copies a Snakefile named Snakefile_multi_Capture_wholegene to the output directory. This Snakefile likely contains the rules and steps for the variant calling workflow.
Executes Snakemake using the copied Snakefile:
snakemake -c4 -s <Snakefile_path>: Runs Snakemake with 4 cores (-c4) using the specified Snakefile (-s).
In essence, this script acts as a wrapper to set up the necessary configuration and trigger a Snakemake workflow for variant calling from a BAM file. The actual variant calling process is defined in the Snakefile_multi_Capture_wholegene, which is not shown in the provided code.
"""
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