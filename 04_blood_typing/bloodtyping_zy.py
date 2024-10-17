#!/usr/bin/env python
# coding=utf-8
# @Content : blood typing software
# @Date    : 2021-09-26
# @Author  : JIANG YUAN
# @Version : bloodtyping_main_v1
# @Require : preprocess.py, bloodtyping_v6.py

import os, re, sys, time, datetime
import argparse




parser = argparse.ArgumentParser(description="blood typing",
                                 prog="bloodtyping_main_v6.py",
                                 usage = "python bloodtyping_main_v6.py -i <input bam file> -o <output directory> -s <system> -g <genome version>")
parser.add_argument("-s", "--system", help = "blood group system", required = True, default = "all")
parser.add_argument("-g", "--genome", help = "human reference genome", required = True)
parser.add_argument("-i", "--inputBam", help = "input the bam file",required = True )
parser.add_argument("-o", "--output", help = "output directory", required = True)
parser.add_argument("-v", "--vcf", help = "variant calling file")
parser.add_argument("-n", "--cnv", help = "copy number analysis file")
parser.add_argument("-rt", "--runtime", help = "get software run time",action="store_true")
args = parser.parse_args()
input = args.inputBam
outputDir = args.output
system = args.system
genome = args.genome
runtime = args.runtime
vcf = args.vcf
cnv = args.cnv
basedir = sys.path[0]
#---------------------------------------------------------------------
# Get software run time
#---------------------------------------------------------------------
# ts = get_run_time(runtime,0)
# get_run_time(runtime,ts)

#---------------------------------------------------------------------
# Run pre-process pipeline script if there is no vcf file
#---------------------------------------------------------------------
if (vcf == None):
    vcf_cmd = "nohup python {pre-process_path} -i {inputBam} -o {outputDir} -g {genome} &"
    vcf_args = {}
    vcf_args["pre-process_path"] = basedir + "/scripts/data_preprocess.py"
    vcf_args["inputBam"] = input
    vcf_args["outputDir"] = outputDir + "/variant_calling_file"
    vcf_args["genome"] = genome
    vcf_cmd.format(**vcf_args)
    os.system(vcf_cmd)
    ## specify generated vcf file
    vcf = outputDir + "/variant_calling_file"
#---------------------------------------------------------------------
# Run software
#---------------------------------------------------------------------

bloodtyping_cmd = ""
