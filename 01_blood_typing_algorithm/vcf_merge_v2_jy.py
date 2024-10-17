#!/usr/bin/env python
# coding=utf-8
# @Function : merge indel.vcf and snp.vcf for QYW2_
# @Date     : 2021-07-25
# @Author   : JIANG YUAN
# @Version  : v2
# @Required : bcftools, bgzip

import os

inputDir = "/oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL"
outputDir =  "/oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL/snp-indel_vcf"
#samples = os.listdir(inputDir)
'''
for sample in samples[1:]:
    path = os.path.abspath(inputDir)
    if (not os.path.exists(path + "/" + sample + ".gz")):
        cmd1 = "bgzip {path}/{sample}".format(path = path, sample = sample)
        os.system(cmd1)
    cmd2 = "bcftools index -t {path}/{sample}.gz".format(path = path, sample = sample)
    os.system(cmd2)
'''
samples = ["1","2","3","4","5","6","7","8","9",
           "14","15","16","17","19",
           "20","21","22","23","24","25","26","27","28","29",
           "30","31","32","33","34","35","36","37","38","39",
           "40","41","42","43","44","45","46","47","48","49",
           "50_1","50_2","51"]
path = os.path.abspath(inputDir)
for sample in samples:
    cmd5 = "bcftools concat -a {path}/QYW2_{sample}_snp.vcf.gz {path}/QYW2_{sample}_indel.vcf.gz -o {outputDir}/QYW2_{sample}.vcf".format(path = path, sample = sample, outputDir = outputDir)
    #print("cmd5:" + cmd5 )
    os.system(cmd5)
