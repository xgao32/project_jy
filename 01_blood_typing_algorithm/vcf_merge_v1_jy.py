#!/usr/bin/env python
# coding=utf-8
# @Function : merge indel.vcf and snp.vcf
# @Date     : 2021-07-25
# @Author   : JIANG YUAN
# @Version  : v1
# @Required : bcftools, bgzip

import os

inputDir = "/oss/biobc/Result/05.GATK"
outputDir =  "/oss/biobc/Result/05.GATK/snp-indel_vcf"
samples = os.listdir(inputDir)

for sample in samples[1:]:
    path = os.path.join(inputDir, sample)
    if (not os.path.exists(path + "/" + sample + "_indel.vcf.gz")):
        cmd1 = "bgzip {path}/{sample}_indel.vcf".format(path = path, sample = sample)
        os.system(cmd1)
    cmd2 = "bcftools index -t {path}/{sample}_indel.vcf.gz".format(path = path, sample = sample)
    os.system(cmd2)
    
    if (not os.path.exists(path + "/" + sample + "_snp.vcf.gz")):
        cmd3 = "bgzip {path}/{sample}_snp.vcf".format(path = path, sample = sample)
        os.system(cmd3)
    cmd4 = "bcftools index -t {path}/{sample}_snp.vcf.gz".format(path = path, sample = sample)
    os.system(cmd4)

    cmd5 = "bcftools concat -a {path}/{sample}_snp.vcf.gz {path}/{sample}_indel.vcf.gz -o {outputDir}/{sample}.vcf".format(path = path, sample = sample, outputDir = outputDir)
    print("cmd5:" + cmd5 )
    os.system(cmd5)
