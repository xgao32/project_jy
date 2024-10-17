#!/usr/bin/env python

# coding=utf-8
import os,re,time

# set system
system = "ABO,P1PK,LU,KEL"

# 第一组样本
path = "/oss/biobc/Result/02.DnaSeqMap_QIHUI/"
files = os.listdir(path)
samples = []

for f in files:
    print(f)
    if (".deduped.bam" not in f):
        continue
    sample = re.search(r"(.*)(.deduped.bam)", f).group(1)
    samples.append(sample)

for sample in samples:
    print(sample)
    cmd = "nohup python bloodtyping_v6.py -c {system} -g hg38 -i /oss/biobc/Result/02.DnaSeqMap_QIHUI/{sample}.deduped.bam -o /root/projects/projects_jy/WGS_blood/04_blood_typing/temp/ -v /root/projects/projects_jy/WGS_blood/04_blood_typing/output/avinput/{sample}.avinput -n /root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43_cnv_result/{sample}.deduped_cnv.cnr &".format(system = system, sample = sample)
    os.system(cmd)


# 第二组样本

path = "/oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL/snp-indel_vcf/"
files = os.listdir(path)
samples = []

for f in files:
    print(f)
    if (".vcf" not in f):
        continue
    sample = re.search(r"(.*)(.vcf)",f).group(1)
    samples.append(sample)
print(samples)

for sample in samples:
    print(sample)
    cmd = "nohup python bloodtyping_v6.py -c {system} -g hg38 -i /oss/biobc/HaploX/HGC20191114003-0001_191130_1202_1205_1206_1208_1211_1212_1214_1220_SHENzhenxueye/deduped_bam/{sample}/deduped.bam -o /root/projects/projects_jy/WGS_blood/04_blood_typing/temp/ -v /root/projects/projects_jy/WGS_blood/04_blood_typing/output/avinput/{sample}.avinput -n /root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43_cnv_result/{sample}.deduped_cnv.cnr &".format(system = system, sample = sample)
    os.system(cmd)
