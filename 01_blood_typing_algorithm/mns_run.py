#!/usr/bin/env python

# coding=utf-8
import os,re,time

path = "/oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL/snp-indel_vcf/"
files = os.listdir(path)
samples = []

for f in files:
    sample = re.search(r"(.*)(.vcf)", f).group(1)
    samples.append(sample)

print(samples)

for sample in samples[2:]:
    if(sample == "QYW2_1"):
        continue
    cmd = "python /root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/mns_bloodtyping_temp_jy.py -c MNS -g hg38 -s {sample} -i /oss/biobc/HaploX/HGC20191114003-0001_191130_1202_1205_1206_1208_1211_1212_1214_1220_SHENzhenxueye/deduped_bam/{sample}/deduped.bam -o /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/ -v /oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL/snp-indel_vcf/{sample}.vcf -n /root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43-7_cnv_result/{sample}.deduped_cnv.cnr".format(sample = sample)
    os.system(cmd)


