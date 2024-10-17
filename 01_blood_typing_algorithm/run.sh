#!/bin/bash

## Run Workflow

# generateBed_v2.py 
#python generateBed_v2.py --db /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_ref/002_MNS_Typer_database_hg38-hg19_v2_jy -o /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/
# generateBed_v4.py
python generateBed_v4.py --db /root/projects/projects_jy/WGS_blood/02_reference_database/db_genomeAnno/039_CTL2_Typer_database_hg38-hg19.txt -o /root/projects/projects_jy/WGS_blood/02_reference_database/db_withNY/

# mns_bloodtyping_v2_jy.py
#python mns_bloodtyping_v2_jy.py -c MNS -g hg38 -v /root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/1.1_snp/QIHUI/IndelApplyVQSR.g.vcf -m /oss/biobc/Result/02.DnaSeqMap_QIHUI/QYW_01.deduped.bam -i /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/ -o /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/ 


# s_v4.out test
#./s_v4.out /oss/biobc/Result/02.DnaSeqMap_QIHUI/QYW_01.deduped.bam /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/temp_bed/MNS_temp_bed /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/Genome_refseq/hg38.fa /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/Genome_refseq/hg38.fa.fai /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/QYW_01_s_v4_bloodRecord.txt

# Annovar
#perl /root/software/api/annovar/convert2annovar.pl --format vcf4 /root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/1.1_snp/QIHUI/IndelApplyVQSR.g.vcf > /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/vcf/QYW_01_vcf.avinput

#perl /root/software/api/annovar/annotate_variation.pl -geneanno -buildver hg38 /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/QYW_01_vcf_bloodRecord /root/software/api/annovar/humandb/ -out /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/annovar

## mns_bloodtyping_v3_jy.py
# sample 1
#python mns_bloodtyping_v3_jy.py -c MNS -g hg38 -i /oss/biobc/Result/02.DnaSeqMap_QIHUI/QYW_16.deduped.bam -o /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/ -v /oss/biobc/Result/05.GATK/snp-indel_vcf/QYW_16.vcf -n /root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43-7_cnv_result/QYW_16.deduped_cnv.cnr 

# sample 2-1
#python mns_bloodtyping_v3_jy.py -c MNS -g hg38 -i /oss/biobc/HaploX/HGC20191114003-0001_191130_1202_1205_1206_1208_1211_1212_1214_1220_SHENzhenxueye/deduped_bam/QYW2_1/QYW2_1.deduped.bam -o /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/ -v /oss/biobc/HaploX/HGC20191114003-0001_20191231/Result/04.Mutation/SNV_INDEL/snp-indel_vcf/QYW2_1.vcf -n /root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.2_cnv/43-7_cnv_result/QYW2_1.deduped_cnv.cnr 
