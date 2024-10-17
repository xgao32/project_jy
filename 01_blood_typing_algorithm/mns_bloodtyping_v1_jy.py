#!/usr/bin/env python
# coding=utf-8
# @Content : MNS blood system typing
# @Date    : 2021-07-10
# @Author  : JIANG YUAN
# @Version : MNS_v1


import os, re, sys, time
import argparse


# Open target file and record information
def txt_to_lines(file_path):
    _f = open(file_path)
    lines = _f.readlines()
    _f.close()
    return(lines)

# setting parameters
all_CLASS = {"ABO":"001", "MNS":"002"}


# type 1: snp/indel analysis based on vcf file
def vcf_maker(dbDir, bedDir, genome, geneRegion_file, vcf_file, bam_file, outputDir, CLASS):
    sample = re.search(r"(\w+)(.deduped.bam)", os.path.basename(bam_file)).group(1)
    # step 1: record gene region information
    with open(geneRegion_file) as fr:
        region = {}
        lines = fr.readlines()
        for line in lines:
            for c in CLASS:
                if c in line:
                    s = line.split("\t")
                    gene = s[1]
                    geneReg = s[6] + ":" + s[2]
                    if (geneReg in region.keys()):
                        region[geneReg].append([int(s[8]), int(s[9])])
                    else:
                        region[geneReg] = [[c],[gene], [int(s[8]),int(s[9])]]
    #print(region)    

    # step 2: record db_bed information, initiate result dict for step 3, create temp bed for step 4
    bed = {}
    result = {}
    temp_bed = {}
    for c in CLASS:
        db_bed = bedDir + all_CLASS[c] + "_" + c + "_bed"
        if c not in bed.keys():
            bed[c] = {}
            result[c] = {}
        temp_bed[c] = outputDir + "temp_bed/" + c + "_temp_bed"
        with open(db_bed) as fb, open(temp_bed[c],"w") as ft: # read bed files of different blood systems
            count = 0
            lines = fb.readlines()
            for line in lines:
                s = line.split("\t")
                result[c][str(count)] = "N"
                if (s[3] != "-"):
                    bed[c][str(count)] = [s[0], s[1], s[2], s[3], s[6].strip()]
                    ft.write("\t".join(s[0:4]) + "\t" + s[6].strip() + "\n")
                count += 1
    
    # step 3: read vcf file and match with bed
    
    with open(vcf_file) as fv, open(outputDir + sample + "_vcf_bloodRecord.txt","w") as fw:
        line = fv.readline()
        while(line):
            if (line.startswith("chr")): # filter annotation line
                s = line.split("\t")
                for r in region.keys():
                    if s[0] in r:
                        start = int(re.search(r"(\w+)(:)(\d+)(-)(\d+)", r).group(3))
                        end = int(re.search(r"(\w+)(:)(\d+)(-)(\d+)", r).group(5))
                        if int(s[1]) in range(start, end): # mutation located in gene
                            c = region[r][0][0] # blood system 
                            match = "F"
                            for count in bed[c].keys():
                                if ((s[0] == bed[c][count][0]) & (s[1] == bed[c][count][1]) & 
                                    (s[3].upper() == bed[c][count][3].upper()) & (bed[c][count][4] in s[4].upper().split(","))): # mutation in reference database
                                    fw.write(line.strip() + "\t" + "match" + "\n")
                                    match = "T"
                                    result[c][count] = "Y"
                                    break
                            if (match == "F"): # mutation not in reference database
                                fw.write(line.strip() + "\t" + "unmatch" + "\n")
            line = fv.readline()                                     

    # step 4: double check with s_v4.out and update result dict
    abpath = sys.path[0]
    s_v4 = abpath + "/s_v4.out"
    genomefai = genome + ".fai"
    fw = open(outputDir + sample + "_vcf_bloodRecord.txt","a") 
    for c in temp_bed.keys():
        cmd = "{s_v4} {bam} {bedfile} {hg} {hgfai} {output}".format(s_v4=s_v4, bam = bam_file, bedfile = temp_bed[c], hg=genome, hgfai=genomefai, output=outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord")
        os.system(cmd) # run s_v4.out
        while(os.path.exists(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt") == "False"):
            time.sleep(0.5)
        with open(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt") as fs:
            for line in fs.readlines(): # update vcf_bloodRecord file for further annotation
                if (line.startswith("chr")):
                    s = line.split("\t")
                    af = re.search(r"(AF=)(.+)", s[3]).group(2)
                    if (af != "0"):
                        chr = re.split(r"\.|\t", line)[0]
                        ori_start = re.split(r"\.|\t", line)[1]
                        ori_end = re.split(r"\.|\t", line)[2]
                        ori_seq = re.split(r"\.|\t", line)[3]
                        sub_seq = re.split(r"\.|\t", line)[4]
                        ad = re.search(r"(AD=)(.+)", s[2]).group(2)
                        dp = re.search(r"(DP=)(.+)", s[1]).group(2)
                        match = "F"
                        for count in bed[c].keys(): # update result dick
                            if ((chr == bed[c][count][0]) & (ori_start == bed[c][count][1]) & (ori_seq.upper() == bed[c][count][3].upper()) & (sub_seq.upper() == bed[c][count][4])):
                                fw.write("\t".join([chr, ori_start, "-", ori_seq, sub_seq, "-", ".", s[3].strip(), "AD:DP"]) + "\t" + ad + ":" + dp + "\t" + "match" + "\n") # write match into vcf format
                                result[c][count] = "Y"
                                match = "T"
                                break
                        if (match == "F"):
                            fw.write("\t".join([chr, ori_start, ".", ori_seq, sub_seq, "-", ".", s[3].strip(), "AD:DP"]) + "\t" + ad + ":" + dp + "\t" + "unmatch" + "\n") # write unmatch into vcf format
        os.remove(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt")
    fw.close()
    

    # step 5: compare with reference databse and find phenotype NNNYYYNNN
    for c in result.keys():
        with open(dbDir + all_CLASS[c] + "_" + c + "_Typer_database_hg38-hg19_v2_jy_NYhg38_NYhg19") as fd, open(outputDir + all_CLASS[c] + "_" + c + "_phenotype_result","w") as fp: 
            fp.write("NY_result:" + "".join(list(result[c].values())) + "\n")
            lines = fd.readlines()
            for line in lines:
                NY_infor = line.split("\t")[-1]
                idx = NY_infor.find("Y")
                matchRecord = ""
                while (idx != -1):
                    if(result[c][str(idx)] == "Y"):
                        matchRecord = matchRecord + "T"
                    else:
                        matchRecord = matchRecord + "F"
                    idx = NY_infor.find("Y", idx +1)
                if(("T" in matchRecord) & ("F" in matchRecord)):
                    fp.write(line.strip()+ "\t" + "part_match\n")
                elif(("T" in matchRecord) & ("F" not in matchRecord)):
                    fp.write(line.strip() + "\t" + "perfect_match\n")










'''
# type 3: cnv analysis
def cnv_maker():

# This is for test
vcf_maker("/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_ref/", "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/hg38_bed/", "/root/projects/projects_jy/WGS_blood/02_reference_database/Genome_refseq/GCF_000001405.25_GRCh37.p13_genomic.fna","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg38_jy.txt","/root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/1.1_snp/QIHUI/IndelApplyVQSR.g.vcf","/root/projects/projects_jy/WGS_blood/02_reference_database/generate_database.py","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/",["MNS"])

'''
# main function
def main(inputDir, outputDir, vcf_file, bam_file, CLASS, SYSTEM):
    # set CLASS: which blood systems are studied
    if (CLASS == "all"): 
        CLASS = list(all_CLASS.keys())
    else:
        classes = CLASS.strip().split(",")
        CLASS = []
        for c in classes:
            if c in all_CLASS.keys():
                CLASS.append(c)
            else:
                print("error -c: {} is not a blood system".format(c))
                exit()
    # set genome, geneRegion_file and bedDir: choose hg19 or hg38 genome for analysis
    if (SYSTEM == "hg38"):
        genome = inputDir + "Genome_refseq/hg38.fa"
        geneRegion_file = inputDir + "db_region/43+7exon_position_hg38_jy.txt"
        bedDir = inputDir + "db_bed/hg38_bed/"
    elif (SYSTEM == "hg19"):
        genome = inputDir + "Genome_refseq/hg19.fa"
        geneRegion_file = inputDir + "db_region/43+7exon_position_hg19_jy.txt"
        bedDir = inputDir +"db_bed/hg19_bed/"
    else:
        print("error -g:please choose hg19 or hg38 as human reference genome")
        exit()
   
   # set dbDir:
    dbDir = inputDir + "db_ref/"
     
    # snp & short indel analysis
    vcf_maker(dbDir, bedDir, genome, geneRegion_file, vcf_file, bam_file, outputDir, CLASS)






parser = argparse.ArgumentParser(description="MNS blood system typing", usage = "python mns_bloodtyping_v1_jy.py")
parser.add_argument("-c", "--classes", help = "blood group system class file", required = True, default = "all")
parser.add_argument("-g", "--genome", help = "human reference genome", required = True, default = "hg38")
parser.add_argument("-v", "--vcf", help = "variant calling file", required = True)
parser.add_argument("-m", "--bam", help = "bam file ", required = True)
parser.add_argument("-i", "--input", help = "input directory containng reference database, gene region file and hg19/hg38 genome",required = True )
parser.add_argument("-o", "--output", help = "output directory", required = True)

args = parser.parse_args()
main(inputDir = args.input, outputDir = args.output, vcf_file = args.vcf, bam_file = args.bam, CLASS = args.classes, SYSTEM = args.genome)



