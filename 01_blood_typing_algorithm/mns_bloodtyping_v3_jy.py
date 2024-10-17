#!/usr/bin/env python
# coding=utf-8
# @Content : MNS blood system typing
# @Date    : 2021-07-23
# @Author  : JIANG YUAN
# @Version : MNS_v3
# @Require : annovar,s_v4.out

import os, re, sys, time
import argparse
from interval import interval


# setting parameters
all_CLASS = {"ABO":"001", "MNS":"002", "P1PK":"003", "RHD":"004", "RHCE":"004", "RHD_negative":"004", "WeakD_Del":"004", 
             "LU":"005", "KEL":"006"}

# make directory if not exists
def mkdir(path):
    if (not os.path.exists(path)):
        os.mkdir(path)


# snp/indel analysis based on vcf file
def vcf_maker(bedDir, genome, geneRegion_file, vcf_file, input, outputDir, CLASS):
    sample = re.search(r"(.*)(.deduped.bam)", os.path.basename(input)).group(1)
    # step 1: record gene region information
    with open(geneRegion_file) as fr:
        region = {}
        lines = fr.readlines()
        for line in lines:
            for c in CLASS:
                if c in line:
                    s = line.strip().split("\t")
                    gene = s[1]
                    geneReg = s[6] + ":" + s[2]
                    if (geneReg in region.keys()):
                        region[geneReg].append([int(s[8]), int(s[9])])
                    else:
                        region[geneReg] = [[c],[gene], [int(s[8]),int(s[9])]] 

    # step 2: record db_bed information, initiate result dict for step 3, create temp bed for step 4
    bed = {}
    result = {}
    temp_bed = {}
    mkdir(outputDir + "temp_bed/")
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
                s = line.strip().split("\t")
                result[c][str(count)] = "N"
                if (s[3] != "*"):
                    bed[c][str(count)] = [s[0], s[1], s[2], s[3], s[6]]
                    ft.write("\t".join(s[0:4]) + "\t" + s[6] + "\n")
                count += 1

    # step 3: convert vcf to avinput
    avinput = "/root/projects/projects_jy/WGS_blood/00_CNV_SNV_data/1.1_snp/HaploX/" + sample +".avinput"
    cmd = "perl /root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/software/convert2annovar.pl --format vcf4 {vcf_file} > {avinput}".format(vcf_file = vcf_file, avinput = avinput)
    if(not os.path.exists(avinput)):
        os.system(cmd)
        while(os.path.exists(avinput) == "False"):
            time.sleep(0.5)
    
    # step 4: read avinput file and record gene mutation
    mkdir(outputDir + "vcf_bloodRecord")
    with open(avinput) as fv, open(outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.txt","w") as fw:
        line = fv.readline()
        while(line): 
            s = line.strip().split("\t")
            for r in region.keys():
                if s[0] in r:
                    start = int(re.search(r"(\w+)(:)(\d+)(-)(\d+)", r).group(3))
                    end = int(re.search(r"(\w+)(:)(\d+)(-)(\d+)", r).group(5))
                    if int(s[1]) in range(start, end): # mutation located in gene
                        c = region[r][0][0] # blood system 
                        match = "F"
                        for count in bed[c].keys():
                            if ((s[0] == bed[c][count][0]) & (s[1] == bed[c][count][1]) & 
                                (s[3].upper() == bed[c][count][3].upper()) & (bed[c][count][4] in s[4].upper().split(";"))): # mutation in reference database
                                fw.write("" + line.strip() + "\t" + "match" + "\n")
                                match = "T"
                                result[c][count] = "Y"
                                break
                        if (match == "F"): # mutation not in reference database
                            fw.write(line.strip() + "\t" + "unmatch" + "\n")
            line = fv.readline()                                     

    # step 5: double check with s_v4.out and update result dict
    abpath = sys.path[0]
    s_v4 = abpath + "/s_v4.out"
    genomefai = genome + ".fai"
    fw = open(outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.txt","a") 
    for c in temp_bed.keys():
        cmd = "{s_v4} {bam} {bedfile} {hg} {hgfai} {output}".format(s_v4=s_v4, bam = input, bedfile = temp_bed[c], hg=genome, 
                                                                    hgfai=genomefai, output=outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord")
        os.system(cmd) # run s_v4.out
        while(os.path.exists(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt") == "False"):
            time.sleep(0.5)
        with open(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt") as fs:
            for line in fs.readlines(): # update vcf_bloodRecord file for further annotation
                if (line.startswith("chr")):
                    s = line.strip().split("\t")
                    af = re.search(r"(AF=)(.+)", s[3]).group(2)
                    if (float(af) > 0.2):
                        chr = re.split(r"\.|\t", line)[0]
                        ori_start = re.split(r"\.|\t", line)[1]
                        ori_end = re.split(r"\.|\t", line)[2]
                        ori_seq = re.split(r"\.|\t", line)[3]
                        sub_seq = re.split(r"\.|\t", line)[4]
                        dp = re.search(r"(DP=)(.+)", s[1]).group(2)
                        if (float(af) >= 0.8 ):
                            haplotype = "hom"
                        elif (float(af) < 0.8):
                            haplotype = "het"
                        match = "F"
                        for count in bed[c].keys(): # update result dick
                            if ((chr == bed[c][count][0]) & (ori_start == bed[c][count][1]) & (ori_seq.upper() == bed[c][count][3].upper()) & (sub_seq.upper() == bed[c][count][4])):
                                fw.write("\t".join([chr, ori_start, ori_end, ori_seq, sub_seq, haplotype, "-", dp]) + "\t" + "match" + "\n") # write match into vcf format
                                result[c][count] = "Y"
                                match = "T"
                                break
                        if (match == "F"):
                            fw.write("\t".join([chr, ori_start, ori_end, ori_seq, sub_seq, haplotype, "-", dp]) + "\t" + "unmatch" + "\n") # write unmatch into vcf format
        os.remove(outputDir + all_CLASS[c] + "_" + c + "_s_v4_bloodRecord_tumor_bmd.vcfDedup_mapping.txt")
    fw.close()
    
    # step 6: annotate vcf_bloodRecord by annovar
    hg = re.search(r"(.*)(.fa)", os.path.basename(genome)).group(1)
    cmd = "perl /root/software/api/annovar/table_annovar.pl -buildver {hg} {outputDir}vcf_bloodRecord/{sample}_vcf_bloodRecord.txt /root/software/api/annovar/humandb/ -out {outputDir}vcf_bloodRecord/{sample} -remove -protocol refGene -operation g -nastring .".format(outputDir = outputDir, sample = sample, hg = hg ) 
    os.system(cmd) # run annovar
    while(not os.path.exists(outputDir + "vcf_bloodRecord/" + sample + ".hg38_multianno.txt")):
        time.sleep(1)
    # integrate snp/indel information 
    with open(outputDir + "vcf_bloodRecord/" + sample + ".hg38_multianno.txt") as fa, open(outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.txt") as fv, open(outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.hg38_multianno.txt","w") as fw: 
        fw.write("chr\tStart\tEnd\tRef\tAlt\thaplotype\tscore\tdp\tmatchInfor\tFunc.refGene\trefGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\n")
        fa_lines = fa.readlines()[1:]
        fv_lines = fv.readlines()
        total = len(fa_lines)
        for i in range(total):
            sa = fa_lines[i].strip().split("\t")
            sv = fv_lines[i].strip().split("\t")
            if (sa[0:5] == sv[0:5]):
                fw.write("\t".join(sv) + "\t" + "\t".join(sa[5:]) + "\n")
    os.remove(outputDir + "vcf_bloodRecord/" + sample + ".hg38_multianno.txt")
    os.remove(outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.txt")

    return(result)

# cnv analysis based on cnv file for MNS system
def cnv_maker(bedDir, cnv_file, CLASS):
    # step 1: record db_bed information, initiate cnv result dict for step 2
    bed = {}
    result = {}
    for c in CLASS:
        db_bed = bedDir + all_CLASS[c] + "_" + c + "_bed"
        if c not in bed.keys():
            bed[c] = {}
            result[c] = {}
            with open(db_bed) as fb: # read bed files of different blood systems
                count = 0
                for line in fb.readlines():
                    s = line.strip().split("\t")
                    result[c][str(count)] = [0,0,0,0]
                    if ((s[6] == "-") or (s[4] != "*")):
                        bed[c][str(count)] = [s[0], s[1], s[2], s[4], s[5], s[6]]
                    count += 1
    # step 2: MNS read cnv_file
    c = "MNS"
    with open(cnv_file) as fc:
        for line in fc.readlines():
            s = line.strip().split("\t")
            cnv_chr = s[0]
            cnv_start = int(s[1])
            cnv_end = int(s[2])
            gene = s[3]
            log2 = float(s[5])
            if ((gene == "GYPA") or (gene == "GYPB") or (gene == "GYPE")):
                # step 3: MNS match with bed
                for count in bed[c].keys():
                    bed_chr = bed[c][count][0]
                    ori_start = int(bed[c][count][1])
                    ori_end = int(bed[c][count][2])
                    if (cnv_chr != bed_chr):
                        continue
                    overlap = interval[cnv_start, cnv_end] & interval[ori_start, ori_end]
                    if (overlap != interval()):
                        result[c][count][1] = result[c][count][1] + overlap[0][1] - overlap[0][0]
                        if (log2 < -1):
                            result[c][count][0] = result[c][count][0] + overlap[0][1] - overlap[0][0]
                    if (bed[c][count][5] == "-"): # deletion
                        result[c][count][3] = -1
                        continue
                    sub_start = int(bed[c][count][3])
                    sub_end = int(bed[c][count][4])
                    overlap = interval[cnv_start, cnv_end] & interval[sub_start, sub_end]
                    if (overlap != interval()):
                        result[c][count][3] = result[c][count][3] + overlap[0][1] - overlap[0][0]
                        if (log2 > 0.585):
                            result[c][count][3] = result[c][count][3] + overlap[0][1] - overlap[0][0]
    # return result dict
    for c in CLASS:
        for count in result[c].keys():
            if ((result[c][count][1] == 0) or (result[c][count][1] == 0)):
                result[c][count] = "N"
            elif((result[c][count][3] == -1) and (result[c][count][0]/result[c][count][1] > 0.5)): # deletion
                result[c][count] = "Y"
            elif((result[c][count][0]/result[c][count][1] > 0.5) and (result[c][count][3]/result[c][count][4] > 0.5)): #conversion
                result[c][count] = "Y"
            else:
                result[c][count] = "N"
    return(result)



# merge result of cnv and snp
def merge(vcf_result, cnv_result, dbDir, outputDir, sample):
    # step 1: merge vcf_result and cnv_result
    result = {}
    for c in vcf_result.keys():
        result[c] = {}
        for count in vcf_result[c].keys():
            if ((vcf_result[c][count] == "Y") or (cnv_result[c][count] == "Y")):
                result[c][count] = "Y"
            else:
                result[c][count] = "N"
    
    # step 2: compare with reference database and find phenotype NNNYYYNNN
    mkdir(outputDir + "43-7_phenotype") # phenotype output directory
    for c in result.keys():
        mkdir(outputDir + "43-7_phenotype/" + all_CLASS[c] + "_" + c +"_phenotype")
        with open(dbDir + all_CLASS[c] + "_" + c + "_Typer_database_hg38-hg19_v2_jy_NYhg38_NYhg19") as fd, open(outputDir + "43-7_phenotype/" + all_CLASS[c] + "_" + c + "_phenotype/"+ sample + "_" + c + "_phenotype_result","w") as fp: 
            fp.write("NY_result:" + "".join(list(result[c].values())) + "\n")
            lines = fd.readlines()
            fp.write(lines[0])
            for line in lines[1:]:
                NY_infor = line.strip().split("\t")[-1]
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
    return(result)










#vcf_maker("/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_ref/", "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/hg38_bed/", "/root/projects/projects_jy/WGS_blood/02_reference_database/Genome_refseq/GCF_000001405.25_GRCh37.p13_genomic.fna","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg38_jy.txt","/root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/1.1_snp/QIHUI/IndelApplyVQSR.g.vcf","/root/projects/projects_jy/WGS_blood/02_reference_database/generate_database.py","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/",["MNS"])




#main function
def main(input, outputDir, vcf_file, cnv_file, CLASS, SYSTEM):
    abpath = sys.path[0]
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
        genome = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/Genome_refseq/hg38.fa"
        geneRegion_file = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg38_jy.txt"
        bedDir = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/hg38_bed/"
    elif (SYSTEM == "hg19"):
        genome = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/Genome_refseq/hg19.fa"
        geneRegion_file = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg19_jy.txt"
        bedDir = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/hg19_bed/"
    else:
        print("error -g:please choose hg19 or hg38 as human reference genome")
        exit()
   
   # set dbDir:
    dbDir = "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_ref/"
    outputDir = os.path.abspath(outputDir) +"/"
    # snp & short indel analysis
    vcf_result = vcf_maker(bedDir, genome, geneRegion_file, vcf_file, input, outputDir, CLASS)
    # cnv analysis
    cnv_result = cnv_maker(bedDir, cnv_file, CLASS)
    # merge cnv and snp/indel result
    sample = re.search(r"(.*)(.deduped.bam)", os.path.basename(input)).group(1)
    result = merge(vcf_result, cnv_result, dbDir, outputDir, sample)





parser = argparse.ArgumentParser(description="MNS blood system typing", 
                                 usage = "python mns_bloodtyping_v1_jy.py -i <input dir> -o <output dir> -v <vcf file> -c <classes>")
parser.add_argument("-c", "--classes", help = "blood group system class file", required = True, default = "all")
parser.add_argument("-g", "--genome", help = "human reference genome", required = True, default = "hg38")
parser.add_argument("-i", "--input", help = "input the bam file",required = True )
parser.add_argument("-o", "--output", help = "output directory", required = True)
parser.add_argument("-v", "--vcf", help = "variant calling file", required = True)
parser.add_argument("-n", "--cnv", help = "copy number analysis file", required = True)

args = parser.parse_args()
main(input = args.input, outputDir = args.output, vcf_file = args.vcf, cnv_file = args.cnv, CLASS = args.classes, SYSTEM = args.genome)



