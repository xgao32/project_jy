#!/usr/bin/env python
# coding=utf-8
# @Content : blood typing software
# @Date    : 2021-11-16
# @Author  : JIANG YUAN
# @Version : all system_v6
# @Require : annovar,s_v4.out

import os, re, sys, time
import argparse
from interval import interval


# setting parameters
all_CLASS = {"ABO":"001", "P1PK":"003", "RH":"004", "LU":"005", \
             "KEL":"006", "LE":"007", "FY":"008", "JK":"009", "DI":"010", \
             "YT":"011", "XG":"012", "SC":"013", "DO":"014", "CO":"015",\
             "LW":"016", "CHRG":"017", "H":"018", "XK":"019", "GE":"020", \
             "CROM":"021", "KN":"022", "IN":"023", "OK":"024", "RAPH":"025",\
             "JMH":"026", "I":"027", "GLOB":"028", "GIL":"029", "RHAG":"030", \
             "FORS":"031", "JR":"032", "LAN":"033", "VEL":"034", "CD59":"035",\
             "AUG":"036", "KANNO":"037", "SID":"038", "CTL2":"039", "PEL":"040", \
             "MAM":"041", "EMM":"042"}

# make directory if not exists
def mkdir(path):
    if (not os.path.exists(path)):
        os.mkdir(path)

# convert blood group genes vcf to avinput
def runavinput(vcf_file, outputDir):
    abpath = sys.path[0]
    sample = re.search(r"(.*).vcf", os.path.basename(vcf_file)).group(1)
    avinput = outputDir + "avinput/" + sample +".avinput"
    cmd = "perl {path}/convert2annovar.pl --format vcf4 {vcf_file} > {avinput}".format(path = abpath + "/supporting_files/supporting_software", vcf_file = vcf_file, avinput = avinput)
    if(not os.path.exists(avinput)):
        os.system(cmd)
        while(os.path.exists(avinput) == "False"):
            time.sleep(0.5)
    return avinput

# mutation target finding
def targetFind(input,bed,genome,target_region,vcf_bloodRecord):
    cmd = "samtools mpileup -r {region} -l {bed} -f {genome} {bam} > {path}".format(region = target_region, bed = bed, genome = genome, bam = input, path = bed + ".mpileup")
    if (not os.path.exists(bed + ".mpileup")):
        os.system(cmd) # run samtools mpileup
    while(os.path.exists(bed + ".mpileup") == "False"):
        time.sleep(0.5)
    target_result = {}
    skipset = {",",".","^","$","]"}
    ntset = {"A","T","C","G","N"}
    with open(bed + ".mpileup", "r") as f, open(vcf_bloodRecord,"a") as fw:
        for line in f.readlines():
            if (not line.startswith("chr")):
                continue
            s = line.split()
            s[4] = s[4].upper()
            if (set(s[4]) <= skipset): # same with reference,skip
                continue
            elif ((set(s[4]) & ntset) != set()): # SNPs/Indels
                while ("+" in s[4]): # insertion
                    loc = s[4].find("+")
                    bp = int(s[4][loc+1])
                    pattern = s[4][loc:(loc+bp+2)]
                    dp = s[4].count(pattern)
                    af = float(dp/int(s[3]))
                    if (af >=0.8):
                        target_result[",".join(s[0:2]) + "," + s[1] + ",-," + pattern[2:]] = "YY"
                        fw.write("\t".join(s[0:2]) + "\t" + s[1] + "\t-\t" + pattern[2:] + "\thom\n")
                    elif (0.2<af<0.8):
                        target_result[",".join(s[0:2]) + "," + s[1] + ",-," + pattern[2:]]= "Y"
                        fw.write("\t".join(s[0:2]) + "\t" + s[1] + "\t-\t" + pattern[2:] + "\thet\n")
                    s[4] = s[4].replace(pattern,"")
                while ("-" in s[4]): # deletion
                    loc = s[4].find("-")
                    bp = int(s[4][loc+1])
                    pattern = s[4][loc:(loc+bp+2)]
                    dp = s[4].count(pattern)
                    af = float(dp/int(s[3]))
                    if (af >=0.8):
                        target_result[s[0] + "," + str(int(s[1])+1) + "," + str(int(s[1])+bp) + \
                                    "," + pattern[2:] + ",-"] = "YY"
                        fw.write(s[0] + "\t" + str(int(s[1])+1) + "\t" + str(int(s[1])+bp) + \
                                    "\t" + pattern[2:] + "\thom\n")
                    elif (0.2<af<0.8):
                        target_result[s[0] + "," + str(int(s[1])+1) + "," + str(int(s[1])+bp) + \
                                    "," + pattern[2:] + ",-"]= "Y"
                        fw.write(s[0] + "\t" + str(int(s[1])+1) + "\t" + str(int(s[1])+bp) + \
                                    "\t" + pattern[2:] + "\thet\n")
                    s[4] = s[4].replace(pattern,"")
                if ((set(s[4]) & ntset) != set()):# SNPs
                    for nt in (set(s[4]) & ntset):
                        dp = s[4].count(nt)
                        af = float(dp/int(s[3]))
                        if (af >=0.8):
                            target_result[",".join(s[0:2]) + "," + ",".join(s[1:3]) + "," + nt] = "YY"
                            fw.write("\t".join(s[0:2]) + "\t" + "\t".join(s[1:3]) + "\t" + nt + "\thom\n")
                        elif (0.2<af<0.8):
                            target_result[",".join(s[0:2]) + "," + ",".join(s[1:3]) + "," + nt] = "Y"
                            fw.write("\t".join(s[0:2]) + "\t" + "\t".join(s[1:3]) + "\t" + nt + "\thet\n")
    return(target_result)

# snp/indel analysis based on vcf file
def vcf_maker(bedDir, genome, geneRegion_file, avinput, input, outputDir, CLASS, SYSTEM):
    sample = re.search(r"(.*)(.avinput)", os.path.basename(avinput)).group(1)
    # step 1: record gene region information
    with open(geneRegion_file) as fr:
        region = {}
        lines = fr.readlines()
        for line in lines:
            for c in CLASS:
                if c in line:
                    s = line.strip().split("\t")
                    gene = s[1].strip()
                    geneRegion = (s[6] + ":" + s[2]).replace(" ","")
                    if (geneRegion in region.keys()):
                        region[geneRegion].append([int(s[7]), int(s[8])])
                    else:
                        region[geneRegion] = [[c],[gene], [int(s[7]),int(s[8])]] 
    # step 2: record db_bed information, initiate result dict for step 3
    bed = {}
    result = {}
    temp_bed = {}
    mkdir(outputDir + "temp_bed/")
    for c in CLASS:
        db_bed = bedDir + all_CLASS[c] + "_" + c + "_" + SYSTEM + "_bed"
        if c not in bed.keys():
            bed[c] = {}
            result[c] = {}
        temp_bed[c] = outputDir + "temp_bed/" + c + "_temp_bed"
        with open(db_bed) as fb: # read bed files of different blood systems
            count = 0
            for line in fb.readlines():
                s = line.strip().split("\t")
                result[c][str(count)] = "N"
                if (s[3] != "*"):
                    bed[c][",".join(s[0:4]) + "," + s[6]] = str(count)
                count += 1
    # step 3: read avinput file and record gene mutation
    mkdir(outputDir + "vcf_bloodRecord")
    vcf_bloodRecord = outputDir + "vcf_bloodRecord/" + sample + "_vcf_bloodRecord.txt"
    with open(avinput) as fv, open(vcf_bloodRecord,"w") as fw:
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
                        for key in bed[c].keys():
                            sk = key.split(",")
                            if ((s[0] == sk[0]) & (s[1] == sk[1]) & 
                                (s[3].upper() == sk[3].upper()) & (sk[4] in s[4].upper().split(";"))): # mutation in reference database
                                fw.write(line.strip() + "\t" + "match" + "\n")
                                match = "T"
                                if (s[5] == "hom"):
                                    result[c][bed[c][key]] = "YY"
                                else:
                                    result[c][bed[c][key]] = "Y"
                                break
                        if (match == "F"): # mutation not in reference database
                            fw.write(line.strip() + "\t" + "unmatch" + "\n")
            line = fv.readline() 
    # step 4: double check with target finding and update result dict
    target_bed = {}
    mkdir(outputDir + "target_bed/") # generate target bed
    for c in bed.keys():
        min = 100000000000
        max = 0
        target_bed[c] = outputDir + "target_bed/" + sample + "_" + c + "_target_bed"    
        with open(target_bed[c], "w") as fw:
            for key in bed[c].keys():
                if (result[c][bed[c][key]] == "N"):
                    sk = key.split(",")
                    if (sk[4] == "-"):
                        pos = int(sk[1]) - 2
                    else:
                        pos = int(sk[1]) - 1
                    if (pos < min):
                        min = pos
                    if (pos > max):
                        max = pos
                    fw.write(sk[0] + "\t" + str(pos) + "\t" + sk[2] + "\n")
        if (min != 0): # target finding
            target_region = sk[0] + ":" + str(min) + "-" + str(max)
            target_result = targetFind(input, target_bed[c], genome, target_region, vcf_bloodRecord)
            for key in target_result.keys(): 
                result[c][bed[c][key]] = target_result[key]
    return(result)

# cnv analysis based on cnv file
def cnv_maker(bedDir, cnv_file, CLASS, SYSTEM):
    # step 1: record db_bed information, initiate cnv result dict for step 2
    bed = {}
    result = {}
    for c in CLASS:
        db_bed = bedDir + all_CLASS[c] + "_" + c + "_" + SYSTEM + "_bed"
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
    '''
    # step 2: MNS read cnv_file
    c = "MNS"
    with open(cnv_file) as fc:
        for line in fc.readlines():
            if ("chromosome" in line):
                continue
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
    '''
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
    homosites = {}
    for c in vcf_result.keys():
        result[c] = {}
        homosites[c] = set()
        for count in vcf_result[c].keys():
            if ((vcf_result[c][count] == "YY") or (cnv_result[c][count] == "YY")):
                result[c][count] = "H"
                homosites[c].add(int(count))
            elif ((vcf_result[c][count] == "Y") or (cnv_result[c][count] == "Y")):
                result[c][count] = "Y"
            else:
                result[c][count] = "N"
    # step 2: compare with reference database and find phenotype NNNYYYNNN
    mkdir(outputDir + "43-7_phenotype") # phenotype output directory
    for c in result.keys():
        mkdir(outputDir + "43-7_phenotype/" + all_CLASS[c] + "_" + c +"_phenotype")
        mutsites = set()
        idx_dict = {}
        with open(dbDir + all_CLASS[c] + "_" + c + "_Typer_database_hg38-hg19") as fd, open(outputDir + "43-7_phenotype/" + all_CLASS[c] + "_" + c + "_phenotype/"+ sample + "_" + c + "_phenotype_result","w") as fp: 
            fp.write("NY_result:" + "".join(list(result[c].values())) + "\n")
            lines = fd.readlines()
            fp.write(lines[0])
            for line in lines[1:]:
                idx_list = []
                NY_infor = line.strip().split("\t")[-1]
                idx = NY_infor.find("Y")
                matchRecord = ""
                while (idx != -1):
                    idx_list.append(idx)
                    if((result[c][str(idx)] == "Y") or (result[c][str(idx)] == "H")):
                        matchRecord = matchRecord + "T"
                    else:
                        matchRecord = matchRecord + "F"
                    idx = NY_infor.find("Y", idx +1)
                if(("T" in matchRecord) & ("F" not in matchRecord)):
                    fp.write(line.strip() + "\t" + "perfect_match\n")
                    idx_dict["(".join(line.strip().split("\t")[1:3]) +")"] = idx_list
                    mutsites = set(idx_list).union(mutsites)
                elif(matchRecord == ""):
                    fp.write(line.strip() + "\t" + "same with hg_reference\n")
                    idx_dict["(".join(line.strip().split("\t")[1:3]) + ")"] = idx_list
        # step 3: Determine diplotype
            if (c != "MNS"):
                fw = open("{outputDir}/{system}_diplotype_all.csv".format(outputDir = outputDir,system = c),"a")
                homoset = homosites[c] & mutsites
                heterset = mutsites.difference(homoset)
                print(homoset)
                print(heterset)
                idx_keylist = list(idx_dict.keys())
                for n in range(len(idx_keylist)):
                    if (homoset <= set(idx_dict[idx_keylist[n]])):
                        for m in range(n, len(idx_keylist)):
                            if (not (homoset <= set(idx_dict[idx_keylist[m]]))):
                                continue
                            interset = set(idx_dict[idx_keylist[n]]) & set(idx_dict[idx_keylist[m]])
                            unionset = set(idx_dict[idx_keylist[n]]).union(set(idx_dict[idx_keylist[m]]))
                            #print(idx_keylist[m] + " " + idx_keylist[n])
                            #print(unionset.difference(interset))
                            if (unionset.difference(interset) == heterset):
                                fp.write("haplotype:" + idx_keylist[n] + "|" + idx_keylist[m] + "\n")
                                fw.write(sample + "," + idx_keylist[n] + "," + idx_keylist[m] + "\n")
    return(result)










#vcf_maker("/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_ref/", "/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_bed/hg38_bed/", "/root/projects/projects_jy/WGS_blood/02_reference_database/Genome_refseq/GCF_000001405.25_GRCh37.p13_genomic.fna","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/input/db_region/43+7exon_position_hg38_jy.txt","/root/projects/projects_jy/WGS_blood/01_blood_typing_algorithm/1.1_snp/QIHUI/IndelApplyVQSR.g.vcf","/root/projects/projects_jy/WGS_blood/02_reference_database/generate_database.py","/root/projects/projects_jy/WGS_blood/03_mns_blood_typing/output/",["MNS"])




#main function
def main(input, outputDir, vcf_file, cnv_file, CLASS, SYSTEM):
    ts = time.time()
    abpath = sys.path[0]
    # set CLASS: which blood systems are studied
    if (CLASS.lower() == "all"): 
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
        genome = abpath + "/supporting_files/Genome_refseq/Homo_sapiens_assembly38.fasta"
        geneRegion_file = abpath + "/supporting_files/db_region/43+7_exon_position_hg38.txt"
        bedDir = abpath +"/supporting_files/db_bed/hg38_bed/"
    elif (SYSTEM == "hg19"):
        genome = abpath + "/supporting_files/Genome_refseq/Homo_sapiens_assembly19.fasta"
        geneRegion_file = abpath + "/supporting_files/db_region/43+7_exon_position_hg19.txt"
        bedDir = abpath +"/supporting_files/db_bed/hg19_bed/"
    else:
        print("error -g:please choose hg19 or hg38 as human reference genome")
        exit()
    
    # set dbDir and outputDir:
    dbDir = abpath + "/supporting_files/db_ref/"
    outputDir = os.path.abspath(outputDir) +"/"
    
    # convert blood group gene vcf to avinput for analysis
    #avinput = runavinput(vcf_file, outputDir)
    avinput = vcf_file
    # snp & short indel analysis
    vcf_result = vcf_maker(bedDir, genome, geneRegion_file, avinput, input, outputDir, CLASS,SYSTEM)
    print("vcf:ok " + str(time.time() -ts))
    # cnv analysis
    cnv_result = cnv_maker(bedDir, cnv_file, CLASS, SYSTEM)
    print("cnv:ok" + str(time.time() -ts))
    # merge cnv and snp/indel result
    sample = re.search(r"(.*)(.avinput)", os.path.basename(avinput)).group(1)
    result = merge(vcf_result, cnv_result, dbDir, outputDir, sample)
    print(sample +" " + str(time.time() -ts))





parser = argparse.ArgumentParser(description="blood typing", 
                                 usage = "python bloodtyping_v4.py -i <input dir> -o <output dir> -v <vcf file> -c <classes>")
parser.add_argument("-c", "--classes", help = "blood group system class file, choose 'all' for all systems", required = True, default = "all")
parser.add_argument("-g", "--genome", help = "human reference genome", required = True, default = "hg38")
parser.add_argument("-i", "--input", help = "input the bam file",required = True )
parser.add_argument("-o", "--output", help = "output directory", required = True)
parser.add_argument("-v", "--vcf", help = "variant calling file", required = True)
parser.add_argument("-n", "--cnv", help = "copy number analysis file", required = True)

args = parser.parse_args()
main(input = args.input, outputDir = args.output, vcf_file = args.vcf, cnv_file = args.cnv, CLASS = args.classes, SYSTEM = args.genome)



