#!/usr/bin/env python
# coding=utf-8
# Version:v0
# Date:2021/06/22
# Author: JIANG YUAN
# Function: convert nt_change to genome change.

import os 
import re
from Bio import Seq
# Set basic parameters
complement={"A":"T", "G":"C", "C":"G", "T":"A",
            "a":"T", "g":"C", "c":"G", "t":"A"}

codon_1to3 = {"A":"Ala", "R":"Arg", "N":"Asn", "D":"Asp", "C":"Cys", 
              "E":"Glu", "Q":"Gln", "G":"Gly", "H":"His", "I":"Ile", 
              "L":"Leu", "K":"Lys", "M":"Met", "F":"Phe", "P":"Pro", 
              "S":"Ser", "T":"Thr", "W":"Trp", "Y":"Tyr", "V":"Val"}

# 01- Match cDNA with genome. Get genome location of each site of cDNA.
def cds_genomeLocation(exonFile, hgFile, hg_version, outdir): 
    # Set parameters
    dict = {}
    siteDict = {}
    # Set exon position file field parameters
    exonSystem_field = 0
    exonGene_field = 1
    cds_start_field = 3 
    cds_end_field = 4 
    cds_order_field = 5
    exonStrand_field = 7
    if(hg_version == "hg38_chr"):  # genome use chr* to annotate, such as chr1 
        exonChr_field = 6
    elif(hg_version == "hg38_NC"): # genome use NC_* to annotate, such as NC_000001.11
        exonChr_field = 11
    elif(hg_version == "hg19_chr"):
        exonChr_field = 6
    elif(hg_version == "hg19_NC"):
        exonChr_field = 11
    
    # Create a dict to record cds location of different genes at genome
    with open(exonFile) as f:
        lines = f.readlines()[1:]
        for line in lines:
            s = line.split()
            system = s[exonSystem_field]
            gene = s[exonGene_field]
            cds_start = s[cds_start_field]
            cds_end = s[cds_end_field]
            chr = s[exonChr_field]
            strand = s[exonStrand_field]
            cds_order = s[cds_order_field]
            if(gene not in dict.keys()):#initialize dict
                if(cds_order != "*"):
                    dict[gene] = {cds_order: [cds_start, cds_end]}
                dict[gene]["strand"] = strand
                dict[gene]["chr"] = chr
                siteDict[gene] = {} #initialize siteDict
            else:
                if(cds_order != "*"):
                    dict[gene][cds_order] = [cds_start, cds_end]
         
    
    # Match cDNA and genome location for each site
    for gene in dict.keys():
        totalNum = len(dict[gene]) - 2 #total cds number
        cmd = "samtools faidx {} ".format(hgFile)
        for num in range(1, totalNum+1):
            g_start = dict[gene][str(num)][0]
            g_end = dict[gene][str(num)][1]
            cmd = cmd + dict[gene]["chr"] + ":" + g_start + "-" + g_end + " "
        cmd = cmd + "> {}{}_{}_cds_{}.fasta".format(outdir, system, gene, hg_version)
        os.system(cmd) # get cds sequence from genome file according to its location
        
        with open("{}{}_{}_cds_{}.fasta".format(outdir, system, gene, hg_version)) as _f:
            # remove "\n" in fasta file
            fa_lines = _f.readlines()
            n = -1
            lines = [] # record new line
            for fa_line in fa_lines:
                if(fa_line[0] == ">"):
                    lines.append("")
                    n = n + 1
                else:
                    lines[n] = lines[n] + fa_line.strip()
            
            # gene located on forward strand
            if(dict[gene]["strand"] == "+"): 
                c_start = 0
                count = 0
                for num in range(1, totalNum+1):
                    line = lines[num - 1]
                    g_start = int(dict[gene][str(num)][0]) # cds start site at genome
                    g_end = int(dict[gene][str(num)][1])   # cds end site at genome
                    for g_site in range(g_start, g_end + 1):
                        count = count + 1 
                        c_site = "c." + str(count + c_start)
                        siteDict[gene][c_site] = [complement[line[count - 1]], 
                                                  ("g." + str(g_site) + line[count - 1])] # match cDNA and genome
                    c_start = c_start + count
            
            # gene located on reverse strand
            elif(dict[gene]["strand"] == "-"):
                c_start = 0
                for num in range(1, totalNum+1):
                    count = 0
                    line = lines[num - 1]
                    g_start = int(dict[gene][str(num)][0]) # cds start site at genome
                    g_end = int(dict[gene][str(num)][1])   # cds end site at genome
                    for g_site in range(g_end, g_start - 1, -1):
                        count = count + 1 
                        c_site = "c." + str(count + c_start)
                        siteDict[gene][c_site] = [complement[line[-count]],
                                                  ("g." + str(g_site) + line[-count])]# match cDNA and genome
                    c_start = c_start + count
    
    return(siteDict) 

# 02- Read reference database, compare with genome location for each mutation.
def snp_genomeLocation(databaseFile, hg19_siteDict, hg38_siteDict, outdir):
    #Set reference database field parameters
    refGene_field = 3
    mutType_filed = 4
    ntChange_field = 5
    aaChange_field = 6
    
    with open(databaseFile) as f:
        outf = open((outdir + os.path.basename(databaseFile)),"w")
        lines = f.readlines()
        for line in lines:
            s = line.split("\t")
            refGene = s[refGene_field]
            mutType = s[mutType_filed]
            ntChange = s[ntChange_field].split(";")
            aaChange = s[aaChange_field].split(";")
            print(mutType)
            # 02.1 SNP type(including SNP in intron)
            if(mutType == "SNP"): 
                hg38_position = []
                hg19_position = []
                #ntChange = s[ntChange_field].split(";")
                #aaChange = s[aaChange_field].split(";")
                for m in range(len(ntChange)):# check if nt_change original site is identical to genome
                    if("+" in ntChange[m]): # SNP in intron
                        c_site = "c." + re.findall(r"\d+", ntChange[m].split("+")[0])
                        add = int(re.findall(r"\d+", ntChange[m].split("+")[1]))
                        symbol = ntChange[m].find(">")
                        cmd = 
                    else: # SNP in exon
                        c_site = "c." + re.findall(r"\d+", ntChange[m])[0]
                        symbol = ntChange[m].find(">")
                        if(ntChange[m][symbol - 1] != hg38_siteDict[refGene][c_site][0]):
                            ntChange[m] = ntChange[m][:symbol-1] + hg38_siteDict[refGene][c_site][0] + ntChange[m][symbol:]
                        hg38_position.append(hg38_siteDict[refGene][c_site][1] + ">" + complement[ntChange[m][symbol + 1]])
                        hg19_position.append(hg19_siteDict[refGene][c_site][1] + ">" + complement[ntChange[m][symbol + 1]])
                for aa in range(len(aaChange)):# check if aaChange original site is identical to genome
                    if(aaChange[aa][0:2] == "p."):
                        aa_site = int(re.findall(r"\d+", aaChange[aa])[0])
                        codon = (hg38_siteDict[refGene]["c." + str(3*aa_site - 2)][0] + 
                                 hg38_siteDict[refGene]["c." + str(3*aa_site - 1)][0] + 
                                 hg38_siteDict[refGene]["c." + str(3*aa_site)][0])
                        aaChange[aa] = "p." + codon_1to3[Seq.translate(codon)] + aaChange[aa][5:]
         
                outf.write("\t".join(s[0:5]) + "\t" + 
                           ";".join(ntChange) + "\t" +
                           ";".join(aaChange) + "\t" +
                           "\t".join(s[7:10]) + "\t" +
                           ";".join(hg38_position) + "\t" +
                           ";".join(hg19_position) + "\n")
                
            # 02.2 CNV type
            elif("CNV" in mutType):
                for n in range(len(ntChange)):
                    infor_start = int(nt_change[n].find("["))
                    infor_end =  int(nt_change[n].find("]")) + 1
                    infor = nt_change[n][infor_start, infor_end]
                    refGene = infor[0:int(infor.find(":"))]






def main(hg38_exonFile, hg19_exonFile, hg38_genome, hg38_version, hg19_genome, hg19_version, databaseFile, outdir):
    hg38_siteDict = cds_genomeLocation(hg38_exonFile, hg38_genome, hg38_version, outdir)
    hg19_siteDict = cds_genomeLocation(hg19_exonFile, hg19_genome, hg19_version, outdir)
    snp_genomeLocation(databaseFile, hg19_siteDict, hg38_siteDict, outdir)


# Set data path  
inputdir = "/root/projects/projects_jy/WGS_blood/02_reference_database/"
outdir = "/root/projects/projects_jy/WGS_blood/02_reference_database/temp/"

hg19_exonFile = inputdir + "43+7exon_position_hg19_jy.txt"
hg38_exonFile = inputdir + "43+7exon_position_hg38_jy.txt"

hg19_genome = "/root/projects/projects_jy/WGS_blood/02_reference_database/Genome_refseq/GCF_000001405.25_GRCh37.p13_genomic.fna"
hg19_version= "hg19_NC" # For hg_version, please select one from "hg19_chr", "hg19_NC". chr means genome use chr* to annotate, such as chr1; NC means genome use NC_* to annotate, such as NC_000001.10

hg38_genome = "/root/projects/projects_jy/WGS_blood/02_reference_database/Genome_refseq/GCF_000001405.39_GRCh38.p13_genomic.fna"
hg38_version = "hg38_NC" # For hg_version, please select one from "hg38_chr", "hg38_NC". chr means genome use chr* to annotate, such as chr1; NC means genome use NC_* to annotate, such as NC_000001.11

databaseFile = inputdir + "db_input/MNS_Typer_database_hg38_v4.txt"


# Run main function
main(hg38_exonFile, hg19_exonFile, hg38_genome, hg38_version, hg19_genome, hg19_version, databaseFile, outdir)
]
