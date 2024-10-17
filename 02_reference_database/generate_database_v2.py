#!/usr/bin/env python
# coding=utf-8
# Version:v2
# Date:2021/09/01
# Author: JIANG YUAN
# Function: convert nt_change to genome change.

import os 
import re
from Bio import Seq
# Set basic parameters
codon_1to3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
              "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
              "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
              "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}

# reverse complement
def complement(sequence):
    trantab = str.maketrans("ACGTacgtRYMKrymkVBHDvbhd","TGCAtgcaYRKMyrkmBVDHbvdh")
    string = sequence[::-1].translate(trantab)
    return string

def searchseq(chr, g_start, g_end):
    cmd = "samtools faidx /root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/Genome_refseq/Homo_sapiens_assembly38.fasta {chr}:{g_start}-{g_end} > temp.fa".format(chr = chr, g_start = g_start, g_end = g_end)
    os.system(cmd)
    with open("temp.fa","r") as f:
        seq = ""
        for line in f.readlines():
            if (line.startswith(">")):
                continue
            else:
                seq = seq + line.strip()
    os.remove("temp.fa")
    return seq


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
    exonStrand_field = 9
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
            s = line.strip().split()
            if (s == []):
                continue
            system = s[exonSystem_field]
            gene = s[exonGene_field]
            cds_start = s[cds_start_field]
            cds_end = s[cds_end_field]
            chr = s[exonChr_field]
            strand = s[exonStrand_field]
            cds_order = s[cds_order_field]
            if(gene not in dict.keys()):#initialize dict
                dict[gene] = {}
                if(cds_order != "*"):
                    dict[gene] = {cds_order: [cds_start, cds_end]}
                dict[gene]["strand"] = strand
                dict[gene]["chr"] = chr
                siteDict[gene] = {} #initialize siteDict
            else:
                if(cds_order != "*"):
                    dict[gene][cds_order] = [cds_start, cds_end]
    print(dict.keys())

    # Match cDNA and genome location for each site
    for gene in dict.keys():
        totalNum = len(dict[gene]) - 2 #total cds number
        cmd = "samtools faidx {} ".format(hgFile)
        for num in range(1, totalNum+1):
            g_start = dict[gene][str(num)][0]
            g_end = dict[gene][str(num)][1]
            cmd = cmd + dict[gene]["chr"] + ":" + g_start + "-" + g_end + " "
        cmd = cmd + "> {}{}_cds_{}.fasta".format(outdir, gene, hg_version)
        os.system(cmd) # get cds sequence from genome file according to its location
        with open("{}{}_cds_{}.fasta".format(outdir, gene, hg_version)) as _f:
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
                siteDict[gene]["strand"] = "+"
                c_start = 0
                for num in range(1, totalNum+1):
                    count = 0
                    line = lines[num - 1]
                    g_start = int(dict[gene][str(num)][0]) # cds start site at genome
                    g_end = int(dict[gene][str(num)][1])   # cds end site at genome
                    for g_site in range(g_start, g_end + 1):
                        count = count + 1
                        c_site = str(count + c_start)
                        siteDict[gene][c_site] = [line[count - 1], str(g_site), line[count - 1]] # match cDNA and genome
                    c_start = c_start + count
            # gene located on reverse strand
            elif(dict[gene]["strand"] == "-"):
                siteDict[gene]["strand"] = "-"
                c_start = 0
                for num in range(1, totalNum+1):
                    count = 0
                    line = lines[num - 1]
                    g_start = int(dict[gene][str(num)][0]) # cds start site at genome
                    g_end = int(dict[gene][str(num)][1])   # cds end site at genome
                    for g_site in range(g_end, g_start - 1, -1):
                        count = count + 1 
                        c_site = str(count + c_start)
                        siteDict[gene][c_site] = [complement(line[-count]),
                                                  str(g_site), line[-count]]# match cDNA and genome
                    c_start = c_start + count
        os.remove("{}{}_cds_{}.fasta".format(outdir, gene, hg_version))
    return(siteDict)

# 02- Read reference database, compare with genome location for each mutation.
def genomeLocation(databaseFile, hg19_siteDict, hg38_siteDict, outdir, chr):
    #Set reference database field parameters
    refGene_field = 3
    mutType_filed = 4
    ntChange_field = 5
    aaChange_field = 6
    with open(databaseFile) as f, open((outdir + os.path.basename(databaseFile) + "_error"),"w") as fw, open((outdir + os.path.basename(databaseFile)),"w") as outf:
        outf.write("Blood_system\tphenotype\tallel_name\tGENE\tType\tnt_change\tAA\texon\tprevalence\tdatabase\thg38_position\thg19_position\tcomments\n")
        lines = f.readlines()[1:]
        n=1
        for line in lines:
            if (line.strip().split() == []):
                continue
            print(str(n) + line)
            n=n+1
            s = line.strip().split("\t")
            refGene = s[refGene_field].split(",")[0]
            mutType = s[mutType_filed]
            ntChange = s[ntChange_field].split(";")
            aaChange = s[aaChange_field]
            hg38_position = []
            hg19_position = []
            if (ntChange[0] == "-"): #same with genome-reference sequence
                outf.write(line)
                fw.write(line)
                continue
            # start analyze different types
            for m in range(len(ntChange)):
                if(ntChange[m].split() == []):
                    continue
                left = ntChange[m].find("[")
                right = ntChange[m].find("]")
                # 1. snp type
                if ((left == -1) & (right == -1)):
                    result = re.findall(r"\d+", ntChange[m])
                    c_site = result[0]
                    symbol = ntChange[m].find(">")
                    if((ntChange[m][symbol - 1] != hg38_siteDict[refGene][c_site][0]) & (len(result) == 1)):
                        ntChange[m] = ntChange[m][:symbol-1] + hg38_siteDict[refGene][c_site][0] + ntChange[m][symbol:]# correct nt if it is not equal to genome corresponding site
                    if (hg38_siteDict[refGene]["strand"] == "-"):
                        if ("-" in ntChange[m]):
                            hg19g_site =  int(hg19_siteDict[refGene][c_site][1]) + int(result[1])
                            hg38g_site = int(hg38_siteDict[refGene][c_site][1]) + int(result[1])
                            hg19nt = hg38nt = complement(ntChange[m][symbol - 1])
                        elif ("+" in ntChange[m]):
                            hg19g_site = int(hg19_siteDict[refGene][c_site][1]) - int(result[1])
                            hg38g_site = int(hg38_siteDict[refGene][c_site][1]) - int(result[1])
                            hg19nt = hg38nt = complement(ntChange[m][symbol - 1])
                        else:
                            hg19nt = hg19_siteDict[refGene][c_site][2]
                            hg38nt = hg38_siteDict[refGene][c_site][2]
                            hg19g_site = hg19_siteDict[refGene][c_site][1]
                            hg38g_site = hg38_siteDict[refGene][c_site][1]
                        hg38_position.append("g." + str(hg38g_site) + hg38nt + ">" + complement(ntChange[m][symbol + 1]))
                        hg19_position.append("g." + str(hg19g_site) + hg19nt + ">" + complement(ntChange[m][symbol + 1]))
                    elif(hg38_siteDict[refGene]["strand"] == "+"):
                        if ("-" in ntChange[m]):
                            hg19g_site =  int(hg19_siteDict[refGene][c_site][1]) - int(result[1])
                            hg38g_site = int(hg38_siteDict[refGene][c_site][1]) - int(result[1])
                            hg19nt = hg38nt = ntChange[m][symbol - 1]
                        elif ("+" in ntChange[m]):
                            hg19g_site = int(hg19_siteDict[refGene][c_site][1]) + int(result[1])
                            hg38g_site = int(hg38_siteDict[refGene][c_site][1]) + int(result[1])
                            hg19nt = hg38nt = ntChange[m][symbol - 1]
                        else:
                            hg19nt = hg19_siteDict[refGene][c_site][2]
                            hg38nt = hg38_siteDict[refGene][c_site][2]
                            hg19g_site = hg19_siteDict[refGene][c_site][1]
                            hg38g_site = hg38_siteDict[refGene][c_site][1]
                        hg38_position.append("g." + str(hg38g_site) + hg38nt + ">" + ntChange[m][symbol + 1])
                        hg19_position.append("g." + str(hg19g_site) + hg19nt + ">" + ntChange[m][symbol + 1])
                # 2. cnv type
                elif ((left != -1) & (right != -1 )):
                    # genome location
                    mut = {"hg19g_site":[], "hg38g_site":[]}
                    infor = ntChange[m][left+1:right]
                    for half in infor.split(">"):
                        refGene = re.search(r"(.*):(.*)", half).group(1)
                        if(":" in half):
                            symbol = half.find(":")
                            half = half[symbol+1:]
                        for h in half.split("_"):
                            result = re.findall(r"\d+", h)
                            c_site = result[0]
                            if (((hg38_siteDict[refGene]["strand"] == "-") & ("-" in h)) or ((hg38_siteDict[refGene]["strand"] == "+") & ("+" in h))):
                                hg19g_site = int(hg19_siteDict[refGene][c_site][1]) + int(result[1])
                                hg38g_site = int(hg38_siteDict[refGene][c_site][1]) + int(result[1])
                            elif (((hg38_siteDict[refGene]["strand"] == "-") & ("+" in h)) or ((hg38_siteDict[refGene]["strand"] == "+") & ("-" in h))):
                                hg19g_site = int(hg19_siteDict[refGene][c_site][1]) - int(result[1])
                                hg38g_site = int(hg38_siteDict[refGene][c_site][1]) - int(result[1])
                            else:
                                hg19g_site = hg19_siteDict[refGene][c_site][1]
                                hg38g_site = hg38_siteDict[refGene][c_site][1]
                            mut["hg19g_site"].append(str(hg19g_site))
                            mut["hg38g_site"].append(str(hg38g_site))

                    # 2.1 delins type
                    if ("delins" in infor):
                        result = re.search(r"(.*):c.(.*)_(.*)delins(\w+)",infor)
                        refGene = result.group(1)
                        g_start = mut["hg38g_site"][0]
                        g_end = mut["hg38g_site"][1]
                        insseq = result.group(4)
                        delseq = ""
                        if (hg38_siteDict[refGene]["strand"] == "-"):  # gene locate on reverse strand
                            if (int(g_start) <= int(g_end) + 15):
                                delseq = searchseq(chr,g_end, g_start)
                            hg38_position.append("g." + g_end + "_" + g_start + "del" + delseq + "ins" + complement(insseq))
                            hg19_position.append("g." + mut["hg19g_site"][1] + "_" +  mut["hg38g_site"][0] +
                                                 "del" + delseq + "ins" + complement(insseq))
                        elif (hg38_siteDict[refGene]["strand"] == "+"):  # gene locate on forward strand
                            if (int(g_end) <= int(g_start) + 15):
                                delseq = searchseq(chr, g_start, g_end)
                            hg38_position.append("g." + g_start + "_" + g_end + "del" + delseq + "ins" + insseq)
                            hg19_position.append("g." + mut["hg19g_site"][0] + "_" + mut["hg19g_site"][1] + "del" + delseq + "ins" + insseq)
                        else:
                            print("error delins: the gene locate on neither forward or reverse strand")
                    # 3. del & dup type
                    elif (("del" in infor) or ("dup" in infor)):
                        if ("del" in infor):
                            type = "del"
                        else:
                            type = "dup"
                        delseq = ""
                        if ("_" not in infor): # ex: c.237del, c.238-1del,c.237dup, c.238-1dup
                            if (("+" in infor) or ("-" in infor)):
                                delseq = searchseq(chr,  mut["hg38g_site"][0], mut["hg38g_site"][0])
                                hg38_position.append("g." + mut["hg38g_site"][0] + type + delseq)
                                hg19_position.append("g." + mut["hg19g_site"][0] + type + delseq)
                            else:
                                refGene = re.search(r"(.*):c.(\d+){}(.*)".format(type),infor).group(1)
                                c_site = re.search(r"(.*):c.(\d+){}(.*)".format(type),infor).group(2)
                                hg38_position.append("g." + hg38_siteDict[refGene][c_site][1] + type + hg38_siteDict[refGene][c_site][2])
                                hg19_position.append("g." + hg19_siteDict[refGene][c_site][1] + type + hg19_siteDict[refGene][c_site][2])
                        else: # ex: c.237_239del, c.237+1_238-1del
                            result = re.search(r"(.*):c.(.*)_(.*){}(.*)".format(type), infor)
                            refGene = result.group(1)
                            g_start = mut["hg38g_site"][0]
                            g_end = mut["hg38g_site"][1]
                            delseq = ""
                            if (hg38_siteDict[refGene]["strand"] == "-"):  # gene locate on reverse strand
                                #print("g_start:"+g_start)
                                if (int(g_start) <= int(g_end) + 15):
                                    delseq = searchseq(chr, g_end, g_start)
                                hg38_position.append("g." + g_end + "_" + g_start + type + delseq)
                                hg19_position.append("g." + mut["hg19g_site"][1] + "_" + mut["hg19g_site"][0] + type + delseq)
                            elif (hg38_siteDict[refGene]["strand"] == "+"):  # gene locate on forward strand
                                if (int(g_end) <= int(g_start) + 15):
                                    delseq = searchseq(chr, g_start, g_end)
                                hg38_position.append("g." + g_start + "_" + g_end + type + delseq)
                                hg19_position.append("g." + mut["hg19g_site"][0] + "_" + mut["hg19g_site"][1] + type + delseq)
                            else:
                                print("error del/dup: the gene locate on neither forward or reverse strand")

                    # 4. ins type
                    elif ("ins" in infor):
                        result = re.search(r"(.*):c.(.*)_(.*)ins(\w+)", infor)
                        refGene = result.group(1)
                        insseq = result.group(4)
                        if (hg38_siteDict[refGene]["strand"] == "-"):  # gene locate on reverse strand
                            hg38_position.append("g." + mut["hg38g_site"][0] + "ins" + complement(insseq))
                            hg19_position.append("g." + mut["hg19g_site"][0] + "ins" + complement(insseq))
                        elif (hg38_siteDict[refGene]["strand"] == "+"):  # gene locate on forward strand
                            hg38_position.append("g." + mut["hg38g_site"][1] + "ins" + insseq)
                            hg19_position.append("g." + mut["hg19g_site"][1] + "ins" + insseq)
                        else:
                            print("error ins: the gene locate on neither forward or reverse strand")
                    # 5. replace type
                    elif(">" in infor):
                        delseq = ""
                        symbol = infor.find(">")
                        if (len(re.findall(":",infor)) == 2):
                            refGene = re.search(r"(.*):(.*)>(.*):(.*)", infor).group(1)
                            refGene1 = re.search(r"(.*):(.*)>(.*):(.*)", infor).group(3)
                            if ("_" not in infor): # ex: GYPA:c.203>GYPB:Ψ.64
                                if (("+" in infor) or ("-" in infor)):
                                    delseq = searchseq(chr, mut["hg38g_site"][0], mut["hg38g_site"][0])
                                    insseq = searchseq(chr, mut["hg38g_site"][1], mut["hg38g_site"][1])
                                    hg38_position.append("g." + mut["hg38g_site"][0] + delseq + ">g." + mut["hg38g_site"][1] + insseq)
                                    hg19_position.append("g." + mut["hg19g_site"][0] + delseq + ">g." + mut["hg19g_site"][1] + insseq)
                                else:
                                    c_site = re.search(r"(\D+)(\d+)(\D+)(\d+)", infor).group(2)
                                    c_site1 = re.search(r"(\D+)(\d+)(\D+)(\d+)", infor).group(4)
                                    hg38_position.append("g." + "".join(hg38_siteDict[refGene][c_site][1:2]) + ">g." + "".join(hg38_siteDict[refGene1][c_site1][1:2]))
                                    hg19_position.append("g." + "".join(hg19_siteDict[refGene][c_site][1:2]) + ">g." + "".join(hg19_siteDict[refGene1][c_site1][1:2]))
                            elif ("_" in infor): # ex:GYPA:c.272_453>GYPB:c.176_276
                                gene = [refGene,refGene1]
                                hg19_record = []
                                hg38_record = []
                                for i in range(2):
                                    g_start = mut["hg38g_site"][2*i+0]
                                    g_end = mut["hg38g_site"][2*i+1]
                                    delseq = ""
                                    if (hg38_siteDict[gene[i]]["strand"] == "-"):  # gene locate on reverse strand
                                        if (int(g_start) <= int(g_end) + 15):
                                            delseq = searchseq(chr, g_end, g_start)
                                        hg19_record.append(mut["hg19g_site"][2*i+1] + "_" +mut["hg19g_site"][2*i+0] + delseq)
                                        hg38_record.append(g_end + "_" + g_start + delseq)
                                    elif (hg38_siteDict[gene[i]]["strand"] == "+"):  # gene locate on forward strand
                                        if (int(g_end) <= int(g_start) + 15):
                                            delseq = searchseq(chr, g_start, g_end)
                                        hg19_record.append(mut["hg19g_site"][2 * i + 0] + "_" + mut["hg19g_site"][2 * i + 1] + delseq)
                                        hg38_record.append(g_start + "_" + g_end + delseq)
                                    else:
                                        print("error del/dup: the gene locate on neither forward or reverse strand")
                                hg19_position.append("g."+ ">".join(hg19_record))
                                hg38_position.append("g." + ">".join(hg38_record))
                        else:
                            fw.write(line)
                            continue
                    else: # write lines that belongs to none of the types
                        fw.write(line)
                        continue
            # AA information
            '''
            for aa in range(len(aaChange)):
                if(aaChange[aa][0:2] == "p."):
                    aa_site = int(re.findall(r"\d+", aaChange[aa])[0])
                    codon = (hg38_siteDict[refGene][str(3*aa_site - 2)][0] + 
                             hg38_siteDict[refGene][str(3*aa_site - 1)][0] + 
                             hg38_siteDict[refGene][str(3*aa_site)][0])
                    aaChange[aa] = "p." + codon_1to3[Seq.translate(codon)] + aaChange[aa][5:]
            '''
            outf.write("\t".join(s[0:5]) + "\t" +
                        ";".join(ntChange) + "\t" +
                        aaChange + "\t" +
                        "\t".join(s[7:10]) + "\t" +
                        ";".join(hg38_position) + "\t" +
                        ";".join(hg19_position) + "\n")





def main(hg38_exonFile, hg19_exonFile, hg38_genome, hg38_version, hg19_genome, hg19_version, databaseFile, outdir):
    hg38_siteDict = cds_genomeLocation(hg38_exonFile, hg38_genome, hg38_version, outdir)
    hg19_siteDict = cds_genomeLocation(hg19_exonFile, hg19_genome, hg19_version, outdir)
    #1
    databaseFile = inputdir + "db_oriInput_old/018_H_Typer_database_hg38-hg19_v5.2_jy.txt"
    print(databaseFile)
    genomeLocation(databaseFile, hg19_siteDict, hg38_siteDict, outdir,"chr4")  
    

    '''
    with open("./hg38_siteDict.txt","w") as f:
        for element in hg38_siteDict.keys():
            f.write(element + "\n")
            for site in hg38_siteDict[element].keys():
                f.write(site + " " + " ".join(hg38_siteDict[element][site]) + "\n")
    with open("./hg19_siteDict.txt","w") as f:
        for element in hg19_siteDict.keys():
            f.write(element + "\n")
            for site in hg19_siteDict[element].keys():
                f.write(site + " " + " ".join(hg19_siteDict[element][site]) + "\n")
    '''


# Set data path  
inputdir = "/root/projects/projects_jy/WGS_blood/02_reference_database/"
outdir = "/root/projects/projects_jy/WGS_blood/02_reference_database/db_genomeAnno_old/"

hg19_exonFile = "/root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/db_region/43+7_exon_position_hg19.txt"
hg38_exonFile = "/root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/db_region/43+7_exon_position_hg38.txt"

hg19_genome = "/root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/Genome_refseq/GCF_000001405.25_GRCh37.p13_genomic.fna"
hg19_version= "hg19_NC" # For hg_version, please select one from "hg19_chr", "hg19_NC". chr means genome use chr* to annotate, such as chr1; NC means genome use NC_* to annotate, such as NC_000001.10

hg38_genome = "/root/projects/projects_jy/WGS_blood/04_blood_typing/supporting_files/Genome_refseq/GCF_000001405.39_GRCh38.p13_genomic.fna"
hg38_version = "hg38_NC" # For hg_version, please select one from "hg38_chr", "hg38_NC". chr means genome use chr* to annotate, such as chr1; NC means genome use NC_* to annotate, such as NC_000001.11

databaseFile = inputdir + "db_oriInput_old/018_H_Typer_database_hg38-hg19_v5.2_jy.txt"



# Run main function
main(hg38_exonFile, hg19_exonFile, hg38_genome, hg38_version, hg19_genome, hg19_version, databaseFile, outdir)

