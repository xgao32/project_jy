#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2021-5-26 9:28
# @Author  : info@genephar.com.cn
# @Link    : http://www.genephar.com.cn//
# @Version : Id
# @Telephone  : 89053423 or 89053788

import os, sys, re, random, time, argparse, datetime, shutil

###======= Pre-processing Functions =======######

# Function-01: Find target file
def find_file(sample,path1,file_suffix):#

	sample_list = os.listdir(path1)
	file_path = ''
	for sample_name in sample_list :
		if sample == sample_name :
			sample_path = path1 + '/' + sample_name
			file_list = os.listdir(sample_path)
			for file_name in file_list :
				if file_suffix in file_name :
					file_path = sample_path + '/' + file_name
					break
	#print (file_path)
	return file_path

# Function-02: Open and read target file
def txt_to_lines(file_path):#用于打开文件
	data = open(file_path)
	data_lines = data.readlines()
	data.close()
	return data_lines

# Function-03: Get sample names
def GetSampleName(INPUT):
	sample_names = []
	files = os.listdir(INPUT)
	for f in files :
		sample_name = f.split('.')[0]
		if sample_name not in sample_names :
			sample_names.append(sample_name)
	return sample_names

# Function-04: Fastq->bam Data mapping and pre-prosessing
def dataprocess(INPUT,OUTPUT,THREAD,sample_name,CLASS):

	absdir = sys.path[0]
	snprefdir = absdir+"/reference/snpeff_reference"
	softdir = absdir+"/software/soft"
	refdir = absdir+"/reference/hg38"
	#refdir = absdir+"/reference/hg19"
	refchr = refdir+"/hg38.fa"
	#refchr = refdir+"/hg19.fa"
	samtoolspro = softdir+"/samtools"
	flashpro = softdir+"/flash"
	bwapro = softdir+"/bwa-0.7.2/bwa"
	fgbiopro = softdir+"/fgbio-0.5.1.jar"
	seqkitpro = softdir+"/seqkit"
	picardpro = softdir+"/picard.jar"
	bbmergepro = softdir+"/bbmap/bbmerge.sh"
	fusepro = softdir+"/bbmap/fuse.sh"
	bampro = softdir+"/bamUtil/bin/bam"
	freebayespro = softdir+"/freebayes"
	tabixpro = softdir+"/tabix"
	bgzipro = softdir+"/bgzip"
	bcftoolspro = softdir+"/bcftools"
	annoconvert = softdir+"/annovar/convert2annovar.pl"
	annotable = softdir+"/annovar/table_annovar.pl"
	humandb = absdir+"/software/annovar_hg38/"
	snpEff = softdir+"/snpEff/snpEff/snpEff.jar"
	transcript = snprefdir+"/transcript.txt"
	path1 = OUTPUT + '/varResult/' + sample_name + '/'
	'''
	#merge read1 read2
	flash_min = 10
	flash_maxratio = 0.25
	cut_qual = 20
	inputfq1 = INPUT +'/'+ sample_name +'_'+'R1_001.fastq.gz'
	inputfq2 = INPUT +'/'+ sample_name +'_'+'R2_001.fastq.gz'
	outputf = path1 + sample_name +'_'+'merged.fastq.gz'
	outputfq1 = path1 + sample_name +'_'+'unmerged_R1.fastq.gz'
	outputfq2 = path1 + sample_name +'_'+'unmerged_R2.fastq.gz'
	cmd = """
		{bbmergepro} minoverlap={flash_min} \
		maxratio={flash_maxratio} \
		trimq={cut_qual} \
		in1={inputfq1} \
		in2={inputfq2} \
		out={outputf} \
		outu1={outputfq1} \
		outu2={outputfq2} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	os.system(cmd)
	os.remove(outputfq1)
	os.remove(outputfq2)
	
	#mapping
	thread = THREAD
	
	inputf = outputf
	outputf = path1 + sample_name + '.sam'
	bwaheader = '\"'+'@RG'+'\\t'+'ID:'+sample_name+'\\t'+'LB:'+sample_name+'\\t'+'PL:illumina'+'\\t'+'SM:'+sample_name+'\\t'+'PU:'+sample_name + '\"'
	cmd = """
		{bwapro} mem \
		-M -t {thread} \
		-R {bwaheader} \
		{refchr} {inputf} > {outputf} 2>/dev/null
		""".format(**dict(**locals()))
	os.system(cmd)

	inputf = outputf
	#thread = THREAD
	#inputf = INPUT +'/'+ sample_name +'.dedup.bam'
	outputf = path1 + sample_name + '.bam'
	cmd = """
		{samtoolspro} view \
		-bS \
		-o {outputf} \
		{inputf} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	os.system(cmd)
	#os.remove(inputf)
	
	thread = THREAD
	#inputf = outputf
	inputf = path1 + sample_name + '.bam'
	outputf = path1 + sample_name + '.sorted.bam'
	cmd = """
		{samtoolspro} sort 	-o {outputf} 		-@ {thread} 		{inputf} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)
	
	cutread1len = 20
	cutread2len = 20
	inputf = outputf
	outputf = path1 + sample_name + '.trimed.bam'
	cmd = """
		{bampro} trimBam \
		{inputf} \
		{outputf} \
		-c -L {cutread1len} \
		-R {cutread2len} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)

	inputf = outputf
	outputf = path1 + sample_name + '.trim.sorted.bam'

	cmd = """
		{samtoolspro} sort -@ {thread} {inputf} -o {outputf} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)

	#inputf = outputf
	thread = THREAD
	inputf = path1 + sample_name + '.sorted.bam'
	cmd = """
		{samtoolspro} index {inputf} -@ {thread} >/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	'''













	coverage = 20
	frequency = 0.04
	qual = 20
	altdepth = 10
	#referencebedfile = absdir + '/input/MNS_v1.bed'##这是引物区间？,全基因组应该不需要这东西
	referencebedfile = absdir + '/input/'+CLASS+'_hg38_freebayes.bed'
	#inputf = path1 + sample_name + '.dedup.bam'
	inputf = INPUT +sample_name + '.deduped.bam'
	outputf = path1 + sample_name + '.raw.vcf'
	cmd = """
		{freebayespro} \
		-f {refchr} \
		--min-coverage {coverage} \
		-F {frequency} \
		-C {altdepth} \
		-q {qual} \
		-t {referencebedfile} \
		{inputf} > {outputf}
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)

	inputf = outputf
	outputf = path1 + sample_name + '.sort.vcf'
	cmd = """
		java -Xmx4g -XX:ParallelGCThreads=8 \
		-jar {picardpro} \
		SortVcf \
		I={inputf} \
		O={outputf} 
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)

	inputf = outputf
	cmd = """
		{bgzipro} -f {inputf}
		""".format(**dict(**locals()))
	os.system(cmd)

	inputf = path1 + sample_name + '.sort.vcf.gz'
	cmd = """
		{tabixpro} {inputf} 2>/dev/null
		""".format(**dict(**locals()))
	os.system(cmd)

	outputf = path1 + sample_name + '.split.vcf'
	cmd = """
		{bcftoolspro} norm \
		-m-both -O v \
		-o {outputf} \
		{inputf} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)

	inputf = outputf
	outputf = path1 + sample_name + '.norm.vcf'
	cmd = """
		{bcftoolspro} norm \
		-f {refchr} \
		-o {outputf} \
		{inputf} 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)

	inputf = outputf
	outputf = path1 + sample_name + '.ann.vcf'
	cmd = """
		java -Xmx4g -XX:ParallelGCThreads=8 -jar \
		{snpEff} \
		-onlyTr {transcript} \
		hg19 \
		{inputf} > {outputf} 2>/dev/null
		""".format(**dict(**locals()))
	os.system(cmd)
	#os.remove(inputf)
	temp1 = absdir + '/snpEff_summary.html'
	temp2 = absdir + '/snpEff_genes.txt'
	#os.remove(temp1)
	#os.remove(temp2)

	inputf = outputf
	outputf = path1 + sample_name + '.avinput'
	cmd = """
		perl {annoconvert} \
		-format vcf4 \
		--includeinfo \
		--allsample \
		--withfreq \
		{inputf} > {outputf} 2>/dev/null
		""".format(**dict(**locals()))
	os.system(cmd)
	#os.remove(inputf)

	inputf = outputf
	outputf = path1 + sample_name
	#perl /home/genephar/workflow/ABO_Typer/software/soft/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene  ./ 下载数据库
	cmd = """
		perl {annotable} {inputf} {humandb} \
		--buildver hg38 \
		--outfile {outputf} \
		--remove \
		--protocol refGene \
		-operation g \
		--argument '-hgvs' \
		--nastring . \
		--otherinfo 1>/dev/null 2>&1
		""".format(**dict(**locals()))
	print(cmd)
	os.system(cmd)
	#os.remove(inputf)


# Function-05: 
def snapshot(outputDir, sample):
	results = open(outputDir + "/varResult/" + sample + "/" + sample + ".Final_variants.txt", "r")
	resultsList = []
	for line in results:
		lines = line.split("\n")[0].split("\t")
		for i in lines[2].split(";"):
			if i != "":
				resultsList.append(i)
		for j in lines[5].split(";"):
			if j != "":
				resultsList.append(j)
	results.close()
	results2 = open(outputDir + "/varResult/" + sample + "/" + sample + ".Other_variants.txt", "r")
	for line in results2:
		lines = line.split("\n")[0]
		if lines == "":
			continue
		else:
			resultsList.append(lines.split("\t")[0])
	results2.close()
	resultsList = list(set(resultsList))
	# resultsListOutput = []

	if not os.path.exists(outputDir + "/varResult/" + sample + "/snapshot"):
		os.makedirs(outputDir + "/varResult/" + sample + "/snapshot")

	for k in resultsList:
		ss.main(outputDir + "/varResult/" + sample + "/" + sample + ".sorted.bam",
			outputDir + "/varResult/" + sample + "/snapshot/" + k + ".png",
			"chr9:" + k)

# Function-06: Generate beds
def genarateBed(outputDir, sample):
	results = open(outputDir + "/varResult/" + sample + "/" + sample + ".Final_variants.txt", "r")
	p = []
	for r in results:
		if ";" in r:
			rs = r.split("\t")
			rss = rs[2] + rs[5]
			p = rss.replace(";\n", "").split(";")
	results.close()

	results2 = open(outputDir + "/varResult/" + sample + "/" + sample + ".Other_variants.txt", "r")
	for r2 in results2:
		if r2 != "\n":
			p.append(r2.split("\t")[0])
	results2.close()
	p = list(set(p))
	bedfile = outputDir + "/varResult/" + sample + "/" + sample + ".mpileup.bed"
	bed = open(bedfile, "w")
	for i in p:
		bed.write("chr9\t" + i + "\n")
	bed.close()
	abpath = sys.path[0]
	samtools = abpath + "/software/soft/samtools"
	hg38 = abpath + "/reference/hg38/hg38.fa"
	cmd = """
		{samtools} mpileup -f {hg19} -l {bedfile} {outputDir}/varResult/{sample}/{sample}.sorted.bam -d 0 \\
			> {outputDir}/varResult/{sample}/{sample}.mpileup
	""".format(samtools=samtools, hg19=hg38, bedfile=bedfile, outputDir=outputDir, sample=sample)
	os.system(cmd)
	pileup = open(outputDir + "/varResult/" + sample + "/" + sample + ".mpileup", "r")
	output = open(outputDir + "/varResult/" + sample + "/" + sample + ".mpileup.results", "w")
	for line in pileup:
		linesplit = line.split("\t")
		chrom = linesplit[0]
		position = linesplit[1]
		ref = linesplit[2]
		A = "A: " + str(linesplit[4].count("A") + linesplit[4].count("a"))
		T = "T: " + str(linesplit[4].count("T") + linesplit[4].count("t"))
		C = "C: " + str(linesplit[4].count("C") + linesplit[4].count("c"))
		G = "G: " + str(linesplit[4].count("G") + linesplit[4].count("g"))

		if ref == "A":
			A = "A:" + str(linesplit[4].count(".") + linesplit[4].count(","))
		elif ref == "T":
			T = "T:" + str(linesplit[4].count(".") + linesplit[4].count(","))
		elif ref == "C":
			C = "C:" + str(linesplit[4].count(".") + linesplit[4].count(","))
		elif ref == "G":
			G = "G:" + str(linesplit[4].count(".") + linesplit[4].count(","))
		else:
			continue
		output.write(chrom + "\t" + position + "\t(" + ", ".join([A, T, C, G]) + ")\n")
	pileup.close()
	output.close()

def GLOB_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GLOB_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GLOB_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GLOB_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GLOB_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	database_path = abpath + '/input/GLOB_Typer_database_hg38_v3.0.01'
	GLOB_Typer_database_l = txt_to_lines(database_path)
	GLOB_Typer_database_l = GLOB_Typer_database_l[1:]
	Haplotypes = GLOB_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GLOB_Type_result.txt'
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	for chr1,chr2 in Haplotypes.items():
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in GLOB_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def MNS_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/MNS_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".MNS_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".MNS_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				#af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
					else:
						if float(af) <= 0.2 :#纯合无突变  chr11 --N--
														#chr21 --N--
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变    chr11 --Y--
														#chr21 --Y--
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.2<float(af)<0.8 :
											 
											  #未修改   杂合突变第一种可能  chr11 --N--     第二种可能chr21--Y--
											  #                             chr12 --Y--              chr22--N--
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
							print(Haplotypes_temp)
				Haplotypes = Haplotypes_temp	
		#print(Haplotypes)
		return Haplotypes					
def MNS_Haplotype_maker2(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	#kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	hg38_cnv = abpath + "/reference/hg38/hg38.fa"
	hg38fai_cnv = abpath + "/reference/hg38/hg38.fa.fai"
	bedfile = abpath + "/input/RHD_bed"
	bed_cnv_file = abpath + "/input/CNV_rhd_rche_gypa_gypb.bed"
	cnvkit = abpath + "/software/soft/cnvkit/cnvkit.py"
	annotate = abpath + "/input/refFlat.txt"

	cmd = """
		{cnvkit} batch {input}{sample}.deduped.bam -n --targets {bed_cnv_file} -f {hg38_cnv} --output-reference {output}/{sample}.cnn --output-dir {output} --annotate {an}
	""".format(cnvkit=cnvkit, hg38_cnv=hg38_cnv, hg38fai_cnv=hg38fai_cnv,bed_cnv_file=bed_cnv_file, sample=sample_name,an=annotate,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_cnv')
	print(cmd)
	os.system(cmd)
	#GYPA-ENST00000641688.2
	#ex7	144,140,718	144,140,599	
	#ex6	144,120,588	144,120,490
	#ex5	144,119,781	144,119,686
	#ex4	144,118,752	144,118,714
	#ex3	144,116,939	144,116,854		
	#ex2	144,114,767	144,114,689
	#ex1	144,111,410	144,109,303
	#GYPB-ENST00000502664.6
	#ex5	144,019,380	144,019,251
	#ex4	144,001,283	144,001,185
	#ex3	143,999,449	143,999,411		
	#ex2	143,997,634	143,997,540
	#ex1	143,996,304	143,996,104
	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_cnv'+"/"+sample_name+".deduped.cnr","r") as f1:
		rhd_cnv_ex2=0.0
		rhd_cnv_dp=0.0
		rhce_cnv_ex2=0.0
		rhce_cnv_dp=0.0
		rh_typing_result=[]
		rhce_typing_result="exon2nochange"
		#global exon2_rhce
		for line in f1.readlines():
			if line.startswith("chromosome"):
				continue
			else:
				ex_start=line.split("\t")[1]
				ex_end = line.split("\t")[2]
				gene_name = line.split("\t")[3]
				dp = float(line.split("\t")[4])
				log2=line.split("\t")[5]
				if gene_name =="RSRP1,RHD":
					rhd_cnv_dp = rhd_cnv_dp+dp
				if gene_name =="RHCE":
					rhce_cnv_dp= rhce_cnv_dp+dp
				#0.585 #-1 #215估计是215个bin位置
				#对于拷贝数的无需要区分纯杂合缺失，只要有的话就爆出来就行了因为至少一个染色体发生了，看能否后期有需要的话加一个频率，此频率和AF计算方法不一样
				if ex_start=="143996104" and ex_end =="143996304":
						if float(log2)>0.585:
							exon1_gypb="gain"
							continue
						elif float(log2)<-1.0:
							exon1_gypb="loss"
							continue
						else:
							exon1_gypb="middle"
							continue
				if ex_start=="143997540" and ex_end =="143997634":
						if float(log2)>0.585:
							exon2_gypb="gain"
							continue
						elif float(log2)<-1.0:
							exon2_gypb="loss"
							continue
						else:
							exon2_rhd="middle"
							continue
				if ex_start=="143999411" and ex_end =="143999449":
						if float(log2)>0.585:
							exon3_gypb="gain"
							continue
						elif float(log2)<-1.0:
							exon3_gypb="loss"
							continue
						else:
							exon3_gypb="middle"
							continue
				if ex_start=="144001185" and ex_end =="144001283":
						if float(log2)>0.585:
							exon4_gypb="gain"
							continue
						elif float(log2)<-1.0:
							exon4_gypb="loss"
							continue
						else:
							exon4_gypb="middle"
							continue
				if ex_start=="144019251" and ex_end =="144019380":
						if float(log2)>0.585:
							exon5_gypb="gain"
						elif float(log2)<-1.0:
							exon5_gypb="loss"					
						else:
							exon5_gypb="middle"
				if ex_start=="144109303" and ex_end =="144111410":
						if float(log2)>0.585:
							exon1_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon1_gypa="loss"	
							continue
						else:
							exon1_gypa="middle"
							continue
				if ex_start=="144114689" and ex_end =="144114767":
						if float(log2)>0.585:
							exon2_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon2_gypa="loss"						
							continue
						else:
							exon2_gypa="middle"
							continue
				if ex_start=="144116854" and ex_end =="144116939":
						if float(log2)>0.585:
							exon3_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon3_gypa="loss"	
							continue
						else:
							exon3_gypa="middle"
							continue
				if ex_start=="144118714" and ex_end =="144118752":
						if float(log2)>0.585:
							exon4_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon4_gypa="loss"	
							continue
						else:
							exon4_gypa="middle"
							continue
				if ex_start=="144119686" and ex_end =="144119781":
						if float(log2)>0.585:
							exon5_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon5_gypa="loss"						
							continue
						else:
							exon5_gypa="middle"
							continue
				if ex_start=="144120490" and ex_end =="144120588":
						if float(log2)>0.585:
							exon6_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon6_gypa="loss"
							continue
						else:
							exon6_gypa="middle"
							continue
				#print("ok1")
				if ex_start=="144140599" and ex_end =="144140718":
						if float(log2)>0.585:
							exon7_gypa="gain"
							continue
						elif float(log2)<-1.0:
							exon7_gypa="loss"
							continue
						else:
							exon7_gypa="middle"
							continue
				#print("ok2")
		if exon3_gypa=="loss":
				rh_typing_result.append("MNS:St(a+)(GYP*Zan)\tGYP*Zan\tdel exon 3\tp.Asp46_Thr77del\texonic\tnonsynonymous SNV\t.\t.\tISBT\texon3")
		if exon4_gypa=="loss" and exon5_gypa=="loss" and exon6_gypa=="loss" and exon4_gypb=="gain" and  exon5_gypb=="gain" and exon6_gypb=="gain":
				rh_typing_result.append("MNS:Hil+(GYP*Hil)\tGYP*Hil\tGYP(A1–232–B233–312)\tGP(A1–77–B78–104)\texonic\tnonsynonymous SNV\t.\t.\tISBT\texon3")
		if exon4_gypa=="loss" and exon5_gypa=="loss" and exon6_gypa=="loss" and exon4_gypb=="gain" and  exon5_gypb=="gain" and exon6_gypb=="gain":
				rh_typing_result.append("MNS:SAT+(GYP*SAT)\tGYP*SAT\tGYP(A1–232–B233–312)\tGP(A1–77–B78–104)\texonic\tnonsynonymous SNV\t.\t.\tISBT\texon3")
		'''
		if  1.6<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2:
				rh_typing_result.append("RHD:D+(RHD*01)\t")
			#RHD*01N.01 (RHD的外显子全部缺失D-)
		if  0.6<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2<1.5:
				rh_typing_result.append("RHD:D+(RHD*01)\tRHD*01\nRHD:D-(RHD*01N.01)\tRHD*01N.01\tRHD deletion\tRHD deletion\texonic\tSV\t.\tRHD deletion\tISBT\tRHD deletion")
			#RHD*01N.01 (RHD的外显子全部缺失D-)
		if  0.0<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2<0.5:
				rh_typing_result.append("RHD:D-(RHD*01N.01)\tRHD*01N.01\tRHD deletion\tRHD deletion\texonic\tSV\t.\tRHD deletion\tISBT\tRHD deletion")
		'''
		return rh_typing_result
def MNS_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/MNS_Typer_database_hg38_v3.0.01'
	MNS_Typer_database_l = txt_to_lines(database_path)
	MNS_Typer_database_l = MNS_Typer_database_l[1:]
	Haplotypes = MNS_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	#print(Haplotypes)
	Type_result = path1 + sample_name + '/' + sample_name + '.MNS_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in MNS_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def LU_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/LU_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LU_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LU_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def LU_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/LU_Typer_database_hg38_v3.0.01'
	LU_Typer_database_l = txt_to_lines(database_path)
	LU_Typer_database_l = LU_Typer_database_l[1:]
	Haplotypes = LU_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.LU_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in LU_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def KEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/KEL_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".KEL_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		


	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".KEL_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def KEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/KEL_Typer_database_hg38_v3.0.01'
	KEL_Typer_database_l = txt_to_lines(database_path)
	KEL_Typer_database_l = KEL_Typer_database_l[1:]
	Haplotypes = KEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.KEL_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in KEL_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def FY_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/FY_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".FY_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".FY_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def FY_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	database_path = abpath + '/input/FY_Typer_database_hg38_v3.0.01'
	FY_Typer_database_l = txt_to_lines(database_path)
	FY_Typer_database_l = FY_Typer_database_l[1:]
	Haplotypes = FY_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.FY_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in FY_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def JK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/JK_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JK_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JK_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def JK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/JK_Typer_database_hg38_v3.0.01'
	JK_Typer_database_l = txt_to_lines(database_path)
	JK_Typer_database_l = JK_Typer_database_l[1:]
	Haplotypes = JK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.JK_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in JK_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def DI_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/DI_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".DI_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 


	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".DI_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def DI_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/DI_Typer_database_hg38_v3.0.01'
	DI_Typer_database_l = txt_to_lines(database_path)
	DI_Typer_database_l = DI_Typer_database_l[1:]
	Haplotypes = DI_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.DI_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in DI_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def YT_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/YT_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".YT_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".YT_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def YT_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/YT_Typer_database_hg38_v3.0.01'
	YT_Typer_database_l = txt_to_lines(database_path)
	YT_Typer_database_l = YT_Typer_database_l[1:]
	Haplotypes = YT_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.YT_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in YT_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def SC_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/SC_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".SC_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".SC_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def SC_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/SC_Typer_database_hg38_v3.0.01'
	SC_Typer_database_l = txt_to_lines(database_path)
	SC_Typer_database_l = SC_Typer_database_l[1:]
	Haplotypes = SC_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.SC_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in SC_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def DO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/DO_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".DO_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)


	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".DO_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def DO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/DO_Typer_database_hg38_v3.0.01'
	DO_Typer_database_l = txt_to_lines(database_path)
	DO_Typer_database_l = DO_Typer_database_l[1:]
	Haplotypes = DO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.DO_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in DO_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def CO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CO_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CO_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CO_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CO_Typer_database_hg38_v3.0.01'
	CO_Typer_database_l = txt_to_lines(database_path)
	CO_Typer_database_l = CO_Typer_database_l[1:]
	Haplotypes = CO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CO_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in CO_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def LW_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/LW_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LW_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LW_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def LW_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/LW_Typer_database_hg38_v3.0.01'
	LW_Typer_database_l = txt_to_lines(database_path)
	LW_Typer_database_l = LW_Typer_database_l[1:]
	Haplotypes = LW_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.LW_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in LW_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def CHRG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CHRG_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CHRG_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CHRG_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				try:
					af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				except:
					af = 0.0
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CHRG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CHRG_Typer_database_hg38_v3.0.01'
	CHRG_Typer_database_l = txt_to_lines(database_path)
	CHRG_Typer_database_l = CHRG_Typer_database_l[1:]
	Haplotypes = CHRG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CHRG_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in CHRG_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def CROM_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CROM_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CROM_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CROM_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CROM_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]
	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CROM_Typer_database_hg38_v3.0.01'
	CROM_Typer_database_l = txt_to_lines(database_path)
	CROM_Typer_database_l = CROM_Typer_database_l[1:]
	Haplotypes = CROM_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CROM_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in CROM_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def RHAG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/RHAG_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHAG_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHAG_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def RHAG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/RHAG_Typer_database_hg38_v3.0.01'
	RHAG_Typer_database_l = txt_to_lines(database_path)
	RHAG_Typer_database_l = RHAG_Typer_database_l[1:]
	Haplotypes = RHAG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.RHAG_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in RHAG_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def KN_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/RHAG_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHAG_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHAG_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def KN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/RHAG_Typer_database_hg38_v3.0.01'
	RHAG_Typer_database_l = txt_to_lines(database_path)
	RHAG_Typer_database_l = RHAG_Typer_database_l[1:]
	Haplotypes = RHAG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.RHAG_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in RHAG_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def IN_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/IN_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".IN_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".IN_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def IN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/IN_Typer_database_hg38_v3.0.01'
	IN_Typer_database_l = txt_to_lines(database_path)
	IN_Typer_database_l = IN_Typer_database_l[1:]
	Haplotypes = RHAG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.IN_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in IN_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def OK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/OK_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".OK_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".OK_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def OK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/OK_Typer_database_hg38_v3.0.01'
	OK_Typer_database_l = txt_to_lines(database_path)
	OK_Typer_database_l = OK_Typer_database_l[1:]
	Haplotypes = OK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.OK_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in OK_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def RAPH_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/RAPH_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RAPH_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RAPH_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def RAPH_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/RAPH_Typer_database_hg38_v3.0.01'
	RAPH_Typer_database_l = txt_to_lines(database_path)
	RAPH_Typer_database_l = RAPH_Typer_database_l[1:]
	Haplotypes = RAPH_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.RAPH_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in RAPH_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def JMH_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/JMH_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JMH_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JMH_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def JMH_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/JMH_Typer_database_hg38_v3.0.01'
	JMH_Typer_database_l = txt_to_lines(database_path)
	JMH_Typer_database_l = JMH_Typer_database_l[1:]
	Haplotypes = JMH_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.JMH_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		print(len(chr1))
		for l3 in JMH_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def I_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/I_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".I_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".I_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def I_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/I_Typer_database_hg38_v3.0.01'
	I_Typer_database_l = txt_to_lines(database_path)
	I_Typer_database_l = I_Typer_database_l[1:]
	Haplotypes = I_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.I_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in I_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def GIL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GIL_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print("oko"+cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GIL_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GIL_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GIL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	#print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GIL_Typer_database_hg38_v3.0.01'
	GIL_Typer_database_l = txt_to_lines(database_path)
	GIL_Typer_database_l = GIL_Typer_database_l[1:]
	Haplotypes = GIL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GIL_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GIL_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def FROS_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/FROS_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".FROS_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".FROS_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def FROS_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/FROS_Typer_database_hg38_v3.0.01'
	FROS_Typer_database_l = txt_to_lines(database_path)
	FROS_Typer_database_l = FROS_Typer_database_l[1:]
	Haplotypes = FROS_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.FROS_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in FROS_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def LAN_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/LAN_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LAN_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LAN_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def LAN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/LAN_Typer_database_hg38_v3.0.01'
	LAN_Typer_database_l = txt_to_lines(database_path)
	LAN_Typer_database_l = LAN_Typer_database_l[1:]
	Haplotypes = LAN_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.LAN_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in LAN_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def VEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/VEL_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".VEL_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".VEL_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def VEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/VEL_Typer_database_hg38_v3.0.01'
	VEL_Typer_database_l = txt_to_lines(database_path)
	VEL_Typer_database_l = VEL_Typer_database_l[1:]
	Haplotypes = VEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.VEL_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in VEL_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def CD59_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CD59_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CD59_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CD59_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CD59_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CD59_Typer_database_hg38_v3.0.01'
	CD59_Typer_database_l = txt_to_lines(database_path)
	CD59_Typer_database_l = CD59_Typer_database_l[1:]
	Haplotypes = CD59_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CD59_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in CD59_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def AUG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/AUG_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".AUG_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".AUG_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def AUG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/AUG_Typer_database_hg38_v3.0.01'
	AUG_Typer_database_l = txt_to_lines(database_path)
	AUG_Typer_database_l = AUG_Typer_database_l[1:]
	Haplotypes = AUG_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.AUG_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in AUG_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def CTL2_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CTL2_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CTL2_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CTL2_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CTL2_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CTL2_Typer_database_hg38_v3.0.01'
	CTL2_Typer_database_l = txt_to_lines(database_path)
	CTL2_Typer_database_l = CTL2_Typer_database_l[1:]
	Haplotypes = CTL2_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CTL2_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in CTL2_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def GE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GE_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GE_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GE_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GE_Typer_database_hg38_v3.0.01'
	GE_Typer_database_l = txt_to_lines(database_path)
	GE_Typer_database_l = GE_Typer_database_l[1:]
	Haplotypes = GE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GE_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GE_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def H_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/H_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".H_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".H_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def H_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/H_Typer_database_hg38_v3.0.01'
	H_Typer_database_l = txt_to_lines(database_path)
	H_Typer_database_l = H_Typer_database_l[1:]
	Haplotypes = H_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.H_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in H_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def KANNO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/KANNO_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".KANNO_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".KANNO_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def KANNO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/KANNO_Typer_database_hg38_v3.0.01'
	KANNO_Typer_database_l = txt_to_lines(database_path)
	KANNO_Typer_database_l = KANNO_Typer_database_l[1:]
	Haplotypes = KANNO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.KANNO_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in KANNO_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def MAM_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/MAM_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".MAM_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".MAM_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def MAM_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/MAM_Typer_database_hg38_v3.0.01'
	MAM_Typer_database_l = txt_to_lines(database_path)
	MAM_Typer_database_l = MAM_Typer_database_l[1:]
	Haplotypes = MAM_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.MAM_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in MAM_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def P1PK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/P1PK_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".P1PK_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".P1PK_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def P1PK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/P1PK_Typer_database_hg38_v3.0.01'
	P1PK_Typer_database_l = txt_to_lines(database_path)
	P1PK_Typer_database_l = P1PK_Typer_database_l[1:]
	Haplotypes = P1PK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.P1PK_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in P1PK_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def SID_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/SID_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".SID_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".SID_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def SID_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/SID_Typer_database_hg38_v3.0.01'
	SID_Typer_database_l = txt_to_lines(database_path)
	SID_Typer_database_l = SID_Typer_database_l[1:]
	Haplotypes = SID_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.SID_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in SID_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def LE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/LE_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LE_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".LE_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def LE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/LE_Typer_database_hg38_v3.0.01'
	LE_Typer_database_l = txt_to_lines(database_path)
	LE_Typer_database_l = LE_Typer_database_l[1:]
	Haplotypes = LE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.LE_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in LE_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def XK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/XK_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".XK_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".XK_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def XK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/XK_Typer_database_hg38_v3.0.01'
	XK_Typer_database_l = txt_to_lines(database_path)
	XK_Typer_database_l = XK_Typer_database_l[1:]
	Haplotypes = XK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.XK_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in XK_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def JR_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/JR_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JR_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".JR_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes			
def JR_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/JR_Typer_database_hg38_v3.0.01'
	JR_Typer_database_l = txt_to_lines(database_path)
	JR_Typer_database_l = JR_Typer_database_l[1:]
	Haplotypes = JR_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.JR_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in JR_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()
def PEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
		if =="no_mutation":
		return "PEL:PEL+"			
def PEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	Haplotypes = PEL_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.PEL_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	result2.write(Haplotypes)
	result2.close()

def ABO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	#bedfile = abpath + "/input/ABO_bed"
	bedfile = abpath + "/input/ABO_hg38_svOUT.bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".ABO_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".ABO_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		af_dic={}
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				pos_alt = re.split(r"\.|\t",line)[4]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print (af)
				af_dic[pos_start+"_"+pos_end+"_"+pos_alt] =af	
        '''
		if float(af_dic["48751247_48751247_G"])>=0.2 or float(af_dic["48751247_48751247_A"])>=0.2 or float(af_dic["48703384_48703384_G"])>=0.2 or float(af_dic["48702960_48702960_G"])>=0.2 or float(af_dic["48702996_48702996_G"])>=0.2 or float(af_dic["48702996_48702996_G"])>=0.2  and float(af_dic["48703069_48703069_T"])>=0.2:
			H="+"
		elif float(af_dic["48702996_48702996_G"])>=0.2  and float(af_dic["48703437_48703437_A"])>=0.2 or float(af_dic["48703335_48703335_T"])>=0.2 or float(af_dic["48703356_48703356_A"])>=0.2 or float(af_dic["48703437_48703437_A"])>=0.2 or float(af_dic["48703621_48703621_A"])>=0.2 or float(af_dic["48703641_48703641_A"])>=0.2 or float(af_dic["48703672_48703672_A"])>=0.2:
			H="+"
		elif float(af_dic["48702996_48702996_G"])>=0.2 or float(af_dic["48703703_48703704_insGTG"])>=0.2 or float(af_dic["48703234_48703234_T"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 and float(af_dic["48703573_48703573_G"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 and float(af_dic["48703797_48703797_T"])>=0.2 or float(af_dic["48703809_48703809_A"])>=0.2:
			H="+"
		else:
			H="-"    
		'''
		if float(af_dic["48750860_48750860_T"])>=0.2 or float(af_dic["48750821_48750821_C"])>=0.2 or float(af_dic["48750820_48750820_T"])>=0.2 or float(af_dic["48750769_48750769_G"])>=0.2 or float(af_dic["48750744_48750744_A"])>=0.2:      
			H="-"
		#elif float(af_dic["48702996_48702996_G"])>=0.2  and float(af_dic["48703437_48703437_A"])>=0.2 or float(af_dic["48703335_48703335_T"])>=0.2 or float(af_dic["48703356_48703356_A"])>=0.2 or float(af_dic["48703437_48703437_A"])>=0.2 or float(af_dic["48703621_48703621_A"])>=0.2 or float(af_dic["48703641_48703641_A"])>=0.2 or float(af_dic["48703672_48703672_A"])>=0.2:
		elif float(af_dic["48703703_48703704_delAG"])>=0.2 or float(af_dic["48702996_48702996_A"])>=0.2 or float(af_dic["48750587_48750587_T"])>=0.2 or float(af_dic["48750557_48750557_C"])>=0.2 or float(af_dic["48750506_48750506_T"])>=0.2 or  float(af_dic["48750497_48750497_T"])>=0.2 and float(af_dic["48750496_48750496_T"])>=0.2 or float(af_dic["48750456_48750456_A"])>=0.2: 
			H="-"
		elif float(af_dic["48750400_48750401_delTT"])>=0.2 or float(af_dic["48750338_48750338_A"])>=0.2 or  float(af_dic["48750334_48750334_C"])>=0.2 or   float(af_dic["48750302_48750302_G"])>=0.2 or float(af_dic["48750235_48750235_G"])>=0.2 or float(af_dic["48750598_48750598_T"])>=0.2 or float(af_dic["48750588_48750588_T"])>=0.2 or float(af_dic["48750514_48750514_delC"])>=0.2 or float(af_dic["48750492_48750492_insG"])>=0.2 or float(af_dic["48750572_48750572_delG"])>=0.2 or float(af_dic["48750828_48750828_delG"])>=0.2:
		#elif float(af_dic["48702996_48702996_G"])>=0.2 or float(af_dic["48703703_48703704_insGTG"])>=0.2 or float(af_dic["48703234_48703234_T"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 and float(af_dic["48703573_48703573_G"])>=0.2 or float(af_dic["48703341_48703341_T"])>=0.2 and float(af_dic["48703797_48703797_T"])>=0.2 or float(af_dic["48703809_48703809_A"])>=0.2:
			H="-"
		elif float(af_dic["48750994_48750994_T"])>=0.2 or float(af_dic["48750859_48750859_T"])>=0.2:
			H="-"
		else:
			H="+"    
	if H!="-":
		#o-（c.261delG(insC)，c.802A）--133257521_133257521_C,133255929_133255929_T
		#a-（c.526C，c.703G,c.796C,c.803G)--133256205_133256205_G,133256028_133256028_C,133255935_133255935_G,133255928_133255928_C
		#b-（c.526G，c.703A,c.796A,c.803C)--133256205_133256205_C,133256028_133256028_T,133255935_133255935_T,133255928_133255928_G
		if 0.8<float(af_dic["133257521_133257521_C"])<=1.0 and float(af_dic["133255929_133255929_T"])<0.2:#ABO基因c.261insC的突变频率是不是大于0.8小于等于1同时c.802A的突变频率小于0.2
			if float(af_dic["133256205_133256205_C"])>0.2 and float(af_dic["133256028_133256028_T"])>0.2 and float(af_dic["133255935_133255935_T"])>0.2 and float(af_dic["133255928_133255928_G"])>0.2:
				if  float(af_dic["133256205_133256205_G"])>0.2 and float(af_dic["133256028_133256028_C"])>0.2 and float(af_dic["133255935_133255935_G"])>0.2 and float(af_dic["133255928_133255928_C"])>0.2:
					return "AB"
				else:	
					if float(af_dic["133256205_133256205_C"])>0.8 and float(af_dic["133256028_133256028_T"])>0.8 and float(af_dic["133255935_133255935_T"])>0.8 and float(af_dic["133255928_133255928_G"])>0.8:	
						return "BB"
					else:
						return "Bb"
			else:
				if  0.8>float(af_dic["133256205_133256205_G"])>0.2 and 0.8>float(af_dic["133256028_133256028_C"])>0.2 and 0.8>float(af_dic["133255935_133255935_G"])>0.2 and 0.8>float(af_dic["133255928_133255928_C"])>0.2:
					return "Aa"
				else:
					if  float(af_dic["133256205_133256205_G"])>0.8 and float(af_dic["133256028_133256028_C"])>0.8 and float(af_dic["133255935_133255935_G"])>0.8 and float(af_dic["133255928_133255928_C"])>0.8: 
						return "AA"
					else:
						return "ab"
		else:
			if float(af_dic["133257521_133257521_C"])<0.2 and 0.8<float(af_dic["133255929_133255929_T"]):
				return "OO"
			else:
				if float(af_dic["133257521_133257521_C"])<0.2 or 0.8<float(af_dic["133255929_133255929_T"]):
					return "OO"
				elif 0.2<float(af_dic["133257521_133257521_C"])<0.8 and 0.2<float(af_dic["133255929_133255929_T"])<0.8:
					return "OO"
				elif 0.2<float(af_dic["133257521_133257521_C"])<0.8 and 0.8<float(af_dic["133255929_133255929_C"]) and 0.8<float(af_dic["133256205_133256205_G"]) and 0.8<float(af_dic["133256028_133256028_C"]) and 0.8<float(af_dic["133255935_133255935_G"]) and 0.8<float(af_dic["133255928_133255928_C"])	:
					return "AO"
				elif 0.2<float(af_dic["133257521_133257521_C"])<0.8 and 0.8<float(af_dic["133255929_133255929_C"]) and 0.2<float(af_dic["133256205_133256205_C"])<0.8 and 0.2<float(af_dic["133256028_133256028_T"])<0.8 and 0.2<float(af_dic["133255935_133255935_T"])<0.8 and 0.2<float(af_dic["133255928_133255928_G"])<0.8:
					return "BO"
				else:
					return "3?"        
	else:
		return "MM"
def ABO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	#database_path = abpath + '/input/ABO_Typer_database_hg38_v3.0.01'
	#ABO_Typer_database_l = txt_to_lines(database_path)
	#ABO_Typer_database_l = ABO_Typer_database_l[1:]
	Haplotypes = ABO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.ABO_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	
	result2 = open(Type_result,'a')
	result2.write(Haplotypes)
	result2.close()

def RHD_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	#kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	hg38_cnv = abpath + "/reference/hg38/hg38.fa"
	hg38fai_cnv = abpath + "/reference/hg38/hg38.fa.fai"
	bedfile = abpath + "/input/RHD_bed"
	bed_cnv_file = abpath + "/input/CNV_rhd_rche.bed"
	cnvkit = abpath + "/software/soft/cnvkit/cnvkit.py"
	annotate = abpath + "/input/refFlat.txt"

	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
	
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHD_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHD_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
					if line.startswith("chr"):
						pos_chr = re.split(r"\.|\t",line)[0]
						pos_start = re.split(r"\.|\t",line)[1]
						pos_end = re.split(r"\.|\t",line)[2]
						dp = re.split(r"\.|\t",line)[5]
						ad = re.split(r"\.|\t",line)[6]
						af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
						#print(line+af)
						Haplotypes_temp = {}
						dic_temp = {}
						for chr1,chr2 in Haplotypes.items():			
							if pos_start in kb500:
								pass
								'''
								if af <= 0.1 :#纯合无突变
									chr11 = chr1 + 'N'
									chr21 = chr2 + 'N'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)							
								if af >= 0.8 :#纯合突变
									chr11 = chr1 + 'Y'
									chr21 = chr2 + 'Y'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)	
								if 0.1<af<0.8 :#杂合突变
									if kb500 in another_ID:
										if another_ID <=0.1
											chr11 = chr1 + 'N'
											chr21 = chr1 + 'Y'
											chr22 = chr2 + 'N'
											chr12 = chr2 + 'Y'
								'''
							else:
								if float(af) <= 0.1 :#纯合无突变
									chr11 = chr1 + 'N'
									chr21 = chr2 + 'N'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)							
								if float(af) >= 0.8 :#纯合突变
									chr11 = chr1 + 'Y'
									chr21 = chr2 + 'Y'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)	
								if 0.1<float(af)<0.8 :#杂合突变
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
									#print(chr11+"ok"+chr12)
									dic_temp = {chr11:chr12,chr21:chr22}
									Haplotypes_temp.update(dic_temp)
						Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes
def RHD_Haplotype_maker2(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	#kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	hg38_cnv = abpath + "/reference/hg38/hg38.fa"
	hg38fai_cnv = abpath + "/reference/hg38/hg38.fa.fai"
	bedfile = abpath + "/input/RHD_bed"
	bed_cnv_file = abpath + "/input/CNV_rhd_rche_gypa_gypb.bed"
	cnvkit = abpath + "/software/soft/cnvkit/cnvkit.py"
	annotate = abpath + "/input/refFlat.txt"

	cmd = """
		{cnvkit} batch {input}{sample}.deduped.bam -n --targets {bed_cnv_file} -f {hg38_cnv} --output-reference {output}/{sample}.cnn --output-dir {output} --annotate {an}
	""".format(cnvkit=cnvkit, hg38_cnv=hg38_cnv, hg38fai_cnv=hg38fai_cnv,bed_cnv_file=bed_cnv_file, sample=sample_name,an=annotate,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_cnv')
	print(cmd)
	os.system(cmd)
	#RHD-ENST00000622561.4
	#exon1	25,272,509	25,272,695
	#exon2	25,284,573	25,284,759
	#exon3	25,290,641	25,290,791
	#exon4	25,300,946	25,301,093
	#exon5	25,301,520	25,301,686
	#exon6	25,303,322	25,303,459
	#exon7	25,306,596	25,306,729
	#exon8	25,317,000	25,317,079
	#exon9	25,321,889	25,321,962
	#exon10	25,328,898	25,330,445

	#RHCE-ENST00000294413.13
	#exon1	25,420,825	25,420,639
	#exon2	25,408,869	25,408,683
	#exon3	25,402,746	25,402,596
	#exon4	25,392,141	25,391,994
	#exon5	25,390,915	25,390,749
	#exon6	25,389,113	25,388,976
	#exon7	25,385,844	25,385,711
	#exon8	25,375,428	25,375,349
	#exon9	25,370,540	25,370,467
	#exon10	25,362,553	25,362,249	 
	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_cnv'+"/"+sample_name+".deduped.cnr","r") as f1:
		rhd_cnv_ex2=0.0
		rhd_cnv_dp=0.0
		rhce_cnv_ex2=0.0
		rhce_cnv_dp=0.0
		rh_typing_result=[]
		rhce_typing_result="exon2nochange"
		global exon2_rhce
		for line in f1.readlines():
			if line.startswith("chromosome"):
				continue
			else:
				ex_start=line.split("\t")[1]
				ex_end = line.split("\t")[2]
				gene_name = line.split("\t")[3]
				dp = float(line.split("\t")[4])
				log2=line.split("\t")[5]
				if gene_name =="RSRP1,RHD":
					rhd_cnv_dp = rhd_cnv_dp+dp
				if gene_name =="RHCE":
					rhce_cnv_dp= rhce_cnv_dp+dp
				#0.585 #-1 #215估计是215个bin位置
				#对于拷贝数的无需要区分纯杂合缺失，只要有的话就爆出来就行了因为至少一个染色体发生了，看能否后期有需要的话加一个频率，此频率和AF计算方法不一样
				if ex_start=="25272509" and ex_end =="25272695":
						if float(log2)>0.585:
							exon1_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon1_rhd="loss"
							continue
						else:
							exon1_rhd="middle"
							continue
				if ex_start=="25284573" and ex_end =="25284759":
						if float(log2)>0.585:
							exon2_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon2_rhd="loss"
							continue
						else:
							exon2_rhd="middle"
							continue
				if ex_start=="25290641" and ex_end =="25290791":
						if float(log2)>0.585:
							exon3_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon3_rhd="loss"
							continue
						else:
							exon3_rhd="middle"
							continue
				if ex_start=="25300946" and ex_end =="25301093":
						if float(log2)>0.585:
							exon4_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon4_rhd="loss"
							continue
						else:
							exon4_rhd="middle"
							continue
				if ex_start=="25301520" and ex_end =="25301686":
						if float(log2)>0.585:
							exon5_rhd="gain"
						elif float(log2)<-1.0:
							exon5_rhd="loss"					
						else:
							exon5_rhd="middle"
				if ex_start=="25303322" and ex_end =="25303459":
						if float(log2)>0.585:
							exon6_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon6_rhd="loss"	
							continue
						else:
							exon6_rhd="middle"
							continue
				if ex_start=="25306596" and ex_end =="25306729":
						if float(log2)>0.585:
							exon7_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon7_rhd="loss"						
							continue
						else:
							exon7_rhd="middle"
							continue
				if ex_start=="25317000" and ex_end =="25317079":
						if float(log2)>0.585:
							exon8_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon8_rhd="loss"	
							continue
						else:
							exon8_rhd="middle"
							continue
				if ex_start=="25321889" and ex_end =="25321962":
						if float(log2)>0.585:
							exon9_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon9_rhd="loss"	
							continue
						else:
							exon9_rhd="middle"
							continue
				if ex_start=="25328898" and ex_end =="25330445":
						if float(log2)>0.585:
							exon10_rhd="gain"
							continue
						elif float(log2)<-1.0:
							exon10_rhd="loss"						
							continue
						else:
							exon10_rhd="middle"
							continue
				if ex_start=="25420639" and ex_end =="25420825":
						if float(log2)>0.585:
							exon1_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon1_rhce="loss"
							continue
						else:
							exon1_rhce="middle"
							continue
				#print("ok1")
				if ex_start=="25408683" and ex_end =="25408869":
						if float(log2)>0.585:
							exon2_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon2_rhce="loss"
							continue
						else:
							exon2_rhce="middle"
							continue
				#print("ok2")
				if ex_start=="25402596" and ex_end =="25402746":
						if float(log2)>0.585:
							exon3_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon3_rhce="loss"
							continue
						else:
							exon3_rhce="middle"
							continue
				if ex_start=="25391994" and ex_end =="25392141" :
						if float(log2)>0.585:
							exon4_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon4_rhce="loss"
							continue
						else:
							exon4_rhce="middle"
							continue
				if ex_start=="25390749" and ex_end =="25390915":
						if float(log2)>0.585:
							exon5_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon5_rhce="loss"
							continue
						else:
							exon5_rhce="middle"
							continue
				if ex_start=="25388976" and ex_end =="25389113" :
						if float(log2)>0.585:
							exon6_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon6_rhce="loss"
							continue
						else:
							exon6_rhce="middle"
							continue
				if ex_start=="25385711"and ex_end =="25385844" :
						if float(log2)>0.585:
							exon7_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon7_rhce="loss"
							continue
						else:
							exon7_rhce="middle"
							continue
				if ex_start=="25375349" and ex_end =="25375428" :
						if float(log2)>0.585:
							exon8_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon8_rhce="loss"
							continue
						else:
							exon8_rhce="middle"
							continue
				if ex_start=="25370467" and ex_end =="25370540":
						if float(log2)>0.585:
							exon9_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon9_rhce="loss"
							continue
						else:
							exon9_rhce="middle"
							continue
				if ex_start=="25362249"and ex_end =="25362553":
						if float(log2)>0.585:
							exon10_rhce="gain"
							continue
						elif float(log2)<-1.0:
							exon10_rhce="loss"
							continue
						else:
							exon10_rhce="middle"
							continue
			#RHD*01
		#if  1.6<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2<2.5:
		if  1.6<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2:
				rh_typing_result.append("RHD:D+(RHD*01)\t")
			#RHD*01N.01 (RHD的外显子全部缺失D-)
		if  0.6<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2<1.5:
				rh_typing_result.append("RHD:D+(RHD*01)\tRHD*01\nRHD:D-(RHD*01N.01)\tRHD*01N.01\tRHD deletion\tRHD deletion\texonic\tSV\t.\tRHD deletion\tISBT\tRHD deletion")
			#RHD*01N.01 (RHD的外显子全部缺失D-)
		if  0.0<((rhd_cnv_dp/10)/(rhce_cnv_dp/10))*2<0.5:
				rh_typing_result.append("RHD:D-(RHD*01N.01)\tRHD*01N.01\tRHD deletion\tRHD deletion\texonic\tSV\t.\tRHD deletion\tISBT\tRHD deletion")
		#RHD*01.01
	#	if  1.6<((rhd_cnv_dp/215)/(rhce_cnv_dp/217))*2<2.5 and g.25272595G>C:
	#			rh_typing_result="RHD:D+(RHD*01)\t"
		print (exon2_rhce)
		#RHD*03.02 (RHCE的2号外显子取代RHD的2号外显子，形成杂合蛋白)
		if exon2_rhce=="gain" and exon2_rhd=="loss":
				rh_typing_result.append("RHD:DIIIb Caucasian(RHD*03.02)\tRHD*03.02\tRHD*D-CE(2)-D\tRHD*D-CE(2)-D\texonic\tCNV\t.\tRHD*D-CE(2)-D\tISBT\tRHD*D-CE(2)-D")
		#RHD*03.03 (RHCE的3号外显子取代RHD的3号外显子，形成杂合蛋白)
		if exon3_rhce=="gain" and exon3_rhd=="loss":
				rh_typing_result.append("RHD:DIIIc(RHD*03.03)\tRHD*03.03\tRHD*D-CE(3)-D\tRHD*D-CE(3)-D\texonic\tCNV\t.\tRHD*D-CE(3)-D\tISBT\tRHD*D-CE(3)-D")
		#RHD*03.07 (RHCE的3号外显子取代RHD的3号外显子+4个SNP，形成杂合蛋白)

		#RHD*04.03 (RHCE的6号7号8号9号外显子，形成杂合蛋白)
		if exon6_rhce=="gain" and exon7_rhce=="gain" and exon8_rhce=="gain" and exon9_rhce=="gain"and exon6_rhd =="loss"and exon7_rhd  =="loss"and exon8_rhd =="loss"and exon9_rhd =="loss":
				rh_typing_result.append("RHD:DIV type 3(RHD*04.03)\tRHD*04.03\tRHD*D-CE(6-9)-D\tRHD*D-CE(6-9)-D\texonic\tCNV\t.\tRHD*D-CE(6-9)-D\tISBT\tRHD*D-CE(6-9)-D")
		#RHD*04.05 (RHCE的7号8号9号外显子，形成杂合蛋白)
		if exon7_rhce=="gain" and exon8_rhce=="gain" and exon9_rhce=="gain"and exon7_rhd  =="loss"and exon8_rhd  =="loss"and exon9_rhd =="loss":
				rh_typing_result.append("RHD:DIV type 5(RHD*04.05)\tRHD*04.05\tRHD*D-CE(7-9)-D\tRHD*D-CE(7-9)-D\texonic\tCNV\t.\tRHD*D-CE(7-9)-D\tISBT\tRHD*D-CE(7-9)-D")

		#RHD*05.02(RHCE的5号外显子取代RHD的5号外显子，形成杂合蛋白)
		if exon5_rhce=="gain" and exon5_rhd=="loss":
				rh_typing_result.append("RHD:DV type 2(RHD*05.02)\tRHD*05.02\tRHD*D-CE(5)-D\tRHD*D-CE(5)-D\texonic\tCNV\t.\tRHD*D-CE(5)-D\tISBT\tRHD*D-CE(5)-D")
				rh_typing_result.append("RHD:DBS1(RHD*13.01)\tRHD*13.01\tRHD*D-CE(5)-D\tRHD*D-CE(5)-D\texonic\tCNV\t.\tRHD*D-CE(5)-D\tISBT\tRHD*D-CE(5)-D")

		#RHD*05.10(RHCE的5号6号外显子，形成杂合蛋白)
		if exon5_rhce=="gain" and  exon6_rhce=="gain"and exon5_rhd=="loss" and exon6_rhd=="loss" :
				rh_typing_result.append("RHD:DIIIc(RHD*05.02)\tRHD*05.02\tRHD*D-CE(5-6)-D\tRHD*D-CE(5-6)-D\texonic\tCNV\t.\tRHD*D-CE(5-6)-D\tISBT\tRHD*D-CE(5-6)-D")
				rh_typing_result.append("RHD:DIIIc(RHD*46)\tRHD*46\tRHD*D-CE(5-6)-D\tRHD*D-CE(5-6)-D\texonic\tCNV\t.\tRHD*D-CE(5-6)-D\tISBT\tRHD*D-CE(5-6)-D")

		#RHD*06.01 (RHCE的4号5号外显子，形成杂合蛋白)
		if exon4_rhce=="gain" and  exon5_rhce=="gain"and exon4_rhd=="loss" and exon5_rhd=="loss" :
				rh_typing_result.append("RHD:DVI type 1(RHD*06.01)\tRHD*06.01\tRHD*D-CE(4-5)-D\tRHD*D-CE(4-5)-D\texonic\tCNV\t.\tRHD*D-CE(4-5)-D\tISBT\tRHD*D-CE(4-5)-D")

		#RHD*06.02  (RHCE的4号5号6号外显子，形成杂合蛋白)
		if exon4_rhce=="gain" and  exon5_rhce=="gain"and exon6_rhce=="gain" and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss" :
				rh_typing_result.append("RHD:DVI type 2(RHD*06.02)\tRHD*06.02\tRHD*D-CE(4-6)-D\tRHD*D-CE(4-6)-D\texonic\tCNV\t.\tRHD*D-CE(4-6)-D\tISBT\tRHD*D-CE(4-6)-D")

		#RHD*06.03  (RHCE的3号4号5号6号外显子，形成杂合蛋白)
		if exon3_rhce =="gain"and exon4_rhce=="gain" and  exon5_rhce=="gain"and exon6_rhce=="gain" and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss" :
				rh_typing_result.append("RHD:DVI type 3(RHD*06.03)\tRHD*06.03\tRHD*D-CE(3-6)-D\tRHD*D-CE(3-6)-D\texonic\tCNV\t.\tRHD*D-CE(3-6)-D\tISBT\tRHD*D-CE(3-6)-D")

		#RHD*06.03.02 (RHCE的3号4号5号6号外显子加c.1195A突变，形成杂合蛋白，单个点的暂未加入)

		#RHD*06.04  (RHCE的3号4号5号外显子，形成杂合蛋白)
		if exon3_rhce =="gain"and exon4_rhce=="gain" and  exon5_rhce=="gain" and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss"  :
				rh_typing_result.append("RHD:DVI type 4(RHD*06.04)\tRHD*06.04\tRHD*D-CE(3-5)-D\tRHD*D-CE(3-5)-D\texonic\tCNV\t.\tRHD*D-CE(3-5)-D\tISBT\tRHD*D-CE(3-5)-D")

		#RHD*13.01  (RHCE的5号外显子取代RHD的5号外显子，形成杂合蛋白,重复了，到时候一起输出)
		
		#RHD*14.01  (RHCE的5号6号7号外显子，形成杂合蛋白)
		if exon5_rhce == "gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon5_rhd=="loss"and exon6_rhd=="loss" and exon7_rhd=="loss"  :
				rh_typing_result.append("RHD:DBT1(RHD*14.01)\tRHD*14.01\tRHD*D-CE(5-7)-D\tRHD*D-CE(5-7)-D\texonic\tCNV\t.\tRHD*D-CE(5-7)-D\tISBT\tRHD*D-CE(5-7)-D")
		#RHD*14.02  (RHCE的5号6号7号8号9号外显子，形成杂合蛋白)
		if exon5_rhce =="gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon8_rhce=="gain" and exon9_rhce=="gain"and exon5_rhd=="loss"and exon6_rhd=="loss" and exon7_rhd=="loss" and exon8_rhd=="loss"  and exon9_rhd=="loss" :
				rh_typing_result.append("RHD:DBT2(RHD*14.02)\tRHD*14.02\tRHD*D-CE(5-9)-D\tRHD*D-CE(5-9)-D\texonic\tCNV\t.\tRHD*D-CE(5-9)-D\tISBT\tRHD*D-CE(5-9)-D")
		#RHD*17.02	(RHCE的4号外显子取代RHD的4号外显子，形成杂合蛋白)
		if exon3_rhce == "gain"and exon4_rhce=="gain" and  exon5_rhce=="gain" and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss"  :
				rh_typing_result.append("RHD:DFR2(RHD*17.02)\tRHD*17.02\tRHD*D-CE(4)-D\tRHD*D-CE(4)-D\texonic\tCNV\t.\tRHD*D-CE(4)-D\tISBT\tRHD*D-CE(4)-D")
		#RHD*17.05  (RHCE的3号4号外显子，形成杂合蛋白)
		if exon3_rhce == "gain"and exon4_rhce=="gain"  and exon3_rhd=="loss"and exon4_rhd=="loss":
				rh_typing_result.append("RHD:DFR5(RHD*17.05)\tRHD*17.05\tRHD*D-CE(3-4)-D\tRHD*D-CE(3-4)-D\texonic\tCNV\t.\tRHD*D-CE(3-4)-D\tISBT\tRHD*D-CE(3-4)-D")
		#RHD*41 (RHCE的5号6号7号外显子，形成杂合蛋白)
		if exon5_rhce == "gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon5_rhd=="loss"and exon6_rhd=="loss" and exon7_rhd=="loss"  :
				rh_typing_result.append("RHD:DBU(RHD*41)\tRHD*41\tRHD*D-CE(5-7)-D\tRHD*D-CE(5-7)-D\texonic\tCNV\t.\tRHD*D-CE(5-7)-D\tISBT\tRHD*D-CE(5-7)-D")
		#RHD*45 (RHCE的2号3号外显子，形成杂合蛋白)
		if exon2_rhce == "gain"and exon3_rhce=="gain" and exon2_rhd=="loss"and exon3_rhd=="loss" :
				rh_typing_result.append("RHD:DKK(RHD*45)\tRHD*45\tRHD*D-CE(2-3)-D\tRHD*D-CE(2-3)-D\texonic\tCNV\t.\tRHD*D-CE(2-3)-D\tISBT\tRHD*D-CE(2-3)-D")
		#RHD*58 (RHCE的7号外显子，形成杂合蛋白)
		if exon7_rhce == "gain" and exon7_rhd=="loss":
				rh_typing_result.append("RHD:RHD*D-CE(7)-D(RHD*58)\tRHD*58\tRHD*D-CE(7)-D\tRHD*D-CE(7)-D\texonic\tCNV\t.\tRHD*D-CE(7)-D\tISBT\tRHD*D-CE(7)-D")
		#RHD*01EL.30 (RHCE的8号外显子的缺失) 
		if exon3_rhce == "gain"and exon4_rhce=="gain" and  exon5_rhce=="gain" and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss"  :
				rh_typing_result.append("RHD:DIIIc(RHD*01EL.30)\tRHD*01EL.30\tRHD*D-CE(8)-D\tRHD*D-CE(8)-D\texonic\tCNV\t.\tRHD*D-CE(8)-D\tISBT\tRHD*D-CE(8)-D")
		#RHD*01EL.44 (RHCE的4号5号6号7号8号9号外显子，形成杂合蛋白)
		if exon4_rhce == "gain"and exon5_rhce=="gain" and  exon6_rhce=="gain" and  exon7_rhce=="gain"and  exon8_rhce=="gain"and  exon9_rhce=="gain"and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss"and exon7_rhd=="loss" and exon8_rhd=="loss" and exon9_rhd == "loss" :
				rh_typing_result.append("RHD:Del44(RHD*01EL.44)\tRHD*01EL.44\tRHD*D-CE(4-9)-D\tRHD*D-CE(4-9)-D\texonic\tCNV\t.\tRHD*D-CE(4-9)-D\tISBT\tRHD*D-CE(4-9)-D")

		#RHD*01N.02 （CE exons 1-9）
		if exon1_rhce == "gain"and exon2_rhce=="gain" and  exon3_rhce=="gain" and  exon4_rhce=="gain"and  exon5_rhce=="gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon8_rhce=="gain" and  exon9_rhce=="gain"and exon1_rhd=="loss"and exon2_rhd=="loss"and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss"and exon7_rhd=="loss" and exon8_rhd=="loss" and exon9_rhd == "loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.02 )\tRHD*01N.02 \tRHD*D-CE(1-9)-D\tRHD*D-CE(1-9)-D\texonic\tCNV\t.\tRHD*D-CE(1-9)-D\tISBT\tRHD*D-CE(1-9)-D")

		#RHD*01N.03 （CE exons 2-9）
		if  exon2_rhce=="gain" and  exon3_rhce=="gain" and  exon4_rhce=="gain"and  exon5_rhce=="gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon8_rhce=="gain"and  exon9_rhce=="gain"and exon2_rhd=="loss"and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss"and exon7_rhd=="loss" and exon8_rhd=="loss" and exon9_rhd == "loss" :
				rh_typing_result.append("RHD:D-(#RHD*01N.03)\t#RHD*01N.03\tRHD*D-CE(2-9)-D\tRHD*D-CE(2-9)-D\texonic\tCNV\t.\tRHD*D-CE(2-9)-D\tISBT\tRHD*D-CE(2-9)-D")
		
		#RHD*01N.04 （CE exons 3-9）
		if  exon3_rhce=="gain" and  exon4_rhce=="gain"and  exon5_rhce=="gain"and exon6_rhce=="gain" and  exon7_rhce=="gain" and exon8_rhce=="gain" and exon9_rhce=="gain"and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss"and exon7_rhd=="loss" and exon8_rhd=="loss" and exon9_rhd == "loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.04 )\tRHD*01N.04 \tRHD*D-CE(3-9)-D\tRHD*D-CE(3-9)-D\texonic\tCNV\t.\tRHD*D-CE(3-9)-D\tISBT\tRHD*D-CE(3-9)-D")

		#RHD*01N.05 （CE exons 2-7）
		if  exon2_rhce=="gain"and exon3_rhce=="gain" and  exon4_rhce=="gain"and  exon5_rhce=="gain"and  exon6_rhce=="gain" and  exon7_rhce=="gain" and exon2_rhd=="loss"and exon3_rhd=="loss"and exon4_rhd=="loss" and exon5_rhd=="loss" and exon6_rhd=="loss"and exon7_rhd=="loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.05)\tRHD*01N.05\tRHD*D-CE(2-7)-D\tRHD*D-CE(2-7)-D\texonic\tCNV\t.\tRHD*D-CE(2-7)-D\tISBT\tRHD*D-CE(2-7)-D")

		#RHD*01N.06 （CE exons 4-7 with 733C>G ，1006G>T）
		#RHD*01N.07 （CE exons 4-7）
		if exon4_rhce=="gain"and  exon5_rhce=="gain"and  exon6_rhce=="gain"and  exon7_rhce=="gain"and exon4_rhd=="loss"and exon5_rhd=="loss" and exon6_rhd=="loss" and exon7_rhd=="loss"and exon7_rhd=="loss" and exon8_rhd=="loss" and exon9_rhd == "loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.07)\tRHD*01N.07\tRHD*D-CE(4-7)-D\tRHD*D-CE(4-7)-D\texonic\tCNV\t.\tRHD*D-CE(4-7)-D\tISBT\tRHD*D-CE(4-7)-D")

		#RHD*01N.42 （CE exons1, 7-10）
		if  exon1_rhce=="gain" and  exon7_rhce=="gain"and  exon10_rhce=="gain"and exon1_rhd=="loss"and exon7_rhd=="loss" and exon10_rhd=="loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.42)\tRHD*01N.42\tRHD*D-CE(1,7-10)-D\tRHD*D-CE(1,7-10)-D\texonic\tCNV\t.\tRHD*D-CE(1,7-10)-D\tISBT\tRHD*D-CE(1,7-10)-D")

		#RHD*01N.43 （CE exons 1-3）
		if  exon1_rhce=="gain" and  exon2_rhce=="gain"and  exon3_rhce=="gain"and   exon1_rhd=="loss"and exon2_rhd=="loss" and exon3_rhd=="loss" :
				rh_typing_result.append("RHD:D-(RHD*01N.43)\tRHD*01N.43\tRHD*D-CE(1-3)-D\tRHD*D-CE(1-3)-D\texonic\tCNV\t.\tRHD*D-CE(1-3)-D\tISBT\tRHD*D-CE(1-3)-D")



	
		
		#if  rhd_cnv_ex2/3>0.585 and rhce_cnv_ex2/2<-1.0:
		#		rhce_typing_result="exon2change"
		#		rh_typing_result=rh_typing_result+"C+"
		#else:
		#	rh_typing_result=rh_typing_result+"C-"

		return rh_typing_result
def RHD_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/RHD_Typer_database_hg38_v3.0.01'
	RHD_Typer_database_l = txt_to_lines(database_path)
	RHD_Typer_database_l = RHD_Typer_database_l[1:]
	Haplotypes = RHD_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Haplotypes2= RHD_Haplotype_maker2(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.RHD_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in RHD_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.write("\n".join(Haplotypes2))
	result2.close()

def RHCE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	#kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	hg38_cnv = abpath + "/reference/hg38/hg38.fa"
	hg38fai_cnv = abpath + "/reference/hg38/hg38.fa.fai"
	bedfile = abpath + "/input/RHCE_bed"
	bed_cnv_file = abpath + "/input/CNV_rhd_rche.bed"
	cnvkit = abpath + "/software/soft/cnvkit/cnvkit.py"
	annotate = abpath + "/input/refFlat.txt"

	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
	
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHCE_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RHCE_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
					if line.startswith("chr"):
						pos_chr = re.split(r"\.|\t",line)[0]
						pos_start = re.split(r"\.|\t",line)[1]
						pos_end = re.split(r"\.|\t",line)[2]
						dp = re.split(r"\.|\t",line)[5]
						ad = re.split(r"\.|\t",line)[6]
						af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
						#print(line+af)
						Haplotypes_temp = {}
						dic_temp = {}
						for chr1,chr2 in Haplotypes.items():			
							if pos_start in kb500:
								pass
								'''
								if af <= 0.1 :#纯合无突变
									chr11 = chr1 + 'N'
									chr21 = chr2 + 'N'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)							
								if af >= 0.8 :#纯合突变
									chr11 = chr1 + 'Y'
									chr21 = chr2 + 'Y'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)	
								if 0.1<af<0.8 :#杂合突变
									if kb500 in another_ID:
										if another_ID <=0.1
											chr11 = chr1 + 'N'
											chr21 = chr1 + 'Y'
											chr22 = chr2 + 'N'
											chr12 = chr2 + 'Y'
								'''
							else:
								if float(af) <= 0.1 :#纯合无突变
									chr11 = chr1 + 'N'
									chr21 = chr2 + 'N'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)							
								if float(af) >= 0.8 :#纯合突变
									chr11 = chr1 + 'Y'
									chr21 = chr2 + 'Y'
									dic_temp = {chr11:chr21}
									Haplotypes_temp.update(dic_temp)	
								if 0.1<float(af)<0.8 :#杂合突变
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
									#print(chr11+"ok"+chr12)
									dic_temp = {chr11:chr12,chr21:chr22}
									Haplotypes_temp.update(dic_temp)
						Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes
def RHCE_Haplotype_maker2(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	#Nou+  Dav+ MAR+ 
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	
	#kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	#class 1（）#mate reads上的顺反式
	#class 2（）#同一个reads上的顺反式

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	hg38_cnv = abpath + "/reference/hg38/hg38.fa"
	hg38fai_cnv = abpath + "/reference/hg38/hg38.fa.fai"
	bedfile = abpath + "/input/RHCE_bed"
	bed_cnv_file = abpath + "/input/CNV_rhd_rche.bed"
	cnvkit = abpath + "/software/soft/cnvkit/cnvkit.py"
	annotate = abpath + "/input/refFlat.txt"
	
	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_cnv'+"/"+sample_name+".deduped.cnr","r") as f1:
			rhd_cnv_ex2=0.0
			rhd_cnv_dp=0.0
			rhce_cnv_ex2=0.0
			rhce_cnv_dp=0.0
			rh_typing_result=""
			rhce_typing_result="exon2nochange"
			for line in f1.readlines():
				if line.startswith("chromosome"):
					continue
				else:
					ex_start=line.split("\t")[1]
					ex_end = line.split("\t")[2]
					gene_name = line.split("\t")[3]
					dp = float(line.split("\t")[4])
					log2=line.split("\t")[5]
					#0.585 #-1 #215估计是215个位置吧
					if ex_start=="25284109" and ex_end =="25284376":
						rhd_cnv_ex2=rhd_cnv_ex2+float(log2)
					if ex_start=="25284376" and ex_end =="25284642":
						rhd_cnv_ex2=rhd_cnv_ex2+float(log2)
					if ex_start=="25284642" and ex_end =="25284908":
						rhd_cnv_ex2=rhd_cnv_ex2+float(log2)
					if ex_start=="25408610" and ex_end =="25408876":
						rhce_cnv_ex2=rhce_cnv_ex2+float(log2)
					if ex_start=="25408876" and ex_end =="25409142":
						rhce_cnv_ex2=rhce_cnv_ex2+float(log2)
					if gene_name =="RSRP1,RHD":
						rhd_cnv_dp = rhd_cnv_dp+dp
					if gene_name =="RHCE":
						rhce_cnv_dp= rhce_cnv_dp+dp
			if  ((rhce_cnv_ex2/2)/(rhce_cnv_dp/217))*2<0.5:#exon2缺失两个
					rh_typing_result=rh_typing_result+"RHCE:C+c-\n"
			if  0.5<=((rhce_cnv_ex2/2)/(rhce_cnv_dp/217))*2<1.4:#exon2缺失一个
					rh_typing_result=rh_typing_result+"RHCE:C+c+\n"
			if  1.4<=((rhce_cnv_ex2/2)/(rhce_cnv_dp/217))*2<1.4:#exon2没有缺失
					rh_typing_result=rh_typing_result+"RHCE:C-c+\n"						
			#if  rhd_cnv_ex2/3>0.585 and rhce_cnv_ex2/2<-1.0:
			#		rhce_typing_result="exon2change"
			#		rh_typing_result=rh_typing_result+"C+"
			#else:
			#	rh_typing_result=rh_typing_result+"C-"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
	
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RH_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".RH_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		af_dic={}
		#c.307C>T(g.25408711G>A) c.676G>C(g.25390874C>G)
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				pos_alt = re.split(r"\.|\t",line)[4]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print (af)
				af_dic[pos_start+"_"+pos_end+"_"+pos_alt]=af
	#	if float(af_dic["25408711_25408711_G"])>0.0 and rhce_typing_result=="exon2nochange":
	#			rh_typing_result=rh_typing_result+"c+"
	#	else:
	#		rh_typing_result=rh_typing_result+"c-"
		if float(af_dic["25390874_25390874_C"])>0.2:
			rh_typing_result=rh_typing_result+"RHCE:E+\n"
		else:
			rh_typing_result=rh_typing_result+"RHCE:E-\n"
		if float(af_dic["25390874_25390874_G"])>0.2:
			rh_typing_result=rh_typing_result+"RHCE:e+\n"	
		else:
			rh_typing_result=rh_typing_result+"RHCE:e-\n"
	#	print(Haplotypes)
		return rh_typing_result	
def RHCE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/RHCE_Typer_database_hg38_v3.0.01'
	RHCE_Typer_database_l = txt_to_lines(database_path)
	RHCE_Typer_database_l = RHCE_Typer_database_l[1:]
	Haplotypes = RHCE_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Haplotypes2 = RHCE_Haplotype_maker2(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.RHCE_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in RHCE_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]+"("+l3_list[1]+")"
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.write(Haplotypes2)
	result2.close()

def GPIIIa_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GPIIIa_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIIIa_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIIIa_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GPIIIa_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GPIIIa_Typer_database_hg38_v3.0.01'
	GPIIIa_Typer_database_l = txt_to_lines(database_path)
	GPIIIa_Typer_database_l = GPIIIa_Typer_database_l[1:]
	Haplotypes = XK_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GPIIIa_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GPIIIa_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def GPIbα_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GPIbα_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIbα_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIbα_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GPIbα_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GPIbα_Typer_database_hg38_v3.0.01'
	GPIbα_Typer_database_l = txt_to_lines(database_path)
	GPIbα_Typer_database_l = GPIbα_Typer_database_l[1:]
	Haplotypes = GPIbα_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GPIbα_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GPIbα_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def GPIIb_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GPIIb_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIIb_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIIb_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GPIIb_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GPIIb_Typer_database_hg38_v3.0.01'
	GPIIb_Typer_database_l = txt_to_lines(database_path)
	GPIIb_Typer_database_l = GPIIb_Typer_database_l[1:]
	Haplotypes = GPIIb_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GPIIb_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GPIIb_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def GPIa_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GPIa_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIa_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIa_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GPIa_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GPIa_Typer_database_hg38_v3.0.01'
	GPIa_Typer_database_l = txt_to_lines(database_path)
	GPIa_Typer_database_l = GPIa_Typer_database_l[1:]
	Haplotypes = GPIa_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GPIa_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GPIa_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def GPIbβ_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/GPIbβ_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIbβ_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".GPIbβ_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def GPIbβ_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/GPIbβ_Typer_database_hg38_v3.0.01'
	GPIbβ_Typer_database_l = txt_to_lines(database_path)
	GPIbβ_Typer_database_l = GPIbβ_Typer_database_l[1:]
	Haplotypes = GPIbβ_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.GPIbβ_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in GPIbβ_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def CD109_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	#all_postions_path = abpath + '/input/all_postions.txt' 
	#all_postions_path = abpath + '/input/GYPA_allposition.txt'
	#all_varriant_lines = txt_to_lines(all_postions_path)
	#varriant_result_lines = txt_to_lines(convert_result_path)
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	#print(BAM_path)
	'''
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]
	#所有位点的纯杂合情况以及500情况先存于一个字典，然后再循环组成NY字符串
	class 1（）#mate reads上的顺反式
	class 2（）#同一个reads上的顺反式
	with open ("QYW_25.deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f:
		for line in f.readlines():
			if line.startswith("chr"):
				comp=line.split("\t")[0]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				if af <= 0.1:
					all_pos[comp]=N
				if af >= 0.8:
					all_pos[comp]=Y
				if 0.1<af<0.8:
					all_pos[comp]=S	
	'''
	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
	bedfile = abpath + "/input/CD109_bed"
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)
		
	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CD109_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd) 

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".CD109_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				#print(line+af)
				Haplotypes_temp = {}
				dic_temp = {}
				for chr1,chr2 in Haplotypes.items():			
					if pos_start in kb500:
						pass
						'''
						if af <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if af >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<af<0.8 :#杂合突变
							if kb500 in another_ID:
								if another_ID <=0.1
									chr11 = chr1 + 'N'
									chr21 = chr1 + 'Y'
									chr22 = chr2 + 'N'
									chr12 = chr2 + 'Y'
						'''
					else:
						if float(af) <= 0.1 :#纯合无突变
							chr11 = chr1 + 'N'
							chr21 = chr2 + 'N'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)							
						if float(af) >= 0.8 :#纯合突变
							chr11 = chr1 + 'Y'
							chr21 = chr2 + 'Y'
							dic_temp = {chr11:chr21}
							Haplotypes_temp.update(dic_temp)	
						if 0.1<float(af)<0.8 :#杂合突变
							chr11 = chr1 + 'N'
							chr21 = chr1 + 'Y'
							chr22 = chr2 + 'N'
							chr12 = chr2 + 'Y'
							#print(chr11+"ok"+chr12)
							dic_temp = {chr11:chr12,chr21:chr22}
							Haplotypes_temp.update(dic_temp)
				Haplotypes = Haplotypes_temp	
		print(Haplotypes)
		return Haplotypes					
def CD109_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	print(sys.path[0])
	abpath = sys.path[0]

	#database_path = abpath + '/input/ABO_Typer_database_v1.0.txt'
	database_path = abpath + '/input/CD109_Typer_database_hg38_v3.0.01'
	CD109_Typer_database_l = txt_to_lines(database_path)
	CD109_Typer_database_l = CD109_Typer_database_l[1:]
	Haplotypes = CD109_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.CD109_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		#chr1 = chr1[2:]#?
		#chr2 = chr2[2:]
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		#print(len(chr1))
		for l3 in CD109_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
			#print("Strings_REF_len"+str(len(Strings_REF)))
			type_name = l3_list[0]
			str_index = 0
			Y_num = 0 
			chr1_Y_match = 0
			chr2_Y_match = 0  
			for S in Strings_REF:#？
				if S == 'Y':
					if chr1[str_index] != 'Y':
						match_test = False
					if chr2[str_index] != 'Y':
						match_test2 = False
				#if str_index in important_alles_index :
				#	if chr1[str_index] != S:
				#		match_test = False
				#	if chr2[str_index] != S:
				#		match_test2 = False
				str_index += 1
			if match_test:
				chr1 = chr1 + '/' + type_name
			if match_test2:
				chr2 = chr2 + '/' + type_name
		Haplotype_type_result[chr1] = chr2
	for S1,S2 in Haplotype_type_result.items():#含有所有可能的分型结果输出
		LL = S1 + '\t' +  S2 + '\n'
		result2.write(LL)
	result2.close()

def typing_summary(OUTPUT,sample_name) :
	Types_result_path = OUTPUT+'varResult/' + sample_name+'/'
	files = os.listdir(Types_result_path)
	with open(OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_typing_results.txt','a')as q:
		q.write('表型结果' + '\n')
	for file in files:
		txt_path=Types_result_path+file
		if re.findall(r"Type_result.txt",file) !=[]:#遍历每一个Type_result.txt
			try:
				blood=re.search(r"(QYW_\d+\.)(\w[^_]+|\w)(_Type_result.txt)",file).group(2)
			except:
				print("what"+file)
			if blood!="GPIa" and blood!="GPIbα" and blood!="GPIbβ" and blood!="GPIIb" and blood!="GPIIIa" and blood!="CD109":
				with open (txt_path,'r')as txt,open(OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_typing_results.txt','a')as q:
						txt=txt.readlines()
						pheno={}
						for line in txt:#读取每个血型分型结果的每行内容
								#if re.findall(r"、+",line) ==[]:
								if line.startswith(" "):#血型结果为空的直接跳过
									q.write(blood+"\t?")
								else:#ABO、RH、MNS、其它，共分成四类
									if blood == "ABO": #or blood == "RH":
										q.write(blood+":"+line+"\n")
									else:
										for i in re.findall(blood+r":[^/|^\t]+",line.strip()):
											if i in pheno.keys():
												continue
											else:
												pheno[i]="1"											
						#if blood == "MNS":
						for i in pheno.keys():
								q.write(i+"\n")		
			else:		
				with open (txt_path,'r')as txt,open(OUTPUT+'varResult/' + sample_name+'/'+sample_name+'_typing_results.txt','a')as q:
						txt=txt.readlines()
						pheno={}
						for line in txt:#读取每个血型分型结果的每行内容
								#if re.findall(r"、+",line) ==[]:
								if line.startswith(" "):#血型结果为空的直接跳过
									q.write(blood+"\t?")
								else:#ABO、RH、MNS、其它，共分成四类
									if blood == "ABO" or blood == "RH":
										q.write(blood+":"+line+"\n")
									else:
										chr1_line = line.split("\t")[0]
										chr2_line = line.split("\t")[1]
										for i in re.findall(r"HPA:[^/|^\t]+",line.strip()):
											if i in pheno.keys():
												continue
											else:
												pheno[i]="1"											
						#if blood == "MNS":
						for i in pheno.keys():
								q.write(i+"\n")								

def main(INPUT, OUTPUT, THREAD,CLASS):
	if CLASS=="ALL":																			
		bloodtyper_start_time = datetime.datetime.now()
		print('BTyper start is ' + str(bloodtyper_start_time))
		SampleNames = GetSampleName(INPUT)#得到输入文件夹内的样本名
		for sample_name in SampleNames :
			result_path = OUTPUT + '/varResult/' + sample_name
			#'''
			try:
				shutil.rmtree(result_path)#递归删除之前的文件夹
			except Exception as e:
				pass
			os.makedirs(result_path)
			#'''
		sample_dataprocess_start_time = datetime.datetime.now()
		print(' ' + sample_name +' dataprocess start :' + str(sample_dataprocess_start_time) )
		#dataprocess(INPUT,OUTPUT,THREAD,sample_name,CLASS)#开始每个样本跑流程
		sample_dataprocess_end_time = datetime.datetime.now()
		print(' ' + sample_name +' dataprocess done :' + str(sample_dataprocess_end_time) )
		path1 = OUTPUT + '/varResult/'
		sample_list = os.listdir(path1)
		for sample_name in sample_list :
			sample_Typing_start_time = datetime.datetime.now()
			print('   ' + sample_name +' Typing start :' + str(sample_Typing_start_time) )
			convert_result_suffix = '.hg38_multianno.txt'
			convert_result_path = find_file(sample_name,path1,convert_result_suffix)
			BAM_path = convert_result_path.replace('.hg38_multianno.txt','.deduped.bam')#替代字符串找到bam的位置 传给get_amplicon_txt函数
			print(BAM_path)
			#get_amplicon_txt(BAM_path)#得到所有比对到目标位置的bam序列
			#convert1(convert_result_path)
			#'''
			try:
				ABO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("ABO-ok")
				RHD_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RHD-ok")			
				RHCE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RHCE-ok")
				GLOB_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("GLOB-ok")
				MNS_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("MNS-ok")
				LU_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("LU-ok")
				KEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("KEL-ok")
				FY_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("FY-ok")
				JK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("JK-ok")
				DI_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("DI-ok")
				YT_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("YT-ok")
				SC_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("SC-ok")
				DO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("DO-ok")
				CO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("CO-ok")
				LW_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("LW-ok")
				CHRG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("CHRG-ok")
				CROM_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("CROM-ok")
				RHAG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RHAG-ok")
				KN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("KN-ok")
				IN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("IN-ok")
				OK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("OK-ok")
				RAPH_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RAPH-ok")
				JMH_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("JMH-ok")
				I_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("I-ok")
				GIL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("GIL-ok")
				FROS_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("FROS-ok")
				LAN_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("LAN-ok")
				VEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("VEL-ok")
				CD59_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("CD59-ok")
				AUG_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("AUG-ok")
				CTL2_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("CTL2-ok")
				GE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("GE-ok")
				H_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("H-ok")
				KANNO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("KANNO-ok")
				MAM_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("MAM-ok")
				P1PK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("P1PK-ok")
				SID_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("SID-ok")
				XK_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)	
				print("XK-ok")
				LE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("LE-ok")
				JR_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("JR-ok")
				PEL_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("PEL-ok")
				GPIIIa_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				GPIbα_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				GPIIb_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				GPIa_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				GPIbβ_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				CD109_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)	
			except Exception as e:
				print(sample_name + "typing error :" + str(e))
				continue
			
			#Type_result = path1 + sample_name + '/' + sample_name +"_"+ CLASS+'.Type_result.txt'
			#问题1:ABO_Typer_database_v1.0.txt的解释-TYPER_STR
			#'''
			sample_Typing_end_time = datetime.datetime.now()
			typing_summary(OUTPUT,sample_name) 
			#snapshot(OUTPUT, sample_name)
			#genarateBed(OUTPUT, sample_name)#得到所有ATCG碱基比例
			print('   ' + sample_name + 'typing done :' + str(sample_Typing_end_time))
		bloodtyper_end_time = datetime.datetime.now()
		Totally_consumed = '%0.2fsec'%(bloodtyper_end_time-bloodtyper_start_time).total_seconds()
		print('BTyper Totally consumed ' + Totally_consumed )

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Btyper",
		prog="BTyper.py",
		usage="python3 BTyper.py -i <input dir> -o <output dir> -t <threads> -c <Class>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.4 20210320")
	parser.add_argument("-i", "--input", type=str,
		help="Input the directory of the fq.gz file")
	parser.add_argument("-o", "--output", type=str,
		help="the output directory")
	parser.add_argument("-t", "--threads", type=int, default=4,
		help="threads, default=4")
	parser.add_argument("-c", "--Class", type=str, 
		help="class of freebays bed")		
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(INPUT=args.input, OUTPUT=args.output, THREAD=args.threads,CLASS=args.Class)

