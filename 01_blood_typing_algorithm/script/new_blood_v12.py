#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2020-12-31 14:19:59
# @Author  : info@genephar.com.cn
# @Link    : http://www.genephar.com.cn//
# @Version : Id
# @Telephone  : 89053423 or 89053788

import os, sys, re, random, time, argparse, datetime, shutil


def find_file(sample,path1,file_suffix):#查找目录文件

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

def txt_to_lines(file_path):#用于打开文件
	data = open(file_path)
	data_lines = data.readlines()
	data.close()
	return data_lines

def GetSampleName(INPUT):#用于得到样本名
	sample_names = []
	files = os.listdir(INPUT)
	for f in files :
		sample_name = f.split('.')[0]
		if sample_name not in sample_names :
			sample_names.append(sample_name)
	return sample_names

#四类算法函数
def Other_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT,Class):
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
	bedfile = abpath + "/input/"+Class+"_bed"
	#print(bedfile)
	cmd = """
		{s_v3} {input}{sample}.deduped.bam {bedfile} {hg38} {hg38fai} {output}>OKO
	""".format(s_v3=s_v3, hg38=hg38, hg38fai=hg38fai,bedfile=bedfile, sample=sample_name,input=INPUT,output=OUTPUT+'varResult/' + sample_name+'/'+sample_name+'.deduped.bam')
	print(cmd)
	os.system(cmd)		

	cmd = "mv %s %s"%(abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".deduped.bam_tumor_bmd.vcfDedup_mapping.txt",abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+"."+Class+"_deduped.bamtumor_bmd.vcfDedup_mapping.txt")
	os.system(cmd)		

	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+"."+Class+"_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
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
def Other_Haplotyp_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT,CLASS):
	print(sys.path[0])
	abpath = sys.path[0]
	database_path = abpath + '/input/'+CLASS+'_Typer_database_hg38_v3.0.01'
	Typer_database_l = txt_to_lines(database_path)
	Typer_database_l = Typer_database_l[1:]
	Haplotypes = Other_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT,CLASS)
	Type_result = path1 + sample_name + '/' + sample_name + '.'+CLASS+'_Type_result.txt'
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
		for l3 in Typer_database_l:
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
	database_path = abpath + '/input/MNS_Typer_database_hg38_v3.0.01'
	MNS_Typer_database_l = txt_to_lines(database_path)
	MNS_Typer_database_l = MNS_Typer_database_l[1:]
	Haplotypes = MNS_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT)
	Type_result = path1 + sample_name + '/' + sample_name + '.MNS_Type_result.txt'#含有所有可能的分型结果
	try :
		os.remove(Type_result)
	except Exception as e:
		pass
	result2 = open(Type_result,'a')
	Haplotype_type_result = {}
	#important_alles_index = [193,134,135,96,78,80,81,56,57,51,52,53,23,201,54] #可能是几个代表血型的位点情况
	for chr1,chr2 in Haplotypes.items():#血型匹配
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		for l3 in MNS_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
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
def ABO_Haplotype_maker(path1,sample_name,convert_result_path,INPUT,OUTPUT):
	abpath = sys.path[0]
	chromosome1 = 'F'
	chromosome2 = 'F' 
	Haplotypes = {chromosome1:chromosome2}
	dic = {}
	all_pos={}
	kb500=["500内的位点"]
	
	kb500c=[["500nei的单倍型组合"],["500nei的单倍型组合"]]

	s_v3 = abpath + "/software/soft/s_v4.out"
	hg38 = abpath + "/reference/hg38/turn_ucsc.hg38.fa"
	hg38fai = abpath + "/reference/hg38/turn_ucsc.hg38.fa.fai"
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
				af_dic[pos_start+"_"+pos_end+"_"+pos_alt] =af	
	with open (abpath+'/'+OUTPUT+'varResult/' + sample_name+'/'+sample_name+".H_deduped.bamtumor_bmd.vcfDedup_mapping.txt","r") as f1:
		H_af_dic={}
		for line in f1.readlines():
			if line.startswith("chr"):
				pos_chr = re.split(r"\.|\t",line)[0]
				pos_start = re.split(r"\.|\t",line)[1]
				pos_end = re.split(r"\.|\t",line)[2]
				pos_ref = re.split(r"\.|\t",line)[3]
				pos_alt = re.split(r"\.|\t",line)[4]
				dp = re.split(r"\.|\t",line)[5]
				ad = re.split(r"\.|\t",line)[6]
				af = re.search(r"(AF=)(.+)",re.split(r"\t",line)[3]).group(2)
				H_af_dic[pos_start+"_"+pos_end+"_"+pos_ref+"_"+pos_alt] =af	
		#仅考虑FUT1基因，FUT2基因全为纯合小h和小se的时候是纯孟买型，暂不考虑类孟买型。
		if float(H_af_dic["48750860_48750860_C_T"])>=0.8 or float(H_af_dic["48750821_48750821_T_C"])>=0.8 or float(H_af_dic["48750820_48750820_G_T"])>=0.8 or float(H_af_dic["48750769_48750769_C_G"])>=0.8 or float(H_af_dic["48750744_48750744_G_A"])>=0.8:      
			H="-"
		elif float(H_af_dic["48750730_48750731_AG_-"])>=0.8 or float(H_af_dic["48750696_48750696_G_A"])>=0.8 or float(H_af_dic["48750587_48750587_C_T"])>=0.8 or float(H_af_dic["48750557_48750557_A_C"])>=0.8 or float(H_af_dic["48750506_48750506_A_T"])>=0.8 or  float(H_af_dic["48750497_48750497_C_T"])>=0.8 and float(H_af_dic["48750496_48750496_G_T"])>=0.8 or float(H_af_dic["48750456_48750456_G_A"])>=0.8: 
			H="-"
		elif float(H_af_dic["48750400_48750401_TT_-"])>=0.8 or float(H_af_dic["48750338_48750338_G_A"])>=0.8 or  float(H_af_dic["48750334_48750334_G_C"])>=0.8 or float(H_af_dic["48750302_48750302_T_G"])>=0.8 or float(H_af_dic["48750235_48750235_C_G"])>=0.8 or float(H_af_dic["48750598_48750598_C_T"])>=0.8 or float(H_af_dic["48750588_48750588_A_G"])>=0.8 or float(H_af_dic["48750514_48750514_C_-"])>=0.8 or float(H_af_dic["48750492_48750492_-_G"])>=0.8 or float(H_af_dic["48750572_48750572_G_-"])>=0.8 or float(H_af_dic["48750828_48750828_G_-"])>=0.8:
			H="-"
		elif float(H_af_dic["48750994_48750994_A_T"])>=0.8:
			H="-" 
		elif float(H_af_dic["48703233_48703233_G_A"])>=0.8  or float(H_af_dic["48703374_48703374_A_T"])>=0.8 and float(H_af_dic["48703417_48703417_G_A"])>=0.8 or float(H_af_dic["48703558_48703558_G_A"])>=0.8 or float(H_af_dic["48703560_48703560_C_T"])>=0.8 or float(H_af_dic["48703617_48703617_C_T"])>=0.8 or float(H_af_dic["48703647_48703647_C_T"])>=0.8 or float(H_af_dic["48703653_48703653_C_T"])>=0.8 or float(H_af_dic["48703674_48703675_GT_-"])>=0.8:
			H="-"
		elif float(H_af_dic["48703677_48703680_GTC_-"])>=0.8 or float(H_af_dic["48703389_48703389_G_A"])>=0.8 and  float(H_af_dic["48703749_48703749_G_A"])>=0.8 or float(H_af_dic["48703767_48703767_C_-"])>=0.8 or float(H_af_dic["48703838_48703838_G_A"])>=0.8 or float(H_af_dic["48703857_48703857_G_A"])>=0.8 or float(H_af_dic["48703939_48703939_C_T"])>=0.8 or float(H_af_dic["48703291_48703291_C_T"])>=0.8 or float(H_af_dic["48703949_48703949_A_G"])>=0.8 or  float(H_af_dic["48703401_48703401_G_A"])>=0.8 or float(H_af_dic["48703807_48703807_C_A"])>=0.8:
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
		chr1 = chr1[1:]
		chr2 = chr2[1:]
		for l3 in RHCE_Typer_database_l:
			match_test = True
			match_test2 = True
			l3_list = l3.split('\t')
			Strings_REF = l3_list[-1].replace('\n','')#找到ABO_Typer_database_v1.0.txt里的最后一列基因组NY序列
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
	result2.write(Haplotypes2)
	result2.close()


#结果汇总
def typing_summary(OUTPUT,sample_name):
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
#主流程
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
		sample_dataprocess_end_time = datetime.datetime.now()
		print(' ' + sample_name +' dataprocess done :' + str(sample_dataprocess_end_time) )
		path1 = OUTPUT + '/varResult/'
		sample_list = os.listdir(path1)
		for sample_name in sample_list :
			sample_Typing_start_time = datetime.datetime.now()
			print('   ' + sample_name +' Typing start :' + str(sample_Typing_start_time) )
			#convert_result_suffix = '.hg38_multianno.txt'
			convert_result_path = find_file(sample_name,path1,convert_result_suffix)
			#BAM_path = convert_result_path.replace('.hg38_multianno.txt','.deduped.bam')#替代字符串找到bam的位置 传给get_amplicon_txt函数
			#print(BAM_path)
			blood_list=["GLOB","LU","KEL","FY","JK","DI","YT","SC","DO","CO","LW","CHRG","CROM","RHAG","KN","IN","OK","RAPH","JMH","I","FROS","GIL","LAN","VEL","CD59","AUG","CTL2","GE","H","KANNO","MAM","SID","LE","JR","GPIIIa","GPIbα","GPIIb","GPIa","GPIbβ","CD109"]			
			#PEL
			try:
				H_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("H-ok")
				ABO_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("ABO-ok")
				RHD_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RHD-ok")
				RHCE_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("RHCE-ok")
				MNS_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT)
				print("MNS-ok")
			except Exception as e:
				print(sample_name + "typing error :" + str(e))
				continue
			for B_Class in blood_list:
				try:	
					Other_Haplotyp_typer(path1,sample_name,convert_result_path,INPUT,OUTPUT,B_Class)
				except Exception as e:
					print(sample_name +B_Class+"typing error :" + str(e))
					continue
			sample_Typing_end_time = datetime.datetime.now()
			typing_summary(OUTPUT,sample_name) 
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
		version="Version 0.4 20200720")
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