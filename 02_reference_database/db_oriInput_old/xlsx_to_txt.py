#!/usr/bin/env python
# coding=utf-8
import pandas as pd

df = pd.read_excel('./005_LU_Typer_database_hg38-hg19_v5.0_01.xlsx', header=0)
print('开始写入txt文件...')
df.to_csv('005_LU_Typer_database_hg38-hg19_v5.0_jy', header= 0, sep='\t', index=False)       # 写入，\t分隔
print('文件写入成功!')
