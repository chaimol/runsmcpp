This is a pipline for run smc++ analysis.
# 该脚本用来进行smc++分析

## Require software(你需要自行安装以下的软件)
- bgzip
- tabix
- smc++
- python2 and python3
- perl
- bwa(只有提供基因组文件时需要)

##Install this pipline
	```
	chmod 757 runsmcpp.bash
	chmod 757 makeMappabilityMask.py
	chmod 757 getvcfchr.py
	```
## usage
	`runsmcpp.bash -i inputvcffile -g G1 -o pwd -n threads -u 6.18e-8`
	`runsmcpp.bash -i inputvcffile -gf Groupfile -o pwd -n threads -u 1.18e-7`
	`runsmcpp.bash -i inputvcffile -g G1 -gf Groupfile -o pwd -n threads -u 2.81e-8`
	`runsmcpp.bash -i inputvcffile -g G1 -gf Groupfile -o pwd -n threads -u 2.81e-8 -G genome.fa`
	`runsmcpp.bash -i inputvcffile -g G1 -gf Groupfile -o pwd -n threads -u 2.81e-8 -m genome.mask.bed.gz`

#### 必须提供的参数
	-i 指定输入的vcf文件，注意输入的vcf必须有vcf完整的头部
	-g 指定分组名称（必须是字目开头的，不能有特殊字符），只能是一个分组名称，不支持多个分组
	-gf 指定分组文件（第一列是样本名称，第二列是分组名称，使用tab分割）
	-o 指定输出的目录
	
### 有默认值的参数
	-n 指定使用的cpu数量，默认是使用32个cpu
	-u 核苷酸突变率u,默认是7e-9
	
### 可以选择是否提供的参数（不提供也可以运行）
	-m 基因组的mask文件  提供此参数的时候，会使用此参数的文件作为mask文件。
	-G或-genome 基因组文件  提供此参数的时候，会重新计算基因组的mask文件，耗时较长，依赖bwa软件
	-m和-G参数只需要提供一个即可，如果同时提供，并不会重新计算mask文件。


# 三种分析模式
- 模式1：只指定-g，则是把vcf所有的样本作为一个分组，此处的字符是输出的分组名称
- 模式2：只指定-gf，则会根据输入的分组文件，提取vcf里每个分组的样本，分别分析。有多少个分组，就会有多少个结果。
- 模式3：同时指定-g和-gf,则只分析指定的分组文件内-g参数指定的分组的结果

#### Author：chaimol@163.com
#### Date:2023.06.29
