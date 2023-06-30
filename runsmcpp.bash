#!/bin/bash
#用于获取脚本所在的路径，保存为变量path1,调用其他脚本都依赖这个路径。
path1="$(cd "$(dirname ${BASH_SOURCE[0]})";pwd)"
#输入vcf文件，一定不能是压缩格式，而且文件必须包含完整的vcf头部信息
if [ $# -lt 2 ];then
	echo "Usage:
	runsmcpp.bash -i inputvcffile -g G1 -o pwd -n threads -u 6.18e-8
	runsmcpp.bash -i inputvcffile -gf Groupfile -o pwd -n threads -u 1.18e-7
	runsmcpp.bash -i inputvcffile -g G1 -gf Groupfile -o pwd -n threads -u 2.81e-8
	
	#该脚本用来进行smc++分析

	####必须提供的参数
	-i 指定输入的vcf文件，注意输入的vcf必须有vcf完整的头部
	-g 指定分组名称（必须是字目开头的，不能有特殊字符），只能是一个分组名称，不支持多个分组
	-gf 指定分组文件（第一列是样本名称，第二列是分组名称，使用tab分割）
	-o 指定输出的目录
	
	###有默认值的参数
	-n 指定使用的cpu数量，默认是使用32个cpu
	-u 核苷酸突变率u,默认是7e-9
	
	###可以选择是否提供的参数（不提供也可以运行）
	-m 基因组的mask文件  提供此参数的时候，会使用此参数的文件作为mask文件。
	-G或-genome 基因组文件  提供此参数的时候，会重新计算基因组的mask文件，耗时较长，依赖bwa软件
	-m和-G参数只需要提供一个即可，如果同时提供，并不会重新计算mask文件。
	
	##三种分析模式
	模式1：只指定-g，则是把vcf所有的样本作为一个分组，此处的字符是输出的分组名称
	模式2：只指定-gf，则会根据输入的分组文件，提取vcf里每个分组的样本，分别分析。有多少个分组，就会有多少个结果。
	模式3：同时指定-g和-gf,则只分析指定的分组文件内-g参数指定的分组的结果
	
	Author：chaimol@163.com
	Date:2023.06.29
	Version:0.1
	
	注意：
	该脚本依赖的软件有
	bgzip
	tabix
	smc++
	
	"
	exit 1
fi

# 遍历所有参数
while [ $# -gt 0 ]; do
    case "$1" in
    -i)
        origin_inputvcf=$2
		shift
        ;;
	-o)
		workdir=$2
		shift
		;;
	-g)
		groupname=$2
		shift
		;;
	-gf)
		groupfile=$2
		shift
		;;
	-n)
		threads=$2
		shift
		;;
	-u)
		rate=$2
		shift
		;;
	-m)
		maskfile=$2
		shift
		;;
	-G|-genome)
		genomefile=$2
		shift
		;;
    *)
        echo "use -h for help!"
        exit 1
		;;
    esac
    shift
done

if [ -z $origin_inputvcf ];then
	"echo input file $origin_inputvcf not found!"
	exit 1
fi 

##处理输出目录，默认是当前路径，如果不存在改路径，则创建该路径
if [ -z $workdir ];then
	workdir="."
fi
if [ -e $workdir ];then
	cd $workdir
else
	mkdir ${workdir}
	cd ${workdir}
fi
if [ $? -ne 0 ];then
	echo "error: outputfile not found or the file access denied!"
	exit 1
fi

##设置默认cpu数量是32
if [ -z $threads ];then
	threads=32
fi

if [ -z $groupname ] && [ -z $groupfile ];then
	echo "You must be provide at least one of -g or -gf."
	exit 1
fi

#设置默认u的值
if [ -z $rate ];then
	rate="7e-9"
fi



##构建索引和压缩的vcf文件
inputvcf=$(basename ${origin_inputvcf}) #获取原始的文件名
bgzip ${origin_inputvcf} -c >${inputvcf}.gz
tabix ${inputvcf}.gz

##获取vcf文件的chr的列表
${path1}/getvcfchr.py ${origin_inputvcf} chr.list

if [ -z $maskfile ];then #不存在mask文件
	if [ -z $genomefile ];then #不存在基因组文件
		maskmode=0 # 此时使用无mask模式运行smc++
	else  #存在基因组文件
		maskmode=1 #此时需要重新计算基因组的mask文件
	fi
else
	maskmode=2 #此时存在mask文件，则运行smc++的mask模式
fi

##准备基因组的mask文件
function getgenomemaskfile(){
#ref_genome="~/Gh_meta/TM_1.genome.fa"
ref_genome=$1
# https://lh3lh3.users.sourceforge.net/snpable.shtml
${path1}/seqbility-20091110/splitfa ${ref_genome} 35 > ref.fq 
bwa index -a bwtsw ${ref_genome}
bwa aln -R 1000000 -O 3 -E 3 -t ${threads} ${ref_genome} ref.fq | bwa samse ${ref_genome} - ref.fq > ref.sam
${path1}/seqbility-20091110/gen_raw_mask.pl ref.sam > ref.sam.mask
${path1}/seqbility-20091110/gen_mask -l 35 -r 0.5 ref.sam.mask > ref.sam.mask.fa
#makeMappabilityMask.py这个是python2的脚本，原始是msmc2提供的，此处是我修改了支持自定义输入和输出文件
${path1}/makeMappabilityMask.py ref.sam.mask.fa ref.{}.mask.bed.gz
rm -rf ref.sam
rm -rf ref.fq
rm -rf ref.sam.mask 
#输出的mask结果文件是
#ref.{}.mask.bed.gz
cat ref.*.mask.bed.gz >ref_genome.mask_bed.gz
}

##运行一个分组的smc++,需要提供参数分组名和 样本名称用逗号分割
function runonesmcpp(){
	groupname=$1
	sampleid=$2
	if ! [ -e ${groupname} ];then
			mkdir ${groupname}
	fi
	cat chr.list|while read line;
	do
		if [ $maskmode -eq 0 ];then #原始的无mask模式
			smc++ vcf2smc ${inputvcf}.gz ${groupname}/${line}.smc.gz $line ${groupname}:${sampleid}
		elif [ $maskmode -eq 2 ];then #有maskfile模式
			if ! [ -e ${line}.mask.bed.gz ];then
				zcat ${maskfile}|awk -v chr=${line} '$1==chr'|gzip >${line}.mask.bed.gz
			fi
			smc++ vcf2smc --mask ${line}.mask.bed.gz ${inputvcf}.gz ${groupname}/${line}.smc.gz $line ${groupname}:${sampleid}
		else #有基因组的文件的模式
			maskfile="ref_genome.mask_bed.gz"
			if ! [ -e ${maskfile} ];then
				getgenomemaskfile ${genomefile} #判断不存在mask文件，则重新获取
			fi
			if ! [ -e ${line}.mask.bed.gz ];then
				zcat ${maskfile}|awk -v chr=${line} '$1==chr'|gzip >${line}.mask.bed.gz 
			fi
			smc++ vcf2smc --mask ${line}.mask.bed.gz ${inputvcf}.gz ${groupname}/${line}.smc.gz $line ${groupname}:${sampleid}
		fi
	done
#计算模型
	smc++ estimate --spline cubic --cores ${threads} -o ${groupname} ${rate} ${groupname}/*.smc.gz
#绘制模型
	smc++ plot ${groupname}.pdf ${groupname}/*.final.json -c
}

#主控制运行
##判断是否存在分组文件
if [ -z ${groupfile} ];then
        if [ -z ${groupname} ];then
                echo "you must be provide one of -g or -gf."
                exit 1
        else
                #不存在分组文件，只分析一个分组
        sampleid=$(head -1000 ${origin_inputvcf}|grep ^#|tail -1|cut -f10-|tr "\t" ",")
        echo "不存在分组文件，只分析一个分组"
        echo ${groupname}:$sampleid
		runonesmcpp $groupname $sampleid
        fi
else
        if [ -z ${groupname} ];then
                #不存在分组名称，此时需要遍历分组文件，每个都得分析
                groupnamelist=$(cut -f2 ${groupfile}|sort -u)
				if [ -z all_group ];then
					mkdir all_group
				fi
                for groupname in ${groupnamelist};
                do
                        sampleid=$(grep ${groupname} ${groupfile}|cut -f1|tr "\n" ","|sed 's/,$/\n/')
                        echo "${groupname}"
                        echo "${sampleid}"
                        runonesmcpp ${groupname} ${sampleid}
                        cp ${groupname}/*.final.json all_group/${groupname}_model.final.json
                done
                smc++ plot All.pdf all_group/*.final.json -c
        else
                #同时有分组名称和分组文件
                sampleid=$(grep $groupname $groupfile|cut -f1|tr "\n" ","|sed 's/,$/\n/')
                runonesmcpp ${groupname} ${sampleid}
                echo ${groupname}:${sampleid}
                echo "存在分组名称和分组文件，只分析一个"
        fi
fi
