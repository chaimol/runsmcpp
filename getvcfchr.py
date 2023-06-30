#!/usr/bin/env python3
#输入 vcf文件和输出文件名
#输出的是vcf的chr的列表

import sys
import gzip

def get_unique_chromosomes_from_vcf(filename):
    chromosomes = set()
    opener = gzip.open if filename.endswith('.gz') else open

    with opener(filename, 'rt') as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                chromosome = line.split('\t')[0]
                chromosomes.add(chromosome)
    return sorted(chromosomes)

def write_chromosomes_to_file(chromosomes, output_file):
    with open(output_file, "w") as file:
        for chromosome in chromosomes:
            file.write(f"{chromosome}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python3 getvcfchr.py inputvcf.gz outputtxt")
    else:
        vcf_filename = sys.argv[1]
        output_filename = sys.argv[2]
        unique_chromosomes = get_unique_chromosomes_from_vcf(vcf_filename)
        write_chromosomes_to_file(unique_chromosomes, output_filename)
        print("Chromosomes written to file:", output_filename)
