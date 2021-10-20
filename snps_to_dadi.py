import re
import os
import pandas as pd
import argparse
import random

parser = argparse.ArgumentParser(description="create pop gen files")
parser.add_argument('--file',
			default=False,
			help="file with sample data")
parser.add_argument('--vcffile',
			default = None,
			help = "input vcf file")
args = parser.parse_args()

snps = args.vcffile
d = pd.read_csv(args.file)
sps = {}
inds_sps = {}
for a, b in zip(d.SAMPLE_ID, d.nDNA_LINEAGE):
	inds_sps[a] = b
	if b not in sps:
		sps[b] = {}
	sps[b][a] = 1

def most_common(lst):
    return max(set(lst), key=lst.count)

def least_common(lst):
    return min(set(lst), key=lst.count)

def parse_snps(sps, snps):
	keep = {}

	f = open(snps, 'r')
	print('Ingroup\tOutgroup\tAllele1\tzastictus\tpallasotus\tAllele2\tzastictus\tpallasotus\tGene\tPosition')
	
	for l in f:
		if re.search('#CHROM', l):
			d = re.split('\t', l.rstrip())
			inds = d[9:]
		if not re.search('^#', l):
			d =	re.split('\t', l.rstrip())
			snps = d[9:]

			good = True
			# get rid of indels & non biallelics
			if len(d[3]) > 1 or len(d[4]) > 1:
				good = False

			if good:
				snps = [re.search('(^\S\S\S)', snp).group(1) for snp in snps]

				a = [re.split('/', snp) for snp in snps]

				# subsample to just focal taxa
				all = []
				cur = {'zastictus': [], 'pallasotus': []}
				for ind, geno in zip(inds, a):
					sp = inds_sps[ind]
					if sp in cur:
						if geno[0] != '.':
							cur[sp] += geno
					else:
						if geno[0] != '.':
							all += geno

				aa = {'0': d[3], '1': d[4]}

				# only keep with low missing
				miss = []
				for sp in cur:
					curmiss = len(cur[sp]) /  (len(sps[sp]) * 2.0)
					miss.append(curmiss)

				if miss[0] > 0.66 and miss[1] > 0.66:
					alleles = []
					for sp in cur:
						alleles += cur[sp]

					mac1 = alleles.count('0')
					mac2 = alleles.count('1')

					if mac1 > 0 and mac2 > 0:

						# get outgroup allele
						if len(all) > 0:
							out = most_common(all)
						else:
							out = most_common(alleles)
						ing = [x for x in alleles if x != out][0]

						info = []
						info.append('-%s-' % aa[ing])
						info.append('-%s-' % aa[out])
						info.append(aa[ing])
						info.append(cur['zastictus'].count(ing))
						info.append(cur['pallasotus'].count(ing))
						info.append(aa[out])
						info.append(cur['zastictus'].count(out))
						info.append(cur['pallasotus'].count(out))
						info.append(d[0])
						info.append(d[1])
						print('\t'.join([str(x) for x in info]))



parse_snps(sps, snps)