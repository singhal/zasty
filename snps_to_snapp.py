import re
import os
import pandas as pd
import argparse
import random

parser = argparse.ArgumentParser(description="create pop gen files")
parser.add_argument('--random',
			default=False,
			action='store_true',
			help="keep just one snp. default is to use all.")
parser.add_argument('--file',
			default=False,
			help="file with sample data")
parser.add_argument('--stem', 
			default=None,
			help="stem for all output files; can just default")
parser.add_argument('--MISS',
			default = None,
			help = "percent missing")
parser.add_argument('--MAC',
			default = None,
			help = "mininum allele count req'd; rec 2")
parser.add_argument('--sample',
			default = None,
			help = "# of snps to randomly sample")
parser.add_argument('--vcffile',
			default = None,
			help = "input vcf file")
parser.add_argument('--outdir',
			default = None,
			help = "output dir")
parser.add_argument('--drop',
			default = None,
			help = "individuals to drop")
args = parser.parse_args()

snps = args.vcffile
outdir = args.outdir
if not os.path.isdir(outdir):
	os.mkdir(outdir)
d = pd.read_csv(args.file)
sps = {}
for a, b in zip(d.SAMPLE_ID, d.nDNA_LINEAGE):
	sps[a] = b
 
def make_code():
	ix = 0
	code = {}
	invcode = {}

	allbp = ['A', 'T', 'C', 'G']

	for i, bp1 in enumerate(allbp):
		if bp1 not in invcode:
			invcode[bp1] = {}
		for bp2 in allbp[i:]:
			if bp2 not in invcode:
				invcode[bp2] = {}
			code[str(ix)] = [bp1, bp2]
			invcode[bp1][bp2] = str(ix)
			invcode[bp2][bp1] = str(ix)
			ix += 1

	invcode['N'] = {}
	invcode['N']['N'] = '-'
	code['-'] = ['N', 'N']
	
	return code, invcode


def get_af(snp, code):
	a = [code[x] for x in snp]
	a = [x for geno in a for x in geno]
	a = [x for x in a if x != 'N']

	uniq = list(set(a))
	cts = [ a.count(bp) for bp in uniq]

	return len(uniq), min(cts)


def get_random(var):
	var2 = {}
	for c in var:
		keep = random.choice(list(var[c].keys()))
		var2[c] = {}
		var2[c][keep] = var[c][keep]

	return var2


def parse_snps(snps, invcode, MISS, MAC):
	keep = {}

	f = open(snps, 'r')
	
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

			snps = [re.search('(^\S\S\S)', snp).group(1) for snp in snps]
			complete = 1 - (snps.count('./.') / float(len(snps)))
			if complete > MISS and good:
				a = [re.split('/', snp) for snp in snps]
				a = [x for b in a for x in b]
				af1 = a.count('0')
				af2 = a.count('1')

				if af1 >= MAC and af2 >= MAC:
					# keep[d[0]][d[1]] = snp
					chr = re.search('_(\d+)', d[0]).group(1)
					if chr not in keep:
						keep[chr] = {}

					snp = []
					alleles = {'0': d[3], '1': d[4], '.': 'N'}
					for geno, ind in zip(snps, inds):
						geno = [alleles[x] for x in re.split('/', geno)]
						snp.append( invcode[geno[0]][geno[1]] )
			
					keep[chr][int(d[1])] = ''.join(snp)
						
	f.close()

	return keep, inds


def make_binary(var, code):
	binary = {}
	for c in var:
		binary[c] = {}
		for pos in var[c]:
			genos = list(var[c][pos])
			genos = [code[geno] for geno in genos]
		
			alleles = [bp for geno in genos for bp in geno]
			alleles = [bp for bp in alleles if bp != 'N']
			alleles = list(set(alleles))
		
			a = {}
			for ix, allele in enumerate(alleles):
				a[allele] = str(ix)
			a['N'] = '-9'

			new = []
			for geno in genos:
				new.append([a[geno[0]], a[geno[1]]])
			binary[c][pos] = new
	return binary


def remove_miss(var, inds):
	tot_len = 0
	indmiss = dict([(ind, 0) for ind in inds])

	for c in var:
		for pos in var[c]:
			tot_len += 1
			for ind, geno in zip(inds, var[c][pos]):
				if geno == '-':
					indmiss[ind] += 1

	to_drop = []
	for ind in inds:
		permiss = indmiss[ind] / float(tot_len)
		if permiss >= 0.5:
			to_drop.append(ind)
	
	return to_drop


def make_snapp(inds, b, stem, drop):
	out = os.path.join(outdir, '%s.snapp.fasta' % stem)
	o = open(out, 'w')

	snp = {'0': {'0': '0', '1': '1'}, 
		   '1': {'0': '1', '1': '2'}, 
		   '-9': {'-9': '-'}}

	for ix, ind in enumerate(inds):
		if ind not in drop:
			gen = ''

			for c in sorted(b.keys()):
				for pos in sorted(b[c].keys()):
					geno = list(b[c][pos][ix])
					gen += snp[geno[0]][geno[1]]
			o.write('>%s_%s\n%s\n' % (sps[ind], ind, gen))

	o.close()


def get_sample(var, sample):
	sample = int(sample)
	var2 = {}
	while sample > 0:
		rc = random.choice(var.keys())
		if rc not in var2:
			var2[rc] = {}
		rcpos = random.choice(var[rc].keys())
		if rcpos not in var2[rc]:
			var2[rc][rcpos] = var[rc][rcpos]
			sample = sample - 1
	return var2 

MISS = float(args.MISS)
MAC = int(args.MAC)

code, invcode = make_code()
var, inds = parse_snps(snps, invcode, MISS, MAC)
stem = args.stem + '.miss%s.MAC%s' % (MISS, MAC)

# randomly subsample one SNP per locus
if args.random:
	var = get_random(var)
	stem = stem + '.thinned'
if args.sample:
	var = get_sample(var, args.sample)
	stem = stem + '.sample%s' % args.sample


drop = remove_miss(var, inds)
if args.drop:
	drop = drop + re.split(',', args.drop)
drop = list(set(drop))
out = os.path.join(outdir, stem + '.drop.out')
of = open(out, 'w')
for l in drop:
	of.write(l + '\tdrop\n')

# get binary sructure
binary = make_binary(var, code)
make_snapp(inds, binary, stem, drop)

of.close()
