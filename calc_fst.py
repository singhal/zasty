import re
import numpy as np

def fst_reich(counts, sample_sizes):
	counts1 = np.array(counts)[0]
	counts2 = np.array(counts)[1]

	sample_sizes1 = np.array(sample_sizes).astype('float')[0]
	sample_sizes2 = np.array(sample_sizes).astype('float')[1]

	h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
	h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))
	
	N = []
	D = []

	for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
		n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
		N.append(n)
		d = n + _h1 + _h2
		D.append(d)

	F = np.sum(N) / np.sum(D)

	return F

# create a hash of
# inds to species
d = r'/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project4.csv'
f = open(d, 'r')
header = f.readline()
inds = {}
sps = {}
for l in f:
	l = l.rstrip()
	dd = re.split(',', l)
	if dd[7] != "NA":
		inds[dd[0]] = dd[7]
		if dd[7] not in sps:
			sps[dd[7]] = {}
		sps[dd[7]][dd[0]] = 1
f.close()


# read in data
d = r'/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/ipyrad_output/atlas_group_outfiles/ctzast1.loci'
f = open(d, 'r')

locnum = 0
loc = {}

loc[locnum] = {}

for l in f:
	l = l.rstrip()
	if re.search('^//', l):
		locnum += 1
		loc[locnum] = {}
	else:
		dd = re.split('\s+', l)

		# seq name
		ind = dd[0]

		loc[locnum][ind] = dd[1].upper()
f.close()

code = {'A': ['A', 'A'],
		'T': ['T', 'T'],
		'C': ['C', 'C'],
		'G': ['G', 'G'],
		'M': ['A', 'C'],
		'R': ['A', 'G'],
		'W': ['A', 'T'],
		'S': ['C', 'G'],
		'Y': ['C', 'T'],
		'K': ['G', 'T']}

out = '/Users/singhal/Desktop/atlas_group.fst.8Sept21.csv'
o = open(out, 'w')
o.write('sp1,sp2,fst,num_snps\n')

splist = list(sps.keys())
for ix, sp1 in enumerate(splist):
	for sp2 in splist[ix + 1:]:

		counts = {sp1: [], sp2: []}
		sizes = {sp1: [], sp2: []}

		for locnum in loc:
			a1 = [loc[locnum][ind] for ind in sps[sp1] if ind in loc[locnum]]
			a2 = [loc[locnum][ind] for ind in sps[sp2] if ind in loc[locnum]]

			if len(a1) > 0 and len(a2) > 0:
				loclen = len(a1[0])
				for i in range(0, loclen):
					a = [ind[i] for ind in a1]
					a = [code[x] for x in a if x in code]
					a = [item for sublist in a for item in sublist]

					b = [ind[i] for ind in a2]
					b = [code[x] for x in b if x in code]
					b = [item for sublist in b for item in sublist]

					# only want two alleles
					alleles = list(set(a + b))
					if len(alleles) == 2:
						allele = alleles[0]

						if len(a) == 0:
							count1 = np.nan
						else:
							count1 = a.count(allele)
						counts[sp1].append(count1)
						sizes[sp1].append(len(a))

						if len(b) == 0:
							count2 = np.nan
						else:
							count2 = b.count(allele)
						counts[sp2].append(count2)
						sizes[sp2].append(len(b))

		alleles = np.array([counts[sp1], counts[sp2]])
		sizes = np.array([sizes[sp1], sizes[sp2]])
		to_mask = np.any(np.isnan(alleles), axis=0)
		alleles = alleles[:, ~to_mask]
		sizes = sizes[:, ~to_mask]

		if len(alleles[0]) > 0:
			fst = fst_reich(alleles, sizes)
			counts = len(alleles[0])
		else:
			fst = np.nan
			counts = 0
		
		o.write('%s,%s,%.5f,%s\n' % (sp1, sp2, fst, counts))
		o.write('%s,%s,%.5f,%s\n' % (sp2, sp1, fst, counts))			
o.close()