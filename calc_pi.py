import re
import numpy as np


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
		sps[dd[7]] = 1
f.close()

# read in data
d = r'/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/ipyrad_output/atlas_group_outfiles/ctzast1.alleles'
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
		ind = re.sub('_\d$', '', ind)
		ind = re.sub('\^', '', ind)

		# no outgroup
		if ind in inds:
			sp = inds[ind]

			if sp not in loc[locnum]:
				loc[locnum][sp] = []
			loc[locnum][sp].append(dd[1].upper())
f.close()

bps = ['A', 'T', 'C', 'G']

out = '/Users/singhal/Desktop/atlas_group.pi.8Sept21.csv'
o = open(out, 'w')
o.write('sp1,pi,length,num_loci\n')

sps = list(sps.keys())
for sp in sps:

	alldiv = []

	for i in loc:
		if sp in loc[i]:

			# all the haplotypes for that species
			haps = [list(x) for x in loc[i][sp]]

			alldenom = []
			alldiff = []

			for ix, hap1 in enumerate(haps):
				haps2 = haps[(ix +1): ]
				for hap2 in haps2:
					denom = 0
					diff = 0
					for a, b in zip(hap1, hap2):
						if a in bps and b in bps:
							denom += 1
							if a != b:
								diff += 1
					alldenom.append(denom)
					alldiff.append(diff)

			divs = [a / float(b) for a, b in zip(alldiff, alldenom)]
			totdiv = sum(divs) / len(alldenom)
			alldiv.append((totdiv, np.mean(alldenom)))

	totlength = 0 
	totdiv = 0
	for x in alldiv:
		totdiv += x[0] * x[1]
		totlength += x[1]
	if totlength > 0:
		totdiv = totdiv / totlength
	else:
		totdiv = np.nan

	o.write('%s,%.5f,%s,%s\n' % (sp, totdiv, totlength, len(alldiv)))		


