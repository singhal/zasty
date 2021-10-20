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

out = '/Users/singhal/Desktop/atlas_group.dxy.8Sept21.csv'
o = open(out, 'w')
o.write('sp1,sp2,dxy,length,num_loci\n')

sps = list(sps.keys())
for ix, sp1 in enumerate(sps):
	for sp2 in sps[ix + 1:]:

		alldiv = []

		for i in loc:
			if sp1 in loc[i] and sp2 in loc[i]:
				hap1 = loc[i][sp1]
				hap2 = loc[i][sp2]

				divs = []

				for h1 in hap1:
					for h2 in hap2:
						totlen = 0
						diff = 0
						for a, b in zip(h1, h2):
							if a in bps and b in bps:
								totlen += 1
								if a != b:
									diff += 1
						div = diff / float(totlen)
						divs.append(div)
				meandiv = np.mean(divs)
				length = len(hap1[0])

				alldiv.append([meandiv, length])
		finaldiv = 0
		totlength = 0
		for i in alldiv:
			totlength += i[1]
			finaldiv += i[1] * i[0]
		if totlength > 0:
			div = finaldiv / float(totlength)
		else:
			div = np.nan
		o.write('%s,%s,%.5f,%s,%s\n' % (sp1, sp2, div, totlength, len(alldiv)))
		o.write('%s,%s,%.5f,%s,%s\n' % (sp2, sp1, div, totlength, len(alldiv)))
o.close()



