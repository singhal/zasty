import re
import numpy as np


# create a hash of
# inds to species
d = r'/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv'
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
d = r'/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/ipyrad_output/duricola_ingroup_outfiles/duricola_ingroup.alleles'
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

totlen = 0
numloc = 0
for locus in loc:
	if 'zastictus' in loc[locus] and 'pallasotus' in loc[locus]:
		sz = len(loc[locus]['zastictus']) / len(sps['zastictus'])
		sp = len(loc[locus]['pallasotus']) / len(sps['pallasotus'])

		if sz > 0.66 and sp > 0.66:
			totlen += len(loc[locus]['zastictus'][0])
			numloc += 1
print(totlen, numloc)