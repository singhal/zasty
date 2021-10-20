import re
import glob

# do admixture
indir= '/home/babs/Desktop/Ct_zast/admixture_sim/'
files = glob.glob("%s/*ped" % indir)
for file in files:
        file = re.sub('.ped', '', file)
        print("~/Desktop/bin/plink --file %s --make-bed --noweb --out %s" % (file, file))
        newfile = file + '.bed'
        for k in range(1, 7):
                print("~/Desktop/bin/admixture_linux-1.3.0/admixture --cv %s %s" % (newfile, k))
