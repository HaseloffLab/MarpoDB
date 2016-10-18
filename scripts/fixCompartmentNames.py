import sys
import subprocess

genomeDBPath = sys.argv[1]
mRNADBPath   = sys.argv[2]

genomeData   = subprocess.Popen(["blastdbcmd", "-db", genomeDBPath, "-entry", "all", "-outfmt", "'%a %t'"], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()[0]
genomeDict   = dict( (row.split()[0][1:], row.split()[1][:-1]) for row in genomeData.split('\n')[:-1] )

mRNAData	 = subprocess.Popen(["blastdbcmd", "-db", mRNADBPath, "-entry", "all", "-outfmt", "'%a %t'"], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()[0]
mRNADict     = dict( (row.split()[0][1:], row.split()[1]) for row in mRNAData.split('\n')[:-1] )

compFile     = open(sys.argv[3])


for line in compFile:
 	if line != "\n":
 		tabs = line.split()
 		mID  = "BL_ORD_ID:{0}".format(tabs[0].split('|')[2])
 		gID  = "BL_ORD_ID:{0}".format(tabs[1].split('|')[2])

 		tabs[0] = mRNADict[mID]
 		tabs[1] = genomeDict[gID]

 		line = '\t'.join(tabs)
 	print line
