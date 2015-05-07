import sys

fname1 = sys.argv[1]
fname2 = sys.argv[2]
with open(fname1) as f1:
    flines1 = f1.readlines()
with open(fname2) as f2:
    flines2 = f2.readlines()
nlen = len(flines1)/2
for i in range(nlen):
    print flines1[2*i].replace(">", "@"), flines1[2*i+1],
    print flines2[2*i].replace(">", "+"), flines2[2*i+1][1:],
