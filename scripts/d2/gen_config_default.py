import sys

cfg_str_format = """#
InFaFile       /work/alurugroup/jnk/cao-runs/datasets/SRR034509/pfastq/SRR034509_1.fq
QFlag           1
OErrFile        %s-output-d%sh%dr%dp%d-out
BatchSize       10000000
KmerLen         12
hd_max          %d
Step            12
ExpectSearch    4
T_ratio         0.5
QThreshold      72
MaxBadQPerKmer  4
Qlb             39
T_expGoodCnt    24
T_card          5
StoreReads      1
WriteOutput     0
WriteSpectrum   0
CacheOptimizedSearch      %d
KmerCacheSize      16
TileCacheSize      8
KmerSpectrumInFile   kmer-spectrum.bin
TileSpectrumInFile   tile-spectrum.bin
RunType            %d
"""

str_format = """#!/bin/bash
#PBS  -o %s_d%sh%dr%dw0t1p%d_out.log
#PBS  -e %s_d%sh%dr%dw0t1p%d_err.log
#PBS -lnodes=%d:ppn=%d:compute,walltime=4:00:00
# Change to directory from which qsub command was issued
cd $PBS_O_WORKDIR
mpirun -np %d /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw0t1p%d.txt
mpirun -np %d /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw0t1p%d.txt
mpirun -np %d /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw0t1p%d.txt

"""

def write_config(oprefix, fsuffix, copt, dx, hd, rt):
    procs = [ 128, 256, 512, 1024, 1280]
    for p in procs:
        x =  cfg_str_format % (oprefix, dx, hd, rt, p,
                               hd, copt, rt)
        fname = oprefix + fsuffix + str(p) + ".txt"
        with open(fname, "w") as of:
            of.write(x + "\n")

def gen_script(pref, fprefix, dx, hd, rt):
    wprocs = [ 128, 256, 512, 1024, 1280]
    wnodes = [   8,  16,  32,   64,   80]
    for pr,nd in zip(wprocs, wnodes):
        ppn = pr/nd
        scx = str_format % (pref, dx, hd, rt, pr, 
                            pref, dx, hd, rt, pr, 
                            nd, ppn,
                            pr, dx, pref, dx, hd, rt, pr, 
                            pr, dx, pref, dx, hd, rt, pr, 
                            pr, dx, pref, dx, hd, rt, pr)
        fname = pref + fprefix + str(pr) + ".sh"
        with open(fname, "w") as of:
            of.write(scx + "\n")


print sys.argv[1:]
if len(sys.argv) != 4:
    print "python gen_config.py <dataset> <hamming dist> <run type>"
    exit()

dx = sys.argv[1]
hd =  int(sys.argv[2])
rt =  int(sys.argv[3])
fsuffix = "_d%sh%dr%dw0t1p" % (dx, hd, rt)
fprefix = "_d%sh%dr%dw0t1p" % (dx, hd, rt)

write_config("dfault", fsuffix, 0, dx, hd, rt)
# write_config("caware", fsuffix, 1, dx, hd, rt)
# write_config("cobliv", fsuffix, 2, dx, hd, rt)
# write_config("sorted", fsuffix, 3, dx, hd, rt)

gen_script("dfault", fprefix, dx, hd, rt)
# gen_script("caware", fprefix, dx, hd, rt)
# gen_script("cobliv", fprefix, dx, hd, rt)
# gen_script("sorted", fprefix, dx, hd, rt)
