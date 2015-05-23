import sys

num_threads = 16
cfg_format = """#
InFaFile       /work/alurugroup/jnk/cao-runs/datasets/drosophila/pfastq/drosophila.fq
QFlag           1
OErrFile        %s-output-d%sh%dr%dp%d-out-
BatchSize       10000000
KmerLen         13
hd_max          %d
Step            12
ExpectSearch    4
T_ratio         0.5
QThreshold      35
MaxBadQPerKmer  4
Qlb             0
T_expGoodCnt    9
T_card          4
StoreReads      1
WriteOutput     0
WriteSpectrum   0
CacheOptimizedSearch      %d
KmerCacheSize      16
TileCacheSize      8
KmerSpectrumInFile   kmer-spectrum.bin
TileSpectrumInFile   tile-spectrum.bin
KmerSpectrumInAbsent 1
RunType            %d
Threads            %d
WorkDistribution    %d
WorkFactor          %d

"""

script_format = """#!/bin/bash
#PBS  -o %s_d%sh%dr%dw%dt%dp%d_out.log
#PBS  -e %s_d%sh%dr%dw%dt%dp%d_err.log
#PBS -lnodes=%d:ppn=%d:compute,walltime=%d:00:00
#PBS -lpmem=125gb
# Change to directory from which qsub command was issued
cd $PBS_O_WORKDIR
mpirun -np %d -f $PBS_NODEFILE -rr /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw%dt%dp%d.txt
mpirun -np %d -f $PBS_NODEFILE -rr /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw%dt%dp%d.txt
mpirun -np %d -f $PBS_NODEFILE -rr /work/alurugroup/jnk/cao-runs/cao-reptile/src/cao-preptile  /work/alurugroup/jnk/cao-runs/d%sruns/%s_d%sh%dr%dw%dt%dp%d.txt

"""


def write_config(oprefix, copt, dx, hd, rt, wd):
    fsuffix = "_d%sh%dr%dw%dt%dp" % (dx, hd, rt, wd, num_threads)
    wprocs = [  8, 16, 32,  64, 80]
    wfacts = [ 80, 40, 20,  10,  8]
    for p, wf in zip(wprocs, wfacts):
        x =  cfg_format % (oprefix, dx, hd, rt, p, hd,
                           copt, rt, num_threads, wd, wf)
        fname = oprefix + fsuffix + str(p * num_threads) + ".txt"
        with open(fname, "w") as of:
            of.write(x + "\n")


def gen_script(pref, dx, hd, rt, wd):
    fprefix = "_d%sh%dr%dw%dt%dp" % (dx, hd, rt, wd, num_threads)
    sprocs = [  8,  16, 32,  64, 80]
    snodes = [  8,  16, 32,  64, 80]
    hours  = [ 10,   8,  4,   4,  4]
    for pr,nd,nhr in zip(sprocs,snodes,hours):
        ppn = 16
        scx = script_format % (pref, dx, hd, rt, wd, num_threads, pr * num_threads, 
                               pref, dx, hd, rt, wd, num_threads, pr * num_threads, 
                               nd, ppn, nhr,
                               pr, dx, pref, dx, hd, rt, wd, num_threads, pr * num_threads, 
                               pr, dx, pref, dx, hd, rt, wd, num_threads, pr * num_threads, 
                               pr, dx, pref, dx, hd, rt, wd, num_threads, pr * num_threads)
        fname = pref + fprefix + str(pr * num_threads) + ".sh"
        with open(fname, "w") as of:
            of.write(scx + "\n")

print sys.argv[1:]
if len(sys.argv) != 4:
    print "python gen_config.py <dataset> <hamming dist> <run type>"
    exit()

dx = sys.argv[1]
hd = int(sys.argv[2])
rt = int(sys.argv[3])
# wd = int(sys.argv[4])

write_config("dfault", 0, dx, hd, rt, 0)
write_config("dfault", 0, dx, hd, rt, 1)
# write_config("caware", 1, dx, hd, rt, 0)
write_config("caware", 1, dx, hd, rt, 1)
# write_config("cobliv", 2, dx, hd, rt, 0)
# write_config("cobliv", 2, dx, hd, rt, 1)
# write_config("sorted", 3, dx, hd, rt, 0)
# write_config("sorted", 3, dx, hd, rt, 1)
# write_config("parcwr", 4, dx, hd, rt, 0)
# write_config("parcwr", 4, dx, hd, rt, 1)

gen_script("dfault", dx, hd, rt, 0)
gen_script("dfault", dx, hd, rt, 1)
# gen_script("caware", dx, hd, rt, 0)
gen_script("caware", dx, hd, rt, 1)
# gen_script("cobliv", dx, hd, rt, 0)
# gen_script("cobliv", dx, hd, rt, 1)
# gen_script("sorted", dx, hd, rt, 0)
# gen_script("sorted", dx, hd, rt, 1)
# gen_script("parcwr", dx, hd, rt, 0)
# gen_script("parcwr", dx, hd, rt, 1)
