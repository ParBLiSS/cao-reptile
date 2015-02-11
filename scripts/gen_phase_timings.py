import sys


def timings(fname):
    felts = fname.split("_")
    assert(len(felts) > 3)
    rtyp = felts[1]
    rcfg = felts[2].replace("reptile", "")
    rnpc = felts[3]
    with open(fname) as f:
        ncounter = 0
        tstate = 0
        for line in f:
            elts = line.split()
            if tstate == 0 and line.startswith("rank"):
                ncounter += 1
                tstate = 1
            if tstate == 0 and line.startswith("nproc"):
                ncounter += 1
                tstate = 1
            if tstate == 1 and line.startswith("Input"):
                tstate = 0
            if tstate == 1 and line.startswith("parameter"):
                tstate = 0
            if tstate == 0:
                continue
            if line.startswith("-"):
                continue
            if len(elts) < 6:
                continue
            phase = "-".join([elts[1], elts[2]])
            pid = "XX" if phase.endswith("global") else elts[0]
            otimes = [str(ncounter), rcfg, rtyp, rnpc, phase,
                      pid, elts[3], elts[4], elts[5]]
            print "\t".join(otimes)


def main(args):
    for fname in args:
        timings(fname)


if __name__ == '__main__':
    main(sys.argv[1:])
