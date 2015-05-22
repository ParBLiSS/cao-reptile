import sys

def parse_cfg(rcfg):
    elts = rcfg.split('h')
    if(len(elts) != 2):
        return "\t".join([rcfg, "", "", ""])
    dset = elts[0].upper()
    rcfg = elts[1]
    elts = rcfg.split('r')
    if(len(elts) != 2):
        return "\t".join([rcfg, "", "", ""])
    hdmx = elts[0] 
    rcfg = elts[1]
    elts = rcfg.split('w')
    if(len(elts) != 2):
        return "\t".join([rcfg, "", "", ""])
    rtyp = elts[0] 
    rcfg = elts[1]
    elts = rcfg.split('t')
    if(len(elts) != 2):
        return "\t".join([rcfg, "", "", ""])
    dist = elts[0] 
    nthr = elts[1]
    return "\t".join([dset, hdmx, dist, nthr])
    
def timings(fname):
    tfname = fname.replace("run_", "").replace("reptile", "").replace("p", "_")
    felts = tfname.split("_")
    if(len(felts) <= 3):
      print fname
    assert(len(felts) > 3)
    rtyp = felts[0]
    rcfg = parse_cfg(felts[1])
    rnpc = felts[2]
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
