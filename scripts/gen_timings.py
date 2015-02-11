import sys

spt_bstr = "K-SPECTRUM CONSTRUCTION TIME"
err_bstr = "ERR CORRECTION TIME"
tot_bstr = "TOTAL TIME"


def timings(fname):
    felts = fname.split("_")
    assert(len(felts) > 3)
    rtyp = felts[1]
    rcfg = felts[2].replace("reptile", "")
    rnpc = felts[3]
    with open(fname) as f:
        ncounter = 0
        for line in f:
            elts = line.split()
            if line.startswith(spt_bstr):
                ncounter += 1
                otimes = [str(ncounter), rcfg, rtyp, rnpc, "constr.", elts[3]]
                print "\t".join(otimes)
            elif line.startswith(err_bstr):
                otimes = [str(ncounter), rcfg, rtyp, rnpc,  "correct", elts[3]]
                print "\t".join(otimes)
            elif line.startswith(tot_bstr):
                otimes = [str(ncounter), rcfg, rtyp, rnpc, "total.t", elts[2]]
                print "\t".join(otimes)


def main(args):
    for fname in args:
        timings(fname)


if __name__ == '__main__':
    main(sys.argv[1:])
