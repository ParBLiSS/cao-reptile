import sys

spt_bstr = "spectrum construction"
err_bstr = "error correction"
tot_bstr = "total\t"


def timings(fname):
    tfname = fname.replace("run_", "").replace("reptile", "").replace("p", "_")
    felts = tfname.split("_")
    if(len(felts) <= 3):
      print fname
    assert(len(felts) > 3)
    rtyp = felts[0]
    rcfg = felts[1]
    rnpc = felts[2]
    with open(fname) as f:
        ncounter = 0
        for line in f:
            elts = line.split()
            if line.startswith(spt_bstr):
                ncounter += 1
                otimes = [str(ncounter), rcfg, rtyp, rnpc, "constr.", elts[2]]
                print "\t".join(otimes)
            elif line.startswith(err_bstr):
                otimes = [str(ncounter), rcfg, rtyp, rnpc,  "correct", elts[2]]
                print "\t".join(otimes)
            elif line.startswith(tot_bstr):
                otimes = [str(ncounter), rcfg, rtyp, rnpc, "total.t", elts[1]]
                print "\t".join(otimes)


def main(args):
    for fname in args:
        timings(fname)


if __name__ == '__main__':
    main(sys.argv[1:])
