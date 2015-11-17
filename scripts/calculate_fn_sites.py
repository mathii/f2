# taking a flat haplotype file (just a 0/1 matrix with nsnp rows and nhaps+1 cols, first col is poition
# where every entry is a doubleton), output a list of three cols where the first col
# is position, and the second and subsequent are the [0-based] indices of everyone who shares that
# variant

from __future__ import division
import gzip, numpy, sys, getopt
import pdb

REP_INT=1e3                     # Report every REP_INT lines

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "haplotype_file":"", "out":"", "n":0, "verbose":False  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:o:n:v:", ["hap","out","n", "verbose"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-h","--hap"]:           options["haplotype_file"] = a
        elif o in ["-o","--out"]:         options["out"] = a
        elif o in ["-n","--n"]:           options["n"] = int(a)
        elif o in ["-v","--ver"]:         options["verbose"] = bool(a)

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):
    """
    main
    """
    if options["haplotype_file"][-3:]==".gz":
        haps_in=gzip.open(options["haplotype_file"], "r")
    else:
        haps_in=open(options["haplotype_file"], "r")

    if options["out"][-3:]!=".gz":
        options["out"]+=".gz"
    out=gzip.open(options["out"], "w")

    out_string="%d"+"".join(["\t%d"]*options["n"])+"\n"

    line_count=1
    for line in haps_in:

        if options["verbose"] and not line_count%REP_INT:
            print "\rLine " + str(line_count),
            sys.stdout.flush()
        line_count+=1

        bits=line.split()
        pos=int(bits[0])

        minor=None
        gt=bits[1:]
        pdb.set_trace()
        if not all([x=="1" or x=="0" for x in gt]): # any missing genotypes
            continue

        indexed_1=[i for i in range(len(gt)) if gt[i]=="1"]
        if len(indexed_1)==options["n"]:
            out.write(out_string % ((pos,)+tuple(indexed_1))) 
        elif len(indexed_1)==len(gt)-options["n"]: # 0 is the minor allele
            indexed_0=[i for i in range(len(gt)) if gt[i]=="0"]
            out.write(out_string % ((pos,)+tuple(indexed_0))) 
        else:
            raise Exception("MAC!=n at pos %d"%(pos,))
 
    out.close()
    haps_in.close()
    
    return

##########################################################################################################


if __name__=="__main__":
    options=parse_options()
    main(options)



