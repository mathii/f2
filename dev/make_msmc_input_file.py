# Make an msmc input file from a file called genotpyes.txt.gz, which is
# a 01 file (ok, actually haplotypes) generated, for example, by macs. 

from __future__ import division
import gzip, sys, getopt, pdb
import numpy as np

##########################################################################################################

def parse_options():
    """
    File paths and indices (comma separated)
    """
    options ={ "macs_data_dir":"", "p0_indices":"", "p1_indices":"", "seq_len":0  , "chr":"chr", "output_suffix":""}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:0:1:s:c:o:", ["dir", "p0", "p1", "len", "chr", "os"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-d","--dir"]:       options["macs_data_dir"] = a
        elif o in ["-0","--p0"]:      options["p0_indices"] = [int(x) for x in a.split(",")]
        elif o in ["-1","--p1"]:      options["p1_indices"] = [int(x) for x in a.split(",")]
        elif o in ["-s","--seq_len"]: options["seq_len"] = int(a)
        elif o in ["-c","--chr"]: options["chr"] = a
        elif o in ["-o","--os"]: options["output_suffix"] = a

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):
    """
    Just load the appropriate indices, do the calculation and output. 
    """

    nhaps=None
    if len(options["p0_indices"])!=len(options["p1_indices"]):
        raise Exception("Must have same number of haplotypes in p0 and p1")
    else:
        nhaps=len(options["p0_indices"])

    pos_file=gzip.open(options["macs_data_dir"]+"/snps.pos.txt.gz")
    pos=pos_file.readlines()
    pos=[int(float(x)*options["seq_len"]) for x in pos[0][:-1].split()]
    pos_file.close()
    
    nsites=len(pos)
    p0_haps=np.zeros((nhaps, nsites), dtype=int)
    p1_haps=np.zeros((nhaps, nsites), dtype=int)

    gt_file=gzip.open(options["macs_data_dir"]+"/genotypes.txt.gz")
    i=0
    i0=0
    i1=0
    for line in gt_file:
        if i in options["p0_indices"]:
            p0_haps[i0]=[int(x) for x in line[:-1]]
            i0+=1
        if i in options["p1_indices"]:
            p1_haps[i1]=[int(x) for x in line[:-1]]
            i1+=1
        i+=1

    seg_sites=p0_haps[0]+p0_haps[1]+p1_haps[0]+p1_haps[1]
    seg_sites=np.array([(x<2*nhaps and x>0) for x in seg_sites])

    pos=np.array(pos)[seg_sites]
    p0_haps=p0_haps[:,seg_sites]
    p1_haps=p1_haps[:,seg_sites]

    out_file=open(options["macs_data_dir"]+"/msmc_input"+options["output_suffix"]+".txt", "w")
    last_pos=0
    for i in xrange(sum(seg_sites)):
        if pos[i]==last_pos:            # only one variant per site
            print "Skipping site at position "+str(pos[i])
        else:
            allele_str= "".join([str(x) for x in p0_haps[:,i]])+"".join([str(x) for x in p1_haps[:,i]])
            out_file.write("\t".join([options["chr"], str(pos[i]), str(pos[i]-last_pos), allele_str]))
            out_file.write("\n")
            last_pos=pos[i]                

##########################################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)
