# Given a flat haplotype file, (i.e. generated with something like
# awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'
# split it up into a directory so that there is one file called pos, and then one for each of the individuals
# with just a vector of genotypes, separated by spaces.

from __future__ import division
import gzip, numpy, sys, getopt,resource

REP_INT=1e3                     # Report every REP_INT lines
# resource.setrlimit(resource.RLIMIT_NOFILE, (2000,-1)) # just to be safe

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "haplotype_file":"", "out":"", "n":0  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:o:s:", ["haplotype_file","out_dir","samples"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-h","--hap"]:           options["haplotype_file"] = a
        elif o in ["-o","--out"]:         options["out_dir"] = a
        elif o in ["-s","--samples"]:     options["samples"] = a

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):
    """
    main
    """
    # load sample names, one per line
    sample_file=open(options["samples"], "r")
    sample_names=sample_file.readlines()
    sample_names=[x[:-1] for x in sample_names] # strip newlines
    sample_file.close()
    nsamples=len(sample_names)

    if options["haplotype_file"][-3:]==".gz":
        haps_in=gzip.open(options["haplotype_file"], "r")
    else:
        haps_in=open(options["haplotype_file"], "r")

    output_filenames=[options["out_dir"]+"/pos.gz"]+[options["out_dir"]+"/"+n+".gt.gz" for n in sample_names]
    output_files=[gzip.open(f,"w") for f in output_filenames]

    line_count=1
    for line in haps_in:
        if not line_count%REP_INT:
            print "\rLine " + str(line_count),
            sys.stdout.flush()

        bits=line.split()
        bits=[int(x) for x in bits]
        if line_count==1:               # Just check first line is the right length
            if len(bits)!=1+2*nsamples:
                raise Exception("length of line 1 is 1+2*length of samples")
            output_files[0].write("%d"%(bits[0],))
            for i in range(nsamples):
                output_files[1+i].write("%d"%(bits[1+2*i]+bits[2+2*i],))
        else:
            output_files[0].write(" %d"%(bits[0],))
            for i in range(nsamples):
                output_files[1+i].write(" %d"%(bits[1+2*i]+bits[2+2*i],))

        line_count+=1

    [f.close() for f in output_files]
    haps_in.close()
        
    return

##########################################################################################################


if __name__=="__main__":
    options=parse_options()
    main(options)



