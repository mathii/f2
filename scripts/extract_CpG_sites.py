#Read a vcf and write a list of the CpG and non CpG sites.
#Only the C bases 
#Only biallelic snps

from __future__ import division, print_function
import gzip, sys, getopt
from pyfaidx import Fasta

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "ref":None, "out":None }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:r:o:", ["vcf", "ref", "out"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:           options["vcf"] = a
        elif o in ["-o","--out"]:         options["out"] = a
        elif o in ["-r","--ref"]:     options["ref"] = a

    print "found options:"
    print options

    return options

##########################################################################################################

options=parse_options()
CpG_out=open(options["out"]+"CpG.txt", "w")
nonCpG_out=open(options["out"]+"nonCpG.txt", "w")

vcf_file=gzip.open(options["vcf"])
reference=Fasta(options["ref"])
for line in vcf_file:
    if line[0]=="#":
        continue

    bits=line.split()
    CHR=bits[0]
    POS=int(bits[1])
    REF,ALT=bits[2:3]

    if len(REF)!=1 or len(ALT)!=1:
        continue

    try:
        ref_context=reference[chr][(pos-1):(pos+1)].seq.upper()
        if ref_context[0]!=REF:
            raise Exception("Reference alelles don't match")
        if ((REF=="C" and ALT=="T") or (REF=="T" and ALT=="C")) and ref_context=="CG":
            CpG_out.write("%s\t%d\n"%(CHR. POS))
        else:
            nonCpG_out.write("%s\t%d\n"%(CHR. POS))
            
    except IndexError:
        continue
        

