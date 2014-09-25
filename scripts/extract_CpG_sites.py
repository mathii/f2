#Read a vcf and write a list of the CpG and non CpG sites.
#Only the C bases 
#Only biallelic snps

from __future__ import division, print_function
import gzip, sys, getopt, pdb
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
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:           options["vcf"] = a
        elif o in ["-o","--out"]:         options["out"] = a
        elif o in ["-r","--ref"]:     options["ref"] = a

    print("found options:")
    print(options)

    return options

##########################################################################################################

options=parse_options()
CpG_out=open(options["out"]+"CpG.txt", "w")
nonCpG_out=open(options["out"]+"nonCpG.txt", "w")

vcf_file=gzip.open(options["vcf"])
reference=Fasta(options["ref"])
i=0
for line in vcf_file:
    if line[0]=="#":
        continue

    if not i%1000000:
        print(str(i))
    i+=1

    bits=line.split()
    CHR=bits[0]
    POS=int(bits[1])
    REF,ALT=bits[3:5]

    if len(REF)!=1 or len(ALT)!=1:
        continue

    try:
        ref_base=reference[CHR][POS-1].seq.upper()
        if ref_base!=REF:
            raise Exception("Reference alelles don't match at "+str(CHR)+" "+str(POS))
        if ((REF=="C" and ALT=="T") or (REF=="T" and ALT=="C")):
            base_after=reference[CHR][POS].seq.upper()
            if base_after=="G":
                CpG_out.write("%s\t%d\n"%(CHR, POS))
            else:
                nonCpG_out.write("%s\t%d\n"%(CHR, POS))

        if ((REF=="G" and ALT=="A") or (REF=="A" and ALT=="G")):
            base_before=reference[CHR][POS-2].seq.upper()
            if base_before=="C":
                CpG_out.write("%s\t%d\n"%(CHR, POS))
            else:
                nonCpG_out.write("%s\t%d\n"%(CHR, POS))
                                                        
    except IndexError:
        continue
        

