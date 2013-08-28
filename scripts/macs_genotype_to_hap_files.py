from __future__ import division
import numpy as np
import gzip, sys, getopt

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "gt_file":"", "pos_file":"", "chr_len":0, "out_root":""  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:p:l:o:", ["gt", "pos", "len", "out"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-g","--gt"]:        options["gt_file"] = a
        elif o in ["-p","--pos"]:     options["pos_file"] = a
        elif o in ["-l","--len"]:     options["chr_len"] = float(a)
        elif o in ["-o","--out"]:     options["out_root"] = a

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):

    gt_file=gzip.open(options["gt_file"], "r")
    pos_file=gzip.open(options["pos_file"], "r")
    out_haps=gzip.open(options["out_root"]+"/haps.gz", "w")
    out_haps_f1=gzip.open(options["out_root"]+"/haps.f1.gz", "w")
    out_haps_f2=gzip.open(options["out_root"]+"/haps.f2.gz", "w")
    out_samples=open(options["out_root"]+"/samples.txt", "w")

    gt=np.genfromtxt(gt_file, delimiter=1)
    pos=np.genfromtxt(pos_file)
    pos=np.floor(pos*options["chr_len"]).astype(int)
    
    gt=gt.transpose().astype(int)
    # This is because on some platforms the np.genfromtxt tries to import the line endings...     
    gt=gt[range(len(pos)),]               
    
    (nsnp,nind)=gt.shape

    ACs=np.sum(gt, axis=1)
    MACs=np.minimum(ACs, nind-ACs)
    for i in range(nsnp):
        out_haps.write(("\t".join(["%d"]*(nind+1))+"\n")%((pos[i],)+tuple(gt[i,])))
        if MACs[i]==1:
            out_haps_f1.write(("\t".join(["%d"]*(nind+1))+"\n")%((pos[i],)+tuple(gt[i,])))
        if MACs[i]==2:
            out_haps_f2.write(("\t".join(["%d"]*(nind+1))+"\n")%((pos[i],)+tuple(gt[i,])))

    for i in range(int(nind/2)):
        out_samples.write("SIM%d\n"%(i+1,))
            
    for fil in [gt_file, pos_file, out_haps, out_haps_f1, out_haps_f2]:
        fil.close()
    
##########################################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)
