from __future__ import division
import numpy as np
import gzip, sys, getopt

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "gt_file":"", "pos_file":"", "chr_len":0, "out_root":"", "from":1, "to":2  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:p:l:o:f:t:", ["gt", "pos", "len", "out", "from", "to"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-g","--gt"]:        options["gt_file"] = a
        elif o in ["-p","--pos"]:     options["pos_file"] = a
        elif o in ["-l","--len"]:     options["chr_len"] = float(a)
        elif o in ["-o","--out"]:     options["out_root"] = a
        elif o in ["-f","--from"]:    options["from"] = int(a)
        elif o in ["-t","--to"]:      options["to"] = int(a)

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):

    freq_range=range(options["from"], options["to"]+1)
    
    gt_file=gzip.open(options["gt_file"], "r")
    pos_file=gzip.open(options["pos_file"], "r")
    out_haps=gzip.open(options["out_root"]+"/haps.gz", "w")
    out_haps_fn=[gzip.open(options["out_root"]+"/haps.f"+str(x)+".gz", "w") for x in freq_range]

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
        if MACs[i]>=options["from"] and MACs[i]<= options["to"]:
            idx=MACs[i]-options["from"]
            out_haps_fn[idx].write(("\t".join(["%d"]*(nind+1))+"\n")%((pos[i],)+tuple(gt[i,])))

    for i in range(int(nind/2)):
        out_samples.write("SIM%d\n"%(i+1,))
            
    for fil in [gt_file, pos_file, out_haps]+out_haps_fn:
        fil.close()
    
##########################################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)
