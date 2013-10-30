# Run the simulations from the schaffner paper, using macs,
# as interpreted by Jared O'Connell. 

import subprocess,sys
from math import log
import random

def bottleneck(pop,coef,gen,N,N0):
    assert coef>0. and coef <1.
    N1 = -1.0 / (2.0 *  log(1.0 - coef))
    print N1
    T1 = str((gen-1)/(4.*N0))
    T2 = str((gen)/(4.*N0))
    prefix = ""
    prefix += " -en " + T1 + " " + str(pop) + " " + str(float(N1)/float(N0))
    prefix +=  " -en " + T2 + " " + str(pop) + " " + str(float(N)/float(N0))
    return prefix

def change_size(pop,gen,N,N0):
    return " -en " + str(float(gen)/(4.*N0)) + " " + str(pop) + " " + str(float(N)/float(N0))

def split(srcpop,dstpop,gen,N0):
    return " -ej " + str(float(gen)/(4.*N0)) + " " + str(dstpop) + " " + str(srcpop)

def migration_rate(srcpop,dstpop,gen,rate,N0):
    m = str(4.*N0*rate)
    T = str(gen/(4.*N0))
    return    " -em "+T+" "+str(dstpop)+" "+str(srcpop)+" "+m

def migration_rate0(srcpop,dstpop,rate,N0):
    m = str(4.*N0*rate)
    return    " -m "+str(dstpop)+" "+str(srcpop)+" "+m

if __name__=="__main__":

    if len(sys.argv)!=7:
        print "Usage python schaffner.py N L macsmap.txt seed output_prefix macs_path"
        quit()

    prefix = sys.argv[5]
    macs_path = sys.argv[6]

    nhp=int(sys.argv[1])
    nsample1=nhp
    nsample2=nhp
    nsample3=nhp
    nsample4=nhp
 
                     
    # 1 africans
    # 2 europeans
    # 3 asians
    # 4 african-americans

    nsample = nsample1+nsample2+nsample3+nsample4
    assert nsample>1

    N0=100000.
    t=4.*N0
    L=int(sys.argv[2])
    mutation_rate=1.25e-8
    theta= 4 * N0 * mutation_rate
    args = " "+str(nsample)+" "+str(L)
    args+= " -t "+str(theta)
    args += " -r " + str(4*N0*1e-8)

    #POPS
    args += " -h 100 -I 4 %d %d %d %d"%(nsample1,nsample2,nsample3,nsample4)

    #africans
    args += change_size(1,17000.,12500.,N0) # "african pop size" 
    args += change_size(1,200.,24000.,N0) # "agriculture - african" 
    args += bottleneck(1,.008,1997,12500.,N0) #"african bottleneck"

    #europeans
    args += split(1,2,3500.,N0) #"out of Africa"
    args += bottleneck(2,.085,3499.,7700,N0)# "OoA bottleneck"
    args += bottleneck(2,.02,1999.,7700,N0)# "european bottleneck"
    args += change_size(2,350,7700,N0)# "agriculture - european"

    #asians
    args +=  split(2,3,2000,N0) #"asian and european split"
    args +=  bottleneck(3,.06,1995,7700,N0) # "asian bottleneck"
    args +=  change_size(3,400,7700,N0)# "agriculture - asian
                   
    #african-americans
    args +=    split(1,4,7,N0)        

    #MIGRATION
    args += migration_rate0(1,2,.000032,N0)#"afr->eur migration"
    args += migration_rate0(2,1,.000032,N0)#eur -> afr
    args += migration_rate0(1,3,.000008,N0)#afr -> as
    args += migration_rate0(3,1,.000008,N0)#as -> afr

    args += migration_rate(1,2,1993,0.0,N0) #"afr->eur migration"
    args += migration_rate(2,1,1992,0.0,N0) #eur -> afr
    args += migration_rate(1,3,1991,0.0,N0)#afr -> as
    args += migration_rate(3,1,1990,0.0,N0)#as -> afr



    args += " -R " +sys.argv[3]
    args += " -s "+sys.argv[4]

    print "THETA =",theta
#    args += " -F ascertainment.txt 1"
#    print("\n\n./macs"+args)
#    print "./macs"+args#+" 2> trees.txt 1> haplotypes.txt"
    print " ".join([macs_path + "/macs " + args+" -T 2> " +prefix+"/trees.txt | " + macs_path + "/msformatter | gzip -cf > " + prefix + "/haplotypes.txt.gz"])
    print 
    # subprocess.call([macs_path + "/macs " + args+" -T 2> " +prefix+".trees.txt 1> "+prefix+".haplotypes.txt"],shell=True)

    subprocess.call([macs_path + "/macs " + args+" -T 2> " +prefix+"/trees.txt | " + macs_path + "/msformatter | gzip -cf > " + prefix + "/haplotypes.txt.gz"],shell=True)

    quit(0)
