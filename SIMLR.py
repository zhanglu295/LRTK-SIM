import sys
from functools import partial
import multiprocessing
import numpy as np
import os
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.plotly as py
import pdb
import time
from collections import defaultdict
#pdb.set_trace()
#Parameter defination
large_droplet=4000000
large_template=1000000
class Molecule(object):
    def __init__(self,length,start,end):
        #self.seq=sequence.upper()
        self.length=length
        self.index_droplet=0
        self.barcode='Null'
        self.shortreads=0
        self.start=start
        self.end=end
class Qual_Substitute(object):
    def __init__(self,phred,change):
        self.phred=phred
        self.substitute=change
class Short_reads_PE(object):
     def __init__(self,seq1,qual1,start1,end1,seq2,qual2,start2,end2):
         self.seq1=seq1
         self.qual1=qual1
         self.start1=start1
         self.end1=end1
         self.seq2=seq2
         self.qual2=qual2
         self.start2=start2
         self.end2=end2
class Short_reads(object):
    def __init__(self,seq,qual):
        self.seq=seq
        self.qual=qual
class parameter(object):
    def __init__(self):
        self.CF=0
        self.CR=0
        self.N_FP=0
        self.Mu_F=0
        self.SR=0
        self.Mu_IS=0
        self.Std_IS=0
        self.Fasta='N'
        self.barcodepool='N'
        self.in_var1='N'
        self.in_var2='N'
        self.hap=1
        self.Seq_error='N'
        self.Seq_qual='N'
        self.Barcode_qual='N'
        self.Fast_mode='Y'
        self.Error_rate=0
#read parameters from file
def input_parameter(argv,parameter_struc):
    deter=1
    f=open(argv,"r")
    for line in f:
        Par=line.split('=')
        if len(Par)==2:
           if Par[0]=='Path_Fasta':
              parameter_struc.Fasta=Par[1].strip('\n')
           elif Par[0]=='Fast_mode':
              parameter_struc.Fast_mode=Par[1].strip('\n')
           elif Par[0]=='Seq_error':
              parameter_struc.Seq_error=Par[1].strip('\n')
           elif Par[0]=='Error_rate':
              parameter_struc.Error_rate=float(Par[1].strip('\n'))
           elif Par[0]=='Path_Seq_qual':
              parameter_struc.Seq_qual=Par[1].strip('\n')
           elif Par[0]=='Path_Barcode_qual':
              parameter_struc.Barcode_qual=Par[1].strip('\n')
           elif Par[0]=='CF':
              parameter_struc.CF=float(Par[1].strip('\n'))
           elif Par[0]=='Mu_IS':
              parameter_struc.Mu_IS=float(Par[1].strip('\n'))
           elif Par[0]=='Std_IS':
              parameter_struc.Std_IS=float(Par[1].strip('\n'))
           elif Par[0]=='CR':
              parameter_struc.CR=float(Par[1].strip('\n'))
           elif Par[0]=='N_FP':
              parameter_struc.N_FP=int(Par[1].strip('\n'))
           elif Par[0]=='Mu_F':
              parameter_struc.Mu_F=float(Par[1].strip('\n'))
           elif Par[0]=='SR':
              parameter_struc.SR=int(Par[1].strip('\n'))
           elif Par[0]=='Path_barcodepool':
              parameter_struc.barcodepool=Par[1].strip('\n')
           elif Par[0]=='Path_VCF1':
              parameter_struc.in_var1=Par[1].strip('\n')
           elif Par[0]=='Path_VCF2':
              parameter_struc.in_var2=Par[1].strip('\n')
           elif Par[0]=='Hap':
              parameter_struc.hap=int(Par[1].strip('\n'))
    f.close()
    if os.path.isfile(parameter_struc.Fasta)==False:
       deter=0
       print('template fasta file (Fasta) does not exist')
    if parameter_struc.Fasta=='N':
       deter=0
       print('Missing template fasta file (Fasta)')
    if parameter_struc.CF==0:
       deter=0
       print('Missing coverage for long fragment (CF)')
    if parameter_struc.CR==0:
       deter=0
       print('Missing coverage of short reads for long fragment (CR)')
    if parameter_struc.N_FP==0:
       deter=0
       print('Missing the average number of molecules for eachh droplet (N_FP)')
    if parameter_struc.Mu_F==0:
       deter=0
       print('Missing the average length for long fragment (Kb) (Mu_F)')
    if parameter_struc.SR==0:
       deter=0
       print('Missing length of short reads (bp) (SR)')
    if parameter_struc.Mu_IS==0:
       deter=0
       print('Missing mean of insert size of short reads (bp) (Mu_IS)')
    if parameter_struc.Std_IS==0:
       deter=0
       print('Missing standard deviation of insert size of short fragment (bp) (Std_IS)')
    if parameter_struc.barcodepool=='N':
       deter=0
       print('Missing barcode list (barcodepool)')
    if os.path.isfile(parameter_struc.barcodepool)==False:
       deter=0
       print('barcode list does not exist')
    if parameter_struc.in_var1!='N':
       if os.path.isfile(parameter_struc.in_var1)==False:
          deter=0
          print('VCF1 does not exist')
    if parameter_struc.in_var2!='N':
       if os.path.isfile(parameter_struc.in_var2)==False:
          deter=0
          print('VCF2 does not exist')
    return deter
#read template sequence
def input_seq(in_path):
    sequence=''
    f=open(in_path,"r")
    line_index=0
    for line in f:
        if line_index>0:
           sequence+=line.strip('\n')
        line_index=line_index+1
    f.close()
    return sequence

#sapmpling long fragments from empirical distribution
def randomlong(Par,seq):
    global N_frag
    N_frag=0
    #calculate required number of molecules
    lensingle=len(seq)
    N_frag=int(lensingle*Par.CF/(Par.Mu_F*1000))
    #randome simulate molecules
    MolSet=[]
    For_hist=[]
    for i in range(N_frag):
        start=int(np.random.uniform(low=0,high=lensingle))
        length=int(np.random.exponential(scale=Par.Mu_F*1000))
        end=start+length-1
        if length==0:
           continue
        if end>lensingle:
           Molseq=seq[start:lensingle]
           lengthnew=lensingle-start
           nPos=Molseq.count('N')
           if nPos>0:
              continue
           NewMol=Molecule(lengthnew,start,lensingle)
           MolSet.append(NewMol)
        else:
           Molseq=seq[start:end]
           nPos=Molseq.count('N')
           if nPos>0:
              continue
           NewMol=Molecule(length-1,start,end)
           MolSet.append(NewMol)
        For_hist.append(length/1000)
    N_frag=len(MolSet)
    #plt.hist(For_hist)
    #plt.xlabel('Molecule length')
    #plt.ylabel('Number of molecules')
    #plt.title('Histogram of molecule length (total '+str(len(For_hist))+' molecules)')
    #plt.show()
    #plt.savefig('Len_Molecule_hist.png')
    #plt.close()
    return MolSet
#assign long fragments to droplet
def deternumdroplet(N_frag,N_FP):
    frag_drop = np.random.poisson(N_FP, large_droplet)
    assign_drop=[]
    totalfrag=0
    for i in range(large_droplet):
        totalfrag=totalfrag+frag_drop[i]
        if totalfrag<=N_frag:
            assign_drop.append(frag_drop[i])
            Figure_len_molecule.append(frag_drop[i])
        else:
            last=N_frag-(totalfrag-frag_drop[i])
            assign_drop.append(last)
            Figure_len_molecule.append(last)
            break
    #plt.hist(assign_drop)
    #plt.xlabel('Number of molecules in droplets')
    #plt.ylabel('Number of molecules')
    #plt.title('Histogram of the number of molecules in droplets (total '+str(len(assign_drop))+' droplet)')
    #plt.show()
    #plt.savefig('Num_Molecule_hist.png')
    #plt.close()
    return assign_drop
#assign barcode to each fragment
def selectbarcode(pool,assign_drop,MolSet,droplet_container):
    #permute index of long fragment for random sampling
    permutnum=np.random.permutation(N_frag)
    #include barcode in the list
    barcode=[]
    N_droplet=len(assign_drop)
    t=1
    f=open(pool,"r")
    for line in f:
        barcode.append(line.strip('\n'))
        t=t+1
        if t>N_droplet:
            break
    f.close()
    start=0
    for i in range(N_droplet):
        num_molecule_per_partition=assign_drop[i]
        index_molecule=permutnum[start:start+num_molecule_per_partition]
        totalseqlen=0
        temp=[]
        start=start+num_molecule_per_partition
        for j in range(num_molecule_per_partition):
            index=index_molecule[j]
            temp.append(index)
            MolSet[index].index_droplet=i
            MolSet[index].barcode=barcode[i]
            totalseqlen=totalseqlen+MolSet[index].length
        droplet_container.append(temp)
    return MolSet
def diploid(Par,deter,lib):
    global Figure_len_molecule
    global Figure_num_molecule
    Figure_len_molecule=[]
    Figure_num_molecule=[]
    os.system('mkdir '+sys.argv[1]+'/lib'+str(lib)+'/hap1')
    os.system('mkdir '+sys.argv[1]+'/lib'+str(lib)+'/hap2')
    if Par.in_var1=='N' or Par.in_var2=='N':
       print('missing variants files, must provide both  VCF1 and VCF2')
    else:#insert variants from vcf
       os.system('tabix -f -p vcf '+Par.in_var1)
       os.system('tabix -f -p vcf '+Par.in_var2)
       os.system('cat '+Par.Fasta+' | vcf-consensus '+Par.in_var1+' > new1.fa')
       os.system('cat '+Par.Fasta+' | vcf-consensus '+Par.in_var2+' > new2.fa')
       for j in range(2):
           droplet_container=[]
           Par.Fasta='new'+str(j+1)+'.fa'
           seq=input_seq(Par.Fasta)
           #recode cut position of long fragment
           print('read template finished (hap'+str(j+1)+')')
           MolSet=randomlong(Par,seq)
           print('molecule sampling finished (hap'+str(j+1)+')')
           #calculate number of droplet
           assign_drop=deternumdroplet(N_frag,Par.N_FP)
           print('assign molecule to droplet finished (hap'+str(j+1)+')')
           #print(assign_drop)
           MolSet=selectbarcode(Par.barcodepool,assign_drop,MolSet,droplet_container)
           print('assign barcodes to molecules finished (hap'+str(j+1)+')')
           SIMSR(MolSet,Par,seq,lib)
           os.system('mv read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib)+'/hap'+str(j+1))
           os.system('mv read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib)+'/hap'+str(j+1))
           os.system('mv Num_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib)+'/hap'+str(j+1))
           os.system('mv Len_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib)+'/hap'+str(j+1))
       if Par.in_var1!='N':
          os.system('rm new1.fa')
       if Par.in_var2!='N':
          os.system('rm new2.fa')
       os.system('zcat '+sys.argv[1]+'/lib'+str(lib)+'/hap1/read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib)+'/hap2/read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq.gz>'+sys.argv[1]+'/lib'+str(lib)+'/read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq.gz')
       os.system('zcat '+sys.argv[1]+'/lib'+str(lib)+'/hap1/read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib)+'/hap2/read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq.gz>'+sys.argv[1]+'/lib'+str(lib)+'/read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq.gz')
       plt.hist(Figure_len_molecule)
       plt.xlabel('Molecule length')
       plt.ylabel('Number of molecules')
       plt.title('Histogram of molecule length (total '+str(len(Figure_len_molecule))+' molecules)')
       plt.show()
       plt.savefig('Len_Molecule_hist_'++'.png')
       plt.close()
       plt.hist(Figure_num_molecule)
       plt.xlabel('Number of molecules in droplets')
       plt.ylabel('Number of molecules')
       plt.title('Histogram of the number of molecules in droplets (total '+str(len(Figure_num_molecule))+' droplet)')
       plt.show()
       plt.savefig('Num_Molecule_hist.png')
       plt.close()
       os.system('mv Len_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib))
       os.system('mv Num_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib))
    print('library'+str(lib)+' simulation completed!')
    return None
def haploid(Par,deter,lib,hap):
    global Figure_len_molecule
    global Figure_num_molecule
    Figure_len_molecule=[]
    Figure_num_molecule=[]
    droplet_container=[]
    if hap=='1':
       if Par.in_var1!='N':
          os.system('tabix -f -p vcf '+Par.in_var1)
          os.system('cat '+Par.Fasta+' | vcf-consensus '+Par.in_var1+' > new_hap'+str(hap)+'_lib'+lib+'_.fa')
          Par.Fasta='new_hap'+str(hap)+'_lib'+lib+'.fa'
    if hap=='2':
       if Par.in_var1!='N':
          os.system('tabix -f -p vcf '+Par.in_var2)
          os.system('cat '+Par.Fasta+' | vcf-consensus '+Par.in_var2+' > new_hap'+str(hap)+'_lib'+lib+'_.fa')
          Par.Fasta='new_hap'+hap+'_lib'+lib+'.fa'
    seq=input_seq(Par.Fasta)
    #recode cut position of long fragment
    print('read template finished (haplotype '+hap+' in library '+lib+')')
    MolSet=randomlong(Par,seq)
    print('molecule sampling finished (haplotype '+hap+' in library '+lib+')')
    #calculate number of droplet
    assign_drop=deternumdroplet(N_frag,Par.N_FP)
    print('assign molecule to droplet finished (haplotype '+hap+' in library '+lib+')')
    MolSet=selectbarcode(Par.barcodepool,assign_drop,MolSet,droplet_container)
    print('assign barcode to molecule finished (haplotype '+hap+' in library '+lib+')')
    SIMSR(MolSet,Par,seq,lib,hap)
    os.system('mv read-RA_si-CCTGGAGA_lib-00'+lib+'-hap-00'+hap+'.fastq.gz '+sys.argv[1]+'/lib'+lib)
    os.system('mv read-I1_si-CCTGGAGA_lib-00'+lib+'-hap-00'+hap+'.fastq.gz '+sys.argv[1]+'/lib'+lib)
    plt.hist(Figure_len_molecule)
    plt.xlabel('Molecule length')
    plt.ylabel('Number of molecules')
    plt.title('Histogram of molecule length (total '+str(len(Figure_len_molecule))+' molecules)')
    plt.show()
    plt.savefig('Len_Molecule_hist_lib'+lib+'_hap'+hap+'.png')
    plt.close()
    plt.hist(Figure_num_molecule)
    plt.xlabel('Number of molecules in droplets')
    plt.ylabel('Number of molecules')
    plt.title('Histogram of the number of molecules in droplets (total '+str(len(Figure_num_molecule))+' droplet)')
    plt.show()
    plt.savefig('Num_Molecule_hist_lib'+lib+'_hap'+hap+'.png')
    plt.close()
    os.system('mv Num_Molecule_hist_lib'+lib+'_hap'+hap+'.png '+sys.argv[1]+'/lib'+lib)
    os.system('mv Len_Molecule_hist_lib'+lib+'_hap'+hap+'.png '+sys.argv[1]+'/lib'+lib)
    if hap==1 and Par.in_var1!='N':
       os.system('rm new1.fa')
    if hap==2 and Par.in_var2!='N':
       os.system('rm new2.fa')
    print('Haplotype '+hap+' in library '+lib+' simulation completed!')
    return None
def reverseq(seq):
    complementary=''
    rev_complementary=''
    for i in range(len(seq)):
        if seq[i]=='A':
           complementary+='T'
        elif seq[i]=='T':
           complementary+='A'
        elif seq[i]=='C':
           complementary+='G'
        elif seq[i]=='G':
           complementary+='C'
    rev_complementary=complementary[::-1]
    return rev_complementary

def Input_BarcodeQual(Par):
    f=open(Par.Barcode_qual,"r")
    line_index=0
    position=0
    Qual_dict=defaultdict(list)
    Prob_dict=defaultdict(list)
    #pos_qual={}
    for line in f:
        if line_index>0:
           linequal=line.strip('\t,\n')
           qualarray=linequal.split('\t')
           Qual_dict[qualarray[0]].append(ord(qualarray[1]))
           Prob_dict[qualarray[0]].append(float(qualarray[2]))
          # if position!=int(qualarray[0]):
          #    Qual_dict[str(position)]=pos_qual
          #    pos_qual={}
          #    position=int(qualarray[0])
           #for i in range(6):
           #    change.append(int(qualarray[i+2]))
          # pos_qual[ord(qualarray[1])]=change
        line_index=line_index+1
    #Qual_dict[str(position)]=pos_qual
    f.close()
    return Qual_dict,Prob_dict
def Input_SeqQual(Par):
    f=open(Par.Seq_qual,"r")
    line_index=0
    position=0
    Qual_dict=defaultdict(list)
    Prob_dict=defaultdict(list)
    Substitute_dict=defaultdict(list)
    for line in f:
        if line_index>0:
           change=[]
           linequal=line.strip('\t,\n')
           qualarray=linequal.split('\t')
          # print(qualarray[0])
          # print(qualarray[1])
           Qual_dict[qualarray[0]].append(ord(qualarray[1]))
           Prob_dict[qualarray[0]].append(float(qualarray[2]))
           Substitute_dict[(qualarray[0],ord(qualarray[1]))]=list(map(float,qualarray[3:]))
        line_index=line_index+1
    f.close()
    return Qual_dict,Prob_dict,Substitute_dict
def SIMSR(MolSet,Par,seq,lib,hap):
    f_reads = open('read-RA_si-CCTGGAGA_lib-00'+lib+'-hap-00'+hap+'.fastq', "w")
    f_sample = open('read-I1_si-CCTGGAGA_lib-00'+lib+'-hap-00'+hap+'.fastq', "w")
    [SeqQual_dict,SeqProb_dict,SeqSubstitute_dict]=Input_SeqQual(Par)
    [BarcodeQual_dict,BarcodeProb_dict]=Input_BarcodeQual(Par)

    #eva_num_reads=[]
    #for i in range(len(MolSet)):
    #    eva_num_reads.append(int(int(MolSet[i].length/(Par.SR*2))*Par.CR))
    #Seq_qual=np.zeros((sum(eva_num_reads)*2,Par.SR),dtype=int)
    #Barcode_qual=np.zeros((sum(eva_num_reads),16),dtype=int)
    #for m in range(16):
    #   Barcode_coll_phred=BarcodeQual_dict[str(m)]
    #   Barcode_coll_pos_qual=[]
    #   Barcode_qual[:,m]=np.random.choice(list(Barcode_coll_phred.keys()),size=(sum(eva_num_reads)))
    #for m in range(Par.SR):
    #   Seq_coll_phred=SeqQual_dict[str(m)]
    #   Seq_coll_pos_qual=[]
    #   Seq_qual[:,m]=np.random.choice(list(Seq_coll_phred.keys()),size=(sum(eva_num_reads)*2))
    last_reads=0
   # print(len(MolSet))
    for i in range(len(MolSet)):
        print('Finish  '+str(i+1)+'/'+str(len(MolSet)),end="\r")
        Seq_rand_qual=[]
        Barcode_rand_qual=[]
        num_reads=int(int(MolSet[i].length/(Par.SR*2))*Par.CR)
        if num_reads==0:
            continue
        All_forward=seq[MolSet[i].start:MolSet[i].end].upper()
        All_reverse=reverseq(All_forward)
        insert_size=np.random.normal(loc=Par.Mu_IS, scale=Par.Std_IS,size=num_reads)
        Totalreads=[]
        new_reads=last_reads+num_reads
        #Barcode_new_qual=Barcode_qual[last_reads:new_reads,:]
        #Seq_new_qual=Seq_qual[last_reads:new_reads,:]
        Seq_new_qual=np.zeros((num_reads*2,Par.SR),dtype=int)
        Barcode_new_qual=np.zeros((num_reads,16),dtype=int)
        for m in range(16):
            Barcode_coll_phred=BarcodeQual_dict[str(m)]
            Barcode_coll_prob=BarcodeProb_dict[str(m)]
            Barcode_new_qual[:,m]=np.random.choice(Barcode_coll_phred,p=Barcode_coll_prob,size=(num_reads))
        for m in range(Par.SR):
            Seq_coll_phred=SeqQual_dict[str(m)]
            Seq_coll_prob=SeqProb_dict[str(m)]
            Seq_new_qual[:,m]=np.random.choice(Seq_coll_phred,p=Seq_coll_prob,size=(num_reads*2))
       # pool = multiprocessing.Pool(5)
       # for j in range(num_reads):
       #     pool.apply_async(pairend,(Par,insert_size,MolSet[i],Barcode_new_qual,Seq_new_qual,All_reverse,j,BarcodeQual_dict,SeqQual_dict,),callback=Totalreads.append)
       # pool.close()
       # pool.join()
        for j in range(num_reads):
            PE_read=pairend(Par,insert_size,MolSet[i],Barcode_new_qual,Seq_new_qual,All_forward,All_reverse,j,SeqSubstitute_dict)
            Totalreads.append(PE_read)
        for j in range(len(Totalreads)):
            read1seq=Totalreads[j].seq1
            read1qual=Totalreads[j].qual1
            read2seq=Totalreads[j].seq2
            read2qual=Totalreads[j].qual2
            readname='@ST-K00126:'+str(i+1)+':H5W53BBXX:'+str(MolSet[i].start)+':'+str(MolSet[i].end)+':'+str(Totalreads[j].start1)+':'+str(Totalreads[j].end1)
            f_reads.write(readname+' 1:N:0\n')
            f_reads.write(read1seq+'\n')
            f_reads.write('+\n')
            f_reads.write(read1qual)
            f_reads.write('\n')
            f_reads.write(readname+' 3:N:0\n')
            f_reads.write(read2seq+'\n')
            f_reads.write('+\n')
            f_reads.write(read2qual+'\n')
            f_sample.write(readname+' 2:N:0\n')
            f_sample.write('CCTGGAGA\n')
            f_sample.write('+\n')
            f_sample.write('AAFFFKKK\n')
    f_reads.close()
    f_sample.close()
    os.system('gzip -f read-RA_si-CCTGGAGA_lib-00'+str(lib)+'-hap-00'+hap+'.fastq')
    os.system('gzip -f read-I1_si-CCTGGAGA_lib-00'+str(lib)+'-hap-00'+hap+'.fastq')
    return None
def pairend(Par,insert_size,MolSetX,Barcode_rand_qual,Seq_rand_qual,All_forward,All_reverse,index,SeqSubstitute_dict):
    is_read=int(np.absolute(insert_size[index]))
    start_for=int(np.random.uniform(low=1,high=MolSetX.length-Par.SR-1))
    if start_for+is_read+1>MolSetX.length:
       is_read=MolSetX.length-start_for-1
    end_for=start_for+int(is_read)
    forward_seq=All_forward[start_for:end_for]
    start_rev=MolSetX.length-end_for
    end_rev=start_rev+int(is_read)
    reverse_seq=All_reverse[start_rev:end_rev]
    read1=forward_seq[23:Par.SR]
    read2=reverse_seq[0:Par.SR]
    read1seq=''
    read2seq=''
    read1qual=''
    read2qual=''
    if Par.Fast_mode=='N':
       #Sim_seq1=SIMSeqQual(SeqQual_dict,read1,Par)
       #Sim_seq2=SIMSeqQual(SeqQual_dict,read2,Par)
       #Sim_barcode=SIMBarcodeQual(BarcodeQual_dict,MolSetX.barcode,Par)
       #read1seq=MolSetX.barcode+'NNNNNNN'+Sim_seq1.seq
       #read2seq=Sim_seq2.seq
       #read1qual=Sim_barcode.qual+'KKKKKKK'+Sim_seq1.qual
       #read2qual=Sim_seq2.qual
       if Par.Seq_error=='N':
          read1seq=MolSetX.barcode+'NNNNNNN'+read1
          read2seq=read2
       if Par.Seq_error=='Y':
          readerror=np.random.choice([0,1],p=[1-Par.Error_rate,Par.Error_rate],size=(Par.SR*2))
          error1=readerror[0:Par.SR].nonzero()
          error2=readerror[Par.SR:Par.SR*2].nonzero()
          read1new=list(read1)
          read2new=list(read2)
          for i in range(len(error1[0])):
               #print((np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])])))
               if read1[i]=='A':
                  rand_nuc=np.random.choice(['C','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][0:4]))
               elif read1[i]=='C':
                  rand_nuc=np.random.choice(['A','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][4:8]))
               elif read1[i]=='G':
                  rand_nuc=np.random.choice(['A','C','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][8:12]))
               elif read1[i]=='T':
                  rand_nuc=np.random.choice(['A','C','G','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][12:16]))
               read1new[i]=rand_nuc[0]

          for i in range(len(error2[0])):
               if read2[i]=='A':
                  rand_nuc=np.random.choice(['C','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][0:4]))
               elif read2[i]=='C':
                  rand_nuc=np.random.choice(['A','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][4:8]))
               elif read2[i]=='G':
                  rand_nuc=np.random.choice(['A','C','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][8:12]))
               elif read2[i]=='T':
                  rand_nuc=np.random.choice(['A','C','G','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual[index,i])][12:16]))
               read2new[i]=rand_nuc[0]
          read1seq=MolSetX.barcode+'NNNNNNN'+''.join(read1new)
          read2seq=''.join(read2new)

       read1qual=''.join(map(chr,Barcode_rand_qual[index,:]))+'KKKKKKK'+''.join(map(chr,Seq_rand_qual[index,23:Par.SR]))
       read2qual=''.join(map(chr,Seq_rand_qual[index,:]))
    else:
       read1seq=MolSetX.barcode+'NNNNNNN'+read1
       read2seq=read2
       read1qual='K'*Par.SR
       read2qual='K'*Par.SR
    return Short_reads_PE(read1seq,read1qual,start_for,end_for,read2seq,read2qual,start_rev,end_rev)
def main():
    pool = multiprocessing.Pool(5)
    os.system('rm -rf '+sys.argv[1]+'/lib*')
    list=os.listdir(sys.argv[1])
    Par=parameter()
    hap=1
    for i in range(len(list)):
        print('processing library '+str(i+1)+' for '+list[i])
        os.system('mkdir '+sys.argv[1]+'/lib'+str(i+1))
        deter=input_parameter(sys.argv[1]+'/'+list[i],Par)
        if deter==1:
           if Par.hap==1:
              pool.apply_async(haploid,(Par,deter,str(i+1),str(1)))
              #haploid(Par,deter,str(i+1),str(Par.hap))
           elif Par.hap==2:
               for j in range(2):
                    pool.apply_async(haploid,(Par,deter,str(i+1),str(j+1)))
              #diploid(Par,deter,i+1)
    pool.close()
    pool.join()
    return None
if __name__=="__main__":
    main()
