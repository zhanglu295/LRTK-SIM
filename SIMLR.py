import sys
import numpy as np
import os
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.plotly as py
import pdb
#pdb.set_trace()
#Parameter defination
large_droplet=4000000
large_template=1000000
class Molecule(object):
    def __init__(self,sequence,length,start,end):
        self.seq=sequence.upper()
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
           NewMol=Molecule(Molseq,lengthnew,start,end)
           MolSet.append(NewMol)
        else:
           Molseq=seq[start:end]
           nPos=Molseq.count('N')
           if nPos>0:
              continue
           NewMol=Molecule(Molseq,length-1,start,end)
           MolSet.append(NewMol)
        For_hist.append(length/1000)
    N_frag=len(MolSet)
    plt.hist(For_hist)
    plt.xlabel('Molecule length')
    plt.ylabel('Number of molecules')
    plt.title('Histogram of molecule length (total '+str(len(For_hist))+' molecules)')
    plt.show()
    plt.savefig('Len_Molecule_hist.png')
    plt.close()
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
    plt.hist(assign_drop)
    plt.xlabel('Number of molecules in droplets')
    plt.ylabel('Number of molecules')
    plt.title('Histogram of the number of molecules in droplets (total '+str(len(assign_drop))+' droplet)')
    plt.show()
    plt.savefig('Num_Molecule_hist.png')
    plt.close()
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
    newinput=''
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
           MolSet=randomlong(Par,seq)
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
       plt.savefig('Len_Molecule_hist.png')
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
def haploid(Par,deter,lib):
    global Figure_len_molecule
    global Figure_num_molecule
    Figure_len_molecule=[]
    Figure_num_molecule=[]
    droplet_container=[]
    newinput=''
    if Par.in_var1=='N':
       newinput=Par.Fasta
    else:#insert variants from vcf
       os.system('tabix -f -p vcf '+Par.in_var1)
       os.system('cat '+Par.Fasta+' | vcf-consensus '+Par.in_var1+' > new.fa')
       Par.Fasta='new.fa'
       print('generate new fasta')

    seq=input_seq(Par.Fasta)
    #recode cut position of long fragment
    print('read template finished')
    MolSet=randomlong(Par,seq)
    print('molecule sampling finished')
    #calculate number of droplet
    assign_drop=deternumdroplet(N_frag,Par.N_FP)
    print('assign molecule to droplet finished')
    MolSet=selectbarcode(Par.barcodepool,assign_drop,MolSet,droplet_container)
    print('assign barcode to molecule finished')
   # print(MolSet[0].barcode)
    SIMSR(MolSet,Par)
    os.system('mv read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib))
    os.system('mv read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq.gz '+sys.argv[1]+'/lib'+str(lib))
    os.system('mv Num_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib))
    os.system('mv Len_Molecule_hist.png '+sys.argv[1]+'/lib'+str(lib))
    if Par.in_var1!='N':
       os.system('rm new.fa')
    print('library'+str(lib)+' simulation completed!')
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
    Qual_dict={}
    pos_qual=[]
    for line in f:
        if line_index>0:
           change=[]
           linequal=line.strip('\t,\n')
           qualarray=linequal.split('\t')
           if position!=int(qualarray[0]):
              Qual_dict[str(position)]=pos_qual
              pos_qual=[]
              position=int(qualarray[0])
           for i in range(5):
               change.append(int(qualarray[i+2]))
           error_profile_line=Qual_Substitute(int(qualarray[1]),change)
           pos_qual.append(error_profile_line)
        line_index=line_index+1
    Qual_dict[str(position)]=pos_qual
    f.close()
    return Qual_dict
def Input_SeqQual(Par):
    f=open(Par.Seq_qual,"r")
    line_index=0
    position=0
    Qual_dict={}
    pos_qual=[]
    for line in f:
        if line_index>0:
           change=[]
           linequal=line.strip('\t,\n')
           qualarray=linequal.split('\t')
           if position!=int(qualarray[0]):
              Qual_dict[str(position)]=pos_qual
              pos_qual=[]
              position=int(qualarray[0])
           for i in range(20):
               change.append(int(qualarray[i+2]))
           error_profile_line=Qual_Substitute(int(qualarray[1]),change)
           pos_qual.append(error_profile_line)
        line_index=line_index+1
    Qual_dict[str(position)]=pos_qual
    f.close()
    return Qual_dict
def deter_nucleo(Seq_Nuc,Num_Nuc):
    Num_Nuc.insert(0,0)
    Seq_Nuc.insert(0,'X')
    length=len(Num_Nuc)
    real_nuc=''
    rand_nuc=np.random.random_integers(low=0, high=Num_Nuc[length-1])
    for i in range(length-1):
        if rand_nuc>= Num_Nuc[i] and rand_nuc<=Num_Nuc[i+1]:
            real_nuc=Seq_Nuc[i+1]
            break
    return real_nuc
def SIMBarcodeQual(Qual_dict,Sequence,Par):
    length=len(Sequence)
    Qual=''
    Sim_fastq=Short_reads(Sequence,Qual)
    Listseq=list(Sim_fastq.seq)
    for i in range(length):
        Pos_qual=Qual_dict.get(str(i))
        validate=0
        rand_qual=0
        while validate==0:
            rand_qual=np.random.random_integers(low=0, high=len(Pos_qual)-1)
            line_select=Pos_qual[rand_qual]
            if Sequence[i]=='A':
               if line_select.substitute[0]!=0:
                  validate=1
            elif Sequence[i]=='T':
               if line_select.substitute[3]!=0:
                       validate=1
            elif Sequence[i]=='C':
               if line_select.substitute[1]!=0:
                       validate=1
            elif Sequence[i]=='G':
               if  line_select.substitute[2]!=0:
                        validate=1
            elif Sequence[i]=='N':
                    validate=1
        Sim_fastq.seq="".join(Listseq)
        Sim_fastq.qual= Sim_fastq.qual+chr(Pos_qual[rand_qual].phred+33)
    return Sim_fastq
def SIMSeqQual(Qual_dict,Sequence,Par):
    length=len(Sequence)
    Qual=''
    Sim_fastq=Short_reads(Sequence,Qual)
    Listseq=list(Sim_fastq.seq)
    if Par.Seq_error=='Y':
        for i in range(length):
            Pos_qual=Qual_dict.get(str(i))
            validate=0
            rand_qual=0
            while validate==0:
                rand_qual=np.random.random_integers(low=0, high=len(Pos_qual)-1)
                line_select=Pos_qual[rand_qual]
                if Sequence[i]=='A':
                    if line_select.substitute[4]!=0:
                       validate=1
                       real_nuc=deter_nucleo(['A','C','G','T','N'],line_select.substitute[0:5])
                       Listseq[i]=real_nuc
                if Sequence[i]=='T':
                    if line_select.substitute[19]!=0:
                       validate=1
                       real_nuc=deter_nucleo(['A','C','G','T','N'],line_select.substitute[15:20])
                       Listseq[i]=real_nuc
                if Sequence[i]=='C':
                    if line_select.substitute[9]!=0:
                       validate=1
                       real_nuc=deter_nucleo(['A','C','G','T','N'],line_select.substitute[5:10])
                       Listseq[i]=real_nuc
                if Sequence[i]=='G':
                    if line_select.substitute[14]!=0:
                       validate=1
                       real_nuc=deter_nucleo(['A','C','G','T','N'],line_select.substitute[10:15])
                       Listseq[i]=real_nuc
                if Sequence[i]=='N':
                    validate=1
            Sim_fastq.qual=Sim_fastq.qual+chr(Pos_qual[rand_qual].phred+33)
    if Par.Seq_error=='N':
        for i in range(length):
            Pos_qual=Qual_dict.get(str(i))
            validate=0
            rand_qual=0
            while validate==0:
                rand_qual=np.random.random_integers(low=0, high=len(Pos_qual)-1)
                line_select=Pos_qual[rand_qual]
                if Sequence[i]=='A':
                    if line_select.substitute[4]!=0:
                       validate=1
                elif Sequence[i]=='T':
                    if line_select.substitute[19]!=0:
                       validate=1
                elif Sequence[i]=='C':
                    if line_select.substitute[9]!=0:
                       validate=1
                elif Sequence[i]=='G':
                    if  line_select.substitute[14]!=0:
                        validate=1
                elif Sequence[i]=='N':
                    validate=1
            Sim_fastq.seq="".join(Listseq)
            Sim_fastq.qual= Sim_fastq.qual+chr(Pos_qual[rand_qual].phred+33)
    return Sim_fastq
def SIMSR(MolSet,Par):
    f_reads = open('read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq', "w")
    f_sample = open('read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq', "w")
    SeqQual_dict=Input_SeqQual(Par)
    BarcodeQual_dict=Input_BarcodeQual(Par)
    for i in range(len(MolSet)):
        print('Finish  '+str(i+1)+'/'+str(len(MolSet)),end="\r")
        num_reads=int(int(MolSet[i].length/(Par.SR*2))*Par.CR)
        All_reverse=reverseq(MolSet[i].seq)
        insert_size=np.random.normal(loc=Par.Mu_IS, scale=Par.Std_IS,size=num_reads)
        for j in range(num_reads):
            is_read=int(np.absolute(insert_size[j]))
            start_for=int(np.random.uniform(low=1,high=MolSet[i].length-Par.SR-1))
            if start_for+is_read+1>MolSet[i].length:
               is_read=MolSet[i].length-start_for-1
            end_for=start_for+int(is_read)
            forward_seq=MolSet[i].seq[start_for:end_for]
            start_rev=MolSet[i].length-end_for
            end_rev=start_rev+int(is_read)
            reverse_seq=All_reverse[start_rev:end_rev]
            read1=forward_seq[23:Par.SR]
            read2=reverse_seq[0:Par.SR]
            read1seq=''
            read2seq=''
            read1qual=''
            read2qual=''
            if Par.Fast_mode=='N':
               Sim_seq1=SIMSeqQual(SeqQual_dict,read1,Par)
               Sim_seq2=SIMSeqQual(SeqQual_dict,read2,Par)
               Sim_barcode=SIMBarcodeQual(BarcodeQual_dict,MolSet[i].barcode,Par)
               read1seq=MolSet[i].barcode+'NNNNNNN'+Sim_seq1.seq
               read2seq=Sim_seq2.seq
               read1qual=Sim_barcode.qual+'KKKKKKK'+Sim_seq1.qual
               read2qual=Sim_seq2.qual
            else:
               read1seq=MolSet[i].barcode+'NNNNNNN'+read1
               read2seq=read2
               read1qual='K'*Par.SR
               read2qual='K'*Par.SR
            readname='@ST-K00126:'+str(i+1)+':H5W53BBXX:'+str(MolSet[i].start)+':'+str(MolSet[i].end)+':'+str(start_for)+':'+str(end_for)
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
    os.system('gzip -f read-RA_si-CCTGGAGA_lane-001-chunk-001.fastq')
    os.system('gzip -f read-I1_si-CCTGGAGA_lane-001-chunk-001.fastq')
    return None

def main():
    os.system('rm -rf '+sys.argv[1]+'/lib*')
    list=os.listdir(sys.argv[1])
    Par=parameter()
    hap=1
    for i in range(len(list)):
        print('processing library '+str(i+1))
        os.system('mkdir '+sys.argv[1]+'/lib'+str(i+1))
        deter=input_parameter(sys.argv[1]+'/'+list[i],Par)
        if deter==1:
           if Par.hap==1:
              haploid(Par,deter,i+1)
           elif Par.hap==2:
              diploid(Par,deter,i+1)
    return None
if __name__=="__main__":
    main()
