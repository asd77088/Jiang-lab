from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from functools import reduce
import os
from multiprocessing.dummy import Pool as MyPool
from time import sleep,ctime

ID2Name_dict={eachline.split(',')[0].strip('"') :  eachline.split(',')[2].strip('"') for eachline in open('ncbi2symbol.csv','r') }


hg19_dict={}
for seq_record in SeqIO.parse("path_to/Hg19.fa", "fasta"):
        hg19_dict[seq_record.id.strip('chr')]=seq_record

def GetSeq(hg19_dict,chr,start) :
	DNA= hg19_dict[chr][int(start)-1]
#	DNA = ''.join([eachnt for eachnt in hg19_dict[chr][int(start)-1:int(end)]])
	return DNA

def Translate(DNA,stand) :
	DNA= Seq(DNA.upper())
	if stand =='-' :
		DNA= DNA.reverse_complement()
	CDS= DNA[DNA.find('ATG'):]
	if len(CDS)== 0 :
		return 'No seq'
	elif len(CDS) == 1:
		return 'No Promoter'
	else :
		AAseq = CDS.translate(to_stop=True)
		if len(AAseq) *3 == len(CDS):
				return 'No end'
		else :
				return ''.join([each for each in AAseq])


def Cal_pos(list1,list2,fuhao='-') :
	if fuhao == '+' :
		for eachlist2 in list2 :
			if int(eachlist2[0]) > int(list1[0][0])  and int(eachlist2[1]) < int(list1[-1][1]) : 
				list1= list1+[eachlist2]
		return list1

	elif fuhao == '-':
		list1 = reduce(lambda x,y:x+y,[list(range(min(int(each[0]),int(each[1])),max(int(each[0]),int(each[1]))+1)) for each in list1])
		list2 = reduce(lambda x,y:x+y,[list(range(min(int(each[0]),int(each[1])),max(int(each[0]),int(each[1]))+1)) for each in list2])
		return sorted(list(set([each for each in list1 if each not in list2])))
	else :
		return sorted(list(set(reduce(lambda x,y:x+y,[list(range(int(each[0]),int(each[1])+1)) for each in list1]))))

def qie(seq,n=9) :
	if len(seq) >= n :
		return set([seq[each:each+n] for each in range(len(seq)-n+1)])
	else :
		return set([])


def cutseq(seq1,seq2,type,n) :
	if seq2 =='NotFound' :
		return []
	if (n and type in ['ES','AA','AD','ME'])  :
		return list(qie(seq1)  - qie(seq2))
	elif (not n and type in ['ES','AA','AD','ME']) :
		return list(qie(seq2) - qie(seq1))
	elif  (n and type in ['RI','AT','AP']) :
		return list(qie(seq2) - qie(seq1))
	else :
		return list(qie(seq1)  - qie(seq2))


def Merge(Gene_pos,Max_Gene_ditc,gene,exon1,type,Tran_pos_dict,Gene_dict) :
	Rid = Max_Gene_ditc[gene]['Rid']
	if type in ['ES','AD','AA'] :
		drop_pos=[Gene_pos[gene][eachexon] for eachexon in exon1.split(':')]
		ref_pos = [Max_Gene_ditc[gene][eachexon][1:] for eachexon in Max_Gene_ditc[gene]['CDS'][1:]]
		stand = Max_Gene_ditc[gene]['CDS'][0]
		chr = Max_Gene_ditc[gene][Max_Gene_ditc[gene]['CDS'][1]][0]
		Ref_dna = ''.join([GetSeq(hg19_dict,chr,eachpos) for eachpos in Cal_pos(ref_pos,[[1,2]],'xx')])
		Mutation_dna = ''.join([GetSeq(hg19_dict,chr,eachpos)  for eachpos in Cal_pos(ref_pos,drop_pos,'-')])
		Ref_aa = Translate(Ref_dna,stand)	
		Mutation_aa = Translate(Mutation_dna,stand)
		if stand == '+' :
			if int(drop_pos[-1][-1] ) <= int(ref_pos[0][0]) :
				A='A5<'
			elif int(drop_pos[0][0]) >= int(ref_pos[-1][-1]) :
				A='A3<'
			else :
				A='B<'
		if stand == '-' :
			if int(drop_pos[0][0]) >= int(ref_pos[-1][-1]) :
				A='A5<'
			elif int(drop_pos[-1][-1]) <= int(ref_pos[0][0]) :
				A='A3<'
			else :
				A='B<'
		return Rid,Rid,Ref_dna,Mutation_dna,Ref_aa,Mutation_aa,A
	elif type =='RI' :
		add_pos = [Gene_pos[gene][eachexon] for eachexon in exon1.split(':')]
		ref_pos = [Max_Gene_ditc[gene][eachexon][1:] for eachexon in Max_Gene_ditc[gene]['CDS'][1:]]
		stand = Max_Gene_ditc[gene]['CDS'][0]
		chr = Max_Gene_ditc[gene][Max_Gene_ditc[gene]['CDS'][1]][0]
		Ref_dna = ''.join([GetSeq(hg19_dict,chr,eachpos) for eachpos in Cal_pos(ref_pos,[[1,2]],'xx')])
		
		Mutation_dna = ''.join([GetSeq(hg19_dict,chr,eachpos)  for eachpos in Cal_pos(Cal_pos(ref_pos,add_pos,'+'),[],'x')])
		Ref_aa = Translate(Ref_dna,stand)
		Mutation_aa = Translate(Mutation_dna,stand)
		if stand == '+' :
			if int(add_pos[-1][-1] ) <= int(ref_pos[0][0]) :
				A='A5>'
			elif int(add_pos[0][0]) >= int(ref_pos[-1][-1]) :
				A='A3>'
			else :
				A='B>'
		if stand == '-' :
			if int(add_pos[0][0]) >= int(ref_pos[-1][-1]) :
				A='A5>'
			elif int(add_pos[-1][-1]) <= int(ref_pos[0][0]) :
				A='A3>'
			else :
				A='B>'
		return Rid,Rid,Ref_dna,Mutation_dna,Ref_aa,Mutation_aa,A
	elif type =='ME':
		#alt_pos = [Gene_pos[gene][eachexon] for eachexon in exon1.split('|')]
		alt_pos1 = [Gene_pos[gene][each11] for each11 in exon1.split('|')[0].split(':')]
		alt_pos2 =  [Gene_pos[gene][each11] for each11 in exon1.split('|')[1].split(':')]
		ref_pos = [Max_Gene_ditc[gene][eachexon][1:] for eachexon in Max_Gene_ditc[gene]['CDS'][1:]]
		stand = Max_Gene_ditc[gene]['CDS'][0]
		chr = Max_Gene_ditc[gene][Max_Gene_ditc[gene]['CDS'][1]][0]
		All_dna = Cal_pos(ref_pos,alt_pos1+alt_pos2,'+')
		Mutation1_dna = ''.join([GetSeq(hg19_dict,chr,eachpos)  for eachpos in  Cal_pos(All_dna,alt_pos1,'-')])
		Mutation2_dna = ''.join([GetSeq(hg19_dict,chr,eachpos)  for eachpos in Cal_pos(All_dna,alt_pos2,'-')])
		Mutation1_aa = Translate(Mutation1_dna,stand)
		Mutation2_aa = Translate(Mutation2_dna,stand)
		if stand == '+' :
			if int(alt_pos1[-1][-1] ) <= int(ref_pos[0][0])  and int(alt_pos2[-1][-1] ) <= int(ref_pos[0][0]) :
				A='A5'
			elif int(alt_pos1[0][0]) >= int(ref_pos[-1][-1]) and  int(alt1_pos2[0][0]) >= int(ref_pos[-1][-1]) :
				A='A3'
			else :
				A='B'
		if stand == '-':
			if int(alt_pos1[0][0]) >= int(ref_pos[-1][-1]) and int(alt_pos2[0][0]) >= int(ref_pos[-1][-1]) :
				A='A5'
			elif int(alt_pos1[-1][-1]) <= int(ref_pos[0][0]) and int(alt_pos2[-1][-1]) <= int(ref_pos[0][0]):
				A='A3'
			else :
				A='B'
		if len(Mutation2_dna) > len(Mutation1_dna) :
			A+='>'
		if len(Mutation2_dna) == len(Mutation1_dna) :
			A+='='
		if len(Mutation2_dna) < len(Mutation1_dna) :
			A+='<'
		return Rid,Rid,Mutation2_dna,Mutation1_dna,Mutation2_aa,Mutation1_aa,A

	elif type in ['AP','AT'] :
		stand = Max_Gene_ditc[gene]['CDS'][0]
		chr = Max_Gene_ditc[gene][Max_Gene_ditc[gene]['CDS'][1]][0]
		if '.' in exon1 :
			exon1list =[each  for each in  Gene_pos[gene] if each.startswith(exon1.split('.')[0]+'.')]
		else :
			exon1list=[exon1]
		for eachexon in exon1list :
			
			if stand =='+' and type == 'AP' :
				thii = Gene_pos[gene][eachexon][1]
				tran_list = [eachtran for eachtran in Tran_pos_dict[gene] if thii in Tran_pos_dict[gene][eachtran]['end']]
				tran_list =[ eachtran for eachtran in tran_list if (Tran_pos_dict[gene][eachtran]['end'].index(thii) == 0)]
			elif stand =='+' and type == 'AT' :
				thii = Gene_pos[gene][eachexon][0]
				tran_list = [eachtran for eachtran in Tran_pos_dict[gene] if thii in Tran_pos_dict[gene][eachtran]['start']]
				tran_list =[ eachtran for eachtran in tran_list if Tran_pos_dict[gene][eachtran]['start'][-1]==thii]
			elif stand =='-' and type == 'AP' :
				thii = Gene_pos[gene][eachexon][1]
				tran_list = [eachtran for eachtran in Tran_pos_dict[gene] if thii in Tran_pos_dict[gene][eachtran]['start']]
				tran_list =[ eachtran for eachtran in tran_list if Tran_pos_dict[gene][eachtran]['start'][-1]==thii]
			else :
				thii = Gene_pos[gene][eachexon][0]
				tran_list = [eachtran for eachtran in Tran_pos_dict[gene] if thii in Tran_pos_dict[gene][eachtran]['end']]
				tran_list =[ eachtran for eachtran in tran_list if (Tran_pos_dict[gene][eachtran]['end'].index(thii) == 0)]
			if tran_list ==[] :
				continue
			else :
				break 
		leng = 0
		if tran_list ==[] :
			ref_pos=[ Max_Gene_ditc[gene][eachexon][1:] for eachexon in Max_Gene_ditc[gene]['CDS'][1:]]
			Ref_dna = ''.join([GetSeq(hg19_dict,chr,eachpos) for eachpos in Cal_pos(ref_pos,[[1,2]],'xx')])
			Ref_aa = Translate(Ref_dna,stand)
			return '','',Ref_dna,'NotFound',Ref_aa,'NotFound',''
		else :
			for eachran  in tran_list :
				AllCDS = Gene_dict[gene][eachtran]['CDS'][1:]
				LenCDS =sum( [int(each2.split('_')[3]) - int(each2.split('_')[2]) for each2 in AllCDS if each2.split('_')[1]=='CDS'])
				if LenCDS > leng :
					leng = LenCDS 
					final_exon = eachran
			ref_pos = [Max_Gene_ditc[gene][eachexon][1:] for eachexon in Max_Gene_ditc[gene]['CDS'][1:]]	
			alt_pos = [Gene_dict[gene][final_exon][eachexon][1:] for eachexon in Gene_dict[gene][final_exon]['CDS'][1:]]
			Ref_dna = ''.join([GetSeq(hg19_dict,chr,eachpos) for eachpos in Cal_pos(ref_pos,[[1,2]],'xx')])
			alt_dna = ''.join([GetSeq(hg19_dict,chr,eachpos) for eachpos in Cal_pos(alt_pos,[[1,2]],'xx')])
			Ref_aa = Translate(Ref_dna,stand)
			alt_aa = Translate(alt_dna,stand)
			if Ref_dna == alt_dna  :
				refid = Max_Gene_ditc[gene]['Rid']
				altid = final_exon
				if type =='AP': 
					A='A5'
					if stand == '+' :
						reflen = abs(int(Gene_dict[gene][refid]['CDS'][1].replace('_codon','codon').split('_')[2])- int( Tran_pos_dict[gene][refid]['start'][0]))
						altlen = abs(int(Gene_dict[gene][altid]['CDS'][1].replace('_codon','codon').split('_')[2])- int( Tran_pos_dict[gene][altid]['start'][0]))
					if stand == '-' :
						reflen = abs(-int(Gene_dict[gene][refid]['CDS'][-1].replace('_codon','codon').split('_')[3])+ int( Tran_pos_dict[gene][refid]['end'][-1]))	
						altlen = abs(-int(Gene_dict[gene][altid]['CDS'][-1].replace('_codon','codon').split('_')[3])+ int( Tran_pos_dict[gene][altid]['end'][-1]))
					if altlen > reflen :
							A+='>'
					if altlen == reflen :
							A+='='
					if altlen < reflen :
							A+='<'
				if  type =='AT':
					A='A3'
					if stand == '+' :
						reflen =abs( -int(Gene_dict[gene][refid]['CDS'][-1].replace('_codon','codon').split('_')[3])+ int( Tran_pos_dict[gene][refid]['end'][-1]))
						altlen = abs(-int(Gene_dict[gene][altid]['CDS'][-1].replace('_codon','codon').split('_')[3])+ int( Tran_pos_dict[gene][altid]['end'][-1]))
					if stand == '-' :
						reflen = abs(int(Gene_dict[gene][refid]['CDS'][1].replace('_codon','codon').split('_')[2])- int( Tran_pos_dict[gene][refid]['start'][0]))
						altlen = abs(int(Gene_dict[gene][altid]['CDS'][1].replace('_codon','codon').split('_')[2])- int( Tran_pos_dict[gene][altid]['start'][0]))
					if altlen > reflen :
						A+='>'
					if altlen == reflen :
						A+='='
					if altlen < reflen :
						A+='<'
			if Ref_dna != alt_dna  :
				if len(Ref_dna) > len(alt_dna) :
					A='B<'
				elif len(Ref_dna) == len(alt_dna) :
					A='B='	
				else :
					A='B>'
			
			return Rid,final_exon,Ref_dna,alt_dna,Ref_aa,alt_aa,A

chr_list=[str(each) for each in range(1,23)]+['X','Y']
type_list=['CDS','stop_codon','start_codon']
Gene_dict ={}


##Get mat transcript in gene
mark =''
for eachline in open('path_to/hg19.gtf','r') :
	eachline = eachline.strip().split('\t')
	chr,a,type,start,end,b,stand,c,id = eachline
	chr = chr.strip('chr')
	if (chr not in chr_list) or (type not in type_list):
		continue
	transcript_id = id.split(';')[1].split(' ')[2].strip('"')
	if transcript_id != mark :
		mark = transcript_id 
		last_n = 0
	
	CDSID = '_'.join([chr,type,start,end,stand])
	GeneName = ID2Name_dict[transcript_id]
	Gene_dict.setdefault(GeneName,{})
	Gene_dict[GeneName].setdefault(transcript_id,{'CDS':[stand]})
	Gene_dict[GeneName][transcript_id]['CDS'].append(CDSID)
	Gene_dict[GeneName][transcript_id][CDSID]= [chr,max(int(start),last_n+1),end]
	last_n = int(end)
Max_Gene_ditc= {}
for eachgene in Gene_dict :
	exon_len= 0
	start=9999999999999999999999999999999999999999999999
	for eachtran in Gene_dict[eachgene] :
		stand = Gene_dict[eachgene][eachtran]['CDS'][0]
		AllCDS = Gene_dict[eachgene][eachtran]['CDS'][1:]
		LenCDS =sum( [int(each2.split('_')[3]) - int(each2.split('_')[2]) for each2 in AllCDS if each2.split('_')[1]=='CDS'])
		if stand == '+' :
			startexon = Gene_dict[eachgene][eachtran]['CDS'][1]
			startpos = Gene_dict[eachgene][eachtran][startexon][1]
		else :
			startexon = Gene_dict[eachgene][eachtran]['CDS'][-1]
			startpos = Gene_dict[eachgene][eachtran][startexon][2]
	
		if int(stand+ str(startpos)) < start :
			start = int(stand+ str(startpos))
			Max_Gene_ditc[eachgene]=Gene_dict[eachgene][eachtran]
			Max_Gene_ditc[eachgene]['Rid']= eachtran
			exon_len = LenCDS
		elif int(stand+ str(startpos)) == start and LenCDS > exon_len :
			exon_len = LenCDS
			Max_Gene_ditc[eachgene]=Gene_dict[eachgene][eachtran]
			Max_Gene_ditc[eachgene]['Rid']= eachtran

##Get thr pos of each tran
type_list =['exon']
Tran_pos_dict = {}
#Tran_pos_dict
for eachline in open('hg19.gtf','r') :
	eachline = eachline.strip().split('\t')
	chr,a,type,start,end,b,stand,c,id = eachline
	chr = chr.strip('chr')
	transcript_id = id.split(';')[1].split(' ')[2].strip('"')
	if (chr not in chr_list) or (type not in type_list) or (transcript_id.startswith('NR')):
		continue
	GeneName = ID2Name_dict[transcript_id]
	Tran_pos_dict.setdefault(GeneName,{})
	Tran_pos_dict[GeneName].setdefault(transcript_id,{})
	Tran_pos_dict[GeneName][transcript_id].setdefault('start',[]).append(start)
	Tran_pos_dict[GeneName][transcript_id].setdefault('end',[]).append(end)
	


##Get the exon pos
Gene_pos = {}
for eachline in open('path_to/TCGA_SpliceSeq_Gene_Structure.txt','r') :
	gene,chr,stand,Exon,start,end = eachline.strip('\n').split('\t')
	Gene_pos.setdefault(gene,{})
	if stand == '+' :
		Gene_pos[gene][Exon]=[start,end.strip('\r')]
	else :
		Gene_pos[gene][Exon]=[end.strip('\r'),start.strip('\r')]
	

#Read list
soft1=[eachlien.split()[0] for eachlien in open('HLATYPE/MHC_pseudo.dat') if eachlien.startswith('HLA')]
soft2=[eachlien.strip().strip('.thr') for eachlien in open('HLATYPE/HLA_other.list') if eachlien.startswith('HLA')]
mylist=[]
for eachline in open('HLATYPE/HLA.list1','r') :
        eachlist= eachline.strip().split('\t')[2]
        A = eachlist.split('*')[0]
        B= eachlist.split('*')[1].split(':')[0:2]
        mylist.append('HLA-'+A+':'.join(B))
mylist =list(set(mylist))
list1 = [each for each in mylist if each in soft1 ]
list2 = [each for each in mylist if each in soft2 ]

#hla_list=['HLA-A01:01','HLA-A02:01']
pwd=os.getcwd()
keephead=['symbol','as-id','splice-type','exons','from-exon','to-exon','FDR','det']
#for eachfile in open('Al.list','r'):
def Run(eachfile) : 
	name =eachfile.strip().split('_')[-1].split('.')[0]
	os.system('mkdir -p %s'%(name))
	fafile = os.path.join(pwd,name,name+'.ForASPN.csv')
	fafile2 = os.path.join(pwd,name,name+'.CDS.csv')
	fafile3 = os.path.join(pwd,name,name+'.ForDomain.fa')
#	hla_result1 = os.path.join(pwd,name,name+'.netMHCpan.myout')
#	hla_result2 = os.path.join(pwd,name,name+'.netCTLpan.myout')
#	shfile= os.path.join(pwd,name,name+'_netMHCpan.shh')
	fa=  open(fafile,'w')
	fa3 =open(fafile3,'w')
	fa2= open(fafile2,'w')
	fa2.write('\t'.join(['AS_event','gene','AS_type','exons','from-exon','to-exon','FDR','det','Unspliced_dna_seq','Spliced_dna_seq','Unspliced_AA_seq','Spliced_AA_seq','Cut_AA_seq','CDS','Cancer'])+'\n')
	ASfile = open(eachfile.strip(),'r')
	head = ASfile.readline().strip('\n').split('\t')
	head_n = [head.index(each) for each in keephead]
	for eachline in ASfile :
		eachline = eachline.strip().split('\t')
		gene,asid,type,exon1,exon2,exon3,FDR,det =[ eachline[each] for each in head_n]
		if float(FDR) > 0.05 or abs(float(det)) < 0.1 or type not in ['ES','AD','AA','RI','ME','AP','AT'] or gene not in Max_Gene_ditc:
			continue
		Rid= Max_Gene_ditc[gene]['Rid']
		tid1,tid2,dna1,dna2,seq1,seq2,A=Merge(Gene_pos,Max_Gene_ditc,gene,exon1,type,Tran_pos_dict,Gene_dict)
#		same='Nosame'
#		if tid2 == Rid :
#			same='same'
#		if len(seq2) < 30 :
#			continue
		lastseq = cutseq(seq1,seq2,type,float(det) > 0)
		CDS1='Different'
		if seq1==seq2 :
			CDS1='Same'
		if seq2 == 'NotFound':
			CDS1='NotFound'
		name1 = '_'.join([gene,type,exon1,exon2,exon3])
		
		fa2.write('\t'.join([gene+'_'+asid,gene,type,exon1,exon2,exon3,FDR,det,dna1,dna2,seq1,seq2,','.join(lastseq),CDS1,name])+'\n')
		if CDS1== 'Different' :
			fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'WT'])+'\n'+seq1 +'\n')
			fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'MT'])+'\n'+seq2 +'\n')
		if len(seq2) >  30 :
			fa.write('\t'.join([name1,str(det),seq1,seq2,','.join(lastseq)])+'\n')
#		if len(lastseq) >= 9 :
#				fa.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3])+'\n'+lastseq +'\n')
#		fa.write('\t'.join([name,str(det),seq1,seq2,','.join(lastseq)])+'\n')
#		if seq1 != seq2 and :
#			fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'WT'])+'\n'+seq1 +'\n')
#			fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'MT'])+'\n'+seq2 +'\n')
#		
#		if len(seq2) < 30:
#			continue
#		fa2.write('\t'.join([gene+'_'+asid,gene,type,exon1,exon2,exon3,FDR,det,dna1,dna2,seq1,seq2,','.join(lastseq),CDS1,name])+'\n')
#		fa.write('\t'.join([name,str(det),seq1,seq2,','.join(lastseq)])+'\n')
#		fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'WT'])+'\n'+seq1 +'\n')
#		fa3.write('>'+'_'.join([gene,Rid,type,exon1,exon2,exon3,'MT'])+'\n'+seq2 +'\n')
	fa2.close()
	fa3.close()
#	open(shfile,'w').write(sh_code)
#	os.system('sh %s'%(shfile))





Par_list = [eachfile.strip() for eachfile in open('Al.list','r') if not eachfile.startswith('#')]

pool=MyPool(8)
pool.map(Run,Par_list)
pool.close()
pool.join()	
