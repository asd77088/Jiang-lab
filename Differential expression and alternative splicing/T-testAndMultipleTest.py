######A t-test was performed to identify different expressions with FPKM value and different alternative splicing with PSI value. P-values were corrected for multiple testing using the Benjamini-Hochberg method. 
from scipy import stats
import math
from scipy.stats import rankdata
import numpy as np
def fdr(p_vals):
	p_vals = np.array(p_vals)
	ranked_p_values = rankdata(p_vals)
	fdr = p_vals * len(p_vals) / ranked_p_values
	fdr[fdr > 1] = 1
	return fdr



for eachfile in open('path_to/example/Input.list','r') :
	if eachfile.startswith('#') :
		continue
	eachfile = eachfile.strip()
	out = open('Cal.'+eachfile,'w')
	type = eachfile.strip('.csv').split('_')[-1]
	if type == 'express' :
		add= 'Log2FC'
	else :
		add='det'
	base = []
	P_list=[]
	add_list =[]
	for eachline in open(eachfile,'r') :
		if eachline.startswith('\t') :
			head = eachline.strip('\n').split('\t')
			head.extend(['P','FDR','Bonferroni',add])
			N_list =[ i  for i,j  in enumerate(head) if j.startswith('TCGA') and j.endswith('Norm')]
			T_list = [i for i,j  in enumerate(head) if j.startswith('TCGA') and not j.endswith('Norm')]
			print(eachfile,N_list,T_list)
			out.write('\t'.join(head)+'\n')
			if len(N_list) == 0 or  len(T_list) == 0:
				break
			continue
		eachlist = eachline.strip().split('\t')
		N_value = [float(eachlist[int(eachn)])  for eachn in N_list]
		T_value = [float(eachlist[int(eachn)])  for eachn in T_list]
		if type == 'express'  and sum(N_value+T_value) < 1 :
			continue
		base.append(eachlist)
		P= stats.ttest_ind(N_value,T_value)[1]
		P_list.append(P)
		if type == 'express' :
			if sum(N_value) ==0 or sum(T_value) == 0 :
				add='NA'
			else :
				add = math.log((len(N_value)*sum(T_value))/(sum(N_value)*len(T_value)),2)
		else :
			add = sum(T_value)/len(T_value)  - sum(N_value)/len(N_value)
		add_list.append(add)
	Bonferroni =  [each*len(P_list) if each*len(P_list) < 1 else 1 for each in P_list ]
	FDR = list(fdr(P_list))
	out.write('\n'.join(['\t'.join(['"'+str(each)+'"' for each in a+[b,c,d,e]]) for a,b,c,d,e in zip(base,P_list,FDR,Bonferroni,add_list)]))
	out.close()


