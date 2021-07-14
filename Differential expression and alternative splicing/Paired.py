#####Screen paired samples
import os
for eachfile in os.listdir() :
	if not eachfile.startswith('path_to/overlap') :
		continue
	out=open('Paired.'+eachfile,'w')
	print(eachfile)
	sample_dict ={}
	for eachline in open(eachfile,'r'):
		if eachline.startswith('"",') :
			keepn=[]
			headlist =[each.strip('"') for each in eachline.strip().split('",')]
			basen = [n for n,each in enumerate(headlist) if not each.startswith('TCGA')]
			for n,each in enumerate(headlist) :
				if each.startswith('TCGA'):
					samplename = each.strip('-Norm')
					sample_dict.setdefault(samplename,[0,0])
					if each.endswith('-Norm'):
						sample_dict[samplename][0] = n
					else :
						sample_dict[samplename][1] = n
			for eachsample in sample_dict :
				if not  0 in sample_dict[eachsample] :
					keepn.extend(sample_dict[eachsample])
			out.write('\t'.join([headlist[each] for each in  basen + keepn])+'\n')
			continue
		eachlist = [each.strip('"') for each in eachline.strip().split(',')]
		out.write('\t'.join([eachlist[each] for each in  basen + keepn])+'\n')
	out.close()
		