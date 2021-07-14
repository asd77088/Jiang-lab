dict={}
for eachfile in open('path_to/example/AS.list','r'):
	eachfile = eachfile.strip()
	cancer = eachfile.split('/')[-2]
	for eachlien in open(eachfile,'r'):
		eachlist = eachlien.strip('\n').split('\t')
		name,det,wt,mt,cut= eachlist
		type= name.split('_')[1]
#		if type not in ['RI','ME'] :
#			continue
		for eachcut in cut.split(',') :	
			dict.setdefault(eachcut,[])
			dict[eachcut].append(cancer+'_'+name)
out1=open('Test.test.fa','w')
out2=open('Test.test.anno','w')
for n,key in enumerate(dict) :
	value = dict[key]
	out1.write('>'+str(n+468891)+'\n'+key+'\n')
	out2.write(str(n+468891)+'\t'+','.join(value)+'\n')
