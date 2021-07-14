from multiprocessing.dummy import Pool as MyPool
from time import sleep,ctime
import os,re

def RunMHCpan(hla) :
	path1=os.getcwd()
	pathsh= os.path.join(path1,'sh_test')
	pathout= os.path.join(path1,'out_test')
	if not os.path.exists(pathsh) :
		os.system('mkdir -p %s %s'%(pathsh,pathout))
	MHC_code='''path_to/netMHCpan-4.1/netMHCpan \\
	Test.test.fa \\
	 -a %s\\
	-s \\
	-BA \\
	-l 9 > %s'''%(hla,os.path.join(pathout,hla+'.netMHCpan.out'))
	shfile=os.path.join(pathsh,hla+'.netMHCpan.sh')
	open(shfile,'w').write(MHC_code)
	os.system('sh %s'%(shfile))


Par_list = open('HLA.list','r').readline().strip().split(',')


pool=MyPool(48)
pool.map(RunMHCpan,Par_list)
pool.close()
pool.join()
