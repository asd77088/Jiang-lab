library(ggsci)

a<-read.csv("path_to/example/Barplot.csv",header = T)

p=ggplot(data = a, mapping = aes(x =  reorder(class,Freq), y = Freq, fill = Var1)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45, hjust=1)) 

p1=p+scale_fill_npg()+theme_bw()+
  ylab("")+xlab("")+ 
  theme_classic()+
  coord_polar(theta = 'x',start=0)+theme_bw() +
  theme(axis.text.x  = element_text(size=20,color="black",face="bold"),axis.title=element_text(size=15,face="bold"),
        axis.text.y  = element_text(size=10,color="black",face="bold"))+
  theme(legend.title = element_text(size=15,color="black",face = "bold"))+ labs(fill="Type")+ theme(legend.text=element_text(size=15,color="black",face="bold"))

p1