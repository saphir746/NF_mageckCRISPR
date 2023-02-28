#library(MAGeCKFlute)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyverse)
library(ggpubr)


# map_meta_data_2<-function(C){
#   C<-C[,colnames(C)[-2]]
#   Cc=merge(C,meta.data,by.x='Label',by.y='Sample_name')
#   #
#   Cc$telomere_status[Cc$telomere_status=="ALT _+"]="ALT_+"
#   Cc$telomere_status <- droplevels(Cc$telomere_status)
#   #
#   Cc<-Cc[order(Cc$Label),]
#   return(Cc)
# }
# counts<-lapply(counts, function(Y) map_meta_data_2(Y))
######################################

# Gini index
Gini_indx<- function(C,What){
 # poool=str_replace(poool,'_', ' ')
  df=C[,c(What,"GiniIndex","POOL")]
  df$x=df[,What]
  df$y=df$GiniIndex
  df$p=df$POOL
  ggplot(data=df, aes(x = x, y = y, fill=p))+
    geom_bar(stat="identity", position=position_dodge())+
    labs(x=NULL, y="Gini index", title=paste0("Evenness of sgRNA reads, by ",What))+
    scale_y_continuous(expand = c(0,0))+
    labs(fill = NULL)+
    theme_bw(base_size = 12)+
    theme(plot.title = element_text(hjust = 0.5, size=16))+
    theme(axis.text.x = element_text(angle = 90,size = 8,hjust = 1))
}
# WHAT<-"cell_line"#"replicate" #"Label"# "Cas9_efficiencysingle_clones"# "telomere_status"#  "Time_point"#, 
# gini.plots<-lapply(names(coutnsummaries), function(Y) Gini_indx(coutnsummaries[[Y]],Y,WHAT) )
# ggarrange(plotlist = gini.plots, ncol=1, nrow=2)


# Missed sgRNAs
miss_count<- function(C,What){
  C$Missed = C$Zerocounts/C$TotalsgRNAs*100
  df=C[,c(What,"Missed","POOL")]
  df$x=df[,What]
  df$y=df$Missed
  df$p=df$POOL
  ggplot(data=df, aes(x = x, y = y, fill=p))+
    geom_bar(stat="identity", position=position_dodge())+
    labs(x=NULL, y="Missed sgRNA, %", title=paste0("Missed sgRNAs, by ",What))+
    scale_y_continuous(expand = c(0,0))+
    labs(fill = NULL)+
    theme(text = element_text(colour="black",size = 14),
          plot.title = element_text(hjust = 0.5, size=18),
          axis.text.x = element_text(angle = 45, hjust=1, vjust = 1,
                                     colour="gray10", face="plain"),
          axis.text.y= element_text(colour="gray10", face="plain"))+
    theme(axis.line = element_line(size=0.5, colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank())+
    theme(plot.title = element_text(hjust = 0.5, size=16))+
    theme(axis.text.x = element_text(angle = 90,size = 8,hjust = 1))
}
# WHAT<-"replicate" #"Label"# "Cas9_efficiencysingle_clones"# "telomere_status"# "cell_line"# "Time_point"#, 
# miss.plots<-lapply(names(coutnsummaries), function(Y) miss_count(coutnsummaries[[Y]],Y,WHAT) )
# ggarrange(plotlist = miss.plots, ncol=1, nrow=2)

# Read mapping
map_rate<-function(C,poool,What){
  # default <-> What="Label"
  poool=str_replace(poool,'_', ' ')
  gg = data.frame(Label=rep(C[, What], 2),
                  read=rep(C[,"Reads"],2),
                  count=c(C[, "Mapped"],
                          C[, "Unmapped"]),
                  category=factor(rep(c("Mapped", "Unmapped"),
                                      c(nrow(C), nrow(C))),
                                  levels = c("Unmapped", "Mapped")))
  gg2<-gg %>% group_by(Label,category) %>% summarise(Reads=mean(read), Counts=mean(count)) %>% as.data.frame()
  if(What=='Label'){
    gg2$percent = NA
  }else{
    gg2$percent = paste0(round(gg2$Counts*100/gg2$Reads, 1), "%")
  }
  gg2$pos = ceiling(gg2$Counts/2)
  gg2$pos[gg2$category=="Unmapped"] = ceiling(gg2$pos[gg2$category=="Unmapped"] + gg2$pos[gg2$category=="Mapped"]*2)
  fill = c("#9BC7E9", "#1C6DAB")
  ggplot(gg2)+
    geom_bar(aes_string(y = "Counts", x = "Label", fill = "category"),
             stat="identity", width=0.8, alpha=0.9)+
    geom_text(aes_string(x = "Label", y = "pos", label = "percent"), size=4)+
    labs(x=NULL, y="Reads", title="Mapping ratio")+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=fill)+
    theme(legend.title = element_blank())+
    theme(text = element_text(colour="black",size = 14),
          plot.title = element_text(hjust = 0.5, size=18),
          axis.text.x = element_text(angle = 45, hjust=1, vjust = 1,
                                     colour="gray10", face="plain"),
          axis.text.y= element_text(colour="gray10", face="plain"))+
    theme(axis.line = element_line(size=0.5, colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank())+
    ggtitle(paste0('Mapping ratios, ',poool,', by ',What))+
    theme(axis.text.x = element_text(angle = 90,size = 8,hjust = 1))
}
# WHAT<-"Label"#"replicate" # "Cas9_efficiencysingle_clones"# "telomere_status"# "cell_line"# "Time_point"#, 
# maprate.plots<-lapply(names(coutnsummaries), function(Y) map_rate(coutnsummaries[[Y]],Y,WHAT) )
# ggarrange(plotlist = maprate.plots, ncol=1, nrow=2)