library("pheatmap")
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(stats)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(tidyr)
library(DT)
library(ggrepel)
library(data.table)

signf_genes<-function(res,c,pthresh){
  keep<-c("Gene","sgRNA")
  keep<-c(keep, colnames(res)[grepl(c,colnames(res))])
  #
  res<-res[,keep]
  colnames(res)[-1:-2]<-c('beta','padj')
  r <-
    res %>%
    na.omit() %>%
   # mutate( padj = stats::p.adjust(padj, method = p.adjust.methods)) %>%
    filter( padj < pthresh ) %>%
    mutate( log10pval = log10(padj)*(-1) ) %>%
    arrange(Gene)
  return(r)
}	

plot_volcano_cutom<-function(res,c){
  #
  keep<-c("Gene","sgRNA")
  keep<-c(keep, colnames(res)[grepl(c,colnames(res))])
  #
  res<-res[,keep]
  colnames(res)[-1:-2]<-c('beta','padj')
  res$log10pval<--log10(res$padj)
  #
  cats<-get_comparisons(c)
  cat.up<-paste0("pos. selected for ",cats[1])
  cat.down<-paste0("lost for ",cats[1])
  res <- res %>%
    mutate(gene_type = case_when(beta > 0 & padj <= P.THRESH~ cat.up,
                                 beta < 0 & padj <= P.THRESH ~ cat.down,
                                 TRUE ~ "ns"))
  #
  cols <- c("#ffad73", "#26b3ff", "grey") 
  names(cols)<-c(cat.up,cat.down,'ns')
  sizes <- c(5,5,1) 
  names(sizes)<-c(cat.up,cat.down,'ns')
  alphas <- c(1,1,0.5)
  names(alphas)<-c(cat.up,cat.down,'ns')
  #
  sig_il_genes <- res %>%
    filter(gene_type != "ns") %>%
    dplyr::select(Gene,log10pval,beta) %>%
    mutate(cart_dist=log10pval^2+beta^2) %>%
    arrange(desc(cart_dist))
  N<-min(45,nrow(sig_il_genes))
  sig_il_gene=sig_il_genes[1:N,'Gene'] # keep top 5
  sig_il_gene<-levels(droplevels(sig_il_gene))
  sig_il_genes <- res %>%
    filter(Gene %in% sig_il_gene)
  #
  volcPlot <- res %>%
    ggplot(aes(
      x = beta,
      y = log10pval,
      fill = gene_type,
      size = gene_type,
      alpha = gene_type
    )) +
    geom_point(shape = 21,
               # Specify shape and colour as fixed local parameters
               colour = "black") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") +
    geom_vline(xintercept = 0,
               linetype = "dashed") +
    geom_vline(xintercept = 1,
               linetype = "dashed") +
    geom_vline(xintercept = -1,
               linetype = "dashed") +
    geom_label_repel(
      data = sig_il_genes,
      # Add labels last to appear as the top layer
      aes(label = Gene),
      color = "white", segment.color = "#292929",#"grey50",
      force = 2,
      nudge_y = 1,
      show.legend = FALSE,  max.overlaps = Inf
    )  +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes, guide ="none") + # Modify point size
    scale_alpha_manual(values = alphas, guide ="none") + # Modify point transparency
    theme_bw() +   theme(
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    labs(title = paste0("Diff. gene KO - ", c), fill = "Genes")
  return(volcPlot)
}

trim_labels<-function(df,labels){
  dfs<-df[df$group==labels,]
  N<-nrow(dfs)
  Nmax=20
  if(N>Nmax){
    dfs<-dfs[order(dfs$padj),]
    labels.new<-levels(droplevels(dfs[1:Nmax,'Gene']))
  }else{
    labels.new<-levels(droplevels(dfs[,'Gene']))
  }
  return(labels.new)
}

get_comparisons<-function(c){
  cats<-str_extract(c,'[A-Za-z0-9]+_vs_(.*)|[A-Za-z0-9]+_x_(.*)')
  cats<-str_split(cats,'_vs_|_x_')[[1]]
  return(cats)
}

rankplot_custom<-function(data,c){
  keep<-c("Gene","sgRNA")
  keep<-c(keep, colnames(data)[grepl(c,colnames(data))])
  #
  data.sub<-data[,keep]
  colnames(data.sub)[-1:-2]<-c('beta','padj')
  #
  data.sub$pvalOrd<-sign(data.sub$beta)*log10(data.sub$padj)*(-1)
  data.sub$Rank = rank(data.sub$pvalOrd)#data.sub$beta
  cats<-get_comparisons(c)
  # cat.up<-paste0("for ",cats[1])
  # cat.down<-paste0("for ",cats[2])
  cat.up<-paste0("pos. selected for ",cats[1])
  cat.down<-paste0("lost for ",cats[1])
  #
  data.sub <- data.sub %>%
    mutate(group = case_when(beta > 0 & padj <= P.THRESH ~ cat.up,
                             beta < 0 & padj <= P.THRESH ~ cat.down,
                             TRUE ~ "ns"))
  #
  idx =data.sub %>% filter(group != 'ns') %>% select_(.dots='Gene')
  idx=levels(droplevels(idx$Gene))
  mycolour = c("#ffad73", "#26b3ff", "grey")
  names(mycolour) <-c(cat.up,cat.down,'ns')
  # trim length of labels
  idx.down<-trim_labels(data.sub,cat.down)
  idx.up<-trim_labels(data.sub,cat.up)
  idx=c(idx.down,idx.up)
  #
  p = ggplot(data.sub) + 
    geom_jitter(aes_string(x = "pvalOrd", y = "Rank", color = "group"), size = 2,show.legend = FALSE)+
    geom_vline(xintercept = 0,
               linetype = "dashed")
  if (length(idx) > 0) {
    p=p + geom_label_repel(
      aes_string(x = "pvalOrd", y = "Rank", fill = "group", label = "Gene"), 
      data = data.sub[data.sub$Gene %in% idx,], 
      fontface = "bold", color = "white", size = 5, 
      force = 2, segment.color = "#292929",
      nudge_y = 1,
      show.legend = FALSE,
      # box.padding = unit(0.4, "lines"), segment.color = "grey50", 
      # point.padding = unit(0.3, "lines"), segment.size = 0.3, 
      max.overlaps = Inf)+
      geom_vline(xintercept = log10(P.THRESH),
                 linetype = "dashed") +
      geom_vline(xintercept = (-1)*log10(P.THRESH),
                 linetype = "dashed") 
  }
  p= p+scale_color_manual(values = mycolour)+ 
    scale_fill_manual(values = mycolour) + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = "black")) +
    theme(text = element_text(colour = "black", size = 14, 
                              family = "Helvetica"), plot.title = element_text(hjust = 0.5, 
                                                                               size = 18), axis.text = element_text(colour = "gray10"))+
    theme(axis.line = element_line(size = 0.5, colour = "black"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(), panel.background = element_blank())+
    labs(x = "abs.log10 p-value (adj)", y = "Rank", 
         title = c, fill = "Genes")# + 
  #theme(legend.position = "none")#left")
  p
  return(p)
}

group_sgRNAs<-function(M){
  m<-M %>% as_tibble() %>%
    group_by(Gene) %>% 
    summarise(across(everything(), mean))
  m <- m %>%
    select(-sgRNA) %>%
    as.data.frame()
  return(m)
}

group_replicates<-function(M){
  m1<-t(as.data.frame(lapply(toupper(colnames(M)[c(-1,-2)]), function(x) str_split(x,'_\\d',n=2))))
  m1<-as.data.frame(m1)
  m1$sample_names<-colnames(M)[c(-1,-2)]
  new.samples<-levels(m1$V1)
  M.bis<-data.frame(M[,c(1,2)])
  for(N in new.samples){
    w<-m1[m1$V1==N,'sample_names']
    W<- c("sgRNA","Gene",w)
    M.sub<-M[,colnames(M) %in% W]
    if(length(w)>1){
      M.sub[,N]<-rowMeans(M[,colnames(M) %in% w])
    }else{
      M.sub[,N]<-M[,w]
    }
    M.bis<-merge(M.bis,M.sub[,c('sgRNA','Gene',N)],by=c('sgRNA','Gene'))
  }
  #
  return(M.bis)
}

make_heatmaps<-function(M,myList,c,pool){
  names<-levels(meta.data$Sample_name)
  names<-names[grep(paste0('_',pool),sort(names))]
  #
  names<-names[!grepl('plasmid',names)]
  M<-M[,!(grepl('plasmid',colnames(M)))]
  #
  myList[,'padj']<-stats::p.adjust(myList[,'padj'], method = p.adjust.methods)
  myList<-myList[order(myList$padj),'Gene']
  nMax<-min(nrow(myList),30)
  myList<-levels(droplevels(myList[1:nMax]))
  tmp.heatmap<-M[M$Gene %in% myList,]
  #rownames(tmp.heatmap)<-tmp.heatmap$sgRNA
  #tmp.heatmap<-group_replicates(tmp.heatmap)
  tmp.heatmap<-group_sgRNAs(tmp.heatmap)
  rownames(tmp.heatmap)<-tmp.heatmap$Gene
  #
  tmp.heatmap<-as.matrix(tmp.heatmap[,colnames(tmp.heatmap)[c(-1)]])
  #tmp.heatmap[,colnames(tmp.heatmap)[c(-1,-2)]])
  # tmp.heatmap<-scale(tmp.heatmap)
  # column order by cell line
  #tmp.heatmap<-tmp.heatmap[,names]
  #
  set.seed(1)
  htmp<-pheatmap(tmp.heatmap, 
                 cluster_rows=TRUE, 
                 cluster_cols=TRUE, #FALSE, 
                 main=paste0(c,', top genes'), scale="none",#"column",#'row',
                 clustering_distance_rows = "correlation")
  htmp
}

boxplot.genes<-function(M,Important.genes){
  M.sub<-M[M$Gene %in% Important.genes,]
  Nvars<-length(colnames(M))
  M.sub.box<-reshape2::melt(M.sub[,colnames(M)[c(-1,-Nvars)]])
  M.sub.box$sample<-unlist(lapply(M.sub.box$variable, function(c) str_split(c,"_\\d",n=2)[[1]][1]))
  p.box<-ggplot(M.sub.box, aes(x=sample, y=value, fill=Gene)) + 
    geom_boxplot()+ facet_grid(. ~ Gene)+theme_minimal() +
    labs(title=" ",x="", y = "Counts") +  theme(
      axis.text.x = element_text(
        angle = 45))+
    scale_y_continuous(trans = "log10")
  return(p.box)
}