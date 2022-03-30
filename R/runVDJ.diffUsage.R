#' The description of this runVDJ.diffUsage
#' @title runVDJ.diffUsage
#' 
#' @description The runVDJ.diffUsage function will calculate the group specific IGHV Usage.
#' 
#' @param object The Seurat S4 object
#' 
#' @param v.call Default "germline_v_call"; The result of change-o has been integrated in Seurat S4 object. The V region information have been save in "germine_v_call" by default
#' 
#' @param sample.label The sample labels in Seurat S4 object.
#'
#' @param group.label The group labels in Seurat S4 object.
#'
#' @param GroupA Vector of samples names belonging to group A
#'
#' @param GroupB Vector of samples names belonging to group B
#'
#' @param reference Vector of samples names belonging to reference
#'
#' @param min.cells Default 10; The effective sample contains min cell numbers
#'
#' @return Complex dotplot of 3 groups' IGHV Usage  

runVDJ.diffUsage<-function(object,v.call="germline_v_call",
                           sample.label="samples",group.label="group",
                           GroupA,GroupB,reference,
                           min.cells=10){
  sample.size=table(object@meta.data[,sample.label] %>% as.vector()) %>% as.data.frame()
  rownames(sample.size)<-sample.size$Var1
  sample.size<-sample.size[c(GroupA,GroupB,reference),]
  keep.samples<-sample.size[sample.size$Freq>=min.cells,]
  keep.samples$group<-NA
  keep.samples[GroupA,]$group<-"GroupA"
  keep.samples[GroupB,]$group<-"GroupB"
  keep.samples[reference,]$group<-"reference"
  
  
  new<-object@meta.data
  new$IGHV<-NA
  i=1
  for(i in 1:dim(new)[1]){
    new[,v.call][i] %>%as.vector() %>%strsplit(split='\\*') ->tmp
    new[i,"IGHV"]<-tmp[[1]][1]
  }
  Usage<-prop.table(table(new$IGHV,new[,sample.label]))
  for(j in 1:dim(Usage)[2]){
    Usage[,j]<-Usage[,j]/sum(Usage[,j])
  }
  Usage <- as.data.frame.array(t(Usage[,rownames(keep.samples)]))
  Usage$group<-keep.samples$group
  Usage[,-dim(Usage)[2]]<-Usage[,-dim(Usage)[2]]*1000
  
  condition <- factor(Usage$group)
  colData <- data.frame(row.names=rownames(Usage), condition)
  Usage<-floor(Usage[,-dim(Usage)[2]])
  Usage<-t(Usage) %>% as.data.frame()
  
  dds <- DESeqDataSetFromMatrix(countData = Usage,
                                colData = colData, 
                                design= ~condition)
  dds2 <- DESeq(dds)
  res1 <- results(dds2, contrast=c("condition","GroupA","reference")) %>% as.data.frame()
  res2 <- results(dds2, contrast=c("condition","GroupB","reference"))%>% as.data.frame()
  
  
  data1<-na.omit(res1)
  data2<-na.omit(res2)
  
  
  m<-Usage
  z<-m %>% apply(2,sum)
  for(i in 1:dim(m)[2]){
    m[,i]<-m[,i] / as.numeric(z[i])
  }
  n<-m %>% apply(1,median) %>% as.data.frame()
  y<-intersect(data1 %>% rownames(),data2%>% rownames())
  k<-n[y,]
  
  data<-data.frame(X1=data1[y,]$log2FoldChange,X2=data2[y,]$log2FoldChange,row.names = y,Size=k)
  data$ID<-rownames(data)
  library(ggrepel)
  data$color<-"0"
  data[data$X1>0 &data$X2>0,]$color<-"1"
  data[data$X1<0 &data$X2>0,]$color<-"2"
  data[data$X1<0 &data$X2<0,]$color<-"3"
  data[data$X1>0 &data$X2<0,]$color<-"4"
  P1<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
    DimPlot_theme+
    geom_text_repel(data=data[data$X1 >0 & data$X2 >0,],aes(x=X1,y=X2,label=ID))+
    geom_hline(aes(yintercept=0))+
    geom_vline(aes(xintercept=0))+scale_color_manual(values = c("black","red","yellow","blue","green"))
  return(P1)
}
