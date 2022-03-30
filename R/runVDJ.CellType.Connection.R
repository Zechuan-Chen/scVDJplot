#' The description of this runVDJ.CellType.Connection
#' @title runVDJ.CellType.Connection
#' 
#' @description The runVDJ.CellType.Connection function will caculate the clonetype connection between different cell types.
#' 
#' @param object The Seurat S4 object
#' 
#' @param cell.type.label The cell type annotation labels in Seurat S4 object
#' 
#' @param clone.id Default "clone_id"; The result of changeo has been intergrated in Seurat S4 object. The clonetypes information have been save in "clone_id" by default.
#' 
#' 
#' @return Heatmap of clone sharing between cell types

runVDJ.CellType.Connection<-function(object,cell.type.label,clone.id="clone_id",cell.type="all"){
  object$clone_id<-object@meta.data[,clone.id]
  object$sub_cell_type2<-object@meta.data[,cell.type.label]
  Idents(object)<-cell.type.label
  levels(object)
  Clonotype<-list()
  for( i in 1:length(levels(object))){
    object_1<-subset(object,idents=levels(object)[i])
    Clonotype[[i]]<-table(object_1$clone_id) %>% as.data.frame()
    Clonotype[[i]]<-Clonotype[[i]][Clonotype[[i]][,"Freq"]>0,]
  }
  
  
  x<-levels(object)
  num<-length(levels(object))-1
  X1<-(1+num)*num/2
  i=1
  group=0
  data1<-list()
  for( i in 1:num){
    data1[[i]]<-data.frame(CellType=c(rep(levels(object)[i],times=c(num+1-i)),levels(object)[c(i:num)+1]),
                           Group=c(rep(group+c(1:c(num+1-i)),times=2)))
    group=group+c(num+1-i)
  }
  data=data1[[1]]
  i=1
  for( i in 1:c(num-1)){
    data<-rbind(data,data1[[i+1]])
  }
  i=1
  data$UMAP_1<-0
  data$UMAP_2<-0
  for( i in 1:length(levels(object))){
    object_1<-subset(object,idents = levels(object)[i])
    UMAP_Table<-object_1@reductions$umap@cell.embeddings %>% as.data.frame()
    UMAP_1<-UMAP_Table$UMAP_1 %>% median()
    UMAP_2<-UMAP_Table$UMAP_2 %>% median()
    data[data$CellType==levels(object)[i],"UMAP_1"]<-UMAP_1
    data[data$CellType==levels(object)[i],"UMAP_2"]<-UMAP_2
  }
  i=1
  meta_total<-object@meta.data
  data$Common_clono_Freq=0
  for(i in 1:max(data$Group)){
    data2<-data[data$Group==i,]
    meta1<-meta_total[meta_total$sub_cell_type2==as.vector(data2$CellType[1]),]
    meta2<-meta_total[meta_total$sub_cell_type2==as.vector(data2$CellType[2]),]
    object_2<-subset(object,idents = as.vector(data2$CellType))
    object_2<-object_2@meta.data
    clono_Freq<-table(object_2$clone_id) %>% as.data.frame()
    rownames(clono_Freq)<-clono_Freq$Var1
    data[data$Group==i,]$Common_clono_Freq<-sum(clono_Freq[intersect(meta1$clone_id,meta2$clone_id),]$Freq)
  }
  data$Common_Freq_scale<-scale(data$Common_clono_Freq,center = F) %>% as.vector()
  data<-data[order(data$Group),]
  
  if(cell.type=="all"){
    P1<-DimPlot(object,label=T)+geom_line(data=data,aes(UMAP_1, UMAP_2, group =Group), size = data$Common_Freq_scale, alpha = 0.8,color = 'red')
    return(P1)
  }else{
    n=match(cell.type,levels(object))
    for(i in data[data$CellType==levels(object)[n],]$Group){
      if(i==data[data$CellType==levels(object)[n],]$Group[1]){
        data_1<-data[data$Group==i,]
      }else{
        data_1<-rbind(data_1,data[data$Group==i,])
      }
    }
    data_1<-data_1[order(data_1$Group),]
    P1<-DimPlot(object,label=T)+geom_line(data=data_1,aes(UMAP_1, UMAP_2, group =Group), size =data_1$Common_Freq_scale, alpha = 0.8,color = 'red')
    return(P1)
  }
}
