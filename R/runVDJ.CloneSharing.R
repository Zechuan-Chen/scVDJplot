#' The description of this runVDJ.CloneSharing
#' @title runVDJ.CloneSharing
#' 
#' @description The runVDJ.CloneSharing function will calculate the clonetype sharing between different cell types.
#' 
#' @param object The Seurat S4 object
#' 
#' @param cell.type.label The cell type annotation labels in Seurat S4 object
#' 
#' @param cell.type The selected cell type in annotation labels will be used to calculate the clonetype connection with other cell types.
#' 
#' @param clone.id Default "clone_id"; The result of change-o has been integrated in Seurat S4 object. The clonetypes information have been save in "clone_id" by default.
#' 
#' @return Heatmap of clone sharing between cell types

runVDJ.CloneSharing<-function(object,cell.type.label,clone.id){
  object$sub_cell_type2<-object@meta.data[,cell.type.label]
  object$clone_id<-object@meta.data[,clone.id]
  
  object$sub_cell_type2_clone_id<-paste0(object$sub_cell_type2,"__",object$clone_id)
  object$count<-1
  data2<-object@meta.data %>% select(c("sub_cell_type2_clone_id","sub_cell_type2","clone_id","count"))%>%group_by(sub_cell_type2,clone_id)%>%summarise_at(c("count"),funs(sum))  %>% as.data.frame()
  data2$sub_cell_type2_clone_id<-paste0(data2$sub_cell_type2,"__",data2$clone_id)
  Idents(object)<-"sub_cell_type2"
  x<-data.frame(row.names = levels(object))
  i=1
  for(i in 1:length(rownames(x))){
    label<-rownames(x)[i]
    tmp<-c()
    for(j in rownames(x)){
      table1<-data2[data2$sub_cell_type2==label,]
      table2<-data2[data2$sub_cell_type2==j,]
      rownames(table1)<-table1$clone_id
      rownames(table2)<-table2$clone_id
      clone_id_list<-intersect(table1$clone_id,table2$clone_id)
      tmp<-c(tmp,c(sum(table1[clone_id_list,]$count)+sum(table2[clone_id_list,]$count))/c(sum(table1$count)+sum(table2$count)))
    }
    x[,i]<-tmp
  }
  
  
  colnames(x)<-rownames(x)
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  x <- get_lower_tri(x)
  x<-x[dim(x)[1]:1,]
  library(data.table)
  melted_cormat <- melt(as.matrix(x), na.rm = TRUE)
  melted_cormat$value<-melted_cormat$value
  
  melted_cormat<-melted_cormat[melted_cormat$value>=0 & melted_cormat$value!=1,]
  melted_cormat$value <-round(melted_cormat$value,4)
  f<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")
  P1<-f+ scale_fill_gradient2(low = '#31888D', high = "#D53E4F", mid = '#FFF9B4', 
                              midpoint = max(melted_cormat$value)/2,
                              limit = c(0,max(melted_cormat$value)),
                              space = "Lab", 
                              name="Pearson\nCorrelation")+theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  P1<-P1+geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+RotatedAxis()
  return(P1)
}
