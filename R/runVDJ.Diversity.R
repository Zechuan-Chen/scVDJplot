#' The description of this runVDJ.Diversity
#' @title runVDJ.Diversity
#' 
#' @description The runVDJ.Diversity function will calculate the clonetype diversity.
#' 
#' @param object The Seurat S4 object
#' 
#' @param cell.type.label The cell type annotation labels in Seurat S4 object
#' 
#' @param clone.id Default "clone_id"; The result of change-o has been integrated in Seurat S4 object. The clonetypes information have been save in "clone_id" by default.
#' 
#' @return Bar plot of clonetype diversity

runVDJ.Diversity<-function(object,cell.type.label,clone.id='clone_id'){

  object$clone_id<-object@meta.data[,clone.id]
  Idents(object)<-cell.type.label
  for(i in 1:length(levels(object))){
    object_1<-subset(object,idents=levels(object)[i])
    if(i==1){
      data<-object_1$clone_id %>% table() %>% as.data.frame()
      Diversity_number<-entropy(data$Freq)
      #Clonality_number<-1-entropy(data$Freq)/log2(dim(data[data$Freq==1,])[1])
    }else{
      data<-object_1$clone_id %>% table() %>% as.data.frame()
      Diversity_number<-c(Diversity_number,entropy(data$Freq))
      #Clonality_number<-c(Clonality_number,1-entropy(data$Freq)/log2(dim(data[data$Freq==1,])[1]))
    }
    
  }
  data<-data.frame(CellType = levels(object),
                   Diversity = Diversity_number)
  P1<-ggplot(data = data, mapping = aes(
    x = CellType, y = Diversity ,fill=CellType))+geom_col()+
    scale_fill_manual(values = color_panal)+theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+RotatedAxis()
  return(P1)
}