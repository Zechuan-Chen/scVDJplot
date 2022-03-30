#' The description of this runVDJ.CloneFreqLevels
#' @title runVDJ.CloneFreqLevels
#' 
#' @description The runVDJ.CloneFreqLevels function will calculate clonetype frequency levels.
#' 
#' @param object The Seurat S4 object
#' 
#' @param clone.id Default "clone_id"; The result of change-o has been integrated in Seurat S4 object. The clonetypes information have been save in "clone_id" by default.
#'
#' @return Seurat S4 object with CloneFreqLevels labels saved   

runVDJ.CloneFreqLevels<-function(object,clone.id="clone_id"){
  new<-object@meta.data
  new$clone_id<-new[,clone.id]
  colnames(new)
  new$clone_id<-as.vector(new$clone_id)
  q<-as.data.frame(table(new$clone_id))
  new$CloneType_Freq<-NA
  for(i in as.vector(q$Var1)){
    new[new$clone_id==i,]$CloneType_Freq<-q[q$Var1==i,]$Freq
  }
  object@meta.data<-new
  object$CloneType_Freq_levels<-NA
  object@meta.data[object$CloneType_Freq>0,]$CloneType_Freq_levels<-"1"
  object@meta.data[object$CloneType_Freq>1,]$CloneType_Freq_levels<-"2"
  object@meta.data[object$CloneType_Freq>2,]$CloneType_Freq_levels<-"3"
  object@meta.data[object$CloneType_Freq>3,]$CloneType_Freq_levels<-"4"
  object@meta.data[object$CloneType_Freq>4,]$CloneType_Freq_levels<-"5"
  object@meta.data[object$CloneType_Freq>5,]$CloneType_Freq_levels<-"6"
  object@meta.data[object$CloneType_Freq>6,]$CloneType_Freq_levels<-"7"
  object@meta.data[object$CloneType_Freq>7,]$CloneType_Freq_levels<-"8"
  object@meta.data[object$CloneType_Freq>8,]$CloneType_Freq_levels<-"9"
  object@meta.data[object$CloneType_Freq>9,]$CloneType_Freq_levels<-">9"
  return(object)
}