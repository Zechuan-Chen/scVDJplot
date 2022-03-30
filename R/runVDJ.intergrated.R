#' The description of this runVDJ.intergrated.
#' @title runVDJ.intergrated
#' 
#' @description The runVDJ.intergrated function will integrated the result of change-o with the Seurat.S4 object.
#' 
#' @param tsv.file The file paths of change-o result. such as "~/*/filtered_heavy_germ-pass.tsv"
#' 
#' @param object The Seurat S4 object
#' 
#' @param runSHM select T/F; Determined to calculate the mutation frequency
#'
#' @return Integrated S4 object

runVDJ.intergrated<-function(tsv.file,Seurat.object,runSHM=T){
  db<-read.table(tsv.file,header=T,sep="\t")
  if(runSHM){
    db<-observedMutations(db,
                          sequenceColumn="sequence_alignment", 
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=NULL, frequency=T, 
                          combine = T, nproc=10)  
  }
  rownames(db)<-db$cell_id
  meta.data<-object@meta.data
  length(as.vector(db$cell_id))
  y<-intersect(as.vector(db$cell_id),rownames(object@meta.data))
  print("Total ",y," cells intergrated")
  rownames(db)<-as.vector(db$cell_id)
  new<-cbind(meta.data[y,],db[y,])
  object<-subset(object,cells = y)
  object@meta.data<-new
  
  
  object$c_call %>% as.vector() %>% substr(1,4) ->object$c_call2
  object@meta.data[object$c_call2=="",]$c_call2<-"Unknow"
  Idents(object)<-"c_call2"
  object<-subset(object,idents = c("IGHA","IGHD","IGHM","IGHG")) 
  levels(object)<-c("IGHD","IGHM","IGHA","IGHG")
  names(x)<-levels(object)
  object<-RenameIdents(object,x)
  object$c_call2<-Idents(object)
  return(object)
}