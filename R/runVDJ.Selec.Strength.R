#' The description of this runVDJ.Selec.Strength
#' @title runVDJ.Selec.Strength
#' 
#' @description The runVDJ.Selec.Strength function will caculate the selection strength in each group.
#' 
#' @param object The Seurat S4 object
#' 
#' @param group.label The group labels in Seurat S4 object
#'
#' @param isotype.label Default "c_call"; The result of change-o has been integrated in Seurat S4 object. The isotype information have been save in "c_call" by default
#' 
#' @param cell.type.label The cell type annotation labels in Seurat S4 object
#' 
#' @param cell.type The selected cell type in annotation labels will be used to calculate the selection strength
#' 
#' @param clone.id Default "clone_id"; The result of change-o has been integrated in Seurat S4 object. The clonetypes information have been save in "clone_id" by default.
#' 
#' @return Complex plot selection strength in cdr3 and fwr  

runVDJ.Selec.Strength<-function(object,group.label='group',
                                isotype.label,isotype="IGHM",
                                cell.type.label,cell.type,clone.id="clone_id"){
  db <- object@meta.data
  db <- db[db[,isotype.label]==isotype & db[,cell.type.label]==cell.type,]
  # -------------------------------------------------------------------------
  
  db <- collapseClones(db, cloneColumn=clone.id,
                       sequenceColumn="sequence_alignment",
                       germlineColumn="germline_alignment_d_mask",
                       method="thresholdedFreq", minimumFrequency=0.6,
                       includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
  
  # Calculate BASELINe
  baseline <- calcBaseline(db, 
                           sequenceColumn="clonal_sequence",
                           germlineColumn="clonal_germline", 
                           testStatistic="focused",
                           regionDefinition=IMGT_V,
                           targetingModel=HH_S5F,
                           nproc=1)
  
  # Grouping the PDFs by the sample and isotype annotations
  grouped <- groupBaseline(baseline, groupBy=c(group.label,isotype.label))
  
  P1<-plotBaselineDensity(grouped,  group.label, isotype.label,colorValues=dittoColors(), 
                          colorElement=group.label, sigmaLimits=c(-1, 1))
  return(P1)
}
