#' The description of this runVDJ.Usage.
#' @title runVDJ.Usage
#' 
#' @description The runVDJ.Usage function will calculate the IGHV fragment Usage.
#' 
#' @param object The Seurat S4 object
#' 
#' @param v.call Default "germline_v_call"; The result of change-o has been integrated in Seurat S4 object. The V region information have been save in "germine_v_call" by default.
#' 
#' @param sample.label The sample labels in Seurat S4 object.
#'
#' @param group.label The group labels in Seurat S4 object.
#'
#' @param min.cells Default 10; The effective sample contains min cell numbers.
#'
#' @return IGHV Usage box plot 

runVDJ.Usage<-function(object,v.call="germline_v_call",sample.label="samples",group.label="group",min.cells=10){
  sample.size=table(object@meta.data[,sample.label] %>% as.vector()) %>% as.data.frame()
  keep.samples<-sample.size[sample.size$Freq>=min.cells,1] %>% as.vector()
  new<-object@meta.data
  new$IGHV<-NA
  i=1
  for(i in 1:dim(new)[1]){
    new[,v.call][i] %>%as.vector() %>%strsplit(split='\\*') ->tmp
    new[i,"IGHV"]<-tmp[[1]][1]
  }
  Usage<-prop.table(table(new$IGHV,new[,sample.label]))
  Usage<-Usage[,keep.samples]
  for(j in 1:dim(Usage)[2]){
    Usage[,j]<-Usage[,j]/sum(Usage[,j])
    Usage_tmp<-data.frame(Samples = colnames(Usage)[j],
                          IGHV=rownames(Usage),
                          Usage=as.vector(Usage[,j]),
                          group= as.vector(object@meta.data[,group.label][match(colnames(Usage)[j],object@meta.data[,sample.label])]))
    if(j==1){Usage_t<-Usage_tmp}else{
      Usage_t<-rbind(Usage_t,Usage_tmp)
    }
  }
  p <- ggplot(Usage_t, aes(x=IGHV, y=Usage, fill=group)) +
    geom_boxplot(width=0.5, outlier.shape = NA) + 
    #labs(title ="Total") +
    labs(title ="Usage") +
    xlab("") + ylab("IGHV Usage") +
    theme_bw() + theme(legend.position = "none") +
    scale_fill_manual(values=dittoColors())
  #coord_cartesian(ylim = c(0,24))
  P1<-p+theme_bw()+theme(aspect.ratio=0.3, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+theme(axis.text.x = element_text(angle = 30, hjust = 1,vjust =1.0))
  return(P1)
}