

# scVDJplot
This repository was used to store the code, which was used in "Single-cell Profiling Reveals Distinct Adaptive Immune Hallmarks in MDA5+ Dermatomyositis with Therapeutic Implications"

---
## **Attention:** 
All the code has been integrated into the R packages "scSensitiveGeneDefine";

This repository will be renamed as "scSensitiveGeneDefine"

## Description:
`scVDJplot` is a R package which can be used in single cell VDJ sequencing caculating and  visualization. 

`scVDJplot` is build based on Seurat >= 3.0.1; data.table >= 1.14.2; entropy >= 1.2.1; ;dplyr >= 1.0.0; shazam >= 1.1.0; alakazam >= 1.2.0; All of these four dependent packages are R package.

`scVDJplot` intend to publish on Nature Communications.

## Installation(in R/Rstudio)
devtools::install_github("Zechuan-Chen/scVDJplot")

## Dependencies
`scVDJplot` requires the following R packages:

 - Seurat (>=3.0.1)
 - data.table （>= 1.14.2）
 - dplyr (>=1.0.0)
 - entropy (>=1.2.1)
 - shazam (>= 1.1.0) 
 - alakazam (>= 1.2.0)
 - NOTE:The version of these depend packages are temporary.

Example code for `scVDJplot`

```
# The count data and meta.data have saved
# We can get data from /data/count.RData  and /data/meta.data

# Step1: intergrated change-o result with Seurat S4 object
object<-runVDJ.intergrated(object,tsv.file="~/changeo-result/filtered_heavy_germ-pass.tsv",Seurat.object=object,runSHM=T) 

# Or you can get a simple data from our package 
load("count.RData")
load("meta.data.RData")
object<-CreateSeuratobject(count=count,meta.data=meta.data)
object %>% NormalizeData() %>% FindVariableFeatures() %>%  # Common Seurat pipline for visualization
  ScaleData(assay="RNA")  %>% 
  RunPCA(npcs = 40,assay="RNA")%>% 
  FindNeighbors(assay="RNA",dims = 1:40)%>%
  FindClusters(assay="RNA",resolution = 0.6)%>%
  RunUMAP(dims=1:40) -> object
```

```
# Step2: Visualization of the IGHV usage in each groups
P1<-runVDJ.Usage(object,v.call="germline_v_call",sample.label="samples",group.label="group",min.cells=10)
```
![image](https://github.com/Zechuan-Chen/scVDJplot/blob/main/picture/1.png)


```
# Step3: Visualization of the IGHV usage difference in 3 groups
P2<-runVDJ.diffUsage(object,v.call="germline_v_call",
                 sample.label="samples",group.label="group",
                 GroupA=c("P05","P06","P07"),GroupB=c("P014","P15","P16"),reference=c("P01","P02","P03"),
                 min.cells=10)
```
![image](https://github.com/Zechuan-Chen/scVDJplot/blob/main/picture/2.png)

```
# Step4: Visualization of the selection strength
P3<-runVDJ.Selec.Strength(object,group.label='group',
                                isotype.label="c_call",isotype="IGHM",
                                cell.type.label="sub_cell_type2",cell.type=" scB8-ASC",clone.id="clone_id")
```
![image](https://github.com/Zechuan-Chen/scVDJplot/blob/main/picture/3.png)
```

# Step5: Caculate and visualizate clonality and diversity
P4<-runVDJ.Diversity(object,cell.type.label="sub_cell_type2",clone.id='clone_id')
P4<-runVDJ.Clonality(object,cell.type.label="sub_cell_type2",clone.id='clone_id')
```

```
# Step6: Caculate and visualizate the clonetype sharing between different cell types.
P5<-runVDJ.CloneSharing(object,cell.type.label = "sub_cell_type2",clone.id="clone_id")
P6<-runVDJ.CellType.Connection(object,cell.type.label = "sub_cell_type2",clone.id="clone_id")
```
![image](https://github.com/Zechuan-Chen/scVDJplot/blob/main/picture/4.png)
![image](https://github.com/Zechuan-Chen/scVDJplot/blob/main/picture/5.png)