

pie_plot <- function(dt){
#colnames must be Freq and type
#####Freq is Freq or percent , type is type or label
library(ggplot2)
library(ggsci)
#dt = dt[order(dt$Freq, decreasing = TRUE),]   
myLabel = paste(dt$type, " ( ", round(dt$Freq/sum(dt$Freq)*100,2), "% )   ", sep = "")   ## 用 round() 对结果保留两位小数  
p = ggplot(dt, aes(x = "", y = Freq, fill = type)) + 
  geom_bar(stat = "identity", width = 1) + scale_fill_jama(name="",label=myLabel)+ 
  coord_polar(theta = "y") + theme_void()+
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank())#, legend.position = "bottom")
return(p)
}



heatmap_house <- function(dt,colData="",title_hp="Heatmap of Different expressed",cluster_rule="no"){

	######suit two class 
	library("pheatmap")
	anno_color <-c("#E0640D", "#228443")
	names(anno_color) <- c("Tumor","Normal") 
	ann_colors = list(Type= anno_color)
  dt=dt[,as.character(rownames(colData))]
	if (cluster_rule=="row") {
		p <- pheatmap(log2(dt+1),cluster_rows=TRUE, cluster_cols=FALSE,,scale = "row",show_rownames=F,show_colnames = F,
         color = colorRampPalette(c(rep('#1C2B6F',3),'black', rep('#E31E26',3)))(50),
         annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
	} else if (cluster_rule=="col"){
		p <- pheatmap(log2(dt+1),cluster_rows=FALSE,cluster_cols=TRUE, scale = "row",show_rownames=F,show_colnames = F,
         color = colorRampPalette(c(rep('#1C2B6F',3),'black', rep('#E31E26',3)))(50),
          annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
	} else if (cluster_rule=="all"){
		p <- pheatmap(log2(dt+1),cluster_rows=TRUE,cluster_cols=TRUE, scale = "row",show_rownames=F,show_colnames = F,
         color = colorRampPalette(c(rep('#1C2B6F',3),'black', rep('#E31E26',3)))(50),
          annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
	} else if (cluster_rule=="no"){
		p <- pheatmap(log2(dt+1),cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",show_rownames=F,show_colnames = F,
         color = colorRampPalette(c(rep('#1C2B6F',3),'black', rep('#E31E26',3)))(50),
         annotation_col=colData,annotation_colors = ann_colors,main = title_hp)

	}

	return(p)
}


box_plot <- function(df,x,y,fill){
library(ggplot2)
library(ggsci)
attach(df)
p <-  ggplot2::ggplot(df,mapping = ggplot2::aes(x = x, y = y)) +scale_fill_jama()+
  ggplot2::geom_boxplot()+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.9,0.85))
detach(df)
  return(p)
}

###################
volcano_plot <- function(res_df,contrast_factor,pval = 0.05,fc=1.5){
  ######res_df is DE result data.frame with all genes (with non-sig gene)
  #contrast_factor=c("normal","tumor")
  i=1
  lfc=signif(log2(fc),2)    
  tab = data.frame(logFC = res_df$logFC, negLogPval = -log10(res_df$FDR))
  nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
  signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
  signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
  gap = max(tab$logFC)/50
  up_count = length(which(signGenes_up))
  down_count = length(which(signGenes_down))
  #plot
  par(mar = c(5, 6, 5, 5))
  plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
  points(tab[nosigGene, ], pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = "grey")
  if (length(unique(signGenes_up)) > 1){
    points(tab[signGenes_up, ], pch = 21, col = "black", bg = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 21, col = "black", bg = "cornflowerblue")
  }
  abline(h = -log10(pval), col = "green3", lty = 2)
  abline(v = c(-lfc, lfc), col = "orange", lty = 2)
  
  mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
  mtext(c(as.character(contrast_factor[2]), as.character(contrast_factor[1])), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
  mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
  legend("top",legend = c("Up regulate","Down regulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  mtext(c(as.character(contrast_factor[2]), as.character(contrast_factor[1])), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
  mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
  legend("top",legend = c("Up regulate","Down regulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  
}

################################
#clussterprofile

GoKegg=function(gene_list,OrgDb,organism,outdir){
  library(clusterProfiler)
  library(DOSE)
  library(ReactomePA)
  library(pathview)
  setwd(outdir)
  if (!dir.exists('GO_KEGG')){
    dir.create('GO_KEGG')
  }
  #options(bitmapType = "cairo")
  gene_id_list=bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene_id_list<-as.character(gene_id_list[,2])
  #KEGG
  kk <- enrichKEGG(gene = gene_id_list,
                   organism ="human",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   minGSSize = 1,
                   use_internal_data =FALSE)
  write.table(as.data.frame(kk@result), file=paste("GO_KEGG/kegg.txt",sep = ""),sep = "\t",quote = F)
  #pdf(paste("GO_KEGG/"kegg.pdf",sep = ""),bg = "transparent")
  pdf("GO_KEGG/kegg.pdf")
  print(barplot(kk,drop=TRUE,showCategory = 10))
  dev.off()
  cat ("cluster: kegg is ok~\n")
  
  pdf("GO_KEGG/GO enrichment analysis.pdf")
  for (myont in c("CC","BP","MF")){
    ego <- enrichGO(gene=gene_id_list,
                    OrgDb=org.Hs.eg.db,
                    ont = myont,
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.1,
                    readable = TRUE)
    write.table(as.data.frame(ego@result), file=paste("GO_KEGG/clust_",myont,".txt",sep = ""),sep = "\t",quote = F)
    print(barplot(ego,drop=TRUE,showCategory = 10)+labs(title=paste("enrichment analysis of GO",myont)))
  }
  
  dev.off()
  cat("cluster: GO is ok~\n")
}

####################correlation plot

  Cor_plot <- function(CoxExpPlotData,gene1,gene2,cormethod="spearman"){
    library("ggplot2")
  	cor_df = data.frame(CoxExpPlotData, 
  		gene1 = CoxExpPlotData[,gene1],
  		gene2 = CoxExpPlotData[,gene2])
  	pcutoff=0.05
  	CoxTest = cor.test(CoxExpPlotData[, gene1], CoxExpPlotData[, gene2], method = cormethod)

    if(cormethod == "pearson"){
      plottitle <- paste("R = ", signif(CoxTest$estimate[[1]], 3), "\nP value = ", signif(CoxTest$p.value, 4), sep = "")
    }else if(cormethod == "spearman"){
      plottitle <- paste("rho = ", signif(CoxTest$estimate[[1]], 3), "\nP value = ", signif(CoxTest$p.value, 4), sep = "")
    }
    
    print(paste("pvalue =", round(CoxTest$p.value, 4)))
    if(CoxTest$p.value < pcutoff){
      CoxExpPoint = ggplot(data = cor_df, aes(x= gene1, y = gene2))+theme_classic()+
        geom_point(size = 2, color = "gray36")+
        ggtitle(plottitle)+
        stat_smooth(method = lm, se = FALSE, colour = "#191970")+
        ylab(paste(gene2, "Exp.", sep = " "))+xlab(paste(gene1, "Exp.",  sep = " "))+
        scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
        theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black", hjust = 0.5, vjust = -2.5),
          axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))}
    return(CoxExpPoint)
  }  

##############ECDF plot

ECDF_plot <- function(df,value_var,group_var,plot_title,ks_test=""){
  p <- ggplot(df,aes(x=value_var,group=group_var,color=group_var))+theme_test()+
    stat_ecdf( size = 1)+theme(legend.position=c(0.9,.1))+
    scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+annotate("text",x=.45,y=.5,label=ks_test)+
  labs(title=plot_title, y="Cumulative fraction")
  return(p)
}


###################edgeR for DE analysis with mode (paired or unpaired)

get_pheno <- function(x){
  #####x is colnames of expr_mat
  temp_T <-cbind(x[grep("T",x)],Type=rep("Tumor",length(x[grep("T",x)])))
  temp_N <-cbind(x[grep("N",x)],Type=rep("Normal",length(x[grep("N",x)])))
  pheno = as.data.frame(rbind(temp_T,temp_N))
  pheno_data = as.data.frame(pheno$Type)
  rownames(pheno_data)=sub("X","S",as.character(pheno[,1]),fixed = T)
  names(pheno_data) <- "Type"
  return(pheno_data)
}

edgeR_test <- function(expre_mat,group_mat="",design_mode = "unpaired",test_method="QLFit",adjust_method = "BH" , pvalue = 0.05, lfc = 0.58){
  if (group_mat==""){
      print( "use default pheno_mat : sample_id should end with T or N")
      group_mat <- get_pheno(colnames(expre_mat))
      }
  #group_mat <- get_pheno(colnames(expre_mat))
  expre_mat <- expre_mat[,as.character(rownames(group_mat))]
  treatment_inFunc <- as.factor(group_mat$Type)
  patient <- as.factor(gsub("[T|N]","",rownames(group_mat)))
  deg_lst <- DGEList(counts = expre_mat,genes = as.character(rownames(expre_mat)))
  #fliter &TMM 
  keep <- rowSums(cpm(deg_lst)>0) >= 2 #a CPM>1 in at least 2 samples
  deg_lst <- deg_lst[keep,]
  deg_lst <- calcNormFactors(deg_lst)
  #design matrix  ###########treatment VS control factor should be place on last column
  if (design_mode=="paired"){
    design_mat <- model.matrix(~patient+treatment_inFunc)
  } else {
    design_mat <- model.matrix(~treatment_inFunc)
  }
  rownames(design_mat)<-colnames(deg_lst)
  
  #disper and test  
  deg_lst<-estimateDisp(deg_lst,design_mat)
  fit_in_func <- glmQLFit(deg_lst,design_mat)
  if ( test_method=="QLFit"){
    mode_res <- glmQLFTest(fit_in_func)
    topTags(mode_res)
  } else if (test_method=="LRT") {
    mode_res <- glmLRT(fit_in_func)
    topTags(mode_res)
  } else {
    cat("please provide test method : QLFit or LRT")
  }
  de_inFunc<-decideTestsDGE(mode_res,adjust.method = "none" , p.value = pvalue, lfc = lfc)
  print("non-adjust")
  print(summary(de_inFunc))
  de_inFunc<-decideTestsDGE(mode_res,adjust.method = adjust_method , p.value = pvalue, lfc = lfc)
  print(adjust_method)
  print(summary(de_inFunc))
  res_inFunction <- data.frame(mode_res$table, FDR = p.adjust(mode_res$table$PValue, method=adjust_method))
  return(res_inFunction)
}


#############PCA
 PCA_plot<- function(DE_list,df,colData){
  #colData <- subset(colData, Type == comb[1,i] | Type == comb[2,i])
  library(ggplot2)
  df <- df[DE_list, rownames(colData)] ########normalized matrix
  df = log2(df+1)
  pcaData <- as.data.frame(prcomp(df[DE_list,])$rotation)
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=colData$Type)) +
    geom_point(size=3) +
    xlab("PC1") +
    ylab("PC2") +
    scale_colour_hue("Type") +
    #  coord_fixed() +
    theme_bw()
  print(pca_plot)
  pca_plot_text <- ggplot(pcaData, aes(PC1, PC2, color=colData$Type)) +
    geom_text(aes(label = row.names(pcaData))) +
    xlab("PC1") +
    ylab("PC2") +
    scale_colour_hue("Type") +
    #  coord_fixed() +
    theme_bw()
  print(pca_plot_text)

}
