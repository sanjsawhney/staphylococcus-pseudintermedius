# staphylococcus-pseudintermedius
This repository contains all code used in Sawhney &amp; Vargas et al. (2023).

# -------------------------Figure 1B, Suppl Figure 3B: Accessory_Genome_PCoA---------------------------------------------------

# load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)

# Read in accessorygenome.csv (.Rtab file)
df_accessorygenome<-read.csv('211207_accessorygenome_pseud_pgap_488.csv',
                            header = T,
                            row.names = 1)

# Calculate jaccard distance. Set to matrix
jaccard_accessorygenome<-vegdist(df_accessorygenome, method='jaccard',binary=TRUE)
jaccard_accessorygenome<-as.matrix(jaccard_accessorygenome)
write.csv(jaccard_accessorygenome,"jaccard_distance_501_pseud_delph.csv")

# Make PCoA (correction is to account for negative eigenvalues) - can use "lingoes" or "cailliez"
# Look through pcoa_accessorygenome_corr on first use to gain an understanding of what $values, $vectors, and $vectors.cor mean
pcoa_accessorygenome_corr<-pcoa(jaccard_accessorygenome, correction = "lingoes")

# Get PCoA vectors to plot ordination as axis x axis. Set to data frame
pcoavectors_accessorygenome_corr<-pcoa_accessorygenome_corr$vectors.cor
pcoavectors_accessorygenome_corr<-as.data.frame(pcoavectors_accessorygenome_corr)

# Add cohort metadata
write.csv(pcoavectors_accessorygenome_corr[,c("Axis.1","Axis.2")],"pcoa_values.csv")
pcoavectors_accessorygenome_corr=read.csv("pcoa_values.csv", sep=",")

# Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 1.0000
rel_eigen_accessorygenome_corr<-pcoa_accessorygenome_corr$values$Rel_corr_eig
rel_eigen_accessorygenome_corr


# Scree plot
barplot(pcoa_accessorygenome_corr$values$Rel_corr_eig[1:10])
biplot.pcoa(pca_accessorygenome_corr, df_accessorygenome)

# Read in metadata.csv
#metadata for Figure 1B
metadata_df<-read.csv('sig_metadata_493_pseud_redefined_only.csv',
                      sep=",",
                      header = T,
                      row.names = 1)

# plot Figure 1B
ggplot(
  pcoavectors_accessorygenome_corr,aes(x=Axis.1, y=Axis.2))+
  geom_point(shape=21, color = "black", size = 1, stroke = 1.2, show.legend = TRUE)+
  theme_bw()+
  theme(axis.title=element_text(size=16,face="bold"), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.text=element_text(face="bold"), legend.position = c(.18,.12), plot.title = element_text(hjust = 0.5,size=20,face="bold"))+
    ggtitle("Accessory Gene Content")+xlab("PCA1 (3.78%)")+ylab("PCA2 (3.41%)")+
  coord_fixed((ratio=1))+
  geom_point(aes(color=metadata_df$Redefined_Cohort))+
  scale_color_manual(values=c("Human Diagnostic"="#C0392B","Animal Diagnostic"="#FF9388", "Household (Col/Env)"="#76D7C4"))
  


# plot Figure S3B
ggplot(
  pcoavectors_accessorygenome_corr,aes(x=Axis.1, y=Axis.2))+
  geom_point(shape=21, color = "black", size = 0.75, stroke = 1.5, show.legend = TRUE)+
  theme_bw()+
  theme(axis.title=element_text(size=16,face="bold"), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(face="bold"), legend.position = c(.165,.15), plot.title = element_text(hjust = 0.5,size=20,face="bold"))+
  ggtitle("Accessory Gene Content")+xlab("PCA1 (3.78%)")+ylab("PCA2 (3.41%)")+
  coord_fixed((ratio=1))+
  geom_point(aes(color=Full_Cohort_Simplified))+
#v1:  scale_color_manual(values=c("Human Colonizing"="#A6CEE3","Human Diagnostic"="#1F78B4","Pet Colonizing"="#FDBF6F","Animal Diagnostic"="#FF7F00","Environment: Bathroom Countertop"="#FFFFFF","Environment: Bathroom Light Switch"="#f2f2f2","Environment: Bathroom Sink Handle"="#e6e6e6","Environment: Bathroom Toilet Handle"="#cccccc","Environment: Bathroom Toilet Seat"="#b3b3b3","Environment: Bathtub"="#a6a6a6","Environment: Bedsheets"="#999999","Environment: Computer, Keyboard, Mouse"="#808080","Environment: Kitchen Fridge Door Handle"="#666666","Environment: Kitchen Sink Handle"="#4d4d4d","Environment: Kitchen Table"="#333333","Environment: Telephone"="#1a1a1a","Environment: TV Remote"="#000000"))
#v2:  scale_color_manual(values=c("Human Colonizing"="#A6CEE3","Human Diagnostic"="#1F78B4","Pet Colonizing"="#FDBF6F","Animal Diagnostic"="#FF7F00","Environment"="white"))
  scale_color_manual(name = "Host type",values=c("Human"="#1F78B4","Animal"="#FF7F00","Environment"="white"))

#Adonis test for significance in clustering (permANOVA)
adonis2(jaccard_accessorygenome~metadata_df$Redefined_Cohort,pcoavectors_accessorygenome_corr[,c(1,2)],permutations=1000)


# ----------------------------------------Figure 1C, Suppl Figure 3C: Pangenome_Boxplot----------------------------------
library(ggplot2)

# read in boxplot csv
sig_pangenome_boxplot<-read.csv('210608_pangenome_boxplot.csv',
                            sep=",",
                            header = T)

sig_pangenome_boxplot$Cohort <- factor(sig_pangenome_boxplot$Cohort,levels = c("AI", "HI", "SF"))

boxplot_colors <- c("#5E5E5E","#E0E0E0","#5E5E5E","#E0E0E0","#5E5E5E","#E0E0E0")
geom_jitter_colors <- c("black","black","black","black","black","black")

ggplot(sig_pangenome_boxplot, aes(x=Cohort, y=Jaccard))+
  geom_boxplot(width=0.75, lwd=0.375, aes(fill=Comparison),position = position_dodge(0.9))+
  scale_fill_manual(values=boxplot_colors)+
  scale_y_continuous(name="Jaccard dissimilarity score", limits = c(0,0.37))+
  #geom_jitter(shape=16, alpha=0.15, size=2, width=0.4, aes(color=Cohort))+
  scale_color_manual(values=geom_jitter_colors)+
  theme_bw()+
  scale_x_discrete(labels=c("Diagnostic -\nAnimal", "Diagnostic -\nHuman", "Household - Colonizing\nor Environmental"))+
  theme(text=element_text(size=15, face="bold"), legend.position="right", panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), legend.title=element_blank())

library(reshape2)
sig_pangenome_boxplot<-read.csv('220921_Jaccard_pangenome.csv',
                                sep=",",
                                header = T)
                                
# Convert matrix to long format
sig_pangenome_boxplot <- melt(sig_pangenome_boxplot)

# Add column names
colnames(sig_pangenome_boxplot) <- c("isolate1","isolate2","jaccard_dissimilarity_score")

# Delete all self-comparisons
sig_pangenome_boxplot<- sig_pangenome_boxplot[!(sig_pangenome_boxplot$isolate1==sig_pangenome_boxplot$isolate2),]

# Remove colnames
names(sig_pangenome_boxplot)<-NULL

# write csv
write.table(sig_pangenome_boxplot,"220921_Jaccard_pangenome_melt.txt",sep = "\t",row.names = FALSE)

## Plot Fig 1C Jaccard boxplot - DIAGNOSTIC VS COLONIZING VS ENVIRONMENTAL
sig_pangenome_boxplot_filtered<-read.csv('220922_Jaccard_pangenome_melt_filtered_boxplot_DIAGNOSTIC.txt',
                                         sep="\t",
                                         header = T)
library(ggplot2)

sig_pangenome_boxplot_filtered$Cohort <- factor(sig_pangenome_boxplot_filtered$Cohort,levels = c("Diagnostic", "Colonizing", "Environment"))
sig_pangenome_boxplot_filtered$Group <- factor(sig_pangenome_boxplot_filtered$Group,levels = c("AD-AD", "HD-HD","AD-HD","Diag-Col/Env","PC-PC","HC-HC","HC-PC","Col-Diag/Env","Env-Env","Env-PC/HC","Env-Diag","empty4"))


#boxplot_colors <- c("#E0E0E0","#5E5E5E","#E0E0E0","#5E5E5E","#E0E0E0","#5E5E5E")
boxplot_colors <- c("#FFFFFF","#FFFFFF","#CCCCCC","#666666","#FFFFFF","#FFFFFF","#CCCCCC","#666666","#FFFFFF","#CCCCCC","#666666","white")

ggplot(sig_pangenome_boxplot_filtered, aes(x=Cohort, y=Jaccard, fill=Group))+
  geom_violin(position=position_dodge(0.8),aes(fill=Group),width=0.75)+
  scale_fill_manual(values=boxplot_colors,labels=c("Animal Diagnostic <->\nAnimal Diagnostic", "Human Diagnostic <->\nHuman Diagnostic","Animal Diagnostic <->\nHuman Diagnostic", "ANY Diagnostic <->\nANY Colonizing OR\nANY Environment","Pet Colonizing <->\nPet Colonizing","Human Colonizing <->\nHuman Colonizing", "Pet Colonizing <->\nHuman Colonizing", "ANY Colonizing <->\nANY Diagnostic OR\nANY Environment","Environment <->\nEnvironment","Environment <->ANY\nColonizing (Household)","Environment <->\nANY Diagnostic"))+
  geom_boxplot(width=0.35, lwd=0.375, position = position_dodge(0.8))+
  scale_y_continuous(name="Jaccard dissimilarity score", limits = c(0,.34))+
  #geom_jitter(shape=16, alpha=0.15, size=2, width=0.4, aes(color=Cohort))+
  scale_color_manual(values=geom_jitter_colors)+
  theme_bw()+
  #scale_x_discrete(labels=c("Animal", "Human", "Environment"))+
  theme(legend.spacing.y = unit(.1, 'cm'),legend.key.size = unit(0.75,"line")) + guides(fill = guide_legend(byrow = TRUE))+
  theme(text=element_text(size=15, face="bold"), legend.position="right", legend.text=element_text(size=6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), legend.title=element_blank())


# Plot Suppl Figure 3C Jaccard boxplot - HOST SPECIES

library(ggplot2)

# read in filtered, long-format data file
sig_pangenome_boxplot_filtered<-read.csv('220921_Jaccard_pangenome_melt_filtered_boxplot_FS1C.txt',
                                sep="\t",
                                header = T)

sig_pangenome_boxplot_filtered$Cohort <- factor(sig_pangenome_boxplot_filtered$Cohort,levels = c("Animal", "Human", "Environment"))
#sig_pangenome_boxplot_filtered$Cohort <- factor(sig_pangenome_boxplot_filtered$Cohort,levels = c("Diagnostic", "Colonizing", "Environment"))
sig_pangenome_boxplot_filtered$Group <- factor(sig_pangenome_boxplot_filtered$Group,levels = c("Within-group", "Between-group"))
#sig_pangenome_boxplot_filtered$Group <- factor(sig_pangenome_boxplot_filtered$Group,levels = c("AD-AD", "PC-PC","AD-PC","Animal-Human/Env","HD-HD","HC-HC","HD-HC","Human-Animal/Env","empty1","Env-Env","Env-Diag","empty4"))


boxplot_colors <- c("#E0E0E0","#5E5E5E","#E0E0E0","#5E5E5E","#E0E0E0","#5E5E5E")
#boxplot_colors <- c("#FFFFFF","#FFFFFF","#CCCCCC","#666666","#FFFFFF","#FFFFFF","#CCCCCC","#666666","white","#FFFFFF","#666666","white")

#geom_jitter_colors <- c("black","black","black","black","black","black")

ggplot(sig_pangenome_boxplot_filtered, aes(x=Cohort, y=Jaccard, fill=Group))+
  geom_violin(position=position_dodge(0.8),aes(fill=Group),width=0.75)+
  scale_fill_manual(values=boxplot_colors)+
  geom_boxplot(width=0.35, lwd=0.375, position = position_dodge(0.8))+
  scale_y_continuous(name="Jaccard dissimilarity score", limits = c(0,.34))+
  #geom_jitter(shape=16, alpha=0.15, size=2, width=0.4, aes(color=Cohort))+
  scale_color_manual(values=geom_jitter_colors)+
  theme_bw()+
 #scale_x_discrete(labels=c("Animal", "Human", "Environment"))+
  theme(text=element_text(size=15, face="bold"), axis.title.x=element_blank(),  legend.position="bottom", legend.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), legend.title=element_blank())


# ----------------------------------------Figure 1D: Pangenome breakdown by Gene carriage----------------------------------------
library(ggplot2)

# Create data frame
df_SIG_pangenome<-read.csv('220904_SIG_pangenome.csv',
                           sep=",",
                           header = T)

df_SIG_pangenome$Gene <- factor(df_SIG_pangenome$Gene, levels = c("Core","Soft_core","Shell","Cloud"))

# Plot Figure 1D
ggplot(data=df_SIG_pangenome, aes(x=Cohort, y=Percent, fill=Gene)) +
  geom_bar(stat="identity",color="black", width=0.5)+
  theme_classic()+
  scale_fill_brewer(palette="Greys", labels=c("Core\n(≥99%)","Soft core\n(95-99%)","Shell\n(15-95%)","Cloud\n(≤15%)"))+
  #scale_fill_manual(values=c(Core="#66CC99",Soft_core="#339966",Shell="#006633",Cloud="#003300"))+
  scale_y_continuous(limits=c(0,101),expand = c(0,0))+
  theme(axis.title = element_text(size=13,face="bold"), axis.text.x = element_blank(), axis.title.x =element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=7))+
  ylab("Percent of Pangenome")


# ----------------------------------------Figure 1E: COG function---------------------------------

bar_cog = read.csv("210608_eggnog.csv",
                   sep=",",
                   header = T)

bar_cog$COG <- factor(bar_cog$COG, levels = c("Unknown","Amino acid metabolism and transport","Translation","Transcription","Inorganic ion transport and metabolism","Replication and repair","Energy production and conversion","Carbohydrate metabolism and transport","Cell wall/membrane/envelop biogenesis","Nucleotide metabolism and transport","Coenzyme metabolism","Post-translational modification, protein turnover, chaperone functions","Lipid metabolism"))
bar_cog$Cohort  <- factor(bar_cog$Cohort, levels = c("Diagnostic – Human","Diagnostic – Animal","Household – Colonizing or Environmental"))
bar_cog$Percentage <- as.numeric(bar_cog$Percentage)

ggplot(bar_cog, aes(COG, Percentage))+
  geom_bar(width=.65, stat = "identity", aes(fill=Cohort), position = 'dodge', color='black')+
  theme_bw()+
  theme(text = element_text(size=11,face = "bold"),axis.text.x = element_text(angle = 330, hjust = 0), panel.grid.minor = element_blank(), panel.border = element_rect(colour = 'black'), panel.grid.major = element_blank())+
  scale_fill_manual(values=c("#C0392B","#FF9388","#76D7C4"))+
  scale_x_discrete(labels=c("Unknown","Amino acid metabolism and transport","Translation","Transcription","Inorganic ion transpor nand metabolism","Replication and repair","Energy production and conversion","Carbohydrate metabolism and transport","Cell wall/membrane/envelop biogenesis","Nucleotide metabolism and transport","Coenzyme metabolism","PTM, protein turnover, chaperone","Lipid metabolism"))+
  scale_y_continuous(limits = c(0,25), expand = c(0,0))


# ----------------------------------------Figure 2D: ARG count histogram distribution----------------------------------------
library(ggplot2)

# Create data frame
df_arg_count_mrsp<-read.csv('r_ARGct_MRSP_isolates.csv',
                           sep=",",
                           header = T)

df_arg_count_mrsp$ARG_ct <- factor(df_arg_count_mrsp$ARG_ct, levels = c("0","1","2","3","4","5","6","7","8","9","10","11"))
df_arg_count_mrsp$Cohort <- factor(df_arg_count_mrsp$Cohort, levels = c("MSSP","MRSP"))


# Plot Figure 2D
ggplot(data=df_arg_count_mrsp, aes(x=ARG_ct, y=Isolate_ct,fill=Cohort)) +
  geom_col(width = 0.45,color="black",size=0.7)+
  geom_vline(xintercept = 9.5,linetype="dashed", color = "black")+
  theme_classic()+
  scale_fill_manual(values=c("white","#996699"))+
  scale_y_continuous(limits=c(0,200),expand = c(0,0))+
  theme(axis.title = element_text(size=13,face="bold"), axis.text = element_text(size=11,face="bold"), axis.line = element_line(size=0.8), legend.title = element_text(face="bold"), legend.text=element_text(face="bold"))+
  xlab("ARG Count")+
  ylab("Number of Isolates")


# ----------------------------------------Suppl Figure 6A: snp-sites-----------------------------------
library(ggplot2)
library(RColorBrewer)

# Read in lineage barplot file
df_core_lineage <- read.csv('211202_snp-sites_0-8500.csv',
                            sep=",",
                            header = T)

sapply(df_core_lineage, class)
df_core_lineage$Range <- factor(df_core_lineage$Range, levels = c("0-50","51-100","101-150","151-200","201-250","251-300","301-350","351-400","401-450","451-500","501-550","551-600","601-650","651-700","701-750","751-800","801-850","851-900","901-950","951-1000","1001-1050","1051-1100","1101-1150","1151-1200","1201-1250","1251-1300","1301-1350","1351-1400","1401-1450","1451-1500","1501-1550","1551-1600","1601-1650","1651-1700","1701-1750","1751-1800","1801-1850","1851-1900","1901-1950","1951-2000","2001-2050","2051-2100","2101-2150","2151-2200","2201-2250","2251-2300","2301-2350","2351-2400","2401-2450","2451-2500","2501-2550","2551-2600","2601-2650","2651-2700","2701-2750","2751-2800","2801-2850","2851-2900","2901-2950","2951-3000","3001-3050","3051-3100","3101-3150","3151-3200","3201-3250","3251-3300","3301-3350","3351-3400","3401-3450","3451-3500","3501-3550","3551-3600","3601-3650","3651-3700","3701-3750","3751-3800","3801-3850","3851-3900","3901-3950","3951-4000","4001-4050","4051-4100","4101-4150","4151-4200","4201-4250","4251-4300","4301-4350","4351-4400","4401-4450","4451-4500","4501-4550","4551-4600","4601-4650","4651-4700","4701-4750","4751-4800","4801-4850","4851-4900","4901-4950","4951-5000","5001-5050","5051-5100","5101-5150","5151-5200","5201-5250","5251-5300","5301-5350","5351-5400","5401-5450","5451-5500","5501-5550","5551-5600","5601-5650","5651-5700","5701-5750","5751-5800","5801-5850","5851-5900","5901-5950","5951-6000","6001-6050","6051-6100","6101-6150","6151-6200","6201-6250","6251-6300","6301-6350","6351-6400","6401-6450","6451-6500","6501-6550","6551-6600","6601-6650","6651-6700","6701-6750","6751-6800","6801-6850","6851-6900","6901-6950","6951-7000","7001-7050","7051-7100","7101-7150","7151-7200","7201-7250","7251-7300","7301-7350","7351-7400","7401-7450","7451-7500","7501-7550","7551-7600","7601-7650","7651-7700","7701-7750","7751-7800","7801-7850","7851-7900","7901-7950","7951-8000","8001-8050","8051-8100","8101-8150","8151-8200","8201-8250","8251-8300","8301-8350","8351-8400","8401-8450","8451-8500"))

#df_core_lineage$Range <- factor(df_core_lineage$Range, levels = c("0-25","26-50","51-75","76-100","101-125","126-150","151-175","176-200","201-225","226-250","251-275","276-300","301-325","326-350","351-375","376-400","401-425","426-450","451-475","476-500","501-525","526-550","551-575","576-600","601-625","626-650","651-675","676-700","701-725","726-750","751-775","776-800","801-825","826-850","851-875","876-900","901-925","926-950","951-975","976-1000","1001-1025","1026-1050","1051-1075","1076-1100","1101-1125","1126-1150","1151-1175","1176-1200","1201-1225","1226-1250","1251-1275","1276-1300","1301-1325","1326-1350","1351-1375","1376-1400","1401-1425","1426-1450","1451-1475","1476-1500"))

df_core_lineage$Count <- as.numeric(df_core_lineage$Count)

# Make lineage barplot
ggplot(data=df_core_lineage, aes(x=Range,y=Count))+
  geom_bar(stat="identity", width=0.75, color="black")+
#scale_y_continuous(name="Number of pairwise comparisons within SNP range", limits=c(0,170), expand=c(0,0))+
  scale_y_continuous(name="Number of pairwise comparisons within SNP range", limits=c(0,7000), expand=c(0,0))+
  theme_classic()+
  theme(axis.title = element_text(face="bold"))+
#theme(text = element_text(),axis.text.x = element_text(size = 10, angle = 270, hjust=0, vjust = 0.5), plot.margin = margin(10, 30, 10, 10))+
  theme(text = element_text(),axis.text.x = element_text(size = 6, angle = 270, hjust=0, vjust = 0.5), plot.margin = margin(10, 30, 10, 10))+
    xlab("Core Genome SNPs (Range)")


# --------------------------------------Figure 3A: Lineage & Strain cluster coverage x ANI--------------------------------------

library(ggplot2)
library(RColorBrewer)
library(ggExtra)

# Read file
df_coverage_ani <- read.csv('221026_clonalcluster_coverage_ani.csv',
                            sep=",",
                            header = T)
df_coverage_ani$Status <- factor(df_coverage_ani$Status, levels=c("Same lineage","Same strain"))

# Plot Fig 3A
#gg_coverage<-ggplot(data=df_coverage_ani, aes(x=Coverage,y=ANI, fill = Status))+
gg_coverage<-ggplot(data=df_coverage_ani, aes(x=Coverage,y=ANI, color = Status))+
#geom_point(shape =21, colour="black")+
  geom_point(size =0.75)+
 #scale_y_continuous(limits=c(0.99949,1.000009))+
 #scale_x_continuous(limits=c(0.919,1.001))+
  geom_hline(yintercept = 0.99999,linetype="dashed", color = "black")+
  geom_vline(xintercept = 0.98, linetype = "dashed", color = "black")+
 #scale_fill_manual(values=c("#FFC20A","#0C7BDC"))+
  scale_color_manual(values=c("#333333","#0C7BDC"))+
  theme_bw()+
  theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size=12,face="bold"), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1.5) )

ggMarginal(gg_coverage, groupColour = TRUE, groupFill = TRUE, type = "histogram", xparams = list(binwidth = 0.0005), yparams = list(binwidth = 0.0000035) )
ggMarginal(gg_coverage, groupColour = TRUE, groupFill = TRUE)

# -----------------------------------------Figure 3B: Cluster Barplot-----------------------------------------
library(ggplot2)

# Read in cohorts barplot file
df_cluster <- read.csv('221026_clusters_barplot.csv',
                            sep=",",
                            header = T)

df_cluster$Isolates <- factor(df_cluster$Isolates, levels = c("2","3","4","5","6","7","8","9","10"))
#df_cluster$Cohort <- factor(df_cluster$Cohort, levels = c("Household","Diagnostic"))

# Plot Fig 3B
ggplot(data=df_cluster, aes(x=Isolates, y=Clusters, fill=Cohort)) +
  geom_bar(stat="identity",width=0.75, color="black", size=0.75)+
  theme_classic()+
 #theme_minimal()+
  scale_y_continuous(limits=c(0,30), expand=c(0,0))+
 #ylim(0,30)+
  scale_fill_manual(values=c("red","#76D7C4"))+
  theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size=12,face="bold"))+
 #xlab("Total Isolates in a Cluster")+
 #ylab("Number of Clusters")+
  theme(legend.position = c(0.75, 0.88), axis.title=element_blank(),  legend.title=element_text(face="bold"),legend.text=element_text(face="bold"))


# --------------------------------------Figure 3C: GWAS heatmap----------------------------------------------------------------
library(reshape2)
library(dendsort)
library(pheatmap)

#gene_pres.abs_df_sig<-read.csv('211205_gene_pres_abs_scoary_diag_col_collapsedclusters_BH_0isolates_removed_annotation.csv', head = T, row.names=1, check.names=F)
gene_pres.abs_df_sig<-read.csv('221026_gene_pres_abs_scoary_pgap_collapsed.csv', head = T, row.names=1, check.names=F)

#metadata_sig<-read.csv ("sig_metadata_493_pseud_redefined_only.csv",
metadata_sig<-read.csv ("221026_sig_metadata_scoary_pgap.csv",
                        header=T,
                        stringsAsFactors = F,
                        comment.char = "",
                        row.names=1,
                        quote = "")

#anno_colors_sig_simplified = list(
#Simplified_Cohort = c(Human_Infection="#CC0000",Animal_Infection="#33FFFF", Pet_colonizing="#0000FF",Human_colonizing="#FF00FF",Environment="#CCCCCC"))

anno_colors_sig_simplified = list(
  Redefined_Cohort = c(Diagnostic_Human="#C0392B",Diagnostic_Animal="#FF9388", Household="#76D7C4"),
  MLST = c(ST2521="#DBC9DF", ST2571="#AE76A3", ST862="#882E72", ST2575="#1965B0", ST2545="#5289C7", ST923="#7BAFDE", ST2550="#4EB265", ST181="#90C987", ST2553="#CEDFB0", ST2606="#F7F056", ST2544="#F6C141", ST2558="#F1932D", ST1207="#E8601C", ST2414="#DC050C", ST759="#777777", Other="#FFFFFF"))

#anno_colors_sig = list(
#Cohort = c(Human_Infection="#CC0000",Animal_Infection="#33FFFF", Pet_colonizing="#0000FF",Human_colonizing="#FF00FF",Environment_bedsheets="#CCCCCC",Environment_KitchenSinkHandle="#CCCCCC",Environment_TVremote="#CCCCCC",Environment_Telephone="#CCCCCC",Environment_KitchenTable="#CCCCCC",Environment_BathroomLightSwitch="#CCCCCC",Environment_BathroomCountertop="#CCCCCC",Environment_BathroomSinkHandle="#CCCCCC",Environment_FridgeDoorHandle="#CCCCCC",Environment_ToiletSeat="#CCCCCC",Environment_Bathtub="#CCCCCC",Environment_ComputerKeyboardMouse="#CCCCCC",Environment_ToiletHandle="#CCCCCC",Delphini_Environment="#FFFF00",Delphini_Pet_Colonizing="#FFFF00"))
#Location = c(Anterior_Nares_NICU="#660000",Blood_Adult_Child ="#CC9999", Blood_NICU="#FFFFFF"),
#MLST = c(ST1="red",ST5="yellow",ST8="#996600",ST15="#FF33FF",ST30="#006600",ST45="purple",ST59="#66CCCC",ST72="#FFCCFF",ST88="#FF9999",ST97="orange",ST398="#0000FF",Other="#CCCCCC"))

#set callback function for sorting dendrogram
callback = function(hc, ...) {
  dendsort(hc)
}

#drows=dist(gene_pres.abs_df_sig, method = 'binary')

pheatmap(gene_pres.abs_df_sig,
         cluster_rows = T,
         cluster_cols = T,
         #clustering_distance_cols=drows,
         clustering_method = 'complete',
         #clustering_callback = callback,
         cellheight = 4.5,
         cellwidth = 3.475,
         color=c('white','#666666'),
         annotation_colors = anno_colors_sig_simplified,
         annotation_col = metadata_sig,
         border_color = 'black',
         fontsize_col = 3.475,
         fontsize_row = 4.5,
         angle_col =270,
         legend = F)


# -----------------------------------------Figure 4A: Scoary CRISPR bubble plot-----------------------------------------
library(ggplot2)
library(scales)

# Read in cohorts barplot file
df_crispr_scoary <- read.csv('220523_scoary_crispr_bubbleplot.csv',
                       sep=",",
                       header = T)

df_crispr_scoary$x <- factor(df_crispr_scoary$x, levels = c("AI_0069","AI_0013","AI_0077","HI_0181","AI_0006","HI_0179","AI_0004","AI_0039","HI_0006","HI_0018","HI_0012","AI_0051","HI_0105","HI_0003","HI_0035","HI_0177","HI_0186"))
df_crispr_scoary$y <- factor(df_crispr_scoary$y, levels = c("HI_0186","HI_0177","HI_0035","HI_0003","HI_0105","AI_0051","HI_0012","HI_0018","HI_0006","AI_0039","AI_0004","HI_0179","AI_0006","HI_0181","AI_0077","AI_0013","AI_0069"))

# Plot Fig 4A
ggplot(df_crispr_scoary, aes(x = x, y = y,size = Shared,fill=Percent_Total))+
  geom_point(shape=21)+
  scale_radius(range = c(0,10),breaks=c(0,5,10,15,20),limits=c(0,25), name="Shared spacers")+
  scale_fill_gradient(low="white", high="#000066", name = "Percent total\nspacers shared",labels=percent)+
  coord_fixed((ratio=1))+
  geom_abline(slope=-1)+
  guides(fill = guide_colorbar(order = 0),size = guide_legend(order = 1))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"))

# ----------------------------------------Figure 4B: Total Cohort - CRISPR bubble plot----------------------------------------
library(reshape2)

# Read in ANI matrix
df_spacer_all_RC <- read.delim('01_all_spacers_clusters_noduplicates_TAB_SEPARATED_distancematrix.txt', header = T)

# Convert matrix to long format
df_spacer_all_RC_long <- melt(df_spacer_all_RC)

# Add column names
colnames(df_spacer_all_RC_long) <- c("spacer1","spacer2","ANI")

# Order rows by ANI
df_spacer_all_RC_long<- df_spacer_all_RC_long[order(df_spacer_all_RC_long$ANI,decreasing = TRUE), ]

# Delete all self-comparisons
df_spacer_all_RC_long<- df_spacer_all_RC_long[!(df_spacer_all_RC_long$spacer1==df_spacer_all_RC_long$spacer2),]

# Make new dataframe with only comparisons of ANI=100%
df_spacer_all_RC_long_ani100<- df_spacer_all_RC_long[(df_spacer_all_RC_long$ANI==100),]

# Remove colnames
names(df_spacer_all_RC_long_ani100)<-NULL

# Write .txt file and remove in Python (y,x) rows when (x,y) already present
write.table(df_spacer_all_RC_long_ani100,"02_220807_df_spacer_all_long_ani100.txt",sep = "\t",row.names = FALSE)



# READ IN OUTPUT FROM s16.9.5_distancematrix_remove_ab_ba.sh
df_spacer_all_ani100_filtered <- read.delim('03_220807_df_spacer_all_long_ani100_filtered.txt', header=F)

# Remove "_RV" from all identifiers
df_spacer_all_ani100_filtered$V1<-gsub("_RV","",as.character(df_spacer_all_ani100_filtered$V1))
df_spacer_all_ani100_filtered$V2<-gsub("_RV","",as.character(df_spacer_all_ani100_filtered$V2))

# Remove duplicate rows
df_spacer_all_ani100_filtered <- df_spacer_all_ani100_filtered[!duplicated(df_spacer_all_ani100_filtered), ]
# Remove colnames
names(df_spacer_all_ani100_filtered)<-NULL
# Write .txt file and run s16.9.6_distancematrix_remove_allconnections.sh
write.table(df_spacer_all_ani100_filtered,"04_220807_df_spacer_all_long_ani100_filtered_cleaned.txt",sep = "\t",row.names = FALSE)



# READ IN CO-OCCURRENCE MATRIX
df_spacer_all_ani100_cooccurrence <- read.delim('06_220807_spacer_all_ani100_cooccurrence_matrix.txt', header = T)

# Convert matrix to long format
df_spacer_all_ani100_cooccurrence <- melt(df_spacer_all_ani100_cooccurrence)

# Add column names
colnames(df_spacer_all_ani100_cooccurrence) <- c("Isolate1","Isolate2","Shared_spacer_ct")

# Write .txt file and run s16.9.6_distancematrix_remove_allconnections.sh
write.table(df_spacer_all_ani100_cooccurrence,"06_220807_spacer_all_ani100_cooccurrence_matrix_long.txt",sep = "\t",row.names = FALSE)



# Read in BUBBLE PLOT CSV

df_crispr_scoary <- read.csv('07_bubbleplot_input_operonspecific.csv',
                             sep=",",
                             header = T)

library(ggplot2)
library(scales)
#library(ggbreak)

#df_crispr_scoary$x <- factor(df_crispr_scoary$x, levels = c("AI_0069","AI_0013","AI_0077","HI_0181","AI_0006","HI_0179","AI_0004","AI_0039","HI_0006","HI_0018","HI_0012","AI_0051","HI_0105","HI_0003","HI_0035","HI_0177","HI_0186"))
#df_crispr_scoary$y <- factor(df_crispr_scoary$y, levels = c("HI_0186","HI_0177","HI_0035","HI_0003","HI_0105","AI_0051","HI_0012","HI_0018","HI_0006","AI_0039","AI_0004","HI_0179","AI_0006","HI_0181","AI_0077","AI_0013","AI_0069"))

a <- ifelse(grepl('^A', df_crispr_scoary$Isolate1), "#FF9388", ifelse(grepl('^b', df_crispr_scoary$Isolate1), "black", ifelse(grepl('^H', df_crispr_scoary$Isolate1), "#C0392B", ifelse(grepl('^i', df_crispr_scoary$Isolate1), "black", ifelse(grepl('^m', df_crispr_scoary$Isolate1), "#1B4F72", ifelse(grepl('^o', df_crispr_scoary$Isolate1), "black", ifelse(grepl('^S', df_crispr_scoary$Isolate1), "#76D7C4", "black")))))))
c <- ifelse(grepl('^S', df_crispr_scoary$Isolate1), (ifelse(grepl('^A', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^H', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^m', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^S', df_crispr_scoary$Isolate2),"#0E6251", "black")))))))),"white")
d <- ifelse(grepl('^m', df_crispr_scoary$Isolate1), (ifelse(grepl('^A', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^H', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^m', df_crispr_scoary$Isolate2),"#1B4F72",(ifelse(grepl('^S', df_crispr_scoary$Isolate2),"black", "black")))))))),c)
e <- ifelse(grepl('^H', df_crispr_scoary$Isolate1), (ifelse(grepl('^A', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^H', df_crispr_scoary$Isolate2),"#C0392B",(ifelse(grepl('^m', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^S', df_crispr_scoary$Isolate2),"black", "black")))))))),d)
f <- ifelse(grepl('^A', df_crispr_scoary$Isolate1), (ifelse(grepl('^A', df_crispr_scoary$Isolate2),"#CC3333",(ifelse(grepl('^H', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^m', df_crispr_scoary$Isolate2),"black",(ifelse(grepl('^S', df_crispr_scoary$Isolate2),"black", "black")))))))),e)

# Plot Fig 4B
ggplot(df_crispr_scoary, aes(x = Isolate1, y = Isolate2, size = Shared_spacer_ct_20cap, fill=Percent_total))+
  geom_point(shape=21,color=f)+
#geom_point(shape=21)+
#scale_radius(range = c(0,5),breaks=c(0,7,14,21,28,35),limits=c(0,36), name="Shared spacers")+
#scale_radius(range = c(0,4),breaks=c(0,1,2,3,4,5),limits=c(0,6), name="Shared spacers (log2)")+
  scale_radius(range = c(0,4),breaks=c(0,5,10,15,20),limits=c(0,20), name="Shared spacers")+
  scale_fill_gradient(low="white", high="#0000CC", name = "Percent total\nspacers shared",labels=percent)+
  coord_fixed((ratio=1))+
  geom_abline(slope=-1)+
  guides(fill = guide_colorbar(order = 0),size = guide_legend(order = 1))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5,size=5,color=a), axis.text.y=element_text(size=5,color=a), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"))


# Reading in files to compare shared spacer cts and %ages with whole genome ANI and Hadamard scores
# melt ANI matrix
ani_matrix=read.csv("ani_matrix.csv", header=T, sep=",")
ani_sig_long<-melt(ani_matrix)

# Add column names
colnames(ani_sig_long) <- c("isolate1","isolate2","ANI")

# Delete all self-comparisons
ani_sig_long<- ani_sig_long[!(ani_sig_long$isolate1==ani_sig_long$isolate2),]
write.csv(ani_sig_long,"sig_ani_long.csv")

p_aln_matrix=read.csv("p_aln_matrix.csv", header=T, sep=",")
aln_sig_long<-melt(p_aln_matrix)

# Add column names
colnames(aln_sig_long) <- c("isolate1","isolate2","aln")

# Delete all self-comparisons
aln_sig_long<- aln_sig_long[!(aln_sig_long$isolate1==aln_sig_long$isolate2),]
write.csv(aln_sig_long,"sig_aln_long.csv")



# Stats testing
# Read csv for shared spacer cts for isolates within cohorts and isolates between cohorts
df_AI_intra_inter_shared_spacer=read.csv("220813_AI_intra-inter_shared_spacer_comparison.csv", header=T, sep=",")
df_AI_intra_inter_spacer_test<-wilcox.test(Shared_spacers~Comparison, data = df_AI_intra_inter_shared_spacer,paired=FALSE, exact=FALSE, conf.int=TRUE,alternative = c("greater"))
df_AI_intra_inter_spacer_test

df_HI_intra_inter_shared_spacer=read.csv("220813_HI_intra-inter_shared_spacer_comparison.csv", header=T, sep=",")
df_HI_intra_inter_spacer_test<-wilcox.test(Shared_spacers~Comparison, data = df_HI_intra_inter_shared_spacer,paired=FALSE, exact=FALSE, conf.int=TRUE,alternative = c("greater"))
df_HI_intra_inter_spacer_test

df_SF_intra_inter_shared_spacer=read.csv("220813_SF_intra-inter_shared_spacer_comparison.csv", header=T, sep=",")
df_SF_intra_inter_spacer_test<-wilcox.test(Shared_spacers~Comparison, data = df_SF_intra_inter_shared_spacer,paired=FALSE, exact=FALSE, conf.int=TRUE,alternative = c("less"))
df_SF_intra_inter_spacer_test

intra_inter_shared_spacer_pvalues<-c(0.0333,0.224,0.0004)

p.adjust(intra_inter_shared_spacer_pvalues, method = "fdr")


# ----------------------------------------Figure 4C: Shared spacer ct x ANI----------------------------------------

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Read in lineage barplot file
df_spacer_ani <- read.csv('220830_ANI_vs_shared_spacer_ct.csv',
                            sep=",",
                            header = T)

# Make lineage barplot
ggplot(data=df_spacer_ani, aes(x=ANI,y=shared_spacer_ct))+
  geom_jitter(width=0.0001,height=.4, size=0.1)+
  scale_y_continuous(limits=c(-2,36),expand = c(0,0))+
 scale_x_continuous(limits=c(0.9895,1.0001),expand = c(0,0))+
  theme_bw()+
  geom_smooth(method = "lm",linetype="dashed", size=0.5)+
  stat_regline_equation(label.x = 0.99, label.y = 33, aes(label = ..rr.label..))+
  ylab("Shared spacers (count)")+
  theme(legend.title=element_blank(), axis.title=element_text(face="bold"), axis.text=element_text(face="bold"))

# -------------Figure 5B, Suppl Figure 7: COG Genes - Non-synonymous vs Synonymous mts Fishers's exact test (BLAST)-------------


# Read in COG data file - NEW (221113)
df_cog <- read.csv('221113_COG_fishersexact_blast.csv',
                   sep=",",
                   header = T, row.names = 1)

# Stats testing for Figure 7

# Fisher's exact test of independence
#p-value = 0.5817
fisher.test(df_cog, simulate.p.value = T)


library(EMT)
# ALL SNS>1 genes vs. SF_0321: Multi COG (ie AB and ABC) categories and "-" category removed
#p.value=1
multinomial.test(c(5,7,5,9,4,6,9,11,1,7,2,1,10,3,26,3,9,7,0),c(0.0176376,0.0464992,0.0240513,0.0481026,0.0865847,0.0887226,0.0609300,0.0571887,0.0096205,0.0710850,0.0192410,0.0138963,0.0603955,0.0171032,0.2271513,0.0277926,0.0732229,0.0502405,0.0005345), useChisq = FALSE, MonteCarlo =TRUE)

# Non-syn genes vs. SF_0321: Multi COG (ie AB and ABC) categories and "-" category removed
#p.value=1
multinomial.test(c(5,3,5,8,4,4,5,8,0,4,2,1,7,2,19,2,7,5,0),c(0.0176376269375,0.0464991982897,0.0240513094602,0.0481026189204,0.0865847140567,0.0887226082309,0.0609299839658,0.0571886691609,0.0096205237841,0.0710849812934,0.0192410475681,0.0138963121325,0.0603955104222,0.0171031533939,0.2271512560128,0.0277926242651,0.0732228754677,0.0502405130946,0.0005344735436), useChisq = FALSE, MonteCarlo =TRUE)

# Pick out specific COG categories to test
# Predict that Defense Mechanisms (V) will undergo more adaptive selection in non-synonymous mt genes.

# Stats testing for Fig 5B
# NSS COG: Defense Mechanisms (V): two-sided = greater = p-value = 0.02278
binom.test(5, 91, 0.0176376269374666, alternative = "greater")

# SSO COG: Defense Mechanisms (V): two-sided = p-value = 1
binom.test(0, 34, 0.0176376269374666, alternative = "two.sided")

# -----------------------------------------Figure 5C: 221113 Permutation test-----------------------------------------

library(dplyr)

# Create empty vector
number_of_genes_mutated_more_than_once <- c()
#seed 27 = 951 (n=1000)
set.seed(27)

# Simulate over 1,000 iterations
for (i in 1:10000) {

# 1,871 genes in the reference isolate assembly that are single COG letter and do not include "partial" in annotation. For each cluster, sample the number of genes that have undergone at least 1 NON-SYNONYMOUS mutation
  cluster01 <- sample(1871, 5, replace = FALSE, prob = NULL)
  cluster02 <- sample(1871, 9, replace = FALSE, prob = NULL)
  cluster03 <- sample(1871, 19, replace = FALSE, prob = NULL)
  cluster04 <- sample(1871, 4, replace = FALSE, prob = NULL)
  cluster05 <- sample(1871, 16, replace = FALSE, prob = NULL)
  cluster11 <- sample(1871, 17, replace = FALSE, prob = NULL)
  cluster13 <- sample(1871, 4, replace = FALSE, prob = NULL)
  cluster14 <- sample(1871, 10, replace = FALSE, prob = NULL)
  cluster28 <- sample(1871, 4, replace = FALSE, prob = NULL)
  cluster34 <- sample(1871, 3, replace = FALSE, prob = NULL)
  
  # Add row names to each cluster for downstream binding
  names(cluster01) <- 1:5
  names(cluster02) <- 1:9
  names(cluster03) <- 1:19
  names(cluster04) <- 1:4
  names(cluster05) <- 1:16
  names(cluster11) <- 1:17
  names(cluster13) <- 1:4
  names(cluster14) <- 1:10
  names(cluster28) <- 1:4
  names(cluster34) <- 1:3
  
  # Bind all clusters into one dataframe
  df_combinations <- as.data.frame(bind_rows(cluster01,cluster02,cluster03,cluster04,cluster05,cluster11,cluster13,cluster14,cluster28,cluster34))
  
  # Combine all columns into one
  df_combinations <- data.frame(x=unlist(df_combinations))
  
  # Sort the column in ascending order and remove NAs
  df_combinations <- sort(df_combinations$x, decreasing = FALSE, na.last = NA)
  
  # Determine frequency of each randomly chosen gene
  df_combinations_freq <- data.frame(table(df_combinations))
  
  # Count the number of genes that are randomly chosen more than once
  number_of_genes_mutated_more_than_once[i] <- length(df_combinations_freq[df_combinations_freq$Freq > 1,][,2])
}

# Sort the list in ascending order
number_of_genes_mutated_more_than_once <- sort(number_of_genes_mutated_more_than_once, decreasing = FALSE)

number_of_genes_mutated_more_than_once[9500]

write.csv(number_of_genes_mutated_more_than_once,"221113_01_mutated_genes_multiple_clusters.csv")

library(ggplot2)

# Read in cohorts barplot file
df_genes_mutated <- read.csv('221113_02_mutated_genes_multiple_clusters_barplot.csv',
                             sep=",",
                             header = T)

df_genes_mutated$Genes <- factor(df_genes_mutated$Genes, levels = c("0","1","2", "3","4","5","6","7","8"))

# plot Fig 5C
ggplot(data=df_genes_mutated, aes(x=Genes, y=Count, fill=Group)) +
  geom_bar(stat="identity", width=0.75, color="black", size=0.75)+
  theme_classic()+
  scale_y_continuous(limits = c(0, 3050),expand = c(0,0))+
  geom_vline(xintercept = 5.16,linetype="longdash", color = "red")+
  scale_fill_manual(values=c("gray","black"))+
  theme(axis.title = element_text(size=12,face="bold"), axis.text = element_text(size=12,face="bold"),legend.position="none")+
  xlab("Number of Genes Accruing Non-\nSynonymous SNSes in Multiple Clusters")+
  ylab("Count (10,000 simulations)")
