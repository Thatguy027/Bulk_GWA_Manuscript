library(xQTLStats)
library(qs)
library(tidyverse)
library(vcfR)
library(AlphaSimR)
library(Rfast)
library(ggpubr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

#pull from github
xQTLSims.dir = '/Users/Stefan/github_repos/xQTLSims/'

source.dir=paste0(xQTLSims.dir, 'R/')
data.dir=paste0(xQTLSims.dir, 'data/')
project.dir=paste0(xQTLSims.dir, 'projects/JU1793_repeat/')

#function to simulate crosses, treatement of X is incomplete/broken
#and additional helper functions
source(paste0(source.dir, 'simWormCrosses.R'))
source(paste0(source.dir, 'helperFxs.R'))
source(paste0(source.dir, 'makeCountTables.R'))

# modified to deal with old data

makeCountTables_old <- function(sample.key, vcf, gt) {
  countdfs=list() 
  for(i in 1:nrow(sample.key) ) { 
    #nrow(sample.key)) {
    p1=sample.key$'parent 1'[i]
    p2=sample.key$'parent 2'[i]
    
    #in this case, mate herm N2 with male CB
    p.names=c(p1, p2) #B4856') #XZ1516')
    mating.matrix=matrix(c(1,2), nrow=1)
    # in this case, given one-way cross, potentially allow N2 herm to self 
    founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F) #c('N2', 'XZ1516'))
    
    genMap=getGenMap(founderPop)
    
    
    sn=sample.key$'sample name'[i]
    scounts=read_tsv(paste0(sample.dir, sn, sample.suffix)) #'.table'))
    scounts.sub=scounts[paste0(scounts$contig, '_', scounts$position) %in% genMap$id,]
    scounts=data.frame(id=paste0(scounts.sub$contig, '_', scounts.sub$position),ref=scounts.sub$refCount, alt=scounts.sub$altCount)
    scounts=left_join(genMap, scounts, by='id') %>%
      dplyr::distinct(id, .keep_all=T)
    #come on GATK aseReadCounter, wtf is up with different length output given different BAM input
    #for know fill out with 0s 
    scounts$ref[is.na(scounts$ref)]=0
    scounts$alt[is.na(scounts$alt)]=0
    names(scounts)[1]='ID'
    
    #phase it 
    countdf=phaseBiparental(scounts, p.names[1], founderPop, genMap)
    
    #note, we need a better structure for keeping track of which parent is which 
    attr(countdf, 'p1')=p.names[1]
    attr(countdf, 'p2')=p.names[2]
    
    snn=paste(sample.key[i,], collapse='_')
    countdfs[[snn]]=countdf
  }
  
  return(countdfs)
}

#unique chromosomes 
uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')
gmap=restructureGeneticMap(gmap.file)

#pretty intensive memory usage here
#include web link to vcf file 
#elegans.isotypes.vcf=paste0(data.dir,'WI.20220216.impute.isotype.vcf.gz')

#filtered vcf as qsave objects
# !!!! find the premade objects here folks !!!! :
# /u/project/kruglyak/jsbloom/xQTL/elegans/ref/
#and place in your data.dir
elegans.isotypes.vcf.qs=paste0(data.dir,'WI.20220216.vcf.qs')

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs=paste0(data.dir,'WI.20220216.vcf.GT.qs')

#run once, Laura skip this ============================================================
#use premade objects to save yourself the memory related headache, but this is how the objects are made
#preprocessVCF(elegans.isotypes.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
#======================================================================================

vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)

#some samples had low depth from "old" stuff below
#this script combined allele counts from multiple sequencing runs
# ~/UCLA/Projects/bulkGWAS/2024_nic_crosses2/scripts/20240430_process_run1_run2.R

sample.key.file=paste0(project.dir, '2040815_JUmig6_Sample_key_xQTL.tsv')


sample.key=read_tsv(sample.key.file) %>%
  dplyr::select(sample:condition) %>%
  dplyr::rename(`sample name` = sample)

rel_samples <- sample.key%>%
  dplyr::filter(`sample name`%in%c("S31","S32","S33","S34"))

sample.dir='aser/'
sample.suffix=".table"
countdfs=makeCountTables(rel_samples, vcf,gt)

# for(cts in 1:length(countdfs)){
#   countdfs_pr[[cts]] <- countdfs_pr[[cts]] 
# }




#calculate the allele frequencies and SEs, params here need some more dialing in 
#should parallelize this .... yawn 
afds=lapply(names(countdfs), function(snn) {
  calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.95, bin.width=10000, eff.length=2000, uchr=uchr)
})
names(afds)=names(countdfs)

tosave <- list(afds,countdfs)
base::save(tosave, file ="af_counts.rda")

plots=lapply(names(afds), function(snn) {
  plotIndividualExperiment(afds[[snn]], snn) 
})
names(plots)=names(countdfs)
#e.g.  to visualize just call
#plots[[1]]

#dump individual plots somewhere ---------------
plot.dir='plots/'
for(snn in names(plots)) {
   ggsave(paste0( plot.dir, snn, '_10000_cuts.png'), plots[[snn]], width=16)
}

lapply(afds, names)

# pruned
comp_df <- data.frame(contrast1=c(2,2,3),
                      contrast2=c(3,4,4))

#I have to restructure this stuff  ---------------------------------------------

for(comp_plots in 1:nrow(comp_df)){
  ind1=comp_df$contrast1[comp_plots]
  ind2=comp_df$contrast2[comp_plots]
  cont1=afds[[ind1]]
  cont2=afds[[ind2]]
  results=calcContrastStats(results=list(cont1, cont2),
                            L=paste0('_', names(afds)[ind1]), R=paste0('_', names(afds)[ind2]) ) #'_high1', R='_unsel1')
  
  interval_df <- results %>%
    dplyr::group_by(chrom)%>%
    # dplyr::filter(chrom=="II", physical.position < 5e6) %>%
    dplyr::filter(LOD > max(LOD)-2) %>%
    dplyr::mutate(lcon = min(physical.position),
                  rcon = max(physical.position)) %>%
    dplyr::arrange(desc(LOD))%>%
    dplyr::distinct(chrom, lcon, rcon, .keep_all = T) %>%
    dplyr::mutate(marker = paste0(chrom,":",lcon,"-",rcon)) %>%
    dplyr::select(marker, physical.position, LOD)
  
  suffix.1=gsub("ref_","",colnames(results)[6])
  suffix.2=gsub("ref_","",colnames(results)[17])
  
  sc=plotContrast(results, suffix1=suffix.1, suffix2=suffix.2) #S3_N2_XZ1516_7_KO')

  s=plotSummary(results, effective.n.tests=2000)
  
  pare1=str_split(suffix.1, pattern = "_")[[1]][2]
  pare2=str_split(suffix.2, pattern = "_")[[1]][3]
  cond1=str_split(suffix.1, pattern = "_")[[1]][5]
  cond2=str_split(suffix.2, pattern = "_")[[1]][5]
  gen1=str_split(suffix.1, pattern = "_")[[1]][4]
  gen2=str_split(suffix.2, pattern = "_")[[1]][4]
  
  write.table(interval_df, file = glue::glue('{plot.dir}{pare1}_{pare2}_F{gen1}-{gen2}_contrast_{cond1}-{cond2}_10000.tsv'), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(results, file = glue::glue('{plot.dir}{pare1}_{pare2}_F{gen1}-{gen2}_contrast_{cond1}-{cond2}_10000_plot_DF.tsv'), col.names = T, row.names = F, quote = F, sep = "\t")
  
  
  
  a=ggarrange(plots[[ind1]],plots[[ind2]], s,  nrow=3)
  ggsave(filename = glue::glue('{plot.dir}{pare1}_{pare2}_F{gen1}-{gen2}_contrast_{cond1}-{cond2}_10000.png'),
         plot = a, height = 16, width = 16)
}


