library(tidyverse)
library(vcfR)
library(memuse)

#Ania Lorenc al16@sanger.ac.uk

#Aim: use genotypes called in samples containing correctly (assumption) identified spike-in to build a consensus genotype to guide vireo of spike-in in other samples
#Input: Vireo-produced vcf files with donor genotypes
#Output: A vcf file to use as a prior for Vireo runs of other samples
#The whole workflow is a bit simplistic and driven by by manual exploration of the genotypes rather than optimal way of doing things

#Read in all vcfs of interest (:good samples") produced by vireo, extract gt (genotypes part only)
all_vcfs_gt <- list()

#adjust for file structure, source of sample names (here cellranger_samples_to_process df)
for( i in 1:length(cellranger_samples_to_process$V1)){
  print(cellranger_samples_to_process[i,])
  filepath=file.path("../out_lite_whole_gnomad", file.path(cellranger_samples_to_process$V1[i],"vireoOutput/GT_donors.vireo.vcf.gz"))
  all_vcfs_gt[[cellranger_samples_to_process$V3[i]]] <- read.vcfR(filepath)%>%extract.gt%>%as_tibble(rownames ="loc")
}

all_vcfs_gt <- all_vcfs_gt%>%
  bind_rows(.id = "sample")

write_tsv(all_vcfs_gt, file="all_vcfs_spiked_samples.csv")


##### Manual checks 
#Summary of identical versus different genotypes: I check for samples with excess of identical sites across genotypes - they are the ones cellsnp+vireo did not work
all_vcfs_gt%>%group_by(sample, id=donor0==donor1)%>%
  summarise(N=n())%>%
  group_by(sample)%>%
  mutate(freq=round(N/sum(N),2))%>%
  ggplot(aes(x=sample, fill=id, y=freq))+ geom_col()+theme_bw()+coord_flip()


#taking into account only the loci which are different. 190 loci consistently different between genotypes in every call
#excluding the chimeric sample FOXP3_G1_BM, as I forced 2 genotypes and there are 3
all_vcfs_gt %>%
  filter(donor0!=donor1, !sample%in%c("FOXP3_G1_BM"))%>%
  group_by(loc)%>%
  summarise(N=n())%>%
  arrange(desc(N))%>%
  .$N%>%table()

#####

#keeping  sites different between spikein and main sample in at least 7 samples only
consistency_cutoff <- 7
samples_not_to_include <-c("FOXP3_G1_BM")
diff_loci <- all_vcfs_gt %>%
  filter(!sample%in%samples_not_to_include)%>%
  filter(donor0!=donor1) %>% #, donor0!="1/0", donor1!=("1/0"))%>%
  group_by(loc)%>%
  summarise(N=n())%>%filter(N>=consistency_cutoff)%>%.$loc


##### Manual check - visualise agreement between genotypes/clsutering of spike-in genotypes

genotypes_matrix <- all_vcfs_gt %>%filter(loc%in%diff_loci)%>%
  filter(!sample%in%samples_not_to_include)%>%
  pivot_longer(names_to="donor",cols=c(donor0, donor1))%>%
  mutate(value=ifelse(value=="0/0",0, ifelse(value=="1/1",2,1)))%>%
  pivot_wider(id_cols=loc,names_from = c(sample,donor) )


as.matrix(genotypes_matrix[,-1])%>%
  t()%>%heatmap(scale = "none")


#####
#####Preparing artificial vcf from genotypes identified as spike-in in spiked samples


consistency_cutoff <- 7
samples_not_to_include <- c("FOXP3_G1_BM","FOXP3_C1_BM")
diff_loci <- all_vcfs_gt %>%
  filter(!sample%in%samples_not_to_include)%>%
  filter(donor0!=donor1) %>% #, donor0!="1/0", donor1!=("1/0"))%>%
  group_by(loc)%>%
  summarise(N=n())%>%filter(N>=consistency_cutoff)%>%.$loc

#now get sites in agreement. First I manually identify which donor is a spike in (could do it with the 'smaller' donor to automate)
batch_donors=c("CTLA4_A1_T_donor1",
               "FOXP3_A1_BM_donor0",
               "CTRL_B1_T_donor1",
               "CTRL_C1_T_donor0",
               "FOXP3_D1_BM_donor0",
               "CTRL_D1_T_donor0",
               "CTLA4_C1_T_donor0",
               "LRBA_F1_T_donor1" )

#subset GT by spike-ins genotypes only, different sites
potential_sites <- all_vcfs_gt %>%
  filter(!sample%in%samples_not_to_include, loc%in%diff_loci)%>%
  pivot_longer(cols=c(donor0, donor1), names_to="donor")%>%
  pivot_wider(names_from=c("sample", "donor"), id_cols="loc",  values_from = "value")%>%
  select(loc, batch_donors)

#identify sites agreeing across spikeins
agreement_sites <- potential_sites%>%
  pivot_longer(cols = 2:9)%>% #!!!correct with column names, as here depends on the number of donors
  group_by(loc, value)%>%
  summarise(N=n())%>%
  pivot_wider(id_cols = loc,names_from = value, values_from=N)%>%
  filter(`0/0`>=consistency_cutoff |`1/0`>=consistency_cutoff |`1/1`>=consistency_cutoff)%>%.$loc


##read in all vcfs FULL 
all_vcfs_raw <- list()

#adjust for file structure, source of sample names (here cellranger_samples_to_process df)
for( i in 1:length(cellranger_samples_to_process$V1)){
  print(cellranger_samples_to_process[i,])
  filepath=file.path("../out_lite_whole_gnomad", file.path(cellranger_samples_to_process$V1[i],"vireoOutput/GT_donors.vireo.vcf.gz"))
  all_vcfs_raw[[cellranger_samples_to_process$V3[i]]] <- read.vcfR(filepath)
}


#I will use one of the real vcfs to subset to the agreement sites &and its spikein sample. Should check whether this one is indeed in full agreement with consensus or overwrite the data so it is
sample_to_use <- setdiff(names(all_vcfs_raw), samples_not_to_include)[[1]]  # chosen randomly from good samples; here batch effect=donor1
one_sample <- all_vcfs_raw[[sample_to_use]]

which_sites <- which((one_sample%>%extract.gt()%>%rownames() )%in%agreement_sites)
one_sample_subset <- one_sample[which_sites, samples="donor1"]

write.vcf(one_sample_subset, file = "one_sample_subset_CTLA4_A1_T_770.vcf.gz")

