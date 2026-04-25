library(ggplot2)
library(stringr)

wdir = "/lustre/scratch123/hgi/teams/hgi/re3/cardinal_plot_discordance/"
datadir = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc_with_GT/Cardinal_45538_Aug_11_2022/results/concordances/'
#datadir = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc_with_GT/ELGH_9th_May_2022/results/concordances/'

panels = c('CRD_CMB13076375', 'CRD_CMB13076377','CRD_CMB13076379','CRD_CMB13076381','CRD_CMB13076383',
'CRD_CMB13076385', 'CRD_CMB13076376'  ,'CRD_CMB13076378'  ,'CRD_CMB13076380','CRD_CMB13076382', 
'CRD_CMB13076384')


#panels = c('CRD_CMB12813646','CRD_CMB12813647','CRD_CMB12813648','CRD_CMB12813649','CRD_CMB12813650',
#'CRD_CMB12813651','CRD_CMB12813652','CRD_CMB12813653','CRD_CMB12813654','CRD_CMB12813655')

mega_df = data.frame()

for (panel in panels){
  print(panel)
  concordance_file = paste(datadir, panel, "/cell_concordance_table.tsv", sep = "")
  df_conc = read.table(concordance_file, header = T, fill = T, sep = "\t")
  swap_file = paste(datadir, panel, '/', panel, '_subsampling_donor_swap_quantification.tsv', sep = "")
  df_swap = read.table(swap_file, header = T, fill = T, sep = "\t")

  panel = rep(panel, length(df_conc$X))
  df_conc$panel = panel
  
  swapcounts = vector(length=length(df_conc$X))
  newdonors = vector(length=length(df_conc$X))
  
  for (i in 1:length(df_conc$X)){
    nswaps = df_swap$Nr.times.becoming.different.donor.in.subsampling[which(df_swap$cell == df_conc$GT.1[i])]
    swapcounts[i] = nswaps
    if (nswaps > 0){
      newdons = df_swap$New.Donor.Identities[which(df_swap$cell == df_conc$GT.1[i])]
      newdonors[i] = newdons
    } else {
      newdonors[i] = '-'
    }
  }
  df_conc$swap_counts = swapcounts
  df_conc$newdonorids = newdonors
  
  mega_df = rbind(mega_df, df_conc)
  
}

mega_df$Nswaps = as.factor(mega_df$swap_counts)

pdf(file=paste(wdir, 'pc_discordant_n_swaps_all_Cardinal_45538_Aug_11_2022.pdf', sep = ''))
p1 = ggplot(mega_df, aes(x=Nswaps, y=Percent_strict_discordant, fill=Nswaps)) + geom_violin() + theme(legend.position = "none") + xlab("Number of swaps") + ylab("Percentage discordant (strict)") + ggtitle("Percentage strict discordance vs number of times cell switches donors", subtitle="Cardinal_45538_Aug_11_2022, all cells")
print(p1)
dev.off()

#remove the rows where there are >0 swaps but the newdonorids does NOT contain '_' as those not containing '_' are not confident assignments
df_no_swaps = mega_df[which(mega_df$swap_counts == 0),]
df_no_unknowns = mega_df[which(str_detect(mega_df$newdonorids, '_') == TRUE ),]
mega_df_filtered = rbind(df_no_swaps, df_no_unknowns)

pdf(file=paste(wdir, 'pc_discordant_n_swaps_filtered_Cardinal_45538_Aug_11_2022.pdf', sep = ''))
p2 = ggplot(mega_df_filtered, aes(x=Nswaps, y=Percent_strict_discordant, fill=Nswaps)) + geom_violin() + theme(legend.position = "none") + xlab("Number of swaps") + ylab("Percentage discordant (strict)") + ggtitle("Percentage strict discordance vs number of times cell switches donors", subtitle="Cardinal_45538_Aug_11_2022,excluding cells switching only to unknown donors")
print(p2)
dev.off()

#number of cells in each category
counts = c(1:10)
unfiltered_counts = vector(length=10)
filtered_counts = vector(length=10)

for (c in counts){
  n = length(which(mega_df$swap_counts == c))
  m = length(which(mega_df_filtered$swap_counts == c))
  unfiltered_counts[c] = n
  filtered_counts[c] = m
}

#pdf(file=paste(wdir, 'counts_unfiltered_ELGH_9th_May_2022.pdf', sep = ''))
#barplot(unfiltered_counts, names = counts, col = 'red', xlab = 'N times becoming different donor', ylab = 'N cells', main = 'Unfiltered data, ELGH_9th_May_2022')
#dev.off()

#pdf(file=paste(wdir, 'counts_filtered_ELGH_9th_May_2022.pdf', sep = ''))
#barplot(filtered_counts, names = counts, col = 'red', xlab = 'N times becoming different donor', ylab = 'N cells', main = 'Filtered to cells switching to known donors, ELGH_9th_May_2022')
#dev.off()

pdf(file=paste(wdir, 'counts_unfiltered_Cardinal_45538_Aug_11_2022.pdf', sep = ''))
barplot(unfiltered_counts, names = counts, col = 'red', xlab = 'N times becoming different donor', ylab = 'N cells', main = 'Unfiltered data, Cardinal_45538_Aug_11_2022')
dev.off()

pdf(file=paste(wdir, 'counts_filtered_Cardinal_45538_Aug_11_2022', sep = ''))
barplot(filtered_counts, names = counts, col = 'red', xlab = 'N times becoming different donor', ylab = 'N cells', main = 'Filtered to cells switching to known donors, Cardinal_45538_Aug_11_2022.pdf')
dev.off()


