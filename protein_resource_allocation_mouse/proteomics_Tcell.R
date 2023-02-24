library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(tidyr)
library(openxlsx)
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file = 'proteomics_Tcell.RData')

load('assign_group_mouse_files.RData')

order.group = names(color.group)

plot.coarse = T
color.group.coarse = c('other'='#D9D9D9','transcription'="#FED9A6",'translation'='#8DD3C7','metabolism'= '#BEBADA','anabolic'='#FCCDE5',
                       'other mitochondria'='#C6DBEF','TCA'='#6BAED6','ox phos'='#4292C6','glycolysis'='#FFED6F')
order.group.coarse = names(color.group.coarse)
anno.coarse <- data.frame(group = order.group, n = 1:length(order.group)) %>% 
  mutate(group.coarse = case_when(n %in% 1:4 ~ 'other',
                                  n %in% 5:6 ~ 'transcription',
                                  n %in% 7:8 ~ 'translation',
                                  n %in% 9:14 ~ 'anabolic',
                                  n %in% 15:16 ~ 'other mitochondria',
                                  n %in% 17:19 ~ 'TCA',
                                  TRUE ~ as.character(group))) %>% select(-n) %>%
  mutate(group = factor(group, levels = order.group),
         group.coarse = factor(group.coarse, levels = order.group.coarse))

theme_bar <- theme_classic() + theme(line = element_line(size = 0.5),
                                     axis.text.y = element_text(size = 24),axis.title.y = element_blank(),
                                     axis.text.x = element_blank(),axis.title.x = element_blank(),
                                     legend.text = element_blank(), legend.title = element_blank(),
                                     plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))

#### Tcell proteomics data 20220916 ####
dt.Tcell.proteome <- read.xlsx('Tcell/protein_abundance_T_cell.xlsx') %>% mutate(tissue = ifelse(grepl('A',sample),'activated','naive')) %>%
  select(Entry,sample, tissue,  conc) %>%
  left_join(anno.group %>% mutate(group = factor(group, levels = names(color.group)))) %>%  ## anno.group from assign_group_mouse.R
  filter(!is.na(Length)) %>%
  group_by(sample, tissue) %>% mutate(f = conc*Length/sum(conc * Length, na.rm = T)) %>% ungroup() %>%
  filter(!is.na(f)) %>%
  mutate(cohort = gsub('[0-9]','',sample)) %>%
  left_join(anno.coarse)

frac.proteome.Tcell <- dt.Tcell.proteome %>%
  group_by(tissue, sample, group) %>%
  summarise(f=sum(f)) %>% ungroup()

frac.proteome.Tcell.ATP <- frac.proteome.Tcell %>%
  filter(group %in% c('glycolysis','ox phos','TCA','FA oxidation','mitochondrial amino acid')) %>%
  mutate(pathway = ifelse(group == 'glycolysis', 'glycolysis', 'respiration')) %>%
  group_by(tissue,sample,pathway) %>% summarise(f = sum(f)) %>% ungroup()

frac.proteome.Tcell.mean <- frac.proteome.Tcell %>%
  group_by(tissue, group) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
  ungroup() %>% group_by(tissue) %>% arrange(desc(group)) %>% mutate(pos = cumsum(f)) %>% ungroup()

frac.proteome.Tcell.ATP.mean <-  frac.proteome.Tcell.ATP %>%
  group_by(tissue, pathway) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
  ungroup()

mass.proteome.Tcell <- frac.proteome.Tcell %>% mutate(protein = ifelse(tissue == 'naive', 25.8, 61.9)) %>%
  mutate(mass = f * protein) %>% mutate(tissue = factor(tissue, levels = c('naive','activated')))
mass.proteome.Tcell.mean <- mass.proteome.Tcell %>%
  group_by(tissue, group) %>% summarise(mass.e = sd(mass)/sqrt(3), mass = mean(mass)) %>% 
  ungroup()
# %>% group_by(tissue) %>% arrange(desc(group)) %>% mutate(pos = cumsum(f)) %>% ungroup()
mass.FC <- mass.proteome.Tcell.mean %>% select(tissue, group, mass) %>%
  pivot_wider(names_from = 'tissue', values_from = 'mass') %>%
  mutate(FC = activated/naive)

## plot

load('proteomics_Tcell.RData')

pdf('Tcell/proteome_fraction.pdf', w = 10, h = 5)
ggplot(frac.proteome.Tcell.mean, aes(x = tissue, y = f, fill = group, group = group)) +
  geom_col(position = position_stack(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
  scale_fill_manual(values = color.group) + scale_x_discrete(limits = c('naive','activated'))+
  scale_y_continuous(limits = c(0,1.02), expand = c(0,0)) +
  egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

if(plot.coarse) {
  frac.proteome.Tcell <- dt.Tcell.proteome %>%
    group_by(tissue, sample, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.proteome.Tcell.mean <- frac.proteome.Tcell %>%
    group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  frac.proteome.Tcell.1 <- frac.proteome.Tcell %>% 
    mutate(group.coarse = ifelse(!group.coarse %in% order.group.coarse[1:3], 'metabolism',as.character(group.coarse))) %>%
    mutate(group.coarse = factor(group.coarse, levels = order.group.coarse)) %>%
    group_by(tissue, sample, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.proteome.Tcell.mean.1 <- frac.proteome.Tcell.1 %>%
    group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  pdf('Tcell/proteome_fraction_coarse_1.pdf', w = 3, h = 4)
  ggplot(frac.proteome.Tcell.mean.1, aes(x = tissue, y = f)) +
    geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.proteome.Tcell.1 %>% 
                  group_by(sample) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) + scale_x_discrete(limits = c('naive','activated'))+
    scale_y_continuous(limits = c(0,1.02), expand = c(0,0), breaks = c(0,0.5,1)) + 
    theme_bar
  dev.off()
  pdf('Tcell/proteome_fraction_coarse_2.pdf', w = 3, h = 4)
  ggplot(frac.proteome.Tcell.mean %>% filter(!group.coarse %in% order.group.coarse[1:3]) %>%
           group_by(tissue) %>% mutate(pos = pos/sum(f), f.e = f.e/sum(f),f = f/sum(f)) %>% ungroup(), 
         aes(x = tissue, y = f)) +
    geom_col(aes( fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.proteome.Tcell %>% filter(!group.coarse %in% order.group.coarse[1:3]) %>% 
                  group_by(sample) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)/sum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) + scale_x_discrete(limits = c('naive','activated'))+
    scale_y_continuous(limits = c(0,1.05), expand = c(0,0), breaks = c(0,0.5,1), position = 'left') + 
    theme_bar
  dev.off()
}

fig <- ggplot(mass.proteome.Tcell.mean , aes(x = group, y = mass, fill = tissue, group = tissue)) +
  geom_col(position = position_dodge(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = mass, ymax = mass+mass.e),position = position_dodge(),width = 0.5) +
  geom_point(data = mass.proteome.Tcell,position = position_dodge2(0.5), size = 0.5)+
  scale_fill_manual(values = c('naive' = 'white', 'activated' = 'grey'), guide = 'none')+
  theme_bar+ theme(line = element_line(size = 0.25),axis.line.x = element_blank(),axis.ticks.x = element_blank()) 
FC.color <- ggplot(mass.FC, aes(x = group, y = 1, fill = log2(FC))) + geom_col(color = 'black') +
  scale_fill_gradientn(colors = c(rev(RColorBrewer::brewer.pal(n=5,'Blues')),RColorBrewer::brewer.pal(n=4,'Reds')),
                       limits = c(log2(0.8),log2(7.2)), breaks = c(0,1.4)) +theme_void()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),text = element_text(size = 14))
anno <- ggplot(anno.coarse %>% group_by(group.coarse) %>% summarise(n = length(group))%>% ungroup() %>%
                 mutate(group.coarse = factor(group.coarse, levels = rev(order.group.coarse))),
               aes(x =1, y = n,fill = group.coarse)) +
  geom_col(color = 'black')+scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color.group.coarse, guide = 'none') + coord_flip()+
  theme_void()
pdf('Tcell/proteome_mass.pdf', w = 9, h = 5)
anno+fig +FC.color +
  patchwork::plot_layout(ncol = 1, nrow = 3, widths = c(1, 1), heights = c(0.5, 4, 0.5))
dev.off()

## plot heatmap
source('../../../R/scripts/plot_heatmap/my_plot_heatmap.R')
dt.heatmap <- dt.Tcell.proteome %>% select(Entry,sample,f) %>%
  pivot_wider(names_from = 'sample', values_from = 'f')

dt.plot <- log2_mean_center(dt.heatmap %>% select(starts_with('N')|starts_with('A')))
row.names(dt.plot) <- dt.heatmap$Entry

figure <- my_plot_heatmap (dt.plot,
                           cohort = colnames(dt.plot),
                           imgName = "Tcell/heatmap",
                           format = "pdf", # pdf, png, or jpeg
                           dpi = 72,
                           height.pixel = NA,
                           width = NA, # define output graph width
                           palette = "RdBu",  # RdBu, gbr, heat, topo, blueyellow, custom
                           cus_colors = "NA",  # accept customized color pallete
                           viewOpt = "overview", # Detail, overview
                           rowV = T, colV = F, # cluster by row/column
                           distance = 'euclidean',
                           border.pixel = F, # border for each pixel
                           grp.ave = F, # group average
                           scale_ub = 4, scale_lb = -4, # heatmap scale
                           # gaps.c = findBreak(str_extract(colnames(dt.plot), '[0-9]+min')),
                           gaps.c = c(3,6),
                           gaps.r = (0))

## volcano plot

dt.anova <- dt.Tcell.proteome %>%
  group_by(Entry,Gene.names,group,Protein.names) %>% do(anova(lm(f~tissue, data = .))) %>% ungroup() %>%
  rename(p = `Pr(>F)`) %>% filter(!is.na(p)) %>%
  mutate(p = p.adjust(p, 'BH'))

dt.anova.plot <- dt.Tcell.proteome %>% 
  group_by(Entry, tissue) %>% summarise(f = mean(f)) %>% ungroup() %>%
  pivot_wider(names_from = 'tissue', values_from = 'f') %>% mutate(FC = activated/naive) %>%
  full_join(dt.anova) %>%
  arrange(desc(FC)) %>% mutate(rank = 1:length(FC)) %>%
  mutate(pathway = case_when(group == 'glycolysis' ~ 'glycolysis',
                             group %in% c('TCA', 'ox phos') ~ 'respiration',
                             TRUE ~ 'other')) %>%
  mutate(pathway = factor(pathway, levels = c('glycolysis','respiration','other'))) %>%
  arrange(desc(pathway))

pdf('volcano plot.pdf', w = 8, h = 4)
ggplot(dt.anova.plot, aes(x = log2(FC), y = -log10(p))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(data = dt.anova.plot, aes(color = pathway), size = 1) +
  scale_color_manual(values = c('other'= 'grey','glycolysis' = 'goldenrod1', 'respiration' = '#4292C6')) +
  # geom_text_repel(data = dt.anova.plot %>% filter(!pathway =='other') %>% filter(abs(FC) > log2(80)),
  #                 aes(label = Gene.names, color = pathway)) +
  egg::theme_presentation() +theme(panel.background = element_blank())
dev.off()

save.image(file = 'proteomics_grouped_Tcell.RData')

## compare to Marchingo et al
tmp <- read.xlsx('../proteomics_other/Myc_mouse_Marchingo_et_al_2020.xlsx',sheet = 'protein mass', startRow = 3, cols = c(1:10,14:19,29))
colnames(tmp)[c(2,4,17)] = c('uniprot','gene','kDa')
colnames(tmp)[5:16] <- paste(c(rep('CD4_naive',3),rep('CD4_TCR',3),rep('CD8_naive',3),rep('CD8_TCR',3)),
                             rep(1:3,4), sep = '_')
anno.group.gene <- anno.group %>% tibble::as_tibble() %>%
  mutate(gene = strsplit(Gene.names, split = ' ')) %>% unnest(cols = 'gene') %>% select(gene, group) %>% distinct()

# this dataset is presented on a gene base
protein.marchingo <- tmp %>% pivot_longer(cols = starts_with('CD'), names_to = 'condition', values_to = 'mass') %>%
  group_by(condition) %>% mutate(f = mass/sum(mass)) %>% ungroup() %>%
  mutate(Entry = unlist(lapply(strsplit(Fasta.headers, split = '\\|'),'[',2))) %>%
  left_join(anno.group %>% select(Entry,group)) %>% mutate(group = ifelse(is.na(group),'other',group))  %>%
  mutate(cell = ifelse(grepl('CD4',condition),'CD4','CD8')) %>%
  mutate(state = ifelse(grepl('naive',condition),'naive','TCR')) %>%
  mutate(group = factor(group, levels = names(color.group)))

frac.proteome.marchingo <- protein.marchingo %>%
  group_by(cell, condition,state, group) %>% summarise(f = sum(f)) %>% ungroup()
frac.proteome.marchingo.ATP <- frac.proteome.marchingo %>%
  filter(group %in% c('glycolysis','ox phos','TCA','FA oxidation','mitochondrial amino acid')) %>%
  mutate(anno = ifelse(group == 'glycolysis', 'glyc', 'resp')) %>%
  group_by(cell,condition,state,anno) %>% summarise(f = sum(f)) %>% ungroup()
frac.proteome.marchingo.ATP.mean <- frac.proteome.marchingo.ATP %>%
  group_by(cell,state,anno) %>% summarise(e = sd(f)/sqrt(length(f)),f = sum(f)) %>% ungroup()

pdf('Tcell/proteome_fraction_marchingo.pdf', w = 12, h = 5)
ggplot(frac.proteome.marchingo, aes(x = condition, y = f, fill = group, group = group)) +
  geom_col(position = position_stack(),width = 0.7, color = 'black') + 
  scale_fill_manual(values = color.group) + # scale_x_discrete(limits = c('naive','activated'))+
  egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
pdf('Tcell/proteome_fraction_marchingo_energy_pathway.pdf', w = 8, h = 5)
ggplot(frac.proteome.marchingo.ATP.mean, aes(x = paste0(cell, '_',state), y = f, fill = anno)) +
  geom_col(position = position_dodge(width = 0.7),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = f-e,ymax = f+e),position = position_dodge(width = 0.7),width = 0.5)+
  scale_fill_manual(values = c('glyc' = '#FFED6F','resp' = '#377EB8')) + 
  egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()


## amino acid composition from proteomics
fasta.mouse <- seqinr:: read.fasta('../fastas/Mouse_proteome_3AUP00000058-2022.10.11.fasta')
tmp <- lapply(fasta.mouse,table)
aa.freq.mouse <- bind_rows(tmp, .id = 'Entry')%>% mutate(Entry = unlist(lapply(strsplit(Entry,split = '\\|'),'[',2)))
save(aa.freq.mouse, file = '../fastas/aa_frequency_mouse.RData')

aa.freq.Tcell <- aa.freq.mouse %>% pivot_longer(cols = -Entry, names_to = 'aa', values_to = 'number') %>% mutate(number = as.numeric(number)) %>%
  right_join(dt.Tcell.proteome) %>%
  group_by(aa, sample, tissue) %>% summarise(conc = sum(conc * number, na.rm = T)) %>% ungroup() %>%
  group_by(sample, tissue) %>% mutate(molar.fraction = conc/sum(conc)) %>% ungroup() %>%
  filter(!aa %in% c('b','u','x','z')) %>%
  mutate(tissue = factor(tissue, levels = c('naive','activated')))
aa.freq.Tcell.mean <- aa.freq.Tcell %>% group_by(tissue,aa) %>%
  summarise(molar.fraction = mean(molar.fraction)) %>% ungroup()

pdf('Tcell/T_cell_amino_acid_molar_composition.pdf', w = 12, h = 5)
ggplot(aa.freq.Tcell, aes(x = aa, y = molar.fraction, fill = tissue, group = tissue)) +
  geom_col(data = aa.freq.Tcell.mean, position = position_dodge2(width = 0.7), width = 0.7, color = 'black') +
  geom_point(position = position_dodge2(width = 0.7)) +
  ylab('molar fraction in protein') + xlab('amino acid') +
  scale_fill_manual(values = c('naive' = 'white','activated' = 'grey')) +
  egg::theme_presentation()
dev.off()

write.xlsx(aa.freq.Tcell.mean, 'Tcell/T_cell_amino_acid_frequency.xlsx')
