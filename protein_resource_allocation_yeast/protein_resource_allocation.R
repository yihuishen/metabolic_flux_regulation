library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)
library(tidyr)
library(egg)
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../../../../R/scripts/plot_heatmap/my_plot_heatmap.R')
source('../../../../R/scripts/misc_functions.R')
options(stringAsFactors = F)

#load(file = 'protein_resource_allocation.RData')
#### theme and annotation ####
theme_bar <- theme_classic() + theme(line = element_line(size = 0.5),
                                     axis.text.y = element_text(size = 24),axis.title.y = element_blank(),
                                     axis.text.x = element_blank(),axis.title.x = element_blank(),
                                     legend.text = element_blank(), legend.title = element_blank(),
                                     plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))

load('assign_group_files.RData') # from gen_mapping_files.R
anno.coarse <- data.frame(group = order.group, n = 1:length(order.group)) %>% 
  mutate(group.coarse = case_when(n %in% 1:2 ~ 'other',
                                  n %in% 3 ~ 'transcription',
                                  n %in% 4 ~ 'translation',
                                  n %in% 5:10 ~ 'anabolic',
                                  n %in% 11:12 ~ 'other mitochondria',
                                  TRUE ~ as.character(group))) %>% select(-n) %>%
  mutate(group.coarse = factor(group.coarse, levels = order.group.coarse))
color.group.coarse = c('other'='#D9D9D9','transcription'="#FED9A6",'translation'='#8DD3C7','metabolism'= '#BEBADA','anabolic'='#FCCDE5',
                       'other mitochondria'='#C6DBEF','TCA'='#6BAED6','ox phos'='#4292C6','glycolysis'='#FFED6F')
order.group.coarse = names(color.group.coarse)

theme_bar <- theme_classic() + theme(line = element_line(size = 0.5),
                                     axis.text.y = element_text(size = 24),axis.title.y = element_blank(),
                                     axis.text.x = element_blank(),axis.title.x = element_blank(),
                                     legend.text = element_blank(), legend.title = element_blank(),
                                     plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))


dir_abs <- '../../abs_proteomics/'
##### S.cerevisiae #####

## absolute concentration
dt.FY4.abs <- read.xlsx(paste0(dir_abs,'abs_quant_SC_FY4_replicates.xlsx')) %>%
  select(f, Entry, source) %>% left_join(map.SC) %>% group_by(Entry, source) %>%
  mutate(f = f/length(geneID)) %>% ungroup() %>% mutate(conc.AA = 869*f) 
## some peptide correspond to multiple genes 

dt.CENPK.abs <- read.xlsx(paste0(dir_abs,'abs_quant_SC_CENPK_replicates.xlsx')) %>%
  select(f, Entry, source) %>% left_join(map.SC) %>% group_by(Entry, source) %>%
  mutate(f = f/length(geneID)) %>% ungroup() %>% mutate(conc.AA = 869*f) 

# 1 OD*L culture = 1.83mL aqueous cell volume = 350mg dry cell weight = 175mg protein weight = 1.59mmol AA
# total AA concentration = 869mM

dt.SC.chemostat <- read.csv('../../SC_proteomics/normalized_chemostats_setB.csv') %>%
  mutate(Entry = unlist(lapply(strsplit(as.character(Entry), split='\\|'),'[',2))) %>%
  select(-Description, -Gene.Symbol) %>%
  pivot_longer(cols = !Entry & !B, names_to = 'condition', values_to = 'conc') %>%
  mutate(B = 1, conc = conc/B) %>%
  pivot_wider(names_from = 'condition', values_from = 'conc') %>%
  pivot_longer(cols = !Entry, names_to = 'condition', values_to = 'FC')

group.SC <-  rbind(dt.CENPK.abs %>% mutate(condition = 'abs_B') %>% mutate(strain = 'CENPK'),
                   dt.FY4.abs %>% mutate(condition = 'abs_B') %>% mutate(strain = 'FY4'),
                   dt.FY4.abs %>%
                     group_by(Entry) %>% summarise(conc.AA = mean(conc.AA)) %>% ungroup() %>% ## average of replicates
                     left_join(dt.FY4.abs %>% select(-conc.AA,-source,-f)%>% distinct()) %>% inner_join(dt.SC.chemostat) %>%
                     mutate(FC = ifelse(is.na(FC), 1, FC)) %>%
                     mutate(conc.AA = conc.AA * FC) %>% select(-FC) %>%
                     group_by(condition) %>% mutate(f = conc.AA/sum(conc.AA), conc.AA = f*869) %>%
                     ungroup() %>% mutate(strain = 'FY4', source = '011421')) %>%
  mutate(dr = case_when(grepl('1',condition) ~ 0.08,
                        grepl('2',condition) ~ 0.16,
                        grepl('3',condition) ~ 0.22,
                        grepl('4',condition) ~ 0.28,
                        strain == 'CENPK' ~ 0.391,
                        TRUE~ 0.255)) %>%
  mutate(nutr = str_extract(condition,'[CNPB]'))

frac.proteome.SC <- group.SC  %>%
  group_by(group, source,condition, dr, nutr,strain) %>% summarise(fraction.proteome = sum(f,na.rm=T), n = length(f)) %>%
  ungroup()

##### I. orientalis #####

dt.IO.abs <- read.xlsx(paste0(dir_abs,'abs_quant_IO_SD108_replicates.xlsx')) %>%
  select(f, Entry, source) %>% mutate(conc.AA = 869*f) %>%
  left_join(map.IO)

dt.IO.chemostat <- read.xlsx('../../IO_proteomics/20200803/normalized.xlsx') %>%
  select(-Description, -Gene.Symbol) %>%
  pivot_longer(cols = !Entry & !B, names_to = 'condition', values_to = 'conc') %>%
  mutate(B = 1, conc = conc/B) %>%
  pivot_wider(names_from = 'condition', values_from = 'conc') %>%
  pivot_longer(cols = !Entry, names_to = 'condition', values_to = 'FC')

group.IO <-  rbind(dt.IO.abs %>% mutate(condition = 'abs_B'),
                   dt.IO.abs %>%
                     group_by(Entry) %>% summarise(conc.AA = mean(conc.AA)) %>% ungroup() %>% ## average of replicates
                     left_join(dt.IO.abs %>% select(-conc.AA,-source,-f) %>% distinct()) %>% inner_join(dt.IO.chemostat) %>%
                     mutate(FC = ifelse(is.na(FC), 1, FC)) %>%
                     mutate(conc.AA = conc.AA * FC) %>% select(-FC) %>%
                     group_by(condition) %>% mutate(f = conc.AA/sum(conc.AA), conc.AA = f*869) %>% ungroup() %>%
                     mutate(source = '080320'))%>%
  mutate(dr = case_when(grepl('1',condition) ~ 0.12,
                        grepl('2',condition) ~ 0.18,
                        grepl('3',condition) ~ 0.23,
                        grepl('4',condition) ~ 0.34,
                        grepl('5',condition) ~ 0.45,
                        TRUE ~ 0.52)) %>%
  mutate(nutr = str_extract(condition,'[CNPB]')) %>% mutate(strain = 'SD108')

frac.proteome.IO <- group.IO  %>%
  group_by(group, condition,source,nutr, dr,strain) %>% summarise(fraction.proteome = sum(f,na.rm=T), n = length(f)) %>%
  ungroup()
  
frac.p = rbind(frac.proteome.IO %>% mutate(yeast = 'IO'),
               frac.proteome.SC %>% mutate(yeast = 'SC')) %>%
  mutate(group = factor(group, levels = order.group)) %>%
  mutate(nutr = factor(nutr, levels = c('C','N','P','B','abs_B'))) %>%
  arrange(nutr) %>%
  mutate(condition = factor(condition, levels = unique(condition)))

##### S. cerevisiae diverse respiratory and fermentative growth environments (RFGE) #####
dt.SC.RFGE <- read.delim('../../SC_proteomics/20220624/20220701_protein_quant_1423_SC.tsv', sep = '\t')[,-2:-7] %>%
  mutate(Entry = unlist(lapply(strsplit(as.character(Protein.Id), split = '\\|'),'[',2))) %>%
  select(-Protein.Id) %>%
  mutate(C = (C1+C2+C3)/3) %>%
  pivot_longer(cols = c(-Entry,-C), names_to = 'sample', values_to = 'FC') %>%
  mutate(FC = FC/C) %>% select(-C) %>%
  mutate(set = gsub('[123]','',sample)) %>%
  left_join(read.xlsx('../../SC_proteomics/20220624/sample.xlsx', sheet = 'conditions', rows = 1:8, cols = c(2,4)))

group.SC.RFGE <-  dt.CENPK.abs %>%
  group_by(Entry) %>% summarise(conc.AA = mean(conc.AA)) %>% ungroup() %>% ## average of replicates
  left_join(dt.CENPK.abs %>% select(-conc.AA,-source,-f) %>% distinct()) %>%
  inner_join(dt.SC.RFGE) %>%
  mutate(FC = ifelse(is.na(FC), 1, FC)) %>%
  mutate(FC = ifelse(FC>100, 100,FC)) %>%
  mutate(conc.AA = conc.AA * FC) %>% select(-FC) %>%
  group_by(sample) %>% mutate(tmp =sum(conc.AA),f = conc.AA/sum(conc.AA), conc.AA = f*869) %>% ungroup()%>%
  mutate(group = factor(group, levels = order.group))

frac.proteome.SC.RFGE <- group.SC.RFGE  %>%
  group_by(group, condition,sample) %>% summarise(fraction.proteome = sum(f,na.rm=T), n = length(f)) %>%
  ungroup()

##### I. orientalis diverse respiratory and fermentative growth environments (RFGE) #####
dt.IO.RFGE <- read.delim('../../IO_proteomics/20220624/20220624_protein_quant_1406_IO.tsv', sep = '\t')[,-2:-7] %>%
  mutate(Entry = unlist(lapply(strsplit(as.character(Protein.Id), split = '\\|'),'[',2))) %>%
  select(-Protein.Id) %>%
  mutate(C = (C1+C2+C3)/3) %>%
  pivot_longer(cols = c(-Entry,-C), names_to = 'sample', values_to = 'FC') %>%
  mutate(FC = FC/C) %>% select(-C) %>%
  mutate(set = gsub('[123]','',sample)) %>%
  inner_join(read.xlsx('../../IO_proteomics/20220624/sample.xlsx', sheet = 'conditions', rows = 1:8, cols = c(2,4)))

group.IO.RFGE <-  dt.IO.abs %>%
  group_by(Entry) %>% summarise(conc.AA = mean(conc.AA)) %>% ungroup() %>% ## average of replicates
  left_join(dt.IO.abs %>% select(-conc.AA,-source,-f) %>% distinct()) %>%
  inner_join(dt.IO.RFGE) %>%
  mutate(FC = ifelse(is.na(FC), 1, FC)) %>%
  mutate(FC = ifelse(FC>100, 100,FC)) %>%
  mutate(conc.AA = conc.AA * FC) %>% select(-FC) %>%
  group_by(sample) %>% mutate(f = conc.AA/sum(conc.AA), conc.AA = f*869) %>% ungroup() %>%
  mutate(group = factor(group, levels = order.group))

frac.proteome.IO.RFGE <- group.IO.RFGE  %>%
  group_by(group, condition,sample) %>% summarise(fraction.proteome = sum(f,na.rm=T), n = length(f)) %>%
  ungroup()

frac.p.RFGE <- rbind(frac.proteome.SC.RFGE %>% mutate(yeast = 'SC', strain = 'CENPK'),
                     frac.proteome.IO.RFGE%>% mutate(yeast = 'IO', strain = 'SD108'))

##

save(file = 'proteomics_grouped.RData', list = c('group.IO','group.SC','frac.p'))
save.image(file = 'protein_resource_allocation.RData')

#### choose to group proteome into less functional sectors ####
load(file = 'protein_resource_allocation.RData')
plot.coarse <- T
if(plot.coarse) {
  frac.p <- frac.p %>% left_join(anno.coarse) %>% select(-group) %>% rename(group = group.coarse) %>%
    mutate(group = factor(group, levels = order.group.coarse)) %>%
    group_by(group, condition, nutr, dr, strain, yeast, source) %>% summarise(fraction.proteome = sum(fraction.proteome),  n = sum(n)) %>% ungroup()
}
if(plot.coarse) {
  frac.proteome.IO.antimycin <- frac.proteome.IO.antimycin %>%
    left_join(anno.coarse) %>% select(-group) %>% rename(group = group.coarse) %>%
    mutate(group = factor(group, levels = order.group.coarse)) %>%
    mutate(group = factor(group, levels = order.group.coarse)) %>%
    group_by(group, condition, time,rep) %>% summarise(fraction.proteome = sum(fraction.proteome),  n = sum(n)) %>% ungroup()
}
if(plot.coarse) {
  frac.p.RFGE <- frac.p.RFGE %>% left_join(anno.coarse) %>% select(-group) %>% rename(group = group.coarse) %>%
    mutate(group = factor(group, levels = order.group.coarse)) %>%
    group_by(group, condition, sample, strain, yeast) %>% summarise(fraction.proteome = sum(fraction.proteome)) %>% ungroup()
}

#### linear response in protein fraction ####
{
  frac.p.lin <- frac.p  %>% filter(!condition == 'abs_B') %>% group_by(yeast, group) %>%
    summarise(slope = lin.coeff(fraction.proteome, dr),
              slope.se = lin.coeff.se(fraction.proteome, dr),
              slope.p = lin.coeff.p(fraction.proteome, dr),
              mean.frac = mean(fraction.proteome),
              effect.size.growth = 1/mean(dr)) %>% ungroup()%>% # 1/mean(dr) is the effect size if flux is fully controlled by protein
    mutate(effect.size = slope/mean.frac ) %>%
    mutate(flux.explained = effect.size /effect.size.growth) %>% 
    arrange(-desc(effect.size)) %>%
    mutate(group = factor(group, levels = unique(group)))
  pdf('protein_resource_allocation_lin_effect_size.pdf', w=6.5, h = 3.5)
  ggplot(frac.p.lin %>% filter(!group == 'unannotated'), aes(x=group, y = flux.explained, fill = group)) +
    geom_hline(yintercept = 0) +
    geom_col(position = position_dodge(), color = 'black', width= 0.8) +
    geom_text(aes(y = ifelse(flux.explained > 0,1,-1) * (abs(flux.explained) + 0.2),
                  label = case_when(slope.p<0.005 ~ '**',
                                    slope.p<0.05 ~ '*',
                                    TRUE ~ '')), position = position_dodge(width = 1)) +
    facet_grid(.~yeast) + ylab('flux explained by proteome allocation') +
    scale_fill_manual(values = color.group,guide = guide_legend(reverse = TRUE)) +
    coord_flip() + theme_classic()
  dev.off()
  
  # growth rate dependence
  # anabolic group
  group.ana <- frac.p.lin %>% filter(slope > 0, slope.p<0.05) %>% select(group) %>% distinct()
  group.central <- data.frame(group = c('TCA','PPP','glycolysis'))
  pdf('protein_resource_allocation_lin_central.pdf', w = 5.5, h = 4)
  ggplot(frac.p %>%
           #right_join(group.central) %>%
           filter(!group == 'unannotated', ! condition == 'abs_B') ,
         aes(x = dr, y = fraction.proteome, group = yeast, color = yeast)) +
    geom_point() + stat_smooth(method = 'lm', fullrange = TRUE, se=F) +
    facet_wrap(.~group, scales = 'free',nrow=2) +
    ylim(0,NA) +xlim(0,NA) + theme_classic()
  dev.off()
}

#### batch culture ####
{
  frac.p.batch <- frac.p %>% filter(condition == 'abs_B', !strain == 'FY4', !source == '051221') %>% rename(f=fraction.proteome)
  
  frac.p.batch.mean <- frac.p.batch %>%
    group_by(group, yeast) %>% summarise(se = sd(f)/sqrt(length(f)),
                                         f = mean(f)) %>% ungroup() %>%
    group_by(yeast) %>% arrange(desc(group)) %>% mutate(anno = cumsum(f)) %>% ungroup()
  
  frac.p.batch.1 <- frac.p.batch %>%
    mutate(group = ifelse(!group %in% order.group.coarse[1:3], 'metabolism',as.character(group)))%>%
    mutate(group= factor(group, levels = order.group.coarse)) %>%
    group_by(group, yeast,source) %>% summarise(f = sum(f)) %>% ungroup()
  frac.p.batch.1.mean <- frac.p.batch.1 %>%
    group_by(group, yeast) %>% summarise(se = sd(f)/sqrt(length(f)),
                                         f = mean(f)) %>% ungroup() %>%
    group_by(yeast) %>% arrange(desc(group)) %>% mutate(anno = cumsum(f)) %>% ungroup()
  
  pdf('proteome_allocation_SC_IO_batch_1.pdf', w = 3, h = 4)
  ggplot(frac.p.batch.1.mean,
         aes(x = yeast, y = f)) +
    geom_col(aes(fill = group),position = position_stack(), color = 'black', width = 0.7) +
    geom_jitter(data = frac.p.batch.1 %>% 
                  group_by(yeast, source) %>% arrange(desc(group)) %>% mutate(f = cumsum(f)) %>% ungroup(),
                width = 0.2, color = 'black', size = 1)+
    geom_errorbar(aes(ymin = anno, ymax = anno + se), width = 0.5,color = 'black')+
    scale_y_continuous(limits = c(0,1.05), expand = c(0,0), breaks = c(0,0.5,1)) +
    scale_x_discrete(limits = c('IO','SC'))+
    scale_fill_manual(values = color.group.coarse)+
    theme_bar
  dev.off()
  
  pdf('proteome_allocation_SC_IO_batch_2.pdf', w = 3, h = 4)
  ggplot(frac.p.batch.mean %>% filter(!group %in% order.group.coarse[1:3]) %>%
           group_by(yeast) %>% mutate(f.tmp = sum(f)) %>% 
           mutate(anno = anno/f.tmp, se = se/f.tmp,f = f/f.tmp) %>% ungroup(),
         aes(x = yeast, y = f)) +
    geom_col(aes(fill = group),position = position_stack(), color = 'black', width = 0.7) +
    geom_jitter(data = frac.p.batch %>% filter(!group %in% order.group.coarse[1:3]) %>%
                  group_by(yeast, source) %>% arrange(desc(group)) %>% mutate(f = cumsum(f)/sum(f)) %>% ungroup(),
                width = 0.2, color = 'black', size = 1)+
    geom_errorbar(aes(ymin = anno, ymax = anno + se), width = 0.5,color = 'black')+
    scale_y_continuous(limits = c(0,1.05), expand = c(0,0), breaks = c(0,0.5,1),position = 'left') +
    scale_x_discrete(limits = c('IO','SC'))+
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_bar
  dev.off()
  
  
  ## mass fold change
  
  dt.yeast.protein.mass <- read.xlsx('../external_fluxes/all_flux_measurement.xlsx', 
                                     sheet = 'biomass composition', startRow = 2, cols = c(1:3,8)) %>%
    filter(nutrient == 'B', !strain == 'FY4') %>% rename(nutr = nutrient) 
  
  mass.yeast <- frac.p %>% filter(condition == 'abs_B', !source == '051221') %>% filter(!group == 'unannotated')%>% 
    rename(f = fraction.proteome) %>%
    right_join(dt.yeast.protein.mass) %>%
    mutate(mass = f * protein)
  mass.yeast.mean <- mass.yeast %>%
    group_by(organism, strain,group) %>%
    summarise(n=length(mass),mass.e = sd(mass)/sqrt(length(mass)),mass =mean(mass)) %>% ungroup()
  mass.yeast.FC <- mass.yeast.mean %>% select(organism, group, mass) %>%
    pivot_wider(names_from = 'organism', values_from = 'mass') %>%
    mutate(FC = SC/IO) 
    
  pdf('protein content.pdf', w = 4, h = 2.5)
  ggplot(dt.yeast.protein.mass, aes(x = organism, y = 100*protein)) +
    geom_col(aes(fill = organism), position = position_dodge(),width = 0.5, color = 'black') + 
    # geom_errorbar(aes(ymin = `%protein`-error, ymax = `%protein`+error),
    #               position = position_dodge(width = 0.7),width = 0.5) +
    scale_y_continuous(limits = c(0,60),breaks = c(0,25,50),expand = c(0,0))+
    scale_fill_manual(values = c('IO' = 'white','SC' = 'grey')) +
    egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_blank())
  dev.off()
  
  fig <- ggplot(mass.yeast.mean , aes(x = group, y = mass, fill = organism, group = organism)) +
    geom_col(position = position_dodge(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = mass, ymax = mass+mass.e),position = position_dodge(),width = 0.5) +
    geom_point(data = mass.yeast ,position = position_dodge2(0.5), size = 0.5)+
    scale_fill_manual(values = c('IO' = 'white','SC' = 'grey'), guide = 'none')+
    theme_bar+ theme(line = element_line(size = 0.25),axis.line.x = element_blank(),axis.ticks.x = element_blank(), strip.text = element_blank()) 
  
  FC.color <- ggplot(mass.yeast.FC , aes(x = group, y = 1, fill = log2(FC))) + geom_col(color = 'black') +
    scale_fill_gradientn(colors = c(rev(RColorBrewer::brewer.pal(n=6,'Blues')),RColorBrewer::brewer.pal(n=6,'Reds')),
                         limits = c(-2,2), breaks = c(-2,0,2)) +theme_void()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),text = element_text(size = 14), strip.text = element_blank())
  anno <- ggplot(anno.coarse %>% filter(!group == 'unannotated') %>%group_by(group.coarse) %>% summarise(n = length(group))%>% ungroup() %>%
                   mutate(group.coarse = factor(group.coarse, levels = rev(order.group.coarse))),
                 aes(x =1, y = n,fill = group.coarse)) +
    geom_col(color = 'black')+scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = color.group.coarse, guide = 'none') + coord_flip()+
    theme_void()
  pdf('proteome_mass_foldchange.pdf', w = 9, h = 5)
  anno+fig +FC.color +
    patchwork::plot_layout(ncol = 1, nrow = 3, widths = c(1, 1), heights = c(0.5, 4, 0.5))
  dev.off()
  
}

#### plot comparison between SC and IO ####
{
  tmp <- frac.p %>% filter(!condition == 'abs_B') %>%
    arrange(desc(yeast),nutr) %>%
    mutate(condition= paste(yeast,condition)) %>% 
    mutate(condition = factor(condition, levels = unique(condition)))
  anno.dr <- ggplot(tmp %>% select(condition, dr, nutr) %>% distinct() %>% mutate(nutr = as.character(nutr)),
                    aes(x= condition, y = dr, fill = nutr)) + geom_col( width = 0.7) + theme_classic()+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_y_continuous(limits = c(0,0.55), breaks = c(0,0.2,0.4), expand = c(0,0))
  fig <- ggplot(tmp,aes(x = condition, y = fraction.proteome, fill = group)) +
    geom_col(position = position_stack(), width = 0.7, color = 'grey40', size = 0.5) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_classic() +
    scale_y_continuous(limits = c(0,1.02), expand = c(0,0)) +
    ylab('fraction in proteome') + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  figure <- ggarrange(anno.dr, fig, ncol = 1, heights = c(0.5,1.8))
  ggsave('protein_resource_allocation.pdf', figure, width = 7, height = 2.4)

  
  pdf('protein_resource_allocation_n_protein.pdf', w = 3.2, h = 3)
  ggplot(frac.p %>% filter(condition == 'abs_B') %>%
           group_by(yeast, group) %>% summarise(n = mean(n)) %>% ungroup() %>% 
           pivot_wider(names_from = 'yeast', values_from = 'n'), aes(x = SC, y = IO, fill = group)) +
    geom_abline(slope = 1, intercept = 0, color = 'grey') +
    geom_point(shape = 21,size = 4) +
    geom_text_repel(aes(label = group),  box.padding = 1, size = 2.5) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + scale_x_log10(limits = c(10,NA))+ scale_y_log10(limits = c(10,NA)) + theme_classic() +
    xlab('number of proteins detected in SC') + ylab('number of proteins detected in IO') + theme(legend.position = 'none')
  dev.off()
  
}

#### plot antimycin response in I.orientalis ####
{
  dt.plot <- frac.proteome.IO.antimycin %>% group_by(time,group) %>%
    summarise(stdev = sd(fraction.proteome), fraction.proteome = mean(fraction.proteome)) %>% ungroup() %>%
    mutate(time = as.numeric(gsub('min','',time)))
  pdf('fraction_proteome_IO_antimycin_bar.pdf', w = 4,h=3)
  ggplot(dt.plot, aes(x = as.character(time), y = fraction.proteome, fill = group)) +
    geom_col(position = position_stack(), color = 'grey40', width = 0.7, size = 0.5) +
    scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
    scale_x_discrete(limits = c('0','10','30','60','120','240')) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_classic() +
    theme(panel.border = element_rect(fill = NA,color = 'black',size = 1)) +
    ylab('fraction in proteome')
  dev.off()
  pdf('fraction_proteome_IO_antimycin.pdf', w = 4,h=2.5)
  ggplot(dt.plot %>% filter(!group %in% c('unannotated','other')), 
         aes(x = time, y = fraction.proteome, color = group, fill = group)) +
    geom_line(color = 'grey40') +
    geom_point(color = 'black',size = 3, pch = 21) +
    geom_errorbar(aes(ymin = fraction.proteome-stdev, ymax = fraction.proteome+stdev), color = 'grey40')+
    scale_y_continuous(limits = c(0,0.3), expand = c(0,0)) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + 
    scale_color_manual(values = color.group, guide = 'none') + 
    theme_classic() +
    theme(panel.border = element_rect(fill = NA,color = 'black',size = 1)) +
    ylab('fraction in proteome')
  dev.off()
  
}


#### plot RFGE result ####
{
  dt.plot <- frac.p.RFGE %>% 
    filter(yeast == 'SC') %>%
    # filter(group %in% c('glycolysis','respiration','translation/ribosome')) %>%
    group_by(yeast,group, condition) %>%
    summarise(stdev = sd(fraction.proteome), fraction.proteome = mean(fraction.proteome)) %>% ungroup() %>% 
    group_by(yeast, condition) %>%
    arrange(desc(group)) %>% 
    mutate(pos = cumsum(fraction.proteome)) %>% ungroup() 
  dt.plot.scatter <- frac.p.RFGE %>% right_join(dt.plot %>% rename(frac.mean = fraction.proteome)) %>% 
    mutate(fraction.proteome = fraction.proteome- frac.mean + pos)
  
  #pdf('fraction_proteome_RFGE_bar_IO.pdf', w = 3.2,h=2.8)
  fig1<- ggplot(dt.plot,
    aes(x = condition, y = fraction.proteome, fill = group)) +
    geom_col(position = position_stack(), color = 'grey40', width = 0.6, size = 0.5) +
    geom_errorbar(aes(ymin = pos, ymax = pos + stdev), color = 'grey40',width = 0.4, size = 0.5) +
    geom_jitter(data = dt.plot.scatter, width = 0.2, pch = 21, height = 0, size = 0.5) +
    # facet_wrap(.~yeast)+
    scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
    scale_x_discrete(limits = c('glucose + O2','glucose - O2', 'glucose + antimycin')) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_classic() +
    theme(#panel.border = element_rect(fill = NA,color = 'black',size = 1),
          #axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank()) +
    ylab('fraction in proteome')
  #dev.off()
  
  # protein fraction in biomass 0.524 (batch culture, I. o), see ATP.R or biomass composition
  fil <- data.frame(group = 'glycolysis', condition = c('glucose + O2','glucose - O2', 'glucose + antimycin'))
  #pdf('fraction_proteome_RFGE_bar_IO_glycolysis.pdf', w = 3.4,h=2.2)
  fig2 <- ggplot(dt.plot %>% right_join(fil),
         aes(x = condition, y = fraction.proteome*0.524, fill = condition)) +
    geom_col(position = position_stack(), color = 'black', width = 0.6, size = 0.5) +
    geom_errorbar(aes(ymin = pos*0.524, ymax = pos*0.524 + stdev*0.524), width = 0.4, size = 0.5) +
    geom_jitter(data = dt.plot.scatter %>% right_join(fil), width = 0.2, size = 1, pch= 21, fill = NA) +
    # facet_wrap(.~yeast)+
    scale_y_continuous(limits = c(0,0.15), expand = c(0,0)) +
    scale_x_discrete(limits = c('glucose + O2','glucose - O2', 'glucose + antimycin')) +
    scale_fill_manual(values = c('glucose + O2' = 'lightskyblue1','glucose - O2' = 'grey','glucose + antimycin' = 'grey27')) + 
    theme_classic() +
    theme(#panel.border = element_rect(fill = NA,color = 'black',size = 1),
      axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('glycolytic protein mass\n[gProtein/gDW]')
  #dev.off()
  figure <- ggarrange(fig1, fig2, ncol = 1, heights = c(3:2))
  ggsave('fraction_proteome_RFGE_IO.pdf',figure,  width = 3.6,h = 4.8)
  
  dt.plot <- dt.IO.RFGE %>% inner_join(map.IO)  %>% mutate(group = ifelse(group %in% c('TCA','ox phos'), 'respiration',group)) %>%
    filter(condition %in% c('glucose + O2','glucose - O2', 'glucose + antimycin'),
                          group %in% c('glycolysis','respiration','translation/ribosome')) %>%
    group_by(Entry, group, condition) %>% summarise(log2FC = mean(log2(FC), na.rm = T)) %>% ungroup()
  ggplot(dt.plot, aes(x = condition, y = log2FC)) +
    geom_point(position = position_dodge2(width = 0.2), alpha = 0.2, size = 1) + 
    # geom_violin() +
    facet_wrap(.~group) +
    scale_x_discrete(limits = c('glucose + O2','glucose - O2', 'glucose + antimycin'))
  
  dt.plot <- frac.p.RFGE %>% 
    group_by(yeast,group, condition) %>%
    summarise(stdev = sd(fraction.proteome), fraction.proteome = mean(fraction.proteome)) %>% ungroup() %>% 
    group_by(yeast, condition) %>%
    arrange(desc(group)) %>% 
    mutate(pos = cumsum(fraction.proteome)) %>% ungroup() %>% 
    mutate(yeast = factor(yeast, levels = c('SC','IO')))
  dt.plot.scatter <- frac.p.RFGE %>% right_join(dt.plot %>% rename(frac.mean = fraction.proteome)) %>% 
    mutate(fraction.proteome = fraction.proteome- frac.mean + pos)%>% 
    mutate(yeast = factor(yeast, levels = c('SC','IO')))
  
  pdf('fraction_proteome_RFGE_bar_SC_IO.pdf', w = 5, h = 5)
  ggplot(dt.plot,
                aes(x = condition, y = fraction.proteome)) +
    geom_col(aes(fill = group),position = position_stack(), color = 'black', width = 0.7) +
    geom_errorbar(aes(ymin = pos, ymax = pos + stdev), color = 'grey40',width = 0.5) +
    geom_jitter(data = dt.plot.scatter, width = 0.2, color = 'black', size = 1) +
    facet_wrap(.~yeast)+
    scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
    scale_x_discrete(limits = c('glucose + O2','glucose - O2', 'glucose + antimycin')) +
    scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_bar +
    theme(#panel.border = element_rect(fill = NA,color = 'black',size = 1),
      #axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank()) +
    ylab('fraction in proteome')
  dev.off()
}

#### plot fraction of proteome for a particular functional group ####
{
  group.name <- 'TCA'
  dt.plot <- group.IO.RFGE %>% filter(grepl(group.name,group)) %>%
    arrange(desc(f)) %>%
    mutate(Protein.names = substr(Protein.names, start = 1, stop =40)) %>%
    mutate(Protein.names = factor(Protein.names, levels = unique(Protein.names))) %>% filter(f > 1E-4)
  
  pdf(paste0('SC_RFGE_',group.name,'_fraction_bargraph.pdf'), w = 20, h = 4.8)
  ggplot(dt.plot,
         aes(x = sample, y = f, fill = Protein.names)) +
    geom_col(position = position_stack(), color = 'grey40') +
    scale_fill_brewer(palette = 'Set3')
  dev.off()
}

#### plot heatmap of a particular functional group of genes ####
{
  {
    tmp <- group.IO.antimycin %>% filter(grepl('amino acid',group)) %>% mutate(condition = paste(time,rep)) %>%
      select(condition, Protein.names, conc.AA, Entry,geneID) %>% distinct()  %>%
      pivot_wider(names_from = 'condition', values_from = 'conc.AA',names_sort = F)
    dt.plot <- log2((tmp %>% select_if(is.numeric))/apply(tmp %>% select(starts_with('0min')), 1, mean))
  }
  {
    tmp <- group.SC.RFGE %>% filter(grepl('TCA',group)|grepl('ox phos',group), set %in% c('C','D','E')) %>%
      mutate(set = factor(set, levels = c('C','E','D'))) %>%
      arrange(set) %>%
      mutate(condition = paste0(condition,'_',str_extract(sample,'[123]'))) %>%
      select(condition, Protein.names, conc.AA, Entry,geneID) %>% distinct()  %>%
      pivot_wider(names_from = 'condition', values_from = 'conc.AA',names_sort = F)
    tmp.2 <- log2_mean_center(tmp %>% select_if(is.numeric))
    dt.plot <- tmp.2-apply(tmp.2 %>% select(starts_with('glucose + O2')), 1, mean)
    write.table(cbind(tmp %>% select(Entry),2^dt.plot),
                file='report_SC_respiration_inhibition.tsv', sep = '\t',quote = F, row.names = F)
  }
  {
    tmp <- group.IO.RFGE %>% filter(grepl('TCA',group)|grepl('ox phos',group), set %in% c('C','D','E')) %>%
      mutate(set = factor(set, levels = c('C','E','D'))) %>%
      arrange(set) %>%
      mutate(condition = paste0(condition,'_',str_extract(sample,'[123]'))) %>%
      select(condition, Protein.names, conc.AA, Entry,geneID) %>% distinct()  %>%
      pivot_wider(names_from = 'condition', values_from = 'conc.AA',names_sort = F)
    tmp.2 <- log2_mean_center(tmp %>% select_if(is.numeric))
    dt.plot <- tmp.2-apply(tmp.2 %>% select(starts_with('glucose + O2')), 1, mean)
    write.table(cbind(tmp %>% select(Entry),2^dt.plot),
                file='report_IO_respiration_inhibition.tsv', sep = '\t',quote = F, row.names = F)
  }
  row.names(dt.plot) <- numDup(substr(tmp$Protein.names, start =1, stop = 45))
  # row.names(dt.plot) <- numDup(paste(substr(tmp$Protein.names, start =1, stop = 30), substr(tmp$rxnname, start = 1,stop = 15)))
  figure <- my_plot_heatmap (dt.plot,
                             cohort = colnames(dt.plot),
                             imgName = "heatmap",
                             format = "pdf", # pdf, png, or jpeg
                             dpi = 72,
                             height.pixel = 2,
                             width = NA, # define output graph width
                             palette = "RdBu",  # RdBu, gbr, heat, topo, blueyellow, custom
                             cus_colors = "NA",  # accept customized color pallete
                             viewOpt = "overview", # Detail, overview
                             rowV = T, colV = F, # cluster by row/column
                             distance = 'euclidean',
                             border.pixel = F, # border for each pixel
                             grp.ave = F, # group average
                             scale_ub = 2, scale_lb = -2, # heatmap scale
                             # gaps.c = findBreak(str_extract(colnames(dt.plot), '[0-9]+min')),
                             gaps.c = findBreak(gsub('\\_[123]$','',colnames(dt.plot))),
                             gaps.r = (0)
  )
}

# line plot of particular gene across conditions

id = c('P04806','Q00055','P00331','P10963')
id = c('P00358','P00359','P00360')
id = 'P30952'
dt.plot <- map.SC %>% # filter(Entry %in% id) %>%
  filter(grepl('Aconitate hydratase',Protein.names)) %>%
  left_join(rbind(dt.SC.chemostat %>% mutate(set = 'chemostat'),
                  dt.SC.RFGE %>% select(sample, FC, Entry) %>% rename(condition = sample) %>% mutate(set = 'RFGE')))
pdf('tmp.pdf', w = 6, h = 30)
ggplot(dt.plot, 
       aes(x = condition, y = log2(FC))) + geom_point() + facet_wrap(Protein.names~set, scales = 'free', ncol = 2) 
dev.off()

#### example data to set up the coarse grained modeling ####

dt.plot <- rbind(frac.p %>% right_join(data.frame(condition = c('C4','B'),
                                                  yeast = 'SC','SC',
                                                  mode = c('subopt','opt')))%>%
                   select(group, yeast, mode,fraction.proteome),
                 frac.proteome.IO.antimycin %>% right_join(data.frame(time = c('0min','240min'),
                                                                      mode = c('opt','subopt'))) %>%
                   group_by(group, mode) %>% summarise(fraction.proteome = mean(fraction.proteome)) %>% ungroup() %>% mutate(yeast = 'IO')) %>%
  filter(group %in% c('glycolysis','respiration'))

pdf('proteome_allocation_SC_IO_resp_ferm_mode.pdf', w = 4, h = 2.2)
ggplot(dt.plot,
       aes(x = mode, y = fraction.proteome, fill = group)) +
  geom_col(position = position_dodge(), color = 'grey40', width = 0.7, size = 0.5) +
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0)) +
  scale_x_discrete(limits = c('opt','subopt'))+
  facet_wrap(.~yeast)+
  scale_fill_manual(values = get(ifelse(plot.coarse,'color.group.coarse','color.group'))) + theme_classic() +
  theme(panel.border = element_rect(fill = NA,color = 'black',size = 1)) +
  ylab('fraction in proteome')
dev.off()


#### report proteomics data ####
# report absolute proteomics
abs.SC.CENPK <- group.SC %>%
  filter(condition == 'abs_B', source =='081921'|grepl('120421',source), strain == 'CENPK')%>%
  select(f, Entry, geneID, source) %>%
  pivot_wider(names_from = 'source', values_from = 'f')

abs.SC.FY4 <- group.SC %>%
  filter(condition == 'abs_B', source =='081921'|grepl('120421',source), strain == 'FY4')%>%
  select(f, Entry, geneID, source) %>%
  pivot_wider(names_from = 'source', values_from = 'f')

abs.IO.SD108 <- group.IO %>%
  filter(condition == 'abs_B', source =='081921'|grepl('120421',source), strain == 'SD108')%>%
  select(f, Entry, geneID, source) %>%
  pivot_wider(names_from = 'source', values_from = 'f')

rel.SC <- dt.SC.chemostat %>% pivot_wider(names_from = 'condition', values_from = 'FC')
rel.IO <- dt.IO.chemostat %>% pivot_wider(names_from = 'condition', values_from = 'FC')
rel.IO.anti <- dt.IO.antimycin %>% mutate(condition = gsub('\\.sum$','',gsub('X[0-9a-z]+\\_','',condition))) %>%
  select(-rep,-time) %>% pivot_wider(names_from = 'condition', values_from = 'FC')

write.xlsx(list(abs.SC.CENPK, abs.SC.FY4,abs.IO.SD108,rel.SC,rel.IO,rel.IO.anti),'report_proteomics.xlsx')

#####
## is enzyme metabolic leverage negatively correlated with its abundance?
metLev.df <- read.delim('../../result/metLev.tsv')
ord <- metLev.df %>% filter(type == 'enzyme') %>% arrange(desc(lvr))

metLev <- metLev.df %>%
  filter(type == 'enzyme') %>% rename(geneID = id) %>%
  left_join(group.IO %>% filter(condition == 'abs_B') %>% select(geneID, f, group)) %>%
  group_by(rxnname, group) %>% summarise(lvr = sum(lvr), f = sum(f)) %>% ungroup() %>%
  mutate(rxnname = factor(rxnname, levels = unique(ord$rxnname))) %>%
  mutate(group = factor(group, levels = levels(frac.p$group))) %>% filter(!is.na(group))

figure <- ggplot(metLev, aes(x=f,y=lvr, fill = group)) +geom_point(size = 3, pch=21, color = 'black')+
  scale_fill_brewer(palette = 'Set3')
pdf('regulation~protein_cost-scatter_plot.pdf', w = 5.5, h = 3.5)
print(figure)
dev.off()

# figure <- ggplot(metLev, aes(x = rxnname, y = f, fill = group)) +
#   geom_col(position = position_stack(), color = 'black') +
#   coord_flip() +
#   scale_fill_brewer(palette = 'Set3')
# 
# pdf('regulation_protein_cost.pdf', w = 8, h = 12)
# print(figure)
# dev.off()

# fisher's test
tmp <- metLev %>% group_by(group) %>% summarise(all = length(group)) %>% ungroup() %>%
  left_join(metLev %>% filter(lvr>0.5) %>% group_by(group) %>% summarise(enzyme = length(group)) %>% ungroup()) %>%
  mutate(all = ifelse(is.na(all),0,all),enzyme = ifelse(is.na(enzyme),0,enzyme))
data <- tmp[,-1]
row.names(data) <- tmp$group
fisher.test(data, simulate.p.value = T)
# not significant :()




# stat_ecdf(geom = "step")