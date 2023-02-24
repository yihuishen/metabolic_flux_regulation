library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(tidyr)
library(openxlsx)
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('mouse.RData')

#### assign group ####
load('assign_group_mouse_files.RData')
load('proteomics_Tcell.RData')
# this will load the theme_bar and color scales

#### tumor proteomics data 20220816 ####
dt.tumor.proteome <- rbind(read.xlsx('tumor/protein_abundance_all.xlsx', sheet = 'pancreas'),
                  read.xlsx('tumor/protein_abundance_all.xlsx', sheet = 'spleen')) %>%
  left_join(read.xlsx('tumor/protein_abundance_all.xlsx', sheet = 'info') %>%
              mutate(tissue = factor(tissue,levels = unique(tissue))) ) %>%
  select(Entry,sample, tissue,  conc) %>%
  left_join(anno.group %>% mutate(group = factor(group, levels = names(color.group)))) %>%  ## anno.group from assign_group_mouse.R
  filter(!is.na(Length)) %>%
  group_by(sample, tissue) %>% mutate(f = conc*Length/sum(conc * Length, na.rm = T)) %>% ungroup()%>%
  filter(!is.na(f)) %>%
  mutate(cohort = gsub('[0-9]','',sample))

#### proteome fraction tumor ####
load('tumor.RData')
frac.proteome.tumor <- dt.tumor.proteome %>%
  group_by(tissue, sample, group) %>%
  summarise(f=sum(f)) %>% ungroup()

frac.proteome.tumor.ATP <- frac.proteome.tumor %>%
  filter(group %in% c('glycolysis','ox phos','TCA','FA oxidation','mitochondrial amino acid')) %>%
  mutate(pathway = ifelse(group == 'glycolysis', 'glycolysis', 'respiration')) %>%
  group_by(tissue,sample,pathway) %>% summarise(f = sum(f)) %>% ungroup()

frac.proteome.tumor.mean <- frac.proteome.tumor %>%
  group_by(tissue, group) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
  ungroup() %>% group_by(tissue) %>% arrange(desc(group)) %>% mutate(pos = cumsum(f)) %>% ungroup()

frac.proteome.tumor.ATP.mean <-  frac.proteome.tumor.ATP %>%
  group_by(tissue, pathway) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
  ungroup()

pdf('tumor/proteome_fraction.pdf', w = 11, h =6)
ggplot(frac.proteome.tumor.mean, aes(x = tissue, y = f, fill = group, group = group)) +
  geom_col(position = position_stack(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
  scale_fill_manual(values = color.group) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0,0))+
  egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# accounting for protein composition in tissue wet weight
dt.tumor.protein.mass <- read.xlsx('tumor/protein_abundance_all.xlsx', sheet = 'info') %>%
  mutate(tissue = factor(tissue, levels = unique(tissue))) %>%
  mutate(tissue.type = factor(ifelse(grepl('spleen',tissue), 'spleen','pancreas'), levels = c('spleen','pancreas')))
dt.tumor.protein.mass.mean <- dt.tumor.protein.mass %>% 
  group_by(tissue, tissue.type) %>% summarise(error = sd(`%protein`)/sqrt(3),`%protein` = mean(`%protein`)) %>% ungroup()

mass.tumor <- frac.proteome.tumor %>%
  left_join(dt.tumor.protein.mass.mean %>% select(-error)) %>%
  mutate(mass = f * `%protein`/100)
mass.tumor.mean <- frac.proteome.tumor.mean %>%
  left_join(dt.tumor.protein.mass.mean %>% select(-error)) %>%
  mutate(mass = f * `%protein`/100,mass.e = f.e* `%protein`/100)
mass.tumor.FC <- mass.tumor.mean %>% select(tissue, group, mass) %>%
  pivot_wider(names_from = 'tissue', values_from = 'mass') %>%
  mutate(FC.spleen = `leukemic spleen`/spleen, FC.pancreas = (`GEMM PDAC` +`flank PDAC`)/2/pancreas) %>% 
  select(group, starts_with('FC')) %>%
  pivot_longer(cols = starts_with('FC'), values_to = 'FC', names_to = 'tissue.type') %>%
  mutate(tissue.type = factor(gsub('FC\\.', '',tissue.type), levels = c('spleen','pancreas')))

color.tissue <- c('spleen' = 'white','leukemic spleen' = 'grey','pancreas' ='white', 'GEMM PDAC' = 'grey','flank PDAC'='grey60')

tissue.plot <- 'pancreas'

pdf('tumor/protein content in tissue wet weight.pdf', w = 5.2, h = 8)
ggplot(dt.tumor.protein.mass.mean, aes(x = tissue, y = `%protein`)) +
  geom_col(aes(fill = tissue), position = position_dodge(),width = 0.5, color = 'black') + 
  geom_errorbar(aes(ymin = `%protein`-error, ymax = `%protein`+error),
                position = position_dodge(width = 0.7),width = 0.5) +
  geom_point(data = dt.tumor.protein.mass, position = position_dodge(width = 0.7)) + 
  scale_y_continuous(limits = c(0,23),expand = c(0,0))+
  scale_fill_manual(values = color.tissue) +
  facet_wrap(.~tissue.type, ncol = 1, scales = 'free')+
  egg::theme_presentation() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_blank())
dev.off()
fig <- ggplot(mass.tumor.mean %>% filter(tissue.type == tissue.plot), aes(x = group, y = mass, fill = tissue, group = tissue)) +
  geom_col(position = position_dodge(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = mass, ymax = mass+mass.e),position = position_dodge(),width = 0.5) +
  geom_point(data = mass.tumor %>% filter(tissue.type == tissue.plot),position = position_dodge2(0.5), size = 0.5)+
  scale_fill_manual(values = color.tissue, guide = 'none')+
  theme_bar+ theme(line = element_line(size = 0.25),axis.line.x = element_blank(),axis.ticks.x = element_blank(), strip.text = element_blank()) 
FC.color <- ggplot(mass.tumor.FC %>% filter(tissue.type == tissue.plot), aes(x = group, y = 1, fill = log2(FC))) + geom_col(color = 'black') +
  scale_fill_gradientn(colors = c(rev(RColorBrewer::brewer.pal(n=6,'Blues')),RColorBrewer::brewer.pal(n=6,'Reds')),
                       limits = c(-2.5,2.5), breaks = c(-2,0,2)) +theme_void()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),text = element_text(size = 14), strip.text = element_blank())
anno <- ggplot(anno.coarse %>% group_by(group.coarse) %>% summarise(n = length(group))%>% ungroup() %>%
                 mutate(group.coarse = factor(group.coarse, levels = rev(order.group.coarse))),
               aes(x =1, y = n,fill = group.coarse)) +
  geom_col(color = 'black')+scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color.group.coarse, guide = 'none') + coord_flip()+
  theme_void()
pdf(paste0('tumor/proteome_mass_',tissue.plot,'.pdf'), w = 9, h = 5)
anno+fig +FC.color +
  patchwork::plot_layout(ncol = 1, nrow = 3, widths = c(1), heights = c(0.5, 5, 0.5))
dev.off()

#### plot coarse proteome fraction ####
if(plot.coarse) {
  dt.tumor.proteome.coarse <- dt.tumor.proteome %>% left_join(anno.coarse) %>% 
    mutate(tissue.type = ifelse(grepl('spleen',tissue), 'spleen','pancreas'))
    
  frac.tumor.coarse <- dt.tumor.proteome.coarse  %>%
    group_by(tissue, tissue.type,sample, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.tumor.coarse.mean <- frac.tumor.coarse %>%
    group_by(tissue, tissue.type,group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  frac.tumor.coarse.1 <- frac.tumor.coarse %>% 
    mutate(group.coarse = ifelse(!group.coarse %in% order.group.coarse[1:3], 'metabolism',as.character(group.coarse))) %>%
    mutate(group.coarse = factor(group.coarse, levels = order.group.coarse)) %>%
    group_by(tissue, tissue.type,sample, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.tumor.coarse.mean.1 <- frac.tumor.coarse.1 %>%
    group_by(tissue,tissue.type, group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue,tissue.type) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  tissue.plot <- 'spleen'
  pdf(paste0('tumor/proteome_fraction_coarse_',tissue.plot,'_1.pdf'), w = 3, h = 4)
  ggplot(frac.tumor.coarse.mean.1 %>% filter(tissue.type== tissue.plot), aes(x = tissue, y = f)) +
    geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.tumor.coarse.1 %>% filter(tissue.type== tissue.plot)%>% 
                  group_by(sample) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) + 
    scale_y_continuous(limits = c(0,1.02), expand = c(0,0), breaks = c(0,0.5,1)) + 
    theme_bar
  dev.off()
  pdf(paste0('tumor/proteome_fraction_coarse_',tissue.plot,'_2.pdf'), w = 3, h = 4)
  ggplot(frac.tumor.coarse.mean %>% filter(tissue.type== tissue.plot)%>% filter(!group.coarse %in% order.group.coarse[1:3]), 
         aes(x = tissue, y = f)) +
    geom_col(aes( fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.tumor.coarse %>% filter(tissue.type== tissue.plot)%>% filter(!group.coarse %in% order.group.coarse[1:3]) %>% 
                  group_by(sample) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) +
    scale_y_continuous(limits = c(0,NA), expand = c(0,0), 
                       breaks = c(0,0.1,0.2), position = 'right') + 
    theme_bar
  dev.off()
}

#### paxdb data ####
# load('mouse.RData')
if(!file.exists('data_and_mapping_mouse.RData')) {
  file_list <- list.files(path="paxdb-abundance-files-v4.1/10090") # 10090 for mouse organism ID
  
  # file_read <- file_list[grepl('integrated', file_list)]
  file_read <- file_list
  
  dt.list <- lapply(lapply(paste0('paxdb-abundance-files-v4.1/10090/',file_read), read.delim, sep = '\t', header = F, comment.char = '#'),'[',c(1:3))
  names(dt.list) <- file_read
  
  description <- read.delim('description.tsv', sep = '\t')
  
  dt <- do.call(rbind, dt.list)
  rm('dt.list')
  dt <- dt %>% mutate(filename = gsub('\\.[0-9]+$','',rownames(dt))) %>%
    left_join(description) %>%
    mutate(string = gsub('10090\\.','',V2)) %>% rename(internalID = V1,ppm = V3) %>% select(-V2)
  uniprot.proteome <- read.delim('mouse_annotation_uniprot_proteome.tab', sep = '\t', stringsAsFactors = F)
  uniprot.full <- read.delim('mouse_annotation_uniprotKB.tsv', sep = '\t', stringsAsFactors = F)
  
  ## mapping string id to uniprot id
  if(!file.exists('mapping_mouse_string_uniprot_paxdb_reduced.tsv')) {
    mapping.tmp <- data.table::fread("full_uniprot_2_string.04_2015.tsv", sep = '\t')
    colnames(mapping.tmp) <- c('species','uniprot_id','string','identity','bit_score')
    mapping.tmp <- mapping.tmp %>% filter(species == '10090') %>%
      mutate(Entry = unlist(lapply(strsplit(as.character(uniprot_id), split = '\\|'), '[',1))) %>%
      mutate(Entry.name = unlist(lapply(strsplit(as.character(uniprot_id), split = '\\|'), '[',2)))
    mapping.reduced <- mapping.tmp  %>% mutate(rank = 1:length(string)) %>%
      mutate(in.uniprot = ifelse(Entry %in% uniprot.proteome$Entry, 1, 0)) %>% # in the list of reviewed proteome from uniprot
      group_by(string) %>% filter(in.uniprot == max(in.uniprot)) %>% filter(identity == max(identity))  %>%
      filter(bit_score == max(bit_score)) %>% filter(rank == min(rank)) %>% ungroup() 
    write.table(mapping.reduced,'mapping_mouse_string_uniprot_paxdb_reduced.tsv', sep= '\t', quote = F, row.names=F)
  } else{mapping.reduced <- read.delim('mapping_mouse_string_uniprot_paxdb_reduced.tsv', sep= '\t')}
  
  mapping.uniprot.add <- read.delim('mapping_from_uniprotKB.tsv', sep = '\t') %>% rename(Entry = From, string = To) %>% 
    anti_join(mapping.reduced %>% select(string)) %>% inner_join(dt %>% select(string))
  # some mappings are absent in the paxdb file, add from uniprot id mapping

  mapping <- rbind(mapping.uniprot.add,
                   mapping.reduced %>% select(Entry, string))
  
  save(list = c('dt','uniprot.proteome','uniprot.full','mapping.reduced','mapping'),file = 'data_and_mapping_mouse.RData')
} else {load('data_and_mapping_mouse.RData')}

#### map string from paxdb to uniprot entry ####
# note that one uniprot id could be mapped to multiple ensembl protein id

mapping.redundancy <- dt %>% select(string) %>% distinct() %>% left_join(mapping) %>%
  group_by(string) %>% mutate(n = length(string)) %>% ungroup() # n, number of Entry mapped to the same string
table(mapping.redundancy %>% select(-Entry) %>% distinct() %>% select(n)) # check how distribution of redundancy
dt.mouse <- dt %>% left_join(mapping.redundancy) %>% mutate(ppm = ppm/n) %>%
  group_by(Entry, dataset,tissue) %>% summarise(f  = sum(ppm)/1E6) %>% ungroup() %>%
  left_join(anno.group)%>%
  filter(!is.na(as.character(tissue)), !as.character(tissue) == '') 

#### proteome fraction paxdb data ####
frac.proteome <- dt.mouse %>%
  filter(!is.na(group)) %>% # some protein strings are not mapped to uniprot entry, but these only account for a small fraction
  group_by(tissue, dataset, group) %>%
  summarise(f = sum(f)) %>% ungroup() %>% filter(!grepl('integrated',tissue)) %>%
  group_by(tissue, dataset) %>% mutate(f = f/sum(f)) %>% ungroup()

frac.proteome.mean <- frac.proteome %>%
  group_by(tissue, group) %>% summarise(f.e = sd(f)/sqrt(length(f)),f = mean(f)) %>% ungroup()
order.tissue = frac.proteome.mean %>% filter(group == 'glycolysis') %>% arrange(desc(f)) %>% select(tissue) %>% distinct()

pdf('tissue/mouse_proteome_fraction_mitochondria_average_across_dataset.pdf', w = 7, h = 6)
ggplot(frac.proteome.mean, aes(x = tissue, y = f, fill = group)) +
  geom_col(position = position_stack(),size = 0.25, color = 'black') + scale_x_discrete(limits = order.tissue$tissue)+
  scale_fill_manual(values = color.group) + coord_flip()+theme(line = element_line(size = 0.25))
dev.off()

#### plot coarse ####
if(plot.coarse) {
  
  frac.mouse.coarse <- frac.proteome %>% left_join(anno.coarse) %>%
    group_by(tissue, dataset, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.mouse.coarse.mean <- frac.mouse.coarse %>%
    group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(length(f)), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  frac.mouse.coarse.1 <- frac.mouse.coarse %>% 
    mutate(group.coarse = ifelse(!group.coarse %in% order.group.coarse[1:3], 'metabolism',as.character(group.coarse))) %>%
    mutate(group.coarse = factor(group.coarse, levels = order.group.coarse)) %>%
    group_by(tissue, dataset, group.coarse) %>%
    summarise(f=sum(f)) %>% ungroup()
  
  frac.mouse.coarse.mean.1 <- frac.mouse.coarse.1 %>%
    group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
    ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()
  
  order.tissue = frac.mouse.coarse.mean %>% filter(group.coarse %in% order.group.coarse[7:9]) %>% 
    group_by(tissue) %>% summarise(f = sum(f)) %>% ungroup() %>% 
    arrange(desc(f)) %>% select(tissue)
  
  pdf('tissue/mouse_proteome_coarse1.pdf', w = 10, h = 5)
  ggplot(frac.mouse.coarse.mean.1,
         aes(x = tissue, y = f)) +
    geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.mouse.coarse.1 %>% 
                  group_by(dataset) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) + 
    scale_x_discrete(limits = order.tissue$tissue) +
    scale_y_continuous(limits = c(0,1.02), expand = c(0,0), breaks = c(0,0.5,1)) + 
    theme_bar + theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1))
  dev.off()
  pdf('tissue/mouse_proteome_coarse2.pdf', w = 10, h = 5)
  ggplot(frac.mouse.coarse.mean %>%  filter(!group.coarse %in% order.group.coarse[1:3]), 
         aes(x = tissue, y = f)) +
    geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
    geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
    geom_jitter(data = frac.mouse.coarse %>% filter(!group.coarse %in% order.group.coarse[1:3]) %>% 
                  group_by(dataset) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
                width = 0.2, color = 'black', size = 1)+
    scale_fill_manual(values = color.group.coarse) +
    scale_x_discrete(limits = order.tissue$tissue) +
    scale_y_continuous(limits = c(0,NA), expand = c(0,0), 
                       breaks = c(0,0.2,0.4,0.6)) + 
    theme_bar+ theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1))
  dev.off()
}

#### calculate proteome efficiency in mouse ####

# using CB measured fluxes
mouse.flux <- rbind(read.xlsx('mouse_flux_revised.xlsx', sheet = 'TCA') %>%  # unit, nmole/min/gTissue (wet tissue weight?)
                      rename(flux = median, lb = lower, ub = upper) %>%
                      select(tissue, flux, lb, ub) %>% mutate(pathway = 'respiration'),
                read.xlsx('mouse_flux_revised.xlsx', sheet = 'glucose') %>%
                  mutate(lb= flux - err, ub = flux + err) %>% select(-err) %>% mutate(pathway = 'glycolysis')) %>%
  mutate(flux.ATP = flux * ifelse(pathway == 'respiration', 14.5, 2)*60/1E6,
         lb = lb * ifelse(pathway == 'respiration', 14.5, 2)*60/1E6,
         ub = ub * ifelse(pathway == 'respiration', 14.5, 2)*60/1E6) # mmole/h/gTissue

# using tissue slice oxygen consumption, and assuming glycolytic flux = OCR/6
# mouse.flux <- rbind(read.xlsx('mouse_flux.xlsx', sheet = 'OCR') %>% mutate(pathway = 'respiration'),
#                 read.xlsx('mouse_flux.xlsx', sheet = 'OCR') %>% mutate(pathway = 'glycolysis') %>% 
#                   mutate(flux = flux/6)) %>% mutate( ub = (flux-err)/6, lb = (flux+err)/6) %>% select(-err) %>%
#   mutate(flux.ATP = flux * ifelse(pathway == 'mitochondria', 2.8*2, 2)) %>% # unit, nmole/min/gTissue (wet tissue weight?)
#   mutate(flux.ATP = flux.ATP *60/1E6) 

tissue.mapping <- read.xlsx('mouse_flux.xlsx', sheet = 'name') 

dt.mouse.ATP <- frac.mouse.coarse.mean %>% 
  right_join(data.frame(group.coarse = c("TCA","ox phos","glycolysis"), 
                        pathway = c('respiration','respiration','glycolysis'))) %>%
  rename(tissue.paxDB = tissue) %>%
  inner_join(tissue.mapping) %>%
  group_by(tissue, pathway) %>% summarise(f = sum(f)) %>% ungroup() %>%
  inner_join(mouse.flux) %>%
  mutate(organism = 'mouse',f.e = NA, `%protein` = 20) %>% rename(note = tissue)

#### combine with tumor ATP data ####

dt.tumor.ATP <- frac.proteome.tumor.ATP.mean %>%
  left_join(dt.tumor.protein.mass.mean %>% select(-error))%>%
  left_join(mouse.flux %>% right_join(data.frame(tissue = c('pancreas','GEMMPDAC','flankPDAC','controlSpleen_forLeukemia','leukemicSpleen'),
                                                 tissue.name = c('pancreas','GEMM PDAC','flank PDAC','spleen','leukemic spleen'))) %>%
              select(-tissue) %>%rename(tissue = tissue.name))

save.image('mouse.RData')
save(list = c('dt.mouse.ATP','dt.tumor.ATP'), file = 'ATP_tissue_tumor.RData')

#### proteome fraction and ATP production, supplementary figure ####
tmp <- frac.mouse.coarse %>% rename(tissue.paxDB = tissue) %>% left_join(tissue.mapping) %>%
  filter(!is.na(tissue)) %>% select(-tissue.paxDB) %>% filter(tissue %in% dt.mouse.ATP$note, !tissue == 'brain')

tmp.mean <- tmp %>%
  group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(length(f)), f = mean(f)) %>% 
  ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()

tmp.1 <- frac.mouse.coarse %>% rename(tissue.paxDB = tissue) %>% left_join(tissue.mapping) %>%
  filter(!is.na(tissue)) %>% select(-tissue.paxDB) %>% filter(tissue %in% dt.mouse.ATP$note, !tissue == 'brain')%>% 
  mutate(group.coarse = ifelse(!group.coarse %in% order.group.coarse[1:3], 'metabolism',as.character(group.coarse))) %>%
  mutate(group.coarse = factor(group.coarse, levels = order.group.coarse)) %>%
  group_by(tissue, dataset, group.coarse) %>%
  summarise(f=sum(f)) %>% ungroup()

tmp.mean.1 <- tmp.1 %>%
  group_by(tissue, group.coarse) %>% summarise(f.e = sd(f)/sqrt(3), f = mean(f)) %>% 
  ungroup() %>% group_by(tissue) %>% arrange(desc(group.coarse)) %>% mutate(pos = cumsum(f)) %>% ungroup()

order.tissue.tmp <- dt.mouse.ATP %>% filter(pathway == 'respiration') %>% arrange(desc(flux.ATP)) %>%
  rename(tissue=note) %>% select(tissue) %>% distinct()


fig1 <- ggplot(tmp.mean.1,
       aes(x = tissue, y = f)) +
  geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
  geom_jitter(data = tmp.1 %>% 
                group_by(dataset) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
              width = 0.2, color = 'black', size = 1)+
  scale_fill_manual(values = color.group.coarse, guide = 'none') + 
  scale_x_discrete(limits = order.tissue.tmp$tissue) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0,0), breaks = c(0,0.5,1)) + 
  theme_bar + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

fig2 <- ggplot(tmp.mean %>%  filter(!group.coarse %in% order.group.coarse[1:3]), 
       aes(x = tissue, y = f)) +
  geom_col(aes(fill = group.coarse), position = position_stack(),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = pos, ymax = pos+f.e),width = 0.5) +
  geom_jitter(data = tmp %>% filter(!group.coarse %in% order.group.coarse[1:3]) %>% 
                group_by(dataset) %>% arrange(desc(group.coarse)) %>% mutate(f = cumsum(f)) %>% ungroup(), 
              width = 0.2, color = 'black', size = 1)+
  scale_fill_manual(values = color.group.coarse, guide = 'none') +
  scale_x_discrete(limits = order.tissue.tmp$tissue) +
  scale_y_continuous(limits = c(0,NA), expand = c(0,0), 
                     breaks = c(0,0.2,0.4,0.6)) + 
  theme_bar+ theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

fig3 <- ggplot(dt.mouse.ATP %>% filter(!note == 'brain') %>%group_by(note) %>% 
                 summarise(ub = sum(ub), flux.ATP = sum(flux.ATP)) %>% ungroup(),
               aes(x = note, y = flux.ATP)) +
  geom_col(position = position_stack(reverse=T),width = 0.7, color = 'black') + 
  geom_errorbar(aes(ymin = flux.ATP, ymax = ub),width = 0.5) +
  scale_x_discrete(limits = order.tissue.tmp$tissue)  +
  coord_cartesian(ylim = c(0,NA),clip = 'on') +scale_y_continuous(expand = c(0,0), breaks = 10*((0:2))) +
  theme_bar+ theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

fig4 <- ggplot(dt.mouse.ATP %>% filter(!note == 'brain') %>%group_by(note) %>% 
                 mutate(f.ATP = flux.ATP/(sum(flux.ATP))) %>% ungroup(),
               aes(x = note, y = f.ATP, fill = pathway)) +
  geom_col(aes(fill = pathway), position = position_stack(reverse=T),width = 0.7, color = 'black') + 
  scale_fill_manual(values = c('glycolysis' = '#FFED6F','respiration' = '#377EB8'), guide = 'none') +
  scale_x_discrete(limits = order.tissue.tmp$tissue)  +
  coord_cartesian(ylim = c(0,1.02),clip = 'on') +scale_y_continuous(expand = c(0,0), breaks = c(0.2,0.5,1)) +
  theme_bar+
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1))
fig5 <- fig4 + coord_cartesian(ylim = c(0,0.05),clip = 'on') +scale_y_continuous(expand = c(0,0), breaks = c(0,0.1))
pdf('tissue/mouse_proteome_flux_selected.pdf', w = 5, h = 20)
fig1+fig2+fig3 + fig4 + fig5 + patchwork::plot_layout(ncol = 1, nrow = 5, widths = c(1), heights = c(1, 1,1,0.8,0.2))
dev.off()


#### volcano plot ####
# tumorous vs healthy tissue (g/gTissue)
dt.tmp <- dt.tumor.proteome %>% filter(!grepl('spleen',tissue),!grepl('GEMM',tissue)) %>% 
  left_join(dt.tumor.protein.mass.mean) %>% mutate(mass = f * `%protein`) %>% 
  mutate(tissue = ifelse(grepl('PDAC',tissue),'tumorous','healthy'))
  #mutate(tissue = ifelse(grepl('leukemic',tissue),'tumorous','healthy'))

dt.anova <- dt.tmp %>% # g/gTissue
  group_by(Entry) %>% do(model = anova(lm(mass ~ tissue, data = .))) %>% broom::tidy(model) %>% 
  filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value, method = 'BH')) %>%
  select(Entry,p.value) %>%
  left_join(dt.tmp %>% group_by(Entry, tissue) %>% summarise(mass = mean(mass)) %>%
              pivot_wider(names_from = 'tissue', values_from = 'mass') %>% mutate(FC = tumorous/healthy)) %>%
  left_join(anno.group)

pdf('volcano_plot_spleen.pdf', w = 10, h = 6)
ggplot(dt.anova %>% filter(!group %in% order.group[17:21]),
       aes(x = log2(FC), y = -log10(p.value))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(color = 'grey',size = 1) + 
  geom_point(data=dt.anova %>% filter(group %in% order.group[17:20]), color = "#4292C6", size = 1) +
  geom_point(data=dt.anova %>% filter(group %in% order.group[21]), pch=21, color = 'orange', fill = "#FFED6F", size = 1) +
  egg::theme_presentation() +theme(panel.background = element_blank())
dev.off()

a <- dt.tmp %>% select(Entry, sample,f) %>% pivot_wider(names_from = 'sample', values_from = 'f') %>%
  left_join(dt.anova %>% select(Entry, FC, p.value))

a=dt.anova %>% filter(group %in% order.group[21]) %>% arrange(desc(FC)) %>% filter(p.value<0.05)
b = a[c(1:5,(nrow(a)-4):nrow(a)),] %>% left_join(anno.group)
