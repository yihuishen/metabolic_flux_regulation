library(dplyr)
library(openxlsx)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = F)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = 'flux_enzyme_correlation_20220128.RData')
# contains all functions and source data files

#####
source('../../../R/scripts/plot_heatmap/my_plot_heatmap.R')
source('../../../R/scripts/misc_functions.R')
options(stringAsFactors = F)

load('protein_resource_allocation/assign_group_files.RData')
load('gpr.RData')

flx.plot <- 'mfa' # or 'pfba'
coarse <- T

#### color annotation ####
subsystem.new <- data.frame(id = c('ACS_c','PDH_m','FACOAE171_c','ASPTAi_m',
                                   'ALCD2i2_c','ALCD2i2_m','PYRDC_c',
                                   'ADPATPt_c_m','PItps_m',
                                   'SBP_c','FBA3_c',
                                   'ME1_m','FBP_c','PC_c','ALDD2y_c','ALDD2y_m','OAADC_c','PPCK_c','PPA_c','ICDHyi_c','LDHDf_c'),
                            new = c('Lipid biosynthesis','Citric acid cycle','Fatty acid degradation','Alanine, aspartate and glutamate metabolism',
                                    'Glycolysis / Gluconeogenesis','Glycolysis / Gluconeogenesis','Glycolysis / Gluconeogenesis',
                                    'Oxidative phosphorylation','Oxidative phosphorylation',
                                    'Amino sugar and nucleotide sugar metabolism','Amino sugar and nucleotide sugar metabolism',
                                    'Unassigned','Unassigned','Unassigned','Unassigned','Unassigned','Unassigned','Unassigned','Unassigned','Unassigned','Unassigned'),
                            stringsAsFactors = F)

anno.subsystem <- read.delim('../abs_proteomics/subsystem_annotation_IO.tsv', sep = '\t') 
# color.subsystem <- c(RColorBrewer::brewer.pal(5,'Set1'), 'gold3',
#                      RColorBrewer::brewer.pal(7,'Set1')[7],'grey10',rep('grey60',9))
# names(color.subsystem) <- c('glycolysis','ox phos','TCA','PPP','transport','biosynthetic','other')

color.subsystem <- c(RColorBrewer::brewer.pal(n=12,'Paired')[c(6,2,4,2,8)],
                     'black',rep('grey60',7))
# color.subsystem <- c(RColorBrewer::brewer.pal(n=12,'Paired')[c(6,2,2,2)],rep('grey50',9))

fill.subsystem <- c(RColorBrewer::brewer.pal(n=12,'Paired')[c(5,1,3,1,7)], rep('grey90',6),NA,NA)
names(color.subsystem) <- c('glycolysis','ox phos','TCA','respiration','PPP',
                            'biosynthetic','folate','sugar','lipid','amino acid','nucleic acid','transport','other')
names(fill.subsystem) <- names(color.subsystem)

save(list=c('anno.subsystem','subsystem.new','color.subsystem','fill.subsystem'),file = 'flux_enzyme_correlation/subsystem.anno.RData')

#### proteomics data ####
load('protein_resource_allocation/proteomics_grouped.RData')
condition.IO <- group.IO %>% select(condition, nutr, dr,strain) %>% distinct() %>% filter(!grepl('abs',condition))
condition.SC <- group.SC %>% select(condition, nutr, dr,strain) %>% distinct() %>%
  mutate(condition = ifelse(strain == 'CENPK','CENPK',condition)) %>% filter(!grepl('abs',condition))
# difference btw abs_B and B:
# abs_b is normalized to all proteins detected in iBAQ of batch culture, B is normalized to all proteins detected in 16-plex

#### flux IO ####
flux.IO <- read.csv('../IO_fluxomics/20220128_flux.csv') %>%
  mutate(id = gsub('q7','qx',id), reaction = gsub('q7','qx',reaction))%>%
  left_join(subsystem.new) %>% mutate(subsystem = ifelse(is.na(new), subsystem, new)) %>% select(-new) %>% left_join(anno.subsystem)

## turn gene in v2 [ioXXXX] to v1 [JL09_gXXXX] so that it can be mapped to protein id ##
mapping.gene <- read.xlsx('../IO_fluxomics/model_genes_map.xlsx', sheet = 'genes') %>% select(id, fwd_hit, rev_hit) %>%
  pivot_longer(!id,names_to = 'hit',values_to = 'geneID.v1') %>% select(-hit) %>% distinct() %>%
  rename(geneID = id) 
gpr.IO <- unnest(tibble(id = flux.IO$id, geneID=str_extract_all(flux.IO$gpr, 'io[0-9]+')),
                 cols = c(geneID)) %>% filter(!is.na(geneID)) %>%
  left_join(mapping.gene) %>% mutate(geneID = geneID.v1) %>% select(-geneID.v1) %>% distinct() 
mod.IO <- flux.IO %>% select(id,name,reaction,subsystem) %>% left_join(gpr.IO)

flux.IO.long <- flux.IO %>%
  select(id, category,ends_with(flx.plot),ends_with('mfaLB'),ends_with('mfaUB')) %>%
  pivot_longer(cols = c(-id,-category), names_to = 'flx.name', values_to = 'flx') %>%
  mutate(condition = str_extract(flx.name, '^[A-Z0-9]+'), type = str_extract(flx.name, '[a-zA-Z]+$')) %>%
  select(-flx.name) %>%
  pivot_wider(names_from = 'type', values_from = 'flx') %>%
  rename(flx = eval(flx.plot), lb = mfaLB, ub = mfaUB) %>%
  left_join(mod.IO) %>%
  mutate(organism = 'IO',strain = 'SD108') %>%
  left_join(condition.IO) %>% select(-condition)

#### flux SC ####

flux.SC <- read.csv('../SC_fluxomics/20220125_flux.csv')%>%
  mutate(id = gsub('q6','qx',id), reaction = gsub('q6','qx',reaction))%>%
  left_join(subsystem.new) %>% mutate(subsystem = ifelse(is.na(new), subsystem, new)) %>% select(-new) %>% left_join(anno.subsystem)

gpr.SC <- unnest(tibble(id = flux.SC$id,geneID=str_extract_all(flux.SC$gpr, 'Y[A-Z0-9]+')),
                 cols = c(geneID)) 

mod.SC <- flux.SC %>% select(id,name,reaction,subsystem) %>%  left_join(gpr.SC)

flux.SC.long <- flux.SC %>%
  select(id, category,ends_with(flx.plot),ends_with('mfaLB'),ends_with('mfaUB')) %>%
  pivot_longer(cols = c(-id,-category), names_to = 'flx.name', values_to = 'flx') %>%
  mutate(condition = str_extract(flx.name, '^[A-Z0-9]+'), type = str_extract(flx.name, '[a-zA-Z]+$')) %>%
  select(-flx.name) %>%
  pivot_wider(names_from = 'type', values_from = 'flx') %>%
  rename(flx = eval(flx.plot), lb = mfaLB, ub = mfaUB)%>%
  left_join(mod.SC) %>%
  mutate(organism = 'SC',strain = ifelse(grepl('CEN',condition), 'CENPK','FY4')) %>%
  left_join(condition.SC) %>% select(-condition)

# report flux
write.xlsx(list(flux.SC, flux.IO), 'report_flux.xlsx')

#### batch culture: combine flux and protein data ####
# batch culture has more proteins measured from deep proteomeomics
# average absolute 
flux.protein.batch <- rbind(flux.SC.long %>% filter(nutr == 'B') %>%
                              left_join(group.SC %>% filter(condition == 'abs_B') %>% group_by(strain, nutr, Entry,geneID) %>%
                                          summarise(f.e = sd(f,na.rm = T)/sqrt(length(f)-1),n=length(f), f = mean(f, na.rm =T)) %>%
                                          ungroup() %>% filter(n>1) %>% select(-n)),# remove protein only measured in one replicate
                            flux.IO.long %>% filter(nutr == 'B') %>%
                              left_join(group.IO %>% filter(condition == 'abs_B') %>% group_by(strain, nutr, Entry,geneID) %>%
                                          summarise(f.e = sd(f,na.rm = T)/sqrt(length(f)-1),n=length(f), f = mean(f, na.rm =T)) %>%
                                          ungroup()%>% filter(n>1) %>% select(-n)))

#### all nutrient conditions: combine flux and protein data ####

flux.protein.all <- rbind(group.IO %>% filter(!condition == 'abs_B') %>% 
                            group_by( strain, nutr, dr, Entry,geneID,group) %>%
                            summarise(f.e = sd(f)/sqrt(length(f)-1), f = mean(f)) %>% ungroup() %>%
                            inner_join(flux.IO.long),
                          group.SC %>% filter(!(condition == 'abs_B' & strain == 'FY4')) %>% 
                            group_by(strain, nutr, dr, Entry,geneID,group) %>%
                            summarise(f.e = sd(f)/sqrt(length(f)-1), f = mean(f)) %>% ungroup() %>%
                            inner_join(flux.SC.long))

save.image(file = 'flux_enzyme_correlation_20220128.RData')

#### batch flux: scatter plot ####
# between-organism comparison, aggregate subunits or isozymes for enzyme concentration
load(file = 'flux_enzyme_correlation_20220128.RData')

dt.plot.agg <- flux.protein.batch %>%
  group_by(organism,strain, dr, id,reaction,name,category,subsystem) %>%
  summarise(flx = mean(flx),lb = mean(lb), ub = mean(ub), f = sum(f,na.rm = T), f.e = sum(f.e,na.rm = T))%>% # aggregate isozyme and subunits
  ungroup() %>%
  mutate(activity = abs((lb+ub))/2/f/30) %>%
  arrange(desc(category))%>%
  mutate(category = ifelse(category %in% c('folate','sugar','lipid','amino acid','nucleic acid'), 'biosynthetic',as.character(category))) %>%
  mutate(category = factor(category, levels = names(color.subsystem)))
# activity [umole/min/mg] = flux [1000umole/60min/gDw] / [0.5*f*1000mg/gDw]
# so activity = flux/f/30, 50% biomass is protein

# plot FY4 or CEN.PK
strain.sc <- 'CENPK' # FY4 was used for chemostat experiment
mu = ifelse(strain.sc == 'FY4',0.255,0.391)

tmp <- dt.plot.agg %>% select(id,flx,lb,ub,category,strain)
min.flx = 0.01
dt.plot.wide <- tmp %>% filter(strain == strain.sc) %>% select(-strain) %>%
  full_join(tmp %>% filter(strain == 'SD108')  %>% select(-strain),
            suffix = c('.SC','.IO'), by = c('id','category')) %>%
  mutate(fil = case_when(id %in% c('NADHqx_m','NADHcplxI_c_m') ~ TRUE,
                         is.na(flx.SC) ~ FALSE,
                         is.na(flx.IO) ~ FALSE,
                         (flx.SC < min.flx) & (lb.IO < ub.SC) ~ FALSE,
                         (flx.IO < min.flx) & (lb.SC < ub.IO) ~ FALSE,
                         (flx.SC < min.flx) & (flx.IO < min.flx) ~ FALSE,
                         TRUE ~ TRUE)) %>% filter(fil) %>% select(-fil) %>%
  # is it okay to leave out reactions included in one yeast but not the other??
  mutate(flx.SC = ifelse(is.na(flx.SC), 0, flx.SC),flx.IO = ifelse(is.na(flx.IO), 0, flx.IO)) %>%
  # mutate(flx.SC = ifelse((flx.SC < 0.025) & (flx.IO > 0.1), 0, flx.SC)) %>%
  # mutate(flx.IO = ifelse((flx.IO < 0.01) & (flx.SC > 0.25), 0, flx.IO)) %>%
  mutate(text.label = case_when(flx.SC==0|flx.IO==0 ~ id,
                                abs(log(flx.SC/flx.IO) - log(mu/0.521))>log(2) ~ id,
                                category == 'transport' ~ '',
                                TRUE ~ ''))

pdf(paste0('compare_flux_batch/flux~organism_',strain.sc,'_',flx.plot,'_wo_anno.pdf'), w=5, h = 3.5)
ggplot(dt.plot.wide %>% filter(!category %in% c('transport')), # %>% filter(!category %in% c('biosynthetic','other')),
       aes(x = ifelse(flx.SC < min.flx, min.flx /4, flx.SC),
           y = ifelse(flx.IO < min.flx, min.flx /4, flx.IO), fill = category, color = category)) +
  geom_abline(intercept = log10(0.52/mu), slope = 1, color = 'grey') +
  # geom_text_repel(aes(label = text.label),size =2) +
  # geom_errorbar(aes(ymin = lb.IO, ymax = ub.IO), size = 1, alpha = 0.2) +
  # geom_errorbarh(aes(xmin = lb.SC, xmax = ub.SC), size = 1, alpha = 0.2) +
  geom_jitter(size = 3, width = 0.02, height = 0.02, pch = 21) + 
  scale_fill_manual(values = fill.subsystem) +scale_color_manual(values = color.subsystem) +
  xlab(paste0('flux in SC ',strain.sc,' (mmole/h/gDW)')) + ylab('flux in IO (mmole/h/gDW)') +
  guides(fill=guide_legend(ncol=1)) +
  scale_x_log10(limits = c(min.flx/5,100),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(min.flx/5,100),labels = scales::trans_format("log10", scales::math_format(10^.x))) + # only consider significant fluxes
  theme_classic() + theme(panel.background = element_rect(color = 'black', size = 1))
dev.off()

#### batch: in vivo enzyme turnover number ####

turnover <- dt.plot.agg %>% filter(!strain == ifelse(strain.sc=='CENPK','FY4','CENPK')) %>% 
  filter(abs(lb+ub)>1E-7) %>%
  select(-strain,-dr,-flx,-lb,-ub,-f,-f.e,-name)  %>% filter(!is.infinite(activity), activity < 1E8)%>% 
  pivot_wider(names_from ='organism',values_from = 'activity') %>%
  arrange(desc(category)) %>%
  mutate(pathway = case_when(category == 'glycolysis' ~ 'glyc',
                             category %in% c('TCA', 'ox phos') ~ 'resp',
                             TRUE ~ 'other'))
pdf(paste0('flux_enzyme_correlation/activity~organism_',strain.sc,'_',flx.plot,'.pdf'), w=5.5, h = 5)
ggplot(turnover %>% filter(pathway == 'other'),
       aes(x = SC, y = IO)) +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  geom_point(size = 2, pch=21, color='grey40') +
  geom_point(data=turnover %>% filter(pathway == 'resp'), size = 2, color='#377EB8') +
  geom_point(data=turnover %>% filter(pathway == 'glyc'), size = 2, pch=21, fill='#FFED6F', color = 'goldenrod') +
  scale_fill_manual(values = color.path, guide = 'none') +
  xlab(paste0('activity in SC ',strain.sc,'\n[umole/min/mg]')) + ylab('activity in IO\n[umole/min/mg]') +
  scale_x_log10(limits= c(3E-4,3E4),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits= c(3E-4,3E4),labels = scales::trans_format("log10", scales::math_format(10^.x))) + # only consider significant fluxes
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()
summary(lm(log10(IO)~log10(SC),turnover))
save(turnover,file='flux_enzyme_correlation/turnover_number.RData')

#### batch: balance of ATP and NADH ####
source('checkbalance.R')
flux.balance <- flux.protein.batch %>% 
  group_by(organism,strain) %>%
  do(a = check.balance(met = 'nadh',mfa = ., ci = F)) %>% tidyr::unnest(a) %>% ungroup() %>%
  left_join(flux.protein.batch %>% select(id,category) %>% distinct()) %>%
  mutate(type = ifelse(flx < 0, 'consumption','production')) %>%
  mutate(flx = abs(flx)) %>%
  arrange(flx, desc(flx)) %>%
  mutate(id =  gsub('\\_[cm]','',as.character(id)),
         id = case_when(id %in% c('AKGDH','MDH','ICDHx','PDH') ~ 'TCA',
                        id %in% c('NADHqx','ALCD2i2','GAPD','NADHcplxI') ~as.character(id),
                        TRUE ~ 'other' )) %>%
  # mutate(id = ifelse(flx > 0.08*max(flx), gsub('\\_[cm]','',as.character(id)), 'other')) %>%
  mutate(id = factor(id, levels = unique(id))) %>%
  group_by(organism, strain,type) %>%
  arrange(desc(id)) %>%
  mutate(anno = flx/2 + c(0, cumsum(flx)[-length(flx)])) %>%
  ungroup() %>% mutate(strain = factor(strain, levels = c('FY4','CENPK','SD108')))


fig2 <-  ggplot(flux.balance, aes(x=type,y=flx,fill = id)) +
  geom_col(alpha = 0.8, size = 0.5, width =0.8) +
  # geom_text(aes(y = anno,
  #               label = ifelse(id == 'other','',as.character(id))), size=2) +
  facet_grid(.~strain)+ scale_y_continuous(expand = c(0,0), limits = c(0,78)) +
  scale_fill_brewer(palette = 'Set2')+
  ylab('flux [mmol/h/gDW]') +theme_classic() +
  theme(# legend.position = 'none',
    axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 0.8, size = 10),
    axis.title.x = element_blank()) #,
# panel.border = element_rect(fill = NA, color = 'black', size = 1))

figure <- ggarrange(fig1,fig2, ncol = 1)

ggsave('flux_balance/atp_nadh.pdf', plot = figure,width=3.6,h=4)

nde <- data.frame(organism = c('SC','SC','IO','IO','IO'),
                  genotype = c('WT','nde-','WT','nde-','complex I -'),
                  gr = c(1,0.37/0.4,1,0.255/0.4,0.129/0.4), 
                  gr.e = c(NA,NA, 0.028/0.4,0.013/0.4,0.013/0.4))%>% mutate(strain = paste(organism, genotype, sep = ' '))
# data for SC nde mutant is from Luttik, ... Pronk et al, JBC, 1998,
# 'The Saccharomyces cerevisiae NDE1 and NDE2 Genes Encode Separate Mitochondrial NADH Dehydrogenases Catalyzing the Oxidation of Cytosolic NADH'
# data for IO nde and complex I is from 2020/8/25

pdf('flux_balance/nadh_dehydrogenase_mutants.pdf', w = 3.6, h = 1.8)
ggplot(nde, aes(x = strain, y = gr, ymin = gr-gr.e, ymax = gr+gr.e, fill = genotype)) +
  geom_col(width = 0.6,  color = 'black', position = position_dodge(width = 0.6))+
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6)) + ylab('relative growth rate')+
  scale_x_discrete(limits = nde$strain)+ scale_y_continuous(limits = c(0,1.2), breaks = c(0,0.5,1), expand = c(0,0)) + theme_classic()+
  scale_fill_manual(values = c('WT' = 'grey','nde-' = 'darkcyan','complex I -' = 'dodgerblue4'))
dev.off()
t.test2(0.255,0.4,0.013,0.028,2,2)
q1 = (0.4-0.255)/(0.013 + 0.028) # Z score of nde deletion
q1 = (0.4-0.129)/(0.013 + 0.028) # Z score of complex I deletion
p = 2*pnorm(q, lower.tail = FALSE)

#### biomass vs energy carbon partitioning ####
rxn.stoi <- rbind(read.xlsx('flux_balance/biomass_energy_rxn_stoi.xlsx') %>%
                    mutate(precursor = factor(precursor, levels = precursor)) %>%
                    mutate(reaction.stoi = paste0('+ ',reaction.stoi)) %>%
                    mutate(id = str_extract_all(reaction.stoi, '[+-]{1} [a-zA-Z0-9_]+\\_[cme]')) %>%
                    unnest(cols = id) %>% mutate(dir = ifelse(grepl('\\+',id), 1,-1),
                                                 id = gsub(' ','',gsub('^[+-]{1}','',id))) %>% select(-reaction.stoi),
                  data.frame(precursor = '3pg/pep/pyr',nC = 1.5,type = 'biomass',id='FDH_c',dir=-1))
# in S.c, some ser overflow to formate, pyr -> ser -> 2 formate + CO2, thus 1 form = 1.5 central carbon

flux.partition.detail <- rbind(flux.IO.long,flux.SC.long) %>% select(organism,strain,nutr,dr,id,flx) %>% distinct() %>%
  right_join(rxn.stoi) %>%
  group_by(organism,strain,nutr,dr,type, precursor, nC) %>%
  summarise(flx = sum(dir *flx)) %>% ungroup() %>% mutate(flxC = flx * nC) 

pdf('flux_balance/partition_biomass_energy_detail.pdf', w = 6, h = 4)
ggplot(flux.partition.detail %>% filter(nutr == 'B', !grepl('FY4',strain)),
       aes(x = type, y = flxC, fill = precursor)) +
  geom_col(position = position_stack(), width = 0.7, color = 'black', size = 0.25) +
  facet_wrap(organism~., scales = 'free') +
  scale_fill_brewer(palette = 'Set3') +
  scale_y_continuous(expand = c(0,0), limits = c(0,150)) +
  egg::theme_presentation() +theme(text = element_text(size = 20))
dev.off()

flux.partition <- flux.partition.detail %>%
  mutate(precursor.type = case_when(precursor == 'co2_m' ~ 'respiration',
                                    precursor == 'etoh' ~ 'glycolysis',
                                    TRUE ~ 'biomass')) %>%
  group_by(organism, strain, nutr, dr, type, precursor.type) %>% summarise(flxC = sum(flxC)) %>% ungroup()

pdf('flux_balance/partition_biomass_energy.pdf', w = 6, h = 2.5)
ggplot(flux.partition %>% filter(nutr == 'B', !grepl('FY4',strain)),
       aes(x = type, y = flxC, fill = precursor.type)) +
  geom_col(position = position_stack(), width = 0.6, color = 'black', size = 0.25, alpha = 0.9) +
  facet_wrap(organism~., scales = 'free') +
  scale_fill_manual(values = c('biomass' = "grey50", 'glycolysis' = "#FF3333",'respiration'="#339933")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,150)) +
  egg::theme_presentation() +
  theme(text = element_text(size = 20), axis.ticks.x = element_blank(),panel.background = element_blank(),
        axis.text.x = element_blank(),axis.title.x = element_blank())
dev.off()

#### Compare S.c and I.o in flux space ####
load(file = 'flux_enzyme_correlation_20220128.RData')
tmp <- flux.protein.all %>% filter(id %in% c('PYK_c','G6PDH2i_c','CS_m')) %>%
  select(organism,strain,nutr,dr,id,flx) %>% distinct() %>%
  pivot_wider(names_from = 'id', values_from = 'flx')


fig1 <- ggplot(tmp, aes(x = PYK_c/dr, y = CS_m/dr, color = organism, shape = nutr)) + 
  geom_point() + #stat_ellipse(size = 0.5) +
  theme_classic() + scale_color_manual(values = c('SC' = 'grey50', 'IO' = 'steelblue2')) + 
  scale_shape_manual(values = c('B' = 1, 'C' = 16, 'N' = 0, 'P' = 2))+
  xlim(0,NA) + ylim(0,NA)
fig2 <- ggplot(tmp, aes(x = PYK_c/dr, y = G6PDH2i_c/dr, color = organism, shape = nutr)) + 
  geom_point() + #stat_ellipse(size = 0.5) +
  theme_classic() + scale_color_manual(values = c('SC' = 'grey50', 'IO' = 'steelblue2')) + 
  scale_shape_manual(values = c('B' = 1, 'C' = 16, 'N' = 0, 'P' = 2))+
  xlim(0,NA) + ylim(0,NA)
figure <- ggarrange(fig1,fig2, ncol = 1)
ggsave('flux_balance/TCA~PPP~glycolysis.pdf', figure,width = 3,height=4)


#### regulation analysis ####
load(file = 'flux_enzyme_correlation_20220128.RData')
yeast = 'SC'
n = ifelse(yeast=='IO',5,4)
flux.protein <- flux.protein.all %>% filter(strain == ifelse(yeast == 'SC','FY4','SD108'))

## plot growth rate
omit0 <- function(x) {gsub('^0\\.','.',x)}
tmp <- flux.protein  %>% left_join(rbind(condition.SC,condition.IO)) %>% select(condition, nutr, dr)%>%
  distinct() %>% mutate(condition = factor(condition, levels = condition))
pdf(paste0('other/annotation_nutrient_conditions_',yeast,'.pdf'), w = 5, h = 1.8)
ggplot(tmp, aes(x = condition, y = dr, fill = nutr)) + geom_col(width = 0.8) +
  ylab(TeX('DR / $h^{-1}$')) +  theme_classic() +
  scale_x_discrete(limits =  c(paste0('C',1:n),paste0('N',1:n),paste0('P',1:n),'B'))+
  scale_y_continuous(limits = c(0,ifelse(yeast == 'SC',0.35,0.55)),expand = c(0,0), breaks = c(0,0.2,0.4),  labels = omit0)+
  theme_presentation() +theme(axis.ticks.x = element_blank())
dev.off()

## how much of flux variation is explained by growth rate?
dt.plot.lm <- flux.protein %>%
  select(id, nutr, dr,flx) %>% distinct() %>% group_by(id) %>%filter(max(flx)>0.0001) %>%
  do(anova(lm(flx ~ dr, data  = .))) 

dt.anova <- dt.plot.lm %>%
  mutate(variance = ifelse(is.na(`F value`), 'unexplained', 'explained')) %>%
  mutate(fraction = `Sum Sq`/sum(`Sum Sq`)) %>% ungroup() %>%
  group_by(variance) %>% summarise(fraction = mean(fraction)) %>% ungroup()

message(paste0('variance in flux change explained by growth rate ',format(dt.anova$fraction[1], digits = 2)))

# pseudo regulation analysis
# flux = [enzyme] * (enzyme occupancy), log(flux) = log[enzyme] + log(enzyme occupancy)
# d(log(enzyme))/d(log(flux)) + d(log(enzyme occupancy))/d(log(flux)) = 1

fil <- flux.protein %>% group_by(category, id, geneID) %>% filter(!flx == 0) %>% filter(!is.na(f)) %>%
  summarise(n=length(flx)) %>% ungroup() %>%
  filter(n>3) %>% # fluxes in at least 3 conditions are determined
  anti_join(dt.plot.lm %>% filter(is.nan(p)) %>% select(id) %>% distinct())

# do not aggregate enzyme groups (to see isozymes)
dt.regulation.plot <- flux.protein  %>%
  right_join(fil %>% select(-n)) %>% filter(!flx == 0,!f==0) %>%
  group_by(category, id, geneID, subsystem)%>%
  mutate(flx.norm = flx/mean(flx), f.norm = f/mean(f)) %>% ungroup()

# calculate regulation coefficient rho_e
dt.regulation <- dt.regulation.plot %>% group_by(category, id, geneID, subsystem) %>%
  summarise(f = mean(f), flx = exp(mean(log(flx))),
            rho_e = lin.coeff(log(f.norm), log(flx.norm)),
            rho_e.p = lin.coeff.p(log(f.norm), log(flx.norm))) %>% ungroup() %>%
  filter(!is.na(rho_e)) %>%
  mutate(category = factor(category, levels = names(color.subsystem))) %>%
  mutate(rho_e.p = p.adjust(rho_e.p, method = 'BH'))

# volcano plot
pdf(paste0('regulation_analysis/volcano_plot_',yeast,'_w_anno.pdf'), w = 6.2,h=2.2)
ggplot(dt.regulation, aes(y = -log10(rho_e.p), x = rho_e, size = f)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(aes(color = category), alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(rho_e > 0.5 & rho_e.p < 0.05,id,''))) +
  scale_color_manual(values = color.subsystem) +
  guides(size = guide_legend(order = 1), color = guide_legend(ncol = 2, order = 2))+
  xlim(-1.5,1.5) + ylab(TeX('-log_{10}p')) + xlab(TeX('enzyme regulation $\\rho_{E}$'))  +
  scale_size_continuous(limits = c(0, 0.05), breaks = c(0.0,0.002,0.01,0.02)) + theme_classic() +
  theme(legend.direction = "vertical", legend.box = "horizontal")
dev.off()

# regulation coefficient of both enzyme and metabolite for the 50 reactions from SIMMER (I.orientalis)
regAnal.df <- read.delim('../../result/regAnalysis.tsv', sep= '\t') %>% rename(id = rxnid) %>%
  left_join(mod.IO %>% select(id, category) %>% distinct()) %>% mutate(category = factor(category, levels = names(color.subsystem))) %>% arrange(desc(category))
# ggplot(regAnal.df %>% pivot_longer(cols = c(rho_e, rho_m), names_to = 'loc', values_to = 'rho'), aes(x = loc, y = rho))+
#   geom_violin()
pdf('regulation_analysis/regulation_coefficient_SIMMER.pdf', w = 4.0, h = 2.2)
ggplot(regAnal.df, aes(x = rho_e, y = rho_m, color = category)) +
  geom_abline(slope = 1, color = 'grey') +  geom_abline(intercept = 1, slope = -1, color = 'grey', linetype = 'dashed') + 
  # geom_hline(yintercept = 0, color = 'black') + geom_vline(xintercept= 0, color = 'black') +
  geom_jitter(pch = 21,size = 2, width = 0.05, height = 0.05, alpha = 0.9) +
  scale_y_continuous(limits = c(-1.2,2.1), expand = c(0,0))+
  scale_x_continuous(limits = c(-1.2,2.1), breaks = c(-1,0,1,2), expand = c(0,0)) +
  scale_color_manual(values = color.subsystem) +
  xlab(TeX('enzyme regulation $\\rho _{E}$')) + ylab(TeX('metabolite regulation $\\rho _{M}$'))+
  #annotate('label', label = TeX(' $\\rho _{E} + \\rho _{M} = 1$'), x = 1.5, y = -0.5, size = 3)+
  theme_classic() + theme(panel.background = element_rect(fill = NA))
dev.off()


#### does enzyme change predict flux change within or cross yeasts? ####
load(file = 'flux_enzyme_correlation_20220128.RData')
dt.all <- flux.protein.all %>%
  group_by(id,reaction,name,category,subsystem, dr, nutr,organism,strain) %>%
  summarise(flx = mean(flx),lb = mean(lb), ub = mean(ub), f = sum(f, na.rm = T)) %>% # aggregate isozyme and subunits
  ungroup()
save.image(file = 'flux_enzyme_correlation_20220128.RData')
## which yeast to plot flux correlation with enzyme across nutrient?
ref.yeast <- 'SC'
avg.pathway <- 'ratio first' # 'ratio first', calculate FC of enzyme and flux, then average within pathway; 'sum first', sum up enzyme first then calculate FC
## comparison with reference to the ref.yeast
dt.plot <- dt.all %>% filter(organism == ref.yeast, strain == ifelse(ref.yeast == 'IO','SD108','FY4')) %>%
  group_by(id) %>% mutate(flxFC = flx/mean(flx),enzFC = f/mean(f)) %>% ungroup() %>%
  mutate(category = factor(category, levels = rev(names(color.subsystem)))) %>%
  filter(!category %in% c('transport','other'),
         flxFC > 0)
if(grepl('ratio',avg.pathway)) {
  dt.mean <- dt.plot %>% group_by(category, organism, dr, nutr) %>%
    summarise(n = sum(!is.na(enzFC)), flxFC = median(flxFC,na.rm=T), 
              enzFC = median(enzFC,na.rm=T),f = sum(f, na.rm = T)) %>% ungroup() 
} else {
  dt.mean <- dt.plot %>% group_by(category, organism, dr, nutr) %>%
    summarise(n = sum(!is.na(enzFC)), flxFC = mean(flx, na.rm=T)/mean(flx.B,na.rm=T), f = sum(f, na.rm = T), enzFC = f/sum(f.B,na.rm=T)) %>% ungroup()
}

if(ref.yeast == 'IO') {fp.within.IO <- dt.mean} else{fp.within.SC <- dt.mean}
# within organism comparison
fit <- dt.mean %>% filter(enzFC>0 & flxFC>0 & is.finite(enzFC) & is.finite(flxFC)) %>%
  summarize(rsq = lin.rsq(log(flxFC),log(enzFC)),
            slope = lin.coeff(log(flxFC),log(enzFC)),
            lin.p = lin.coeff.p(log(flxFC),log(enzFC),intercept = 0),
            P.R = cor(log(flxFC),log(enzFC), method = 'pearson'))
# variance explained (assuming y = x model)
var(log(dt.mean$enzFC))/var(log(dt.mean$flxFC))
cov(log(dt.mean$enzFC),log(dt.mean$flxFC))/var(log(dt.mean$flxFC))
cov(log(dt.mean$enzFC),log(dt.mean$flxFC))/sqrt(var(log(dt.mean$flxFC))*var(log(dt.mean$enzFC)))
1-var(log(dt.mean$enzFC)-log(dt.mean$flxFC))/var(log(dt.mean$flxFC))

pdf(paste0('flux_enzyme_correlation/enzymeFC~fluxFC_scatter_within_organism_',
           ref.yeast,'_',avg.pathway,'_w_anno.pdf'), w=4, h = 4)

ggplot(dt.mean%>% arrange(category), aes(x = enzFC, y = flxFC, color = category)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed')+
  geom_jitter(pch = 21, height = 0.1, size = 3) + 
  scale_color_manual(values = color.subsystem) +
  scale_x_log10(limits = c(0.05,10), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(limits = c(0.05,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
dev.off()

dens.x <- ggplot(dt.plot, aes(x = enzFC, fill = case_when(category == 'glycolysis' ~ 'glycolysis',
                                                          category %in% c('TCA','ox phos') ~ 'respiration',
                                                          TRUE ~ 'biosynthetic'))) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + scale_fill_manual(values = color.subsystem)+
  scale_x_log10(limits = c(4E-2,10), breaks = c(0.1,0.3,1,3))

dens.y <- ggplot(dt.plot, aes(x = flxFC, fill =  case_when(category == 'glycolysis' ~ 'glycolysis',
                                                           category %in% c('TCA','ox phos') ~ 'respiration',
                                                           TRUE ~ 'biosynthetic'))) + 
  geom_density(alpha = 0.4) + 
  scale_x_log10(limits = c(4E-2,10), breaks = c(0.1,0.3,1,3))+
  theme_void() + 
  theme(legend.position = "none") + scale_fill_manual(values = color.subsystem) + coord_flip()

dens.x + patchwork::plot_spacer() + fig + dens.y + 
  patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

# between organism comparison
avg.pathway <- 'ratio first'
dt.plot <- rbind(dt.all %>% filter(nutr == 'B' | dr==0.23, organism == 'IO') %>% mutate(dr = 0.22) ,
                 dt.all %>% filter(nutr == 'B' | dr==0.22, organism == 'SC')) %>%
  group_by(id,category, reaction, name, nutr, organism) %>% summarise(flx = mean(flx), f = mean(f)) %>% ungroup() %>%
  # take the average of two S.c strains from batch culture
  group_by(id, nutr) %>% mutate(flxFC = flx/mean(flx), enzFC = f/mean(f)) %>% ungroup() %>%
  mutate(category = factor(category, levels = rev(names(color.subsystem)))) %>%
  filter(!category %in% c('transport','other'),flxFC>0)

if(grepl('ratio',avg.pathway)) {
  dt.mean <- dt.plot %>% group_by(category, organism, nutr) %>%
    summarise(n = sum(!is.na(enzFC)), flxFC = median(flxFC,na.rm = T), enzFC = median(enzFC,na.rm = T), 
              f = sum(f, na.rm = T)) %>% ungroup()
}else {
  dt.mean <- dt.plot %>% group_by(category, organism,nutr) %>%
    summarise(n = sum(!is.na(enzFC)), flxFC = mean(flx,na.rm = T)/mean(flx.IO,na.rm = T), 
              f = sum(f, na.rm = T), enzFC = f/sum(f.IO,na.rm=T)) %>% ungroup()
}
fp.btw <- dt.mean
fit <- dt.mean %>% filter(enzFC>0 & flxFC>0 & is.finite(enzFC) & is.finite(flxFC)) %>%
  summarize(rsq = lin.rsq(log(flxFC),log(enzFC)),
            slope = lin.coeff(log(flxFC),log(enzFC)),
            lin.p = lin.coeff.p(log(flxFC),log(enzFC)),
            P.R = cor(log(flxFC),log(enzFC), method = 'pearson'))

# variance explained (assuming y = x model)
var(log(dt.mean$enzFC))/var(log(dt.mean$flxFC))
var(log(dt.mean$enzFC))/var(log(dt.mean$flxFC))
cov(log(dt.mean$enzFC),log(dt.mean$flxFC))/var(log(dt.mean$flxFC))
cov(log(dt.mean$enzFC),log(dt.mean$flxFC))/sqrt(var(log(dt.mean$flxFC))*var(log(dt.mean$enzFC)))
1-var(log(dt.mean$enzFC)-log(dt.mean$flxFC))/var(log(dt.mean$flxFC))

pdf(paste0('flux_enzyme_correlation/enzymeFC~fluxFC_scatter_cross_organism','_',avg.pathway,'.pdf'), w = 4, h = 4)
ggplot(dt.mean,
       aes(x = enzFC, y = flxFC, color = category, shape = organism)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed')+
  geom_jitter( height = 0.1, size = 3) + 
  scale_color_manual(values = color.subsystem) +
  scale_x_log10(limits = c(0.05,10), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(limits = c(0.05,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  scale_shape_manual(values = c('SC' = 21, 'IO'=24)) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
dev.off()
save.image(file = 'flux_enzyme_correlation_20220128.RData')

#### what about metabolites? ####
load(file = 'flux_enzyme_correlation_20220128.RData')

mod <- mod.IO %>% select(-geneID) %>% distinct() %>%
  mutate(subs = unlist(lapply(strsplit(reaction, split = '>'),'[',1))) 

mod.subs <- unnest(tibble(id=mod$id, subs = str_extract_all(mod$subs, '[a-zA-Z0-9_]+\\_[cem]')), cols = c(subs))%>%
  left_join(mod %>% select(-subs))%>% mutate(BiGG = gsub('\\_[cme]','',subs))
# write.table(mod.subs %>% filter(!grepl('[pP]suedo',subsystem)) %>% select(BiGG, subsystem,category) %>% distinct(),
#             '../../../../DATA_MS/id_conversion/met_pathway.tsv',sep = '\t',quote = F,row.names = F)
# parse if all the metabolites are in genome scale model
# mod_met <- read.xlsx('../IO_fluxomics/20220128_ioMFA_GAM216/imIsor_GEM.xlsx', sheet = 'Metabolites')
# mod.subs$subs[which(!mod.subs$subs %in% mod_met$metid)]

met.IO <- read.csv('../IO_metabolomics/data_met_rel_20200724.csv') %>% filter(!is.na(BiGG)) %>% select(-chebi,-source,-metname) %>%
  select(contains('mean')|contains('BiGG')) %>% 
  mutate(B=B_mean) %>%
  pivot_longer(cols = contains('mean'), names_to = 'condition', values_to = 'log2FC') %>%
  mutate(log2FC = log2FC-B) %>% select(-B) %>%
  mutate(strain = 'SD108',yeast = 'IO', condition = gsub('\\_mean','',condition)) %>% left_join(condition.IO) %>% filter(!is.na(log2FC))

met.SC <- read.delim('../SC_metabolomics/SC_chemostat_idmatched.tsv',sep = '\t') %>% filter(!is.na(BiGG)) %>%
  filter(!duplicated(compound)) %>%
  select(-compound,-KEGG) %>%
  select(contains('mean')|contains('BiGG')) %>% 
  pivot_longer(cols = contains('mean'), names_to = 'condition', values_to = 'log2FC') %>%
  mutate(strain = 'FY4',yeast = 'SC', condition = gsub('\\_mean','',condition)) %>% left_join(condition.SC) %>% filter(!is.na(log2FC))

#### does change in metabolites predict flux change within same yeast? ####
yeast = 'SC'

met.within <-  mod.subs %>% left_join(anno.subsystem) %>% select(BiGG, category) %>% distinct() %>%
  inner_join(get(paste0('met.',yeast))%>%group_by(BiGG) %>% mutate(log2FC = log2FC - mean(log2FC)) %>% ungroup() )

fpm.within <- met.within %>% group_by(dr, nutr,category) %>%
  summarise(metFC = 2^median(log2FC)) %>% ungroup() %>%
  left_join(get(paste0('fp.within.',yeast))) %>% filter(!category %in% c('transport','other')) %>%
  mutate(category = factor(category, levels = rev(names(color.subsystem))))

fit <- fpm.within %>% filter(is.finite(metFC) & is.finite(flxFC)) %>%
  summarize(rsq = lin.rsq(log(flxFC),log(metFC)),
            slope = lin.coeff(log(flxFC),log(metFC)),
            lin.p = lin.coeff.p(log(flxFC),log(metFC)),
            P.R = cor(log(flxFC),log(metFC), method = 'pearson'))
var(log(fpm.within$metFC))/var(log(fpm.within$flxFC))
cov(log(fpm.within$metFC),log(fpm.within$flxFC))/var(log(fpm.within$flxFC))
cov(log(fpm.within$metFC),log(fpm.within$flxFC))/sqrt(var(log(fpm.within$flxFC))*var(log(fpm.within$metFC)))
1-var(log(fpm.within$metFC)-log(fpm.within$flxFC))/var(log(fpm.within$flxFC))

pdf(paste0('flux_enzyme_correlation/metFC~fluxFC_scatter_within_organism_',yeast,'_',avg.pathway,'.pdf'), 
    w=4, h = 4)

ggplot(fpm.within %>% arrange(category), aes(x = metFC, y = flxFC, color = category)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed')+
  geom_jitter(pch = 21, height = 0.1, size = 3) +
  scale_color_manual(values = color.subsystem) +
  scale_x_log10(limits = c(0.05,10), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(0.05,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
dev.off()

#### does change in metabolites predict flux change cross different yeast? ####
met.ratio.B <- read.delim('../SC_metabolomics/20210407-all_met_conc_log2ratio_IO_SC_id_matched.tsv',sep = '\t') %>% filter(!is.na(BiGG)) %>%
  inner_join(mod.subs %>% left_join(anno.subsystem) %>% select(BiGG, category) %>% distinct()) %>%
  mutate(dr = 0.255, nutr = 'B',log2FC = -logr) %>% select(-logr) # logr was calculated as log2(IO/SC)

met.ratio.B.CENPK <- read.delim('../SC_metabolomics/20220216-all_met_conc_ratio_IO_SC_CENPK.tsv',sep = '\t') %>% filter(!is.na(BiGG)) %>%
  inner_join(mod.subs%>% left_join(anno.subsystem)  %>% select(BiGG, category) %>% distinct()) %>%
  mutate(dr = 0.391, nutr = 'B',log2FC = -logr) %>% select(-logr) # logr was calculated as log2(IO/SC)

# first calculate log2FC to IO_B
met.btw <- rbind(met.ratio.B.CENPK %>% select(BiGG, category, nutr, log2FC,dr),
                 met.ratio.B %>% rename(B = log2FC) %>% select(BiGG,B,category) %>%
                   inner_join(met.SC %>% filter(dr == 0.22|nutr == 'B') %>% select(BiGG,nutr,log2FC) %>% rename(SC = log2FC)) %>% # FC respect to SC.FY4.B
                   inner_join(met.IO %>% filter(dr == 0.23|nutr == 'B') %>% select(BiGG,nutr,log2FC) %>% rename(IO = log2FC)) %>%
                   mutate(log2FC = SC-IO-B) %>% select(-SC,-IO,-B) %>%
                   pivot_wider(names_from = 'nutr',values_from = 'log2FC') %>%
                   pivot_longer(cols = c(B,C,N,P), names_to = 'nutr', values_to = 'log2FC') %>% mutate(dr = ifelse(nutr == 'B',0.255,0.22)))

met.btw <- rbind(met.ratio.B.CENPK %>% select(BiGG, category, nutr, log2FC,dr) %>% mutate(organism = 'SC'), # log2FC(SC_CENPK_B/IO_B)
                 met.ratio.B %>% rename(B = log2FC) %>% select(BiGG,B,category) %>%
                   inner_join(met.SC %>% filter(dr == 0.22|nutr == 'B') %>% select(BiGG,nutr,log2FC,dr)) %>% # log2FC(SC_x/SC_FY_B)
                   mutate(log2FC = log2FC+B) %>% select(-B)%>% mutate(organism = 'SC'), # log2FC(SC_x/SC_FY_B) + log2FC(SC_FY_B/IO_B)
                 met.ratio.B %>% select(BiGG, category) %>% distinct() %>%
                   inner_join(met.IO %>% filter(dr == 0.23|nutr == 'B') %>% select(BiGG, nutr,log2FC,dr)) %>% # log2FC(IO_x/IO_B)
                   mutate(organism = 'IO')) %>%
  group_by(BiGG, category,nutr, organism) %>% summarise(log2FC = mean(log2FC)) %>% ungroup() %>% 
  # take the average of two S.c strains from batch culture
  group_by(BiGG, category,nutr) %>% mutate(log2FC = log2FC-mean(log2FC)) %>% ungroup()


fpm.btw <- met.btw %>% group_by(nutr,category,organism) %>%
  summarise(metFC = 2^median(log2FC)) %>% ungroup() %>%
  left_join(fp.btw) %>% filter(!category %in% c('transport','other')) %>%
  mutate(category = factor(category, levels = rev(names(color.subsystem))))

fit <- fpm.btw %>% filter(is.finite(enzFC) & is.finite(flxFC)) %>%
  summarize(rsq = lin.rsq(log(flxFC),log(metFC)),
            slope = lin.coeff(log(flxFC),log(metFC)),
            lin.p = lin.coeff.p(log(flxFC),log(metFC)),
            P.R = cor(log(flxFC),log(metFC), method = 'pearson'))
var(log(fpm.btw$metFC))/var(log(fpm.btw$flxFC))
cov(log(fpm.btw$metFC),log(fpm.btw$flxFC))/var(log(fpm.btw$flxFC))
cov(log(fpm.btw$metFC),log(fpm.btw$flxFC))/sqrt(var(log(fpm.btw$flxFC))*var(log(fpm.btw$metFC)))
1-var(log(fpm.btw$metFC)-log(fpm.btw$flxFC))/var(log(fpm.btw$flxFC))
####### multivariate analysis? #######
{
  library(nlmrt)
  tmp <- fpm.within %>% filter(enzFC>0 & flxFC>0 & is.finite(enzFC) & is.finite(flxFC))
  fm <- nlxb(log(flxFC) ~ (b1 * log(enzFC) + b2 * log(metFC))/(b1+b2) + b3,
             data = tmp,
             upper = c(1,500,500),
             lower = c(0,0,-1),
             start = list(b1 = 1, b2 = 1, b3 =0))
}

pdf(paste0('flux_enzyme_correlation/metFC~fluxFC_scatter_cross_organism_',avg.pathway,'.pdf'),
    w=4, h = 4)
ggplot(fpm.btw, # %>% filter(!dr == '0.22'),
       aes(x = metFC, y = flxFC, color = category, shape = organism)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed')+
  geom_jitter( height = 0.1, size = 3) + 
  scale_color_manual(values = color.subsystem) +
  scale_x_log10(limits = c(0.05,10), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(limits = c(0.05,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  scale_shape_manual(values = c('SC' = 21, 'IO'=24)) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
dev.off()


#### reserve capacity in individual enzymes in S.c ####
# which yeast to plot flux correlation with enzyme across nutrient?
ref.yeast <- 'SC'
## comparison with reference to the ref.yeast
dt.plot <- dt.all %>% 
  inner_join(dt.all %>% filter(organism == ref.yeast, strain == 'CENPK', nutr == 'B') %>%
               rename(flx.B = flx, f.B = f) %>% select(id, flx.B,f.B)) %>%
  mutate(category = factor(category, levels = rev(names(color.subsystem)))) %>%
  # mutate(flxFC = flx/flx.B) 
  group_by(category, organism, strain,dr, nutr) %>%
  summarise(flxFC = median(flx/flx.B,na.rm = T), 
            f = sum(f,na.rm=T),f.B = sum(f.B,na.rm = T),
            enzFC = f/f.B) %>% ungroup()

dt.plot.category <- dt.plot %>%
  filter(dr > 0.24 & !dr==0.391, category == 'glycolysis')

pdf('flux_enzyme_correlation/reserve_capacity_in_glycolysis_SC.pdf', w = 2.8, h = 2)
ggplot(dt.plot.category, aes(x = flxFC, y = f, color = organism)) + 
  geom_point(pch=21)+
  xlim(0,1.5) + ylim(0,0.3)+
  geom_abline(slope = 0.25) +
  theme_classic()
dev.off()


#### does a metabolite level correlate with flux through a reaction? ####
tmp.met <- met.ratio.B %>% rename(ref_B = log2FC) %>% select(BiGG,ref_B,category) %>%
  inner_join(rbind(met.SC %>% select(BiGG,nutr,dr, log2FC,yeast),met.IO %>% select(BiGG,nutr,dr,log2FC,yeast))) %>%
  mutate(log2FC = ifelse(yeast == 'SC', log2FC - ref_B, log2FC)) %>%
  filter(BiGG == 'pyr')

tmp.flx <- flux.protein.all %>% filter(grepl('PYRDC',id)) %>% rename(yeast = organism) %>%
  select(yeast, flx, nutr,dr) %>% distinct()

tmp <- tmp.met %>% full_join(tmp.flx)

ggplot(tmp, aes(x = log2FC, y = log2(flx), color = nutr, shape = yeast)) + geom_point()
