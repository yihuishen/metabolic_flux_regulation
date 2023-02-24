library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(tidyr)
library(openxlsx)
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = 'proteome_efficiency_ATP_pathway.RData')

# load NCI-60 data
load('../NCI-60/NCI60.RData')
dt.NCI60.median <- dt.NCI60 %>% mutate(pathway = ifelse(pathway == 'mitochondria', 'respiration', pathway)) %>% rename(flux.ATP = flux) %>%
  select(pathway, f, flux.ATP, anno) %>%
  mutate(organism = 'NCI60', note = 'average')%>%
  group_by(note, anno,organism,pathway) %>%
  summarise(f.e = sd(f), f = median(f),  flux.ATP = median(flux.ATP)) %>% ungroup() %>%
  mutate(lb = quantile(flux.ATP,0.025), ub = quantile(flux.ATP,0.975))
dt.NCI60.by.tissue <- dt.NCI60 %>% group_by(tumor.type,anno,pathway) %>% summarise(f = mean(f), flux.ATP = mean(flux)) %>% ungroup() %>%
  rename(note = tumor.type) %>% mutate(organism = 'cancer cell lines')

# load yeast data
load('../manuscript/ATP/ATP_all_conditions.RData')
dt.yeast <- dt.ATP %>% # filter((nutr == 'B' &( ! strain == 'FY4'))| (nutr == 'C' & dr == 0.28)) %>%
  distinct() %>% rename(flux.ATP = flux) %>% mutate(organism = paste0('yeast_',organism)) %>%
  mutate(note = paste0(strain,'_',nutr,'_',dr))%>%
  select(organism, pathway, note,f, flux.ATP, anno, f.e, lb,ub,protein)

# load E.coli data
dt.ecoli <- read.xlsx('../Ecoli/integrated.xlsx', sheet = 'summary', cols = 1:5, rows = 20:24) %>%
  pivot_longer(cols = c(value, stdev), names_to = 'stat', values_to = 'value') %>%
  mutate(parameter = paste0(parameter,ifelse(stat == 'value','','.e'))) %>%
  select(-stat) %>%
  pivot_wider(values_from = 'value', names_from = 'parameter')%>%
  mutate(organism = 'E.coli',note = 'average') %>%
  mutate(lb = flux.ATP-1.96*flux.ATP.e, ub = flux.ATP+1.96*flux.ATP.e) %>% select(-flux.ATP.e)

# load T cell data
load('../manuscript/flux_enzyme_correlation_T_cell/Tcell_ATP.RData')
dt.Tcell <- dt.Tcell.ATP %>% mutate(organism = 'T cell', anno = substr(pathway, 1,4)) %>%
  rename(protein = `%protein`, note = tissue)

# load mouse data
load('ATP_tissue_tumor.RData')
dt.mouse <- dt.mouse.ATP %>% mutate(organism = 'mouse_paxdb', anno = substr(pathway, 1,4)) %>%
  rename(protein = `%protein`) %>% select(-flux) %>% mutate(protein = protein/100)

dt.tumor <- dt.tumor.ATP %>% mutate(organism = paste0('mouse_',tissue.type), anno = substr(pathway, 1,4)) %>%
  rename(protein = `%protein`, note = tissue) %>% select(-flux, -tissue.type) %>%
  mutate(protein = protein/100)

# combine
# order.tissue.new <- order.tissue %>% rename(tissue.paxDB = tissue) %>%
#   inner_join(tissue.mapping) %>%filter(!is.na(tissue)) %>%
#   select(tissue) %>% distinct()

proteome.efficiency <- rbind(dt.yeast,
                             dt.Tcell,
                             dt.tumor,
                             dt.mouse,
                             rbind(dt.NCI60 %>% mutate(lb=NA,ub=NA,f.e = NA, organism = 'cancer cell lines') %>%
                                     rename(flux.ATP = flux, note = cell) %>% select(-enzyme,-tumor.type),
                                   dt.ecoli) %>% mutate(protein = 0.5))%>% 
  mutate(efficiency = flux.ATP/f/protein, efficiency.lb = lb/protein/(f+1.96*f.e), efficiency.ub = ub/protein/(f-1.96*f.e)) 

write.xlsx(proteome.efficiency,'proteome_efficiency_across_organisms_20230216.xlsx')
rm(list = c('dt.list','mapping','mapping.tmp','dt','dt.anno','annotation'))
save.image(file = 'proteome_efficiency_ATP_pathway.RData')

##### plot yeast, T cell, tumor #####
load(file = 'proteome_efficiency_ATP_pathway.RData')
color.path <- c('glyc' = '#FFED6F', 'mito' = '#BC80BD','resp' = '#377EB8')
theme_bar <- theme_classic() + theme(line = element_line(size = 0.5),
                                     axis.text.y = element_text(size = 20),axis.title.y = element_text(size=24),
                                     axis.text.x = element_blank(),axis.title.x = element_blank(),
                                     legend.text = element_blank(), legend.title = element_blank(),
                                     plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))

organism.plot <- 'yeast' # yeast, T cell, mouse_spleen, mouse_pancreas
anno.plot <- 'resp' # resp or whole mito
eff.plot <- proteome.efficiency %>% filter(grepl(organism.plot,organism)) %>%
  filter(anno %in% c('glyc',anno.plot))
if (organism.plot == 'yeast') {
  eff.plot <- eff.plot %>%
    filter(grepl('B',note), !grepl('FY4',note))
  unit.f <- '[mmol/h/gDW]'
  unit.p <- '[gProtein/gDW]'
  unit.e <- '[mmol/h/gProtein]'
  x.limit <- c('SD108_B_0.52','CENPK_B_0.391')
} else if (organism.plot == 'T cell') {
  unit.f <- '[mmol/h/gDW]'
  unit.p <- '[gProtein/gDW]'
  unit.e <- '[mmol/h/gProtein]'
  x.limit <- c('naive','activated')
} else {
  unit.f <- '[mmol/h/gTissue]'
  unit.p <- '[gProtein/gTissue]'
  unit.e <- '[mmol/h/gTissue]'
  if (organism.plot == 'mouse_spleen') {
    x.limit <- c('spleen','leukemic spleen')
  } else { x.limit <- c('pancreas','GEMM PDAC','flank PDAC')}
  
}

fig1y <- ggplot(eff.plot, aes(x = note, y = flux.ATP, fill = anno)) +
  #geom_col(position = position_dodge(),width = 0.7, color = 'black')+
  geom_col(position = position_stack(),width = 0.7, color = 'black')+
  #geom_errorbar(aes(ymin = lb, ymax = ub),position = position_dodge(width = 0.7), width = 0.5)+
  ylab(paste0('flux\n',unit.f)) +
  scale_y_continuous(expand = c(0,0))+ scale_x_discrete(limits = x.limit)+
  scale_fill_manual(values = color.path, guide = 'none') +
  theme_bar
fig2s <- ggplot(eff.plot, aes(x = note, y = f*protein, fill = anno)) +
  geom_col(position = position_dodge(),width = 0.7, color = 'black')+
  geom_errorbar(aes(ymin = (f-f.e)*protein, ymax = (f+f.e)*protein),position = position_dodge(width = 0.7), width = 0.5)+
  ylab(paste0('protein\n',unit.p)) +scale_x_discrete(limits = x.limit)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = color.path, guide = 'none') +
  theme_bar
# fig3s <- ggplot(eff.plot, aes(x = note, y = efficiency, fill = anno)) + 
#   geom_col(position = position_dodge(),width = 0.7, color = 'black')+
#   geom_errorbar(aes(ymin =  efficiency.lb, ymax = efficiency.ub),position = position_dodge(width = 0.7), width = 0.5)+
#   ylab(paste0('efficiency\n',unit.e)) + scale_x_discrete(limits = x.limit)+
#   scale_y_continuous(expand = c(0,0))+
#   scale_fill_manual(values = color.path) +
#   theme_bar

fig4y <- ggplot(eff.plot, aes(x = anno, y = efficiency, fill = anno)) + 
  geom_col(position = position_dodge(),width = 0.7, color = 'black')+
  geom_errorbar(aes(ymin =  efficiency.lb, ymax = efficiency.ub),position = position_dodge(width = 0.7), width = 0.5)+
  ylab(paste0('efficiency\n',unit.e)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,2500), breaks = c(0,1000,2000))+
  scale_fill_manual(values = color.path, guide = 'none') +
  facet_wrap(.~rev(note), scales = 'free')+
  theme_bar + theme(axis.ticks = element_blank(),strip.text = element_blank())
fig4t <- ggplot(eff.plot, aes(x = anno, y = efficiency, fill = anno)) + 
  geom_col(position = position_dodge(),width = 0.7, color = 'black')+
  geom_errorbar(aes(ymin =  ifelse(efficiency.lb<0,0,efficiency.lb), ymax = efficiency.ub),position = position_dodge(width = 0.7), width = 0.5)+
  ylab(paste0('efficiency\n',unit.e)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,500), breaks = c(0,200,400))+
  scale_fill_manual(values = color.path, guide = 'none') +
  facet_wrap(.~rev(note), scales = 'free')+
  theme_bar + theme(axis.ticks = element_blank(),strip.text = element_blank()) +coord_cartesian(clip = 'on')

#pdf(paste0('proteome_efficiency/',organism.plot,'.pdf'), w = 12, h = 3.5)

pdf('proteome_efficiency_all.pdf', w = 16, h = 14)
fig1y + fig2y + fig3y +
  fig1t + fig2t + fig3t +
  fig1p + fig2p + fig3p +
  fig1s + fig2s + fig3s + 
  patchwork::plot_layout(ncol = 3, nrow = 4, widths = c(1,1, 1), heights = c(1, 1,1,1))
dev.off()

pdf('ATP_flux_yeast_Tcell.pdf', w = 6, h = 4)
fig1y + fig1t + 
  patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(1,1), heights = c(1))
dev.off()

pdf('proteome_efficiency_cells_diff_scale.pdf', w =12, h = 3)
fig4y+fig4t+
  patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(1,1), heights = c(1))
dev.off()

save.image(file = 'proteome_efficiency_ATP_pathway.RData')


##### plot our data with database or other studies ##### 
pdf('ATP_proteome_fraction.pdf', w = 4,h=4)
ggplot(proteome.efficiency %>% select(-flux.ATP,-lb,-ub,-f.e,-pathway) %>%
         pivot_wider(values_from = 'f', names_from = 'anno') %>%
         mutate(mito = mito - resp) %>%
         pivot_longer(cols = c(mito, resp,glyc), values_to = 'f', names_to = 'anno') %>%
         mutate(anno = factor(anno, levels = c('mito','resp','glyc'))),
       aes(x=paste0(organism,' ',note), y = f, fill = anno)) +
  geom_col(position = position_stack(), color = 'black') +
  annotate('text', label = 'mitochondrial fraction\nexclude respiratory enzymes', x = 14, y = 0.44, size = 2.5)+
  scale_fill_manual(values = c('glyc' = '#FFED6F', 'mito' = '#BC80BD','resp' = '#377EB8')) +theme_classic() + 
  ylab('proteome fraction') + coord_flip()
dev.off()

fit <- summary(lm(log10(flux.ATP)~ log10(f) + yeast,
                  proteome.efficiency %>% filter(flux.ATP > 0) %>% mutate(yeast = grepl('yeast',organism))))

pdf('ATP_fraction_flux.pdf', w = 3.6,h=2.2)
ggplot(proteome.efficiency %>% filter(!anno == 'mito') %>% filter((!grepl('yeast',organism))|grepl('B',note)),
       aes(x=f, y = flux.ATP, fill = anno, shape = organism)) +
  geom_abline(slope = 1, intercept = log10(486), color = 'grey')+
  # geom_errorbar(aes(ymin = lb, ymax = ub)) + geom_errorbarh(aes(xmin = f-f.e, xmax = f+f.e))+
  geom_point(aes(size = (organism=='mouse'))) + 
  theme_classic() +
  scale_y_log10(limits= c(5E-3,NA)) +scale_x_log10() +
  # scale_y_log10(limits = c(1E-3,1E2)) +scale_x_log10(limits = c(1E-5,1)) +
  scale_fill_manual(values = c('glyc' = '#FFED6F', 'mito' = '#BC80BD','resp' = '#377EB8')) +
  scale_shape_manual(values = c('yeast_SC' = 24, 'yeast_IO' = 24, 'T cell' = 20,
                                'mouse_pancreas' = 25,'mouse_spleen' = 25, 'mouse_paxdb' =21, 
                                'cancer cell lines' = 23,'E.coli' = 22))+
  # scale_size_manual(values = c('TRUE' = 2, 'FALSE' = 3), guide = guide_legend('none'))+
  xlab('proteome fraction') + ylab('ATP flux [mmole/h/gDW]') 
dev.off()

# order.organism <- c('yeast_IO','yeast_SC','E.coli','cancer cell lines','mouse')
# color.organism <- c('#377EB8','mediumseagreen','#FB9A99','grey','#A6CEE3')
# names(color.organism) <- order.organism
order.type <- c('naive T cell','activated T cell','yeast_IO','yeast_SC C limit','yeast_SC',
                    #'spleen','leukemic spleen','pancreas','GEMM PDAC','flank PDAC',
                    'E.coli','cancer cell lines','mouse','mouse tissues','mouse tumors')
color.type<- c("#B2DF8A","#33A02C","#6BAED6","#FDBF6F","#FF7F00",
               #"#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#FF7F00",
               "#6A3D9A",'grey',"#CAB2D6","#CAB2D6","#6A3D9A")
names(color.type) <- order.type
anno.org <- proteome.efficiency %>% select(organism, note) %>% distinct() %>%
  mutate(type = case_when(organism == 'yeast_SC' ~ paste0(organism,ifelse(grepl('\\_C',note),' C limit','')),
                          organism == 'T cell' ~ paste(note,organism,sep = ' '),
                          note == 'leukemic spleen'|grepl('PDAC',note) ~ 'mouse tumors',
                          grepl('mouse',organism) ~ 'mouse tissues',
                          TRUE ~ as.character(organism))) %>%
  mutate(org = case_when(grepl('yeast', organism) ~ organism,
                         grepl('mouse',organism) ~ 'mouse tissues',
                         TRUE ~ as.character(organism))) %>%
  mutate(type=factor(type, levels = order.type)) %>% arrange(desc(type))
dt.plot <- proteome.efficiency  %>% 
  select(organism,starts_with('efficiency'),note,anno) %>%
  pivot_longer(cols = starts_with('efficiency'), names_to = 'stat', values_to = 'value') %>%
  mutate(stat = paste0(anno, gsub('efficiency','',stat))) %>% select(-anno) %>%
  pivot_wider(names_from = 'stat', values_from = 'value') %>%
  # mutate(type = case_when(grepl('yeast', organism) ~ organism,
  #                         organism == 'T cell' ~ paste(note,organism,sep = ' '),
  #                         organism == 'mouse_paxdb' ~ 'mouse',
  #                         organism %in% c('mouse_spleen','mouse_pancreas') ~ note,
  #                         TRUE ~ as.character(organism))) %>%
  # mutate(type=factor(type, levels = order.type)) %>% arrange(desc(type))
  left_join(anno.org) %>% arrange(desc(type))

pdf('ATP_proteome_efficiency_respiration~glycolysis.pdf', w = 9,h=6)
ggplot(dt.plot %>% filter(!organism == 'E.coli'),
       aes(x=glyc, y = resp, color = type)) +
  geom_abline(slope = 1, intercept =log10(2), color = 'grey', linetype = 'dashed')+
  geom_abline(slope = 1, intercept =log10(0.5), color = 'grey', linetype = 'dashed')+
  geom_abline(slope = 1, intercept =0, color = 'grey')+
  # geom_errorbar(aes(ymin = ifelse(resp.lb<=0,resp*0.1,resp.lb), ymax = resp.ub), alpha = 0.7) +
  # geom_errorbarh(aes(xmin = ifelse(glyc.lb<=0,glyc*0.1,glyc.lb), xmax = glyc.ub), alpha = 0.7)+
  geom_point(size = 3) +
  scale_color_manual(values= color.type)+
  theme_classic() +
  scale_y_log10(limits = c(2E-1,3000),breaks = c(1,10,100,1000), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(limits = c(2E-1,3000),breaks = c(1,10,100,1000), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab('glycolytic\nproteome efficiency') + ylab('respiratory\nproteome efficiency') +
  theme_classic() + theme(line = element_line(size = 0.5),
                          axis.text = element_text(size = 20),axis.title = element_text(size=24),
                          legend.text = element_text(size = 16), legend.title = element_blank(),
                          plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))
dev.off()

plot.1 <- rbind(proteome.efficiency %>% select(note, anno, f,organism) %>% filter(!grepl('yeast',organism)),
                 dt.yeast %>%  group_by(organism,note,anno) %>% summarise(f = mean(f)) %>% ungroup()
                 #dt.NCI60 %>% group_by(cell,anno) %>% summarise(f = mean(f)) %>% ungroup() %>% rename(note = cell) %>% mutate(organism = 'cancer cell lines')
                 ) %>% 
  pivot_wider(names_from = 'anno', values_from = 'f') %>% mutate(r = glyc/resp) %>%
  left_join(anno.org)

pdf('proteome_allocation_across_organisms.pdf', w = 8,h=6)
ggplot(plot.1 %>% filter(!type == 'E.coli'), aes(x = type, y = r,color = type)) +
  stat_boxplot(aes(group = type), width = 0.6,position = position_identity(), outlier.shape = NA) +
  geom_jitter(width = 0.15, pch=21, size = 2) +
  theme_classic() + 
  scale_y_log10(limits = c(2^-9, 2^3.5),breaks = 2^c(-8,-6,-4,-2,0,2), labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  # scale_y_log10(breaks = c(0.5,1,2,4,8,16),limits = c(0.2,NA)) +
  scale_color_manual(values = color.type) +
  #scale_x_discrete(limits = c('yeast','T cell','cancer cell lines','mouse tissues'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('glycolytic / respiratory\nproteome allocation')+
  geom_hline(yintercept =1, linetype= 'dashed') + theme(text = element_text(size = 24))
dev.off()

plot.2 <- rbind(proteome.efficiency %>% select(note, anno, flux.ATP,organism) %>% filter(!grepl('yeast',organism)),
                 dt.yeast %>% select(organism,note,anno, flux.ATP) %>% distinct()
                 # dt.NCI60 %>% select(cell,anno,flux) %>% distinct() %>% 
                 #   rename(note = cell, flux.ATP = flux) %>% mutate(organism = 'cancer cell lines')
                 ) %>% 
  # filter((!organism == 'mouse') | note == 'liver') %>%
  pivot_wider(names_from = 'anno', values_from = 'flux.ATP') %>% mutate(r = glyc/resp) %>%
  left_join(anno.org)

pdf('pathway_ATP_contribution_across_organisms.pdf', w = 8,h=6)
ggplot(plot.2 %>% filter(!type == 'E.coli'), aes(x = type, y = r,color = type)) +
  geom_hline(yintercept =1, linetype= 'dashed')+
  stat_boxplot(aes(group = type), width = 0.6,position = position_identity(), outlier.shape = NA) +
  geom_point(width = 0.15, pch=21, size=2) +
  theme_classic() + 
  scale_y_log10(limits = c(2^-9, 2^3.5),breaks = 2^c(-8,-6,-4,-2,0,2), labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  # scale_y_log10(breaks = c(0.5,1,2,4,8,16),limits = c(0.2,NA)) +
  scale_color_manual(values = color.type) +
  # scale_x_discrete(limits = c('yeast','T cell','cancer cell lines','mouse tissues'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('glycolytic / respiratory\nATP flux')+ theme(text = element_text(size = 24))
dev.off()

plot.3 <- plot.1 %>% rename(r.f = r)%>% select(note,organism,r.f,type,org) %>% 
  full_join(plot.2 %>% rename(r.flx = r) %>% select(note,organism,r.flx,type,org))

pdf('energy preference~proteome allocation.pdf', w = 9,h=6)
ggplot(plot.3 %>% filter(!organism == 'E.coli') %>% arrange(org),
       aes(x=r.f, y = r.flx, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values= color.type)+
  theme_classic() +
  scale_y_log10(limits = c(1E-4,1E1),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(limits = c(1E-4,1E1),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  xlab('proteome allocation\n[glyc/resp]') + ylab('ATP contribution\n[glyc/resp]') +
  theme_classic() + theme(line = element_line(size = 0.5),
                          axis.text = element_text(size = 20),axis.title = element_text(size=24),
                          legend.text = element_text(size = 16), legend.title = element_blank(),
                          plot.margin =unit(c(0.05,0.05,0.05,0.05), 'npc'))
dev.off()

summary(lm(log(r.flx)~log(r.f), plot.3))
