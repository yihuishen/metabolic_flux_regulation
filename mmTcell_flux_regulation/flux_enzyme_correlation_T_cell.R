library(dplyr)
library(openxlsx)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = F)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file = 'flux_enzyme_correlation_T_cell.RData')

load('flux_enzyme_correlation_T_cell/subsystem.anno.RData')

protein.mass <- data.frame(state = c('naive','activated'),
                           protein = c(25.8, 61.9), # ug protein / E6 cells
                           cell.weight = c(29.33,70.4)) # ug(protein + DNA + RNA) / E6 cells

## load 13C mfa flux data, in nmol/h/gCell

load(file='../Tcell_fluxomics/flux_T_cell.RData')

mod.mouse <- flux.long %>% select(id,name,reaction,subsystem,gpr) %>% distinct() %>%
  as_tibble() %>% mutate(geneID=str_extract_all(gpr, '[A-Z][a-z0-9]+')) %>%
  unnest(cols = geneID) %>% 
  left_join(subsystem.new) %>% mutate(subsystem = ifelse(is.na(new), subsystem, new)) %>%
  select(-new) %>% left_join(anno.subsystem)

mod <- mod.mouse %>% select(-geneID) %>% distinct() %>%
  mutate(subs = unlist(lapply(strsplit(reaction, split = '>'),'[',1))) 

mod.subs <- mod %>% as_tibble() %>% mutate(subs = str_extract_all(subs, '[a-zA-Z0-9_]+\\_[cem]')) %>%
  unnest(cols = c(subs)) %>% mutate(subs = gsub('\\_[cme]', '', subs))

flux <- flux.long %>% # flux in nmol/h/gCell
  select(-subsystem) %>%
  full_join(mod.mouse %>% select(id, name, reaction,subsystem, category) %>% distinct()) %>%
  mutate(category = ifelse(is.na(category),'other',category))

ggplot(flux, aes(x = state, y = log10(abs(flux)))) + geom_point()
# the smallest non-zero flux is abut 5E-4 in activated and 2E-4 in naive
min(flux$flux[flux$flux > 2E-6], na.rm=T)
min.flux <- 5E-5

flux.min.adjust <- flux %>% right_join(expand.grid(id = unique(flux$id), state = c('naive','activated')))  %>%
  mutate(flux = ifelse((lb < min.flux)|(abs(ub)<min.flux), (lb+ub)/2, flux))%>%
  mutate(flux = ifelse((!is.finite(flux)), min.flux, flux))
# overview of flux change
tmp <- flux.min.adjust %>% select(id,name,flux,state,category,subsystem) %>% 
  pivot_wider(names_from = 'state', values_from = 'flux') %>%
  mutate(FC =  activated/naive) %>% 
  mutate(category.plot = case_when(category %in% c('amino acid','folate','lipid','nucleic acid','sugar')~'biosynthetic',
                                   category %in% c('PPP','glycolysis','TCA','ox phos') ~ category,
                                   TRUE ~ 'other'))
order.tmp <- c('glycolysis','TCA','ox phos','PPP','biosynthetic','other')
color.tmp <- c("#E31A1C","#33A02C", "#1F78B4","#FF7F00","#BEBADA", 'grey60')
names(color.tmp) <- order.tmp
pdf('flux_enzyme_correlation_T_cell/overview_flux_FC.pdf', w= 8, h = 4)
ggplot(tmp, aes(x = category.plot, y= abs(FC), color = category.plot)) + 
  geom_jitter(width = 0.2, alpha = 0.8)  + coord_flip() +
  # scale_y_continuous(limits = c(0,16), expand = c(0,0))+
  scale_y_log10(limits = c(0.5,2^11),
                labels = scales::trans_format("log2", scales::math_format(2^.x)),
                expand = c(0,0),breaks = 2^c(0,4,8)) + 
  scale_x_discrete(limits = order.tmp) + scale_color_manual(values = color.tmp) +
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()
## load proteomics data
load('../Tcell_proteomics/proteomics_grouped_Tcell.RData')
mapping.gene <- dt.Tcell.proteome %>% select(Entry, Gene.names) %>% distinct() %>%
  as_tibble() %>% mutate(geneID = strsplit(Gene.names, split = ' ')) %>% unnest(cols = geneID)

# calculate protein fraction per reaction, averaged between replicates
protein.gene <- dt.Tcell.proteome %>% group_by(Entry, Gene.names,group,tissue) %>%
  summarise(f = mean(f)) %>% ungroup() %>% rename(state = tissue) %>%
  inner_join(mapping.gene) %>% 
  inner_join(mod.mouse)
protein.rxn <- protein.gene %>% group_by(id, name, state) %>% summarise(f = sum(f)) %>% ungroup()

flux.protein <-  flux.min.adjust %>%
  # only for plotting fold change, when the confidence interval spans 0, use the (lb+ub)/2 to represent flux
  left_join(protein.rxn) %>%
  left_join(protein.mass) %>% mutate(enzyme = f* protein/cell.weight) 

flux.protein.Tcell <- flux.protein %>% select(-f,-protein,-cell.weight) %>% 
  group_by(id) %>% mutate(flxFC = flux/2^mean(log2(flux)), enzFC = enzyme/2^mean(log2(enzyme))) %>% ungroup() %>%
  arrange(flux)

tmp <- flux.protein.Tcell %>% filter(is.finite(log(enzFC)),is.finite(log(flxFC))) # %>% 
  filter(category %in% c('glycolysis','PPP','TCA','ox phos'))
cov(log(tmp$enzFC),log(tmp$flxFC))/sqrt(var(log(tmp$flxFC))*var(log(tmp$enzFC)))
cov(log(tmp$enzFC),log(tmp$flxFC))/var(log(tmp$flxFC))

plot.lb = 0.03
plot.ub = 50

fig <- ggplot(flux.protein.Tcell %>% filter(is.finite(enzFC),is.finite(flxFC)), aes( x = enzFC,y = flxFC, color = category, shape = state)) +
  geom_abline(slope = 1, color = 'black', linetype= 'dashed')+
  geom_jitter( height = 0.02, size = 3) + 
  # geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  scale_color_manual(values = color.subsystem) +
  scale_shape_manual(values = c('naive' = 21, 'activated' = 24)) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  # geom_text_repel(aes(label = ifelse(enzFC > flxFC, geneID, '')))+
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
pdf('flux_enzyme_correlation_T_cell/flx~enz.pdf', w=4, h = 4)
print(fig)
dev.off()

dens.x <- ggplot(flux.protein.Tcell,
                 aes(x = enzFC, 
                     fill = ifelse(category %in% names(color.subsystem)[1:4],as.character(category),'other'))) + 
  geom_density(alpha = 0.4, color = NA) + 
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values = color.subsystem) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme_void() + 
  theme(legend.position = "none") 
  

dens.y <- ggplot(flux.protein.Tcell, 
                 aes(x = flxFC, 
                     fill = ifelse(category %in% names(color.subsystem)[1:4],as.character(category),'other'))) + 
  geom_density(alpha = 0.4, color = NA) + 
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values = color.subsystem) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme_void() + 
  theme(legend.position = "none") +  coord_flip()

pdf('flux_enzyme_correlation_T_cell/flx~enz_w_density.pdf', w=5.5, h = 5.5)
dens.x + patchwork::plot_spacer() + fig + dens.y + 
  patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4, 1.5), heights = c(1.5, 4))
dev.off()

met.Tcell <- read.xlsx('../Tcell_metabolomics/13C_poolsize_anno_report.xlsx',  
                       startRow = 2, cols = c(3,12,13,27:38)) %>% # ion count normalized to cell number
  rename(subs = BiGG) %>%
  pivot_longer(cols = starts_with('A')|starts_with('N'), 
               values_to = 'ion.count', names_to = 'sample') %>%
  mutate(state = ifelse(grepl('A',sample),'activated','naive')) %>%
  group_by(subs, BiGG.name, Compound, state) %>%
  summarise(ion.count = mean(ion.count)) %>% ungroup() %>%
  mutate(ion.count = ion.count/ifelse(state == 'naive',29.3,70.4)) %>%  # normalized to cell weight
  group_by(Compound) %>% mutate(metFC = ion.count/2^mean(log2(ion.count))) %>% ungroup()

fpm.Tcell <-mod.subs %>% inner_join(met.Tcell) %>%
  group_by(id, reaction,state) %>% summarise(metFC = 2^mean(log2(metFC), na.rm = T)) %>% ungroup() %>%
  full_join(flux.protein.Tcell) 

tmp <- fpm.Tcell %>% filter(is.finite(log(metFC)),is.finite(log(flxFC))) #%>% 
  filter(category %in% c('glycolysis','PPP','TCA','ox phos'))
cov(log(tmp$metFC),log(tmp$flxFC))/sqrt(var(log(tmp$flxFC))*var(log(tmp$metFC)))
cov(log(tmp$metFC),log(tmp$flxFC))/var(log(tmp$flxFC))


fig <- ggplot(fpm.Tcell %>% filter(is.finite(flxFC)), aes( x = metFC,y = flxFC, color = category, shape = state)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed')+
  geom_jitter(height = 0.02, size = 3) + 
  scale_color_manual(values = color.subsystem) +
  scale_shape_manual(values = c('naive' = 21, 'activated' = 24)) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(limits = c(plot.lb,plot.ub),labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  theme_classic() + theme(legend.position = "none", text = element_text(size = 24))
pdf('flux_enzyme_correlation_T_cell/flx~met.pdf', w=4, h = 4)
print(fig)
dev.off()

dens.x <- ggplot(fpm.Tcell,
                 aes(x = metFC, 
                     fill = ifelse(category %in% names(color.subsystem)[1:4],as.character(category),'other'))) + 
  geom_density(alpha = 0.4, color = NA) + 
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values = color.subsystem) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme_void() + 
  theme(legend.position = "none") 


dens.y <- ggplot(fpm.Tcell, 
                 aes(x = flxFC, 
                     fill = ifelse(category %in% names(color.subsystem)[1:4],as.character(category),'other'))) + 
  geom_density(alpha = 0.4, color = NA) + 
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values = color.subsystem) +
  scale_x_log10(limits = c(plot.lb,plot.ub), labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  theme_void() + 
  theme(legend.position = "none") +  coord_flip()


pdf('flux_enzyme_correlation_T_cell/met~flx_w_density.pdf', w = 5.5, h = 5.5)
dens.x + patchwork::plot_spacer() + fig + dens.y + 
  patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4, 1.5), heights = c(1.5, 4))
dev.off()

pdf('flux_enzyme_correlation_T_cell/met_enzyme_regulation.pdf', w = 3.2, h = 3)
ggplot(fpm.Tcell %>% filter(state == 'activated'),
       aes( x = log2(enzFC)/log2(abs(flxFC)),y = log2(metFC)/log2(abs(flxFC)), color = category)) +
  geom_abline(slope = -1, intercept = 1, linetype = 'dashed') +
  geom_jitter(pch = 21, size = 2) + 
  scale_color_manual(values = color.subsystem) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #scale_x_log10(limits = c(1,500)) + scale_y_log10(limits = c(1,500)) + 
  scale_x_continuous(limits = c(-0.6,3), expand = c(0,0), breaks = c(0,1,2)) +
  scale_y_continuous(limits = c(-0.6,3), expand = c(0,0), breaks = c(0,1,2)) + 
  theme_classic() + theme(legend.position = "none", text = element_text(size = 16),
                          axis.line = element_blank())
dev.off()

tmp <- fpm.Tcell %>% mutate(a = metFC/enzFC) %>% filter(state == 'activated', is.finite(a)) %>%
  select(id, enzFC, metFC, flxFC, category) %>%
  pivot_longer(cols = c(enzFC, metFC), values_to='FC',names_to = 'component') %>%
  mutate(rho = log(FC)/log(abs(flxFC))) %>% arrange(category)
tmp.ttest <- tmp %>% select(id,rho,category,component) %>%
  pivot_wider(names_from = 'component',values_from = 'rho') %>%
  filter(category %in% names(color.subsystem)[1:5])
t.test(tmp.ttest$enzFC,tmp.ttest$metFC, 'less')

# paired one-way t.test rho.m vs rho.e
# p = 4E-6 for central metabolism, mean difference = 1.08
# p = 0.03 for central metabolism, mean difference = 0.16

pdf('flux_enzyme_correlation_T_cell/regulation_coefficient_central.pdf', w = 2.6, h = 3)
ggplot(tmp%>% filter(category %in% names(color.subsystem)[1:5]),
       aes( x = component,y = rho, color = category, group = id)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point() + geom_line(alpha = 0.7) + ylim(-1,4) +
  scale_color_manual(values = color.subsystem) +
  egg::theme_presentation() +
  theme(legend.position = "none", panel.background = element_blank(), axis.title.x = element_blank())
dev.off()

corr.coeff <- data.frame(scope = rep(c('between yeast','S. cere across conditions','I. orie across conditions','T cell activation'),2),
                         component = c(rep('S',4),rep('E',4)),
                         cov = c(-0.03, 0.2, 0.44,0.42,0.44, 0.11, -0.04, 0.27))

pdf('flux_enzyme_correlation_T_cell/covariance_yeast_Tcell.pdf', w = 5.5, h = 6)
ggplot(corr.coeff, aes(x = scope, y = cov, fill = component)) +
  geom_hline(yintercept = 0) +
  #geom_col(position = position_dodge2(width = 0.5), width = 0.7, alpha = 0.7, color = 'black') +
  geom_col(position = position_stack(), width = 0.6, alpha = 0.7, color = 'black') +
  scale_x_discrete(limits = unique(corr.coeff$scope)) +
  egg::theme_presentation() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c('E' = 'grey', 'S' = '#3DAAD6'))
dev.off()

## plot balance
source('../manuscript/checkbalance.R')
met.check <- 'pyr_m'

flux.balance <- flux %>% select(id, reaction, state, flux, lb, ub) %>% distinct() %>%
  rename(flx=flux) %>%
  group_by(state) %>%
  do(a = check.balance(met = met.check,mfa = ., ci = T)) %>% tidyr::unnest(a) %>% ungroup() %>%
  mutate(type = ifelse(flx < 0, 'consumption','production')) %>%
  mutate(flx = abs(flx), lb.tmp = ifelse(abs(lb)<abs(ub), abs(lb),abs(ub)),ub.tmp = ifelse(abs(lb)<abs(ub), abs(ub),abs(lb))) %>%
  select(-lb,-ub) %>% rename(lb = lb.tmp, ub = ub.tmp) %>%
  group_by(state) %>%
  arrange(flx, desc(flx)) %>%
  mutate(tmp = gsub('\\_[cm]','',as.character(id))) %>%
  mutate(id.label = ifelse(flx > 0.05*max(flx), tmp, 'other')) %>% # do not show minor fluxes
  mutate(id.label = ifelse(tmp %in% id.label, tmp, 'other')) %>% select(-tmp) %>%
  ungroup() %>%
  mutate(id.label = factor(id.label, levels = unique(id.label))) %>%
  group_by(state, type) %>%
  arrange(desc(id.label)) %>%
  mutate(anno = flx/2 + c(0, cumsum(flx)[-length(flx)])) %>%
  ungroup()

pdf(paste0('flux_balance_T_cell/flux_balance_',gsub('\\\\','',met.check),'.pdf'), w = 6, h = 4)
ggplot(flux.balance %>% mutate(state = factor(state,levels = c('naive','activated'))),
       aes(x = type, y = flx, fill = id.label)) +
  geom_col(alpha = 0.8,  width = 0.8, color = 'white', size = 0.25) +
  geom_text(aes(y = anno,
                label = ifelse(id.label == 'other','',as.character(id.label))), size=4) +
  facet_wrap(.~state, scales = 'free')+ 
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0,NA)) +
  scale_fill_brewer(palette = 'Set3', guide = 'none')+
  ylab('flux [mmol/h/gDW]') +egg::theme_presentation() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
dev.off()

## 4C balance of TCA
met.check.list <- c('cit','acon_C','icit','asp__L','oaa','mal__L','fum','succ','succoa','argsuc',
                    'akg','glu__L','glu5p')

flux.balance.list <- list()

for (met.check in c(paste0(met.check.list,'_m'),paste0(met.check.list,'_c'))) {
  flux.balance.list[[met.check]] <- flux %>% rename(flx=flux) %>% select(id, reaction, state, flx, lb, ub) %>% distinct() %>%
    group_by(state) %>%
    do(a = check.balance(met = met.check,mfa = ., ci = T)) %>% tidyr::unnest(a) %>% ungroup()
}
flux.balance <- data.table::rbindlist(flux.balance.list) %>%
  group_by(state,id,reaction) %>% summarise(flx = sum(flx)) %>% filter(!abs(flx) <1E-7)%>%
  left_join(flux.long %>% select(id,subsystem) %>% distinct()) %>%
  mutate(tmp = gsub('\\_[cm]','',as.character(id))) %>%
  mutate(id.label = case_when(grepl('BIOMASS',tmp) ~ 'BIOMASS',
                              tmp == 'G5SDy' ~ 'BIOMASS',
                              id == 'PC_m' ~ id,
                              subsystem %in% c('Pyrimidine metabolism','Purine metabolism') ~ subsystem,
                              TRUE ~ tmp))%>% select(-tmp)

flux.balance.plot <- flux.balance %>% # summarise biomass
  mutate(id.label = ifelse(id.label %in% c('Pyrimidine metabolism','Purine metabolism'),'BIOMASS',id.label)) %>%
  group_by(state,id.label) %>% summarise(flx=sum(flx)) %>% ungroup() %>%
  mutate(type = ifelse(flx < 0, 'consumption','production')) %>%
  mutate(flx = abs(flx)) %>%  
  arrange(flx, desc(flx)) %>%
  mutate(id.label = factor(id.label, levels = unique(id.label))) %>%
  group_by(state, type) %>%
  arrange(desc(id.label)) %>%
  mutate(anno = flx/2 + c(0, cumsum(flx)[-length(flx)])) %>%
  ungroup()  %>% mutate(state = factor(state, levels = c('naive','activated')))

# glu5p --(G5SDy_m)--> glu5sa = arginine and proline sink to protein

pdf('flux_balance_T_cell/flux_balance_TCA_4C.pdf', w = 6, h = 4)
ggplot(flux.balance.plot, aes(x = type, y = flx, fill = id.label)) +
  geom_col(alpha = 0.8,  width = 0.8,color = 'white', size = 0.25) +
  geom_text(aes(y = anno,
                label = ifelse(id.label == 'other','',as.character(id.label))), size=4) +
  facet_wrap(.~state, scales = 'free')+ 
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0,NA)) +
  scale_fill_brewer(palette = 'Set3', guide = 'none')+
  ylab('flux [nmol/h/gDW]') +egg::theme_presentation() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
dev.off()

save.image(file = 'flux_enzyme_correlation_T_cell.RData')

#### calculate proteome efficiency in T cells ####
met.check <- 'atp_c'

flux.balance <- flux %>% select(id, reaction, state, flux, lb, ub) %>% distinct() %>%
  rename(flx=flux) %>%
  group_by(state) %>%
  do(a = check.balance(met = met.check,mfa = ., ci = T)) %>% tidyr::unnest(a) %>% ungroup() %>%
  left_join(flux %>% select(id,category) %>% distinct()) %>%
  filter(category == c('glycolysis') | id == 'ADPATPt_c_m')

ATP.flux <- flux.balance %>% group_by(category, state) %>%
  summarise(flux.ATP = sum(flx),
            lb = sum(lb), ub=sum(ub)) %>% rename(tissue=state) %>% ungroup() %>%
  mutate(pathway = ifelse(category == 'glycolysis','glycolysis','respiration')) %>% 
  select(-category)

dt.Tcell.ATP <- frac.proteome.Tcell.ATP.mean %>%
  left_join(protein.mass %>% rename(tissue = state)) %>%
  left_join(ATP.flux) %>% mutate(tissue = factor(tissue, levels = c('naive','activated'))) %>%
  mutate(`%protein` = protein/cell.weight) %>% select(-cell.weight, -protein)

save(dt.Tcell.ATP, file = 'flux_enzyme_correlation_T_cell/Tcell_ATP.RData')


#### compare to yeast in vivo turnover number
load(file='flux_enzyme_correlation/turnover_number.RData')
# enzyme activity in umole/min/mg
turnover.yeast <- turnover
turnover <- flux.protein %>% 
  mutate(id = gsub('q9','qx',id)) %>% mutate(activity = abs(flux)/(enzyme)/60) %>%
  select(id, activity,state) %>% filter(!activity<1E-8) %>%
  pivot_wider(names_from = 'state', values_from = 'activity') %>% full_join(turnover.yeast)

pdf('flux_enzyme_correlation_T_cell/activity~actTcell_vs_naiTcell.pdf', w=5.5, h = 5)
ggplot(turnover %>% filter(pathway == 'other'),
       aes(x = naive, y = activated)) +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  geom_point(size = 2, pch=21, color='grey40') +
  geom_point(data=turnover %>% filter(pathway == 'resp'), size = 2, color='#377EB8') +
  geom_point(data=turnover %>% filter(pathway == 'glyc'), size = 2, pch=21, fill='#FFED6F', color = 'goldenrod') +
  # scale_fill_manual(values = color.path, guide = 'none') +
  xlab('activity in naive\n[umole/min/mg]') + ylab('activity in activated\n[umole/min/mg]') +
  scale_x_log10(limits= c(3E-4,3E4),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits= c(3E-4,3E4),labels = scales::trans_format("log10", scales::math_format(10^.x))) + # only consider significant fluxes
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()
summary(lm(log10(activated)~log10(IO),turnover))
summary(lm(log10(activated)~log10(naive),turnover))
summary(lm(log10(IO)~log10(SC),turnover))
save(turnover,file='flux_enzyme_correlation_T_cell/turnover_number.RData')

tmp <- turnover %>% mutate(log2FC = log2(activated/naive))
tmp.resp <- tmp %>% filter(pathway == 'resp')
tmp.glyc <- tmp %>% filter(pathway == 'glyc')
pdf('flux_enzyme_correlation_T_cell/activity_log2FC_actT_naiT_pathway.pdf', w=6, h = 4)
ggplot(tmp %>% mutate(pathway = 'all'), aes(x = pathway, y = log2FC, fill = pathway)) +
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_violin(size = 0.5)+ geom_violin(data=tmp %>% filter(pathway %in% c('glyc','resp')), size = 0.5)+
  geom_jitter(width = 0.2,size = 1, color = 'black')+
  geom_jitter(data = tmp %>% filter(pathway %in% c('glyc','resp')),width = 0.2,size = 1, color = 'black')+
  # geom_jitter(data=tmp.resp, width = 0.3,size = 2, color='#377EB8') +
  # geom_jitter(data=tmp.glyc, width = 0.3,size = 2, pch=21, fill='#FFED6F', color = 'goldenrod') +
  scale_fill_manual(values = c('glyc' = '#FFED6F','resp' = '#377EB8','all' = 'grey'))+
  scale_x_discrete(limits =c('all','glyc','resp'))+
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()
ks.test(tmp$log2FC, tmp.glyc$log2FC)
