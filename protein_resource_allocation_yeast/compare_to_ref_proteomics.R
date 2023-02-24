rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# need group.SC
load('protein_resource_allocation/proteomics_grouped.RData')
load('protein_resource_allocation/assign_group_files.RData')

#### reference protein abundance data ####
# Ho et al. 2018, Cell Systems (https://doi.org/10.1016/j.cels.2017.12.004)
# with FDB: Bartolomeo et al. 2020, PNAS (www.pnas.org/cgi/doi/10.1073/pnas.1918216117)
ref.anno <- read.xlsx('protein_abundance_SGD/table1.xlsx',startRow = 2)
ref.filter <- rbind(ref.anno %>% filter(#grepl('minimal',Medium),
  grepl('log',Growth.Phase)) %>%
    rename(source = Abbreviation, method = Type.of.Study) %>% select(source, method),
  data.frame(source = c('FDB','OURS'), method ='mass spectrometry'))

ref.SGD <- read.xlsx('protein_abundance_SGD/tableS4.xlsx', startRow = 2)[,c(-2:-4,-6)] %>%
  rename(geneID = Systematic.Name, median.copy = Median.molecules.per.cell) %>% 
  filter(!is.na(median.copy)) %>%
  left_join(map.SC %>% select(geneID, Length))

tAA.SGD <- sum(ref.SGD$median.copy*ref.SGD$Length, na.rm=T)
# total number of AA (SGD median)
# all data will be normalized to this value
# some numbers for sanity check
# 4.2E-14L per cell, 2E-11 gDW per cell, 50% DW is protein, MW in kD, so about 1E-11g protein is expected per cell
tAA.SGD/6.02E23*113 # protein weight per cell
tAA.SGD/6.02E23/4.2E-14*1E3 # total conc.AA in mmole/L

## our measurement on SC
mySC <- group.SC %>% filter(condition == 'abs_B', grepl('120421',source)|source == '081921') %>%
  group_by(geneID, Protein.names, group, Length) %>%
  summarise(OURS.f = mean(f)) %>% ungroup() %>%
  mutate(OURS = OURS.f * tAA.SGD/Length) %>% # normalize to SGD reference total protein mass
  mutate(geneID = str_extract(geneID,'[RYQ][^;]+')) # some protein is mapped to multiple genes, select the first one
mySC.rank <- mySC %>% rename(copy.mySC = OURS) %>% arrange(desc(copy.mySC)) %>% mutate(rank.mySC = 1:length(copy.mySC))

## reference proteomics from Jens Nilsen 2020 PNAS
tmp <- read.xlsx('protein_abundance_Bartolomeo/JNielsen2020PNAS_S2_whole_cell.xlsx', startRow = 2)
ref.FDB <- data.frame(geneID = tmp$GeneNameOrdered, f = rowMeans(tmp[,5:7])) %>% # f is represented as g/gCellDW
  mutate(geneID = str_extract(geneID,'[RYQ][^;]+')) %>%
  left_join(map.SC %>% select(geneID, Length)) %>%
  mutate(f = f/sum(f)) %>% # normalize to whole proteome
  mutate(FDB = f*tAA.SGD/Length) %>% select(-f)

ref.all <- ref.SGD %>%
  full_join(ref.FDB) %>% 
  full_join(map.SC %>% select(geneID, Protein.names, group, Length)) %>%
  full_join(mySC %>% select(geneID, OURS)) %>%
  arrange(desc(median.copy)) %>%
  mutate(rank = 1:length(median.copy))

ref.long <- ref.all %>%
  pivot_longer(cols = where(is.numeric) & (!rank) & (!median.copy) & (!Length),
               names_to = 'source',values_to = 'copy') %>%
  filter(!is.na(copy)) %>%
  right_join(ref.filter) %>%
  left_join(mySC.rank %>% select(geneID, copy.mySC, rank.mySC))%>% select(-Length) %>%
  left_join(map.SC %>% select(geneID, Length)) %>% mutate(method = gsub('TAP-','',method)) %>%
  arrange(desc(method)) %>% mutate(source=factor(source, levels = unique(source)))

save.image('compare_to_ref_proteomics.RData')
load('compare_to_ref_proteomics.RData')

# plot whole proteome

pdf('compare_ref_proteome/all_proteome_compare.pdf', w= 7, h = 8)
ggplot(ref.long , aes(x = rank.mySC, y = copy, color = method)) +
  geom_point(size = 0.2) +
  geom_line(aes(x = rank.mySC, y = copy.mySC, group = source), color = 'black', size = 0.5) +
  facet_wrap(method~source,nrow=4)+
  scale_y_log10(breaks = c(1,1E2,1E4,1E6)) +theme_classic() +theme(legend.position= 'none')+
  ggtitle('black line indicate our measurement')
dev.off()

# omit those with narrow dynamic range: mean(min 10 proteins) > 1E-3 * mean(max 10 proteins)
dynamic.range <- ref.long %>% group_by(source) %>%
  arrange(desc(copy)) %>% mutate(rank.tmp = 1:length(copy)) %>%
  filter((rank.tmp < 11)|(rank.tmp > max(rank.tmp)-10)) %>%
  mutate(minmax = ifelse(rank.tmp < 11,'max', 'min')) %>% ungroup() %>%
  group_by(minmax, source, method) %>% summarise(copy.mean = mean(copy)) %>% ungroup() %>%
  pivot_wider(names_from = 'minmax', values_from = 'copy.mean') %>%
  mutate(DR = max/min) %>%
  full_join(ref.long %>% group_by(source) %>% summarise(n = sum(!is.na(copy))) %>% ungroup())

pdf('compare_ref_proteome/all_proteome_compare_truth_plot.pdf', w= 24, h = 7)
ggplot(ref.long %>% filter(!source == 'OURS') ,
       aes(x = copy.mySC, y = copy, color = method)) + geom_abline(slope = 1, color = 'black')+
  geom_point(size = 0.2) +
  geom_text(data = dynamic.range%>% filter(!source == 'OURS') , aes(label = paste0('n=',n), x = 8E3, y = 1E1), size = 8, color = 'black') +
  scale_x_log10(breaks = c(1E0,1E3,1E6),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(1E0,1E3,1E6),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(legend.position= 'none')+ 
  ylab('median copy number of SGD reference') + xlab('')+ scale_color_brewer(palette = 'Set1')+
  ggtitle('red line indicate our measurement')+ facet_wrap(.~source, nrow = 2)+
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()

pdf('compare_ref_proteome/compare_dynamic_range_n.pdf', w= 8, h = 4)
ggplot(dynamic.range, aes(x = n/1000, y = DR,  color = method, size = (source == 'OURS'))) +
  geom_point(alpha = 0.7) +scale_y_log10(breaks = c(1E0,1E3,1E6), limits = c(1,1E7),
                              labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(limits = c(0,5.1), expand = expansion(mult = c(0, .1)))+
  scale_color_brewer(palette = 'Set1')+
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()
# compare to median of datasets with good dynamic range
ref.median.filtered <- ref.long %>% filter(!source == 'OURS') %>% right_join(dynamic.range %>% filter(DR > 0) %>% select(source)) %>%
  group_by(rank.mySC, copy.mySC, geneID) %>% summarise(median.copy = median(copy, na.rm = T)) %>% ungroup()

pdf('compare_ref_proteome/median_proteome_compare_dynamic_range_greater_than_2000.pdf', w= 3.5, h = 3.5)
ggplot(ref.median.filtered,
       aes(x = rank.mySC, y = median.copy)) +
  geom_point(size = 0.2) +
  geom_line(aes(x = rank.mySC, y = copy.mySC), color = 'red', size = 0.5) +
  scale_y_log10(breaks = c(1,1E2,1E4,1E6)) +theme_classic() +theme(legend.position= 'none')+
  ylab('median copy number of SGD reference') + xlab('')+
  ggtitle('red line indicate our measurement')+
  egg::theme_presentation() + theme(panel.background = element_blank())
dev.off()

# plot subgroup
pdf('compare_ref_proteome/glycolysis.pdf', w= 7, h = 8)
ggplot(ref.long %>% filter(group == 'glycolysis')%>% filter(!source == 'OURS'),
       aes(x = as.factor(rank.mySC), y = copy, color = method)) +
  geom_point(size = 0.4) +
  geom_point(data = ref.long %>% filter(group == 'glycolysis')%>% filter(!source == 'OURS') %>%
               group_by(method, geneID, rank.mySC) %>% summarise(copy = mean(copy)) %>% ungroup(),
             size = 2, pch = 22, fill = NA) +
  geom_line(aes(x = as.factor(rank.mySC), y = copy.mySC, group = method), color = 'black') +
  facet_wrap(method~.,nrow=2, scales = 'free_x')+
  scale_y_log10() +
  egg::theme_presentation() +
  theme(axis.text.x= element_blank(),axis.ticks.x = element_blank(),legend.position= 'none', panel.background = element_blank()) +
  xlab('glycolytic protein') +ggtitle('black line: our measurement')
dev.off()

ref.glyc <- ref.long %>% filter(group == 'glycolysis', !source == 'OURS')
pdf('compare_ref_proteome/glycolysis_scatter.pdf', w= 3.5, h = 3.5)
ggplot(ref.glyc, aes(x = copy.mySC, y = copy)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_point(data = ref.glyc %>% select(median.copy, copy.mySC) %>% distinct(),
             aes(y = median.copy), size = 1.5, alpha = 0.8, color = 'coral')+
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10(breaks = c(1E3,1E5,1E7), limits = c(1E2,1E7)) +
  scale_x_log10(breaks = c(1E3,1E5,1E7), limits = c(1E2,1E7)) +
  theme_light() + annotate('label', label = 'y=x', x = 10^6.6, y = 10^6.6)+
  xlab('copy number from this study') +
  ylab('copy number in reference proteomics') + 
  ggtitle('glycolytic proteins\nred points indicate median of SGD')
dev.off()

# ref.median <- ref.long %>% filter(!is.na(median.copy)) %>%
#   group_by(source) %>% summarize(n.total = length(!is.na(copy)), median = round(median(copy)),
#                                  median.ref = round(median(median.copy)), median.ours = round(median(copy.mySC, na.rm = T))) %>% ungroup()
# Considering median normalization for proteins detected

dt.plot.frac.median <- ref.long %>% filter(source == 'OURS') %>% select(-source) %>%
  pivot_longer(cols = c(copy,median.copy), names_to = 'source', values_to = 'copy') %>%
  mutate(source = ifelse(source == 'copy', 'OURS', 'SGD')) %>%
  filter(!is.na(copy), !is.na(Length)) %>%
  group_by(source) %>% mutate(f = copy *Length/sum(copy*Length)) %>% ungroup() %>%
  group_by(group, source) %>%
  summarise(fraction.proteome = sum(f), n = length(f)) %>% ungroup() %>%
  mutate(group = factor(group, levels = order.group))
n.glyc <- ref.long %>% filter(group == 'glycolysis') %>% group_by(source) %>%
  summarise(n = length(source)) %>% ungroup()
n.all <- ref.long %>% group_by(source) %>%
  summarise(n = length(source)) %>% ungroup()

dt.plot.frac <- ref.long %>% # right_join(n.glyc %>% filter(n>20) %>% select(source)) %>%
  filter(!is.na(copy)) %>%
  group_by(source) %>%
  mutate(f = copy *Length/sum(copy*Length)) %>% ungroup() %>%
  group_by(group, source, method) %>%
  summarise(fraction.proteome = sum(f), n = length(f)) %>% ungroup() %>%
  mutate(group = factor(group, levels = order.group))

pdf('compare_ref_proteome/reference_SGD_median_protein_resource_allocation.pdf', w = 5.8, h = 4)
ggplot(dt.plot.frac.median, aes(x = source, y = fraction.proteome, fill = group)) +
  geom_col(position = position_stack(), width = 0.7, color = 'grey40') +
  scale_fill_manual(values = color.group)+ xlim('SGD','OURS')+
  egg::theme_presentation()+theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5))
dev.off()

pdf('compare_ref_proteome/reference_SGD_protein_resource_allocation.pdf', w = 10, h = 4)
ggplot(dt.plot.frac %>% filter(!source == 'OURS') %>% arrange(method) %>% 
         mutate(source = factor(source, levels = unique(source))),
       aes(x = source, y = fraction.proteome, fill = group)) +
  geom_col(position = position_stack(), width = 0.7, color = 'grey40') +
  scale_fill_manual(values = color.group) +
  egg::theme_presentation()+theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5))
dev.off()

ggplot(dt.plot.frac %>% filter(group %in% c('glycolysis','mitochondria')),
                 aes(x = source, y = fraction.proteome)) +
  geom_col(width = 0.7, position = position_dodge())+ facet_grid(group~., scales = 'free')

# compare proteome allocation from our measurement to database
dt.plot <- rbind(frac.p %>% filter(strain == 'CENPK') %>% select(fraction.proteome,group,n) %>%
                   mutate(source = 'OURS',method = 'mass spectrometry'),
                 dt.plot.frac %>% filter(!source == 'OURS'))
# significantly different?
tmp.anova <- dt.plot %>% mutate(ours = (source == 'OURS')) %>% group_by(group) %>%
  do(model = anova(lm(fraction.proteome ~ ours, data = .))) %>% broom::tidy(model)

# proteome allocation derived from median protein copy number of 19 references vs mean of our replicates
dt.plot.median <- frac.p %>% filter(strain == 'CENPK') %>% 
  group_by(group) %>% summarise(f.e = sd(fraction.proteome)/sqrt(length(group))) %>% ungroup() %>% mutate(source = 'OURS') %>%
  full_join(dt.plot.frac.median)

pdf('compare_ref_proteome/reference_SGD_protein_resource_allocation_bar_scatter.pdf', w = 14, h = 8)
ggplot(dt.plot,
       aes(x = group, y = fraction.proteome, fill = (source == 'OURS'), group = (source == 'OURS'))) +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_col(data = dt.plot.median, width = 0.6, color = 'black', position= position_dodge(width = 0.6)) +
  geom_errorbar(data = dt.plot.median, aes(ymin = fraction.proteome-f.e,ymax = fraction.proteome+f.e), 
                position = position_dodge(width = 0.6),width = 0.4)+
  geom_point(pch = 21, fill = NA, aes(color = grepl('mass',method)),position = position_dodge(width = 0.6), 
             alpha = 0.8, size= 1.5) +
  scale_fill_manual(values =  c('TRUE' =  'grey','FALSE'='white')) +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'red')) +
  ggtitle(paste0('fraction of proteome, normalized within each dataset,\n',
                 'large circles represent median of all the SGD references')) +
  ylab('proteome fraction') +
  scale_y_continuous(breaks = 0.1*(0:6))+
  scale_x_discrete(limits = rev(order.group))+
  egg::theme_presentation()+theme(panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle= 90, hjust =1, vjust = 0.5))

dev.off()

#### reference ribosome profiling ####
# http://sysbio.gzzoc.com/rpfdb/download.html
dt.rbp <- read.delim('../ref_Sacchromyces_RPF/GSE56622.txt', sep = '\t') %>% 
  select(Gene_ID, Gene_Name,SRX543755,SRX513393) %>% mutate(RBP = (SRX543755+SRX513393)/2) %>% 
  # read.delim('../ref_Sacchromyces_RPF/GSE55400.txt', sep = '\t') %>% 
  # select(Gene_ID, Gene_Name,SRX476345) %>% mutate(RBP = SRX476345) %>% 
  rename(geneID = Gene_ID) %>%
  full_join(mySC.rank) %>%
  mutate(AbsProt = OURS.f*tAA.SGD/Length)

ggplot(dt.rbp , aes(x = log10(AbsProt), y = log10(RBP))) +
  geom_point(alpha = 0.6, color = 'grey') +
  geom_abline(slope = 1, intercept = 0) +
  # geom_smooth(method = 'lm', color = 'black', linetype = 'dashed', size = 0.5) +
  geom_point(data = dt.rbp %>% filter(group == 'glycolysis'),pch=21, fill = '#FFED6F', color = 'gold3') +
  geom_point(data = dt.rbp %>% filter(group %in% c('TCA','ox phos')),pch=21, fill = '#6BAED6', color = '#4292C6') +
  scale_y_continuous(breaks = c(0,2,4,6), limits = c(-1,7)) +
  scale_x_continuous(breaks = c(0,2,4,6), limits = c(-1,7)) +
  egg::theme_presentation()

summary(lm(log10(RBP)~log10(AbsProt), dt.rbp %>% filter(!is.na(RBP),!is.na(AbsProt),is.finite(AbsProt/RBP))))
summary(lm(log10(RBP)~log10(AbsProt), 
           dt.rbp %>% filter(!is.na(RBP),!is.na(AbsProt),is.finite(AbsProt/RBP)) %>%
             filter(group %in% c('glycolysis','TCA','ox phos'))))


dt.rbp.frac <- dt.rbp %>% mutate(RBP.f = RBP*Length/sum(RBP*Length, na.rm=T), OURS.f = OURS.f/sum(OURS.f,na.rm=T)) %>%
  select(geneID, group, RBP.f, OURS.f) %>%
  pivot_longer(cols = c(RBP.f, OURS.f), names_to = 'source', values_to = 'f') %>%
  group_by(source, group) %>% summarise(f = sum(f, na.rm = T)) %>% ungroup() %>%
  mutate(group = factor(group, levels = order.group))

ggplot(dt.rbp.frac, aes(x = source, y = f, fill = group)) +
  geom_col(position = position_stack(), width = 0.7, color = 'grey40') +
  scale_fill_manual(values = color.group)+ 
  egg::theme_presentation()+theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5))
