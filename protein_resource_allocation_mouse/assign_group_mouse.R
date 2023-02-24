rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('assign_group_files/from_mouse_proteomics.RData')
# load('data_and_mapping_mouse.RData')
annotation <- read.delim('mouse_annotation_uniprotKB.tsv', sep = '\t', stringsAsFactors = F) %>%
  rename(Gene.names = Gene.Names, Gene.ontology..GO. = Gene.Ontology..GO.)

dt.tmp <- dt.tumor.proteome %>% select(Entry,f,sample) %>% distinct() %>% 
  group_by(Entry) %>% summarise(f = sum(f)/15)  %>% ungroup()

KO.mapping.mouse <- read.delim( '../proteomap/KO_mouse_mapping.tsv',sep = '\t',stringsAsFactors = F)
mod.rxn <- read.delim('../models/Mouse-GEM/mouse_reactions_exported.txt', sep = '\t', stringsAsFactors = F) %>%
  as_tibble() %>% mutate(gene = strsplit(grRules, split = ' ')) %>%
  unnest(cols = 'gene') %>% filter(!gene %in% c('and','or')) %>% mutate(gene = gsub('[()]','',gene)) %>% select(-grRules)
anno.gene.long <- annotation %>% 
  select(Entry, Protein.names, Length, Gene.names,Gene.ontology..GO.,Pathway) %>% distinct() %>%
  as_tibble() %>% mutate(gene = strsplit(as.character(Gene.names), split = ' ')) %>% unnest(cols = 'gene') 

# mito proteins
mito.carta <- readxl::read_excel('Mouse.MitoCarta2.0.xls',sheet = 'A Mouse MitoCarta2.0')
mito.Entry <- anno.gene.long %>% filter(gene %in% mito.carta$Symbol) %>% select(Entry) %>% distinct()

# enzymes
enzyme.Entry <- anno.gene.long %>% right_join(mod.rxn) %>% select(Entry) %>% distinct()

anno.group <- anno.gene.long %>%
  left_join(mod.rxn %>% filter(!is.na(gene)))%>%
  left_join(KO.mapping.mouse %>% filter(!is.na(gene))) %>%
  select(-gene) %>% distinct() %>%
  group_by(Entry, Protein.names, Length, Gene.names,Gene.ontology..GO.,Pathway) %>%
  summarise(subSystems = paste(ifelse(is.na(subSystems),'',subSystems), collapse = ';'),
            rxnNames = paste(ifelse(is.na(rxnNames),'',rxnNames), collapse = ';'),
            proteomap = paste(ifelse(is.na(lv3),'',lv3), collapse = ';')) %>% ungroup() %>%
  mutate(subSystems = gsub('[;]+',';',subSystems),
         rxnNames = gsub('[;]+',';',rxnNames),
         proteomap = gsub('[;]+',';',proteomap)) %>%
  mutate(group = case_when(grepl('glycolysis',Pathway)|grepl('pyruvate fermentation',Pathway)|
                             (grepl('glycolytic process',Gene.ontology..GO.) & grepl('Glycolysis',subSystems))|
                             grepl('[gG]lucose transporter', Protein.names)~ 'glycolysis',
                           grepl('tricarboxylic', Pathway)|grepl('tricarboxylic acid cycle',Gene.ontology..GO.)|grepl('Tricarboxylic acid cycle',subSystems) ~ 'TCA',
                           (grepl('oxidative phosphorylation', Pathway)|grepl('mitochondrial electron transport',Gene.ontology..GO.)|
                              grepl('ATP synthase', Gene.ontology..GO.))|grepl('respiratory chain',Gene.ontology..GO.)|grepl('Oxidative phosphorylation',subSystems) ~ 'ox phos',
                           grepl('Fatty acid oxidation',subSystems)|grepl('Beta oxidation',subSystems) ~ 'FA oxidation',
                           grepl('tetrahydrofolate interconversion',Gene.ontology..GO.)|grepl('one-carbon metabolic process',Gene.ontology..GO.) ~ 'one carbon',
                           grepl('pentose-phosphate',Gene.ontology..GO.)|grepl('ribose-phosphate',Gene.ontology..GO.) | (grepl('pentose phosphate pathway',Pathway)) ~ 'PPP',
                           grepl('[lL]ipid',Pathway)|grepl('steroid',Gene.ontology..GO.)|proteomap == 'Lipid and steroid metabolism' ~ 'lipid' ,
                           grepl('[aA]mino-acid',Pathway)|grepl('[aA]mino acid',Gene.ontology..GO.)|grepl('Amino acid',proteomap) ~ 'amino acid',
                           grepl('[pP]urine',Pathway)|grepl('[pP]yrimidine',Pathway)|grepl('[nN]ucleo',Pathway)|
                             grepl('[ATCUG]TP biosynthetic process',Gene.ontology..GO.)|grepl('Nucleotide metabolism',subSystems) ~ 'nucleic acid',
                           (Entry %in% enzyme.Entry$Entry) &(Entry %in% mito.Entry$Entry) ~'mitochondrial other enzymes',
                           Entry %in% enzyme.Entry$Entry ~'other enzymes',
                           Entry %in% mito.Entry$Entry ~ 'mitochondria',
                           grepl('[rR]ibosom',Gene.ontology..GO.)|grepl('[tT]ranslation',Gene.ontology..GO.)|grepl('tRNA',Gene.ontology..GO.) ~ 'translation/ribosome',
                           grepl('protein stabilization',Gene.ontology..GO.)|grepl('chaperon',Gene.ontology..GO.) ~'chaperon',
                           grepl('[tT]ranscription',Protein.names)|grepl('[tT]ranscription',Gene.ontology..GO.)|grepl('mRNA splicing',Gene.ontology..GO.) ~ 'transcription',
                           grepl('chromosome',Gene.ontology..GO.)|grepl('nucleosome',Gene.ontology..GO.) ~ 'chromosome',
                           grepl('[eE]ndoplasmic reticulum',Gene.ontology..GO.) ~ 'ER',
                           grepl('[gG]olgi',Gene.ontology..GO.) ~ 'Golgi',
                           grepl('cytoskeleton',Gene.ontology..GO.) ~ 'cytoskeleton',
                           TRUE ~ 'other'
  )) %>%
  mutate(group = ifelse((Entry %in% mito.Entry$Entry)&(group == 'amino acid'), paste0('mitochondrial ',group),group))

order.group <- c('other','cytoskeleton','ER','Golgi','chromosome','transcription','chaperon','translation/ribosome',
                 'other enzymes','amino acid','nucleic acid','lipid','one carbon','PPP',
                 'mitochondria','mitochondrial other enzymes','mitochondrial amino acid',
                 'FA oxidation','TCA','ox phos','glycolysis')
color.group <- c(RColorBrewer::brewer.pal(12, 'Set3')[c(9,2,2,2)],
                 "#FED9A6","#FED9A6",
                 RColorBrewer::brewer.pal(12, 'Set3')[c(1,1,3,3,3,3,8,8)],
                 '#DEEBF7','#C6DBEF','#9ECAE1',
                 '#9ECAE1','#6BAED6','#4292C6','#FFED6F')
# color.group <- c(RColorBrewer::brewer.pal(12, 'Set3')[c(9,11,1,8,8,8,8,8,8,8)],'#C6DBEF','#6BAED6','#4292C6','#FFED6F')
names(color.group) <- order.group
anno.group <- anno.group %>% mutate(group = factor(group, levels = order.group))

# plot color scale
# ggplot(data.frame(g=factor(order.group,levels = order.group), f = 1), aes(x = 1, y = f,fill = g)) +
#   geom_col(position = position_stack(), color = 'black')+scale_fill_manual(values = color.group)

save(list = c('color.group','anno.group'),file = 'assign_group_mouse_files.RData')

assign.tmp <- dt.tmp %>% left_join(anno.group %>% select(Entry,Protein.names, group) %>% distinct()) %>% filter(!is.na(group))

dt.plot <- assign.tmp %>% group_by(group) %>% summarise(f = sum(f)) %>% ungroup()

ggplot(dt.plot %>% mutate(group = factor(group, levels = order.group)),
       aes(x = 1, y = f, fill = group)) + geom_col(position = position_stack(), color = 'black') +
  scale_fill_manual(values = color.group)

assign.tumor <- dt.tumor.proteome %>% select(sample,Entry,Protein.names, f) %>%
  left_join(anno.group %>% select(Entry,Protein.names, group) ) %>% distinct()

dt.plot <- assign.tumor %>% group_by(group,sample) %>% summarise(f = sum(f)) %>% ungroup()
dt.plot.mean <- dt.plot %>% group_by()
ggplot(dt.plot %>% mutate(group = factor(group, levels = order.group)),
       aes(x = sample, y = f, fill = group)) + geom_col(position = position_stack(), color = 'black') +
  scale_fill_manual(values = color.group)



## pathway view
pth <- 'glycolysis'
list.pth <- anno.group %>% select(Entry,Protein.names, group) %>% distinct() %>% 
  filter(group == pth) %>%
  left_join(anno.gene.long %>% inner_join(mod.rxn) %>% 
              group_by(Entry, gene) %>% select(Entry, gene, rxnNames) %>% distinct() %>%
              summarise(rxnNames = paste(rxnNames,collapse = ';')) %>% ungroup())
write.table(list.pth, paste0('assign_group_files/manual_annotation_yeast_reactions_',pth,'.tsv'), sep = '\t', quote = F, row.names = F)

pth.anno2yeast <- read.xlsx(paste0('assign_group_files/manual_annotation_yeast_reactions_',pth,'.xlsx'))
pth.order <- read.xlsx(paste0('assign_group_files/manual_annotation_yeast_reactions_',pth,'.xlsx'), sheet = 2)

dt.pth <- dt.tmp %>% right_join(pth.anno2yeast %>% select(Entry,rxn.Yeast) %>% distinct()) %>%
  group_by(rxn.Yeast) %>% arrange(desc(f)) %>% mutate(isozyme = as.character(1: length(f))) %>% ungroup()

# load yeast data
# load('../manuscript/flux_enzyme_correlation_20220128.RData') # flux.protein.all
# save(flux.protein.all, file = 'flux_enzyme_20220128_CENPK.RData')
load('flux_enzyme_20220128_CENPK.RData')

dt.pth.yeast <- flux.protein.all %>% filter(strain == 'CENPK') %>% select(Entry,id,f) %>% distinct() %>%
  filter(id %in% pth.anno2yeast$rxn.Yeast) %>%
  group_by(id) %>% arrange(desc(f)) %>% mutate(isozyme = as.character(1: length(f))) %>% ungroup()

dt.pth.plot <- rbind(dt.pth.yeast %>% select(Entry, isozyme, id, f) %>% rename(rxn.Yeast = id) %>% mutate(organism = 'S. cerevisiae'),
                     dt.pth %>% mutate(organism = 'mouse'))
pdf(paste0('compare_mouse_to_yeast/',pth,'.pdf'), w = 6, h = 5)
ggplot(dt.pth.plot, aes(x = rxn.Yeast, y = f, fill = isozyme)) +
  geom_col(position = position_stack(reverse = T), color = 'black')+
  facet_grid(organism~., scales = 'free_y') +
  scale_x_discrete(limits = pth.order$id) +
  scale_fill_brewer(palette = 'Blues') +
  egg::theme_presentation() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 12))
dev.off()

# Brenda info
par.brenda <- read.delim('../../code/integration/brenda/all_km_archive.tsv', sep = '\t') %>% filter(organism %in% c('Homo sapiens','Saccharomyces cerevisiae','Mus musculus'))
ec.fil <- '4.2.1.11' # ENO
sub.fil <- c('2-phospho-D-glycerate','2-phosphoglycerate')
ec.fil <- '4.1.2.13' # FBA
sub.fil <- 'D-fructose 1,6-bisphosphate'
ec.fil <- '2.7.1.40' # PYK
sub.fil <- 'phosphoenolpyruvate'
ec.fil <- '5.3.1.1' # TPI
sub.fil <- 'D-glyceraldehyde 3-phosphate'
ec.fil <- '2.7.1.1' # HEX
sub.fil <-'ATP'

tmp <- par.brenda %>% filter(ecNumber == ec.fil, substrate %in% sub.fil) %>% filter(kmValue > 0)
ggplot(tmp, aes(x = organism, y = kmValue)) + geom_jitter(width =0.2)

# the most abundant one in the unmapped
a <- dt.tmp %>% left_join(anno.group %>% select(Entry,Protein.names, Gene.ontology..GO., group,lv1,lv2,lv3) %>% distinct())  %>%
  filter(group == 'other') %>% filter(f == max(f)) %>% select(Entry, Protein.names, f,Gene.ontology..GO.,lv1,lv2,lv3)
a$Gene.ontology..GO.
a$Protein.names

b <- uniprot2string <- unnest(tibble::tibble(Entry = anno.group$Entry, 
                                             GO = strsplit(as.character(anno.group$Gene.ontology..GO.), split = '; ')),
                              cols = c(GO)) %>%
  distinct() %>% inner_join(KO.mapping.mouse %>% select(Entry, Gene.names,lv1,lv2,lv3) %>% filter(!is.na(Entry),!is.na(lv1))) %>%
  group_by(lv1,lv2,lv3,GO) %>% summarise(n = length(GO)) %>% ungroup() %>% filter(lv1 == 'Metabolism')
