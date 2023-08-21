## for generating mapping files to annotate protein Entry with
# GO, geneID, Length, Protein.names (based on uniprot annotation)
# pathway (based on gpr in GEM)
# group (based on GO, Protein.names, pathway)

library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)
library(tidyr)
rm(list=ls())
options(stringAsFactors = F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('assign_group.R')

#### improve annotation using orthologs ####
orthologs <- read.delim('../../fastas/ortholog_IO_SC.tsv', sep = '\t')

##### S.cerevisiae #####
anno.subsystem.SC <- read.delim('../../abs_proteomics/subsystem_annotation_SC.tsv', sep = '\t')
tmp <- read.xlsx('../../models/yeastGEM.xlsx') %>%
  rename(EC = `EC-NUMBER`, geneID = GENE.ASSOCIATION, subsystem = SUBSYSTEM) %>%
  select(geneID, subsystem, EC)
gpr.SC<- tibble(geneID = str_extract_all(tmp$gene, '[A-Z][0-9A-Z]+'),EC= tmp$EC, subsystem = strsplit(tmp$subsystem, split = ';')) %>%
  unnest(cols = geneID) %>%
  unnest(cols = subsystem) %>%
  left_join(anno.subsystem.SC) %>%
  mutate(category = ifelse(is.na(category), 'other', as.character(category))) %>%
  select(-subsystem) %>% distinct() %>%
  group_by(geneID) %>% summarise(pathway = paste(category, collapse = ';')) %>% ungroup()

tmp <- read.delim('../../abs_proteomics/annotation_SC.tab', sep = '\t') %>%
  rename(geneID = Gene.names...ordered.locus..) %>% select(-Gene.ontology.IDs) %>%
  pivot_longer(cols = starts_with('Gene.ontology'), names_to = 'GO.category',values_to = 'GO') %>%
  mutate(Entry = as.character(Entry),GO = as.character(GO)) %>% 
  select(Entry, geneID, Length, GO,Entry.name, Protein.names) 

anno.SC <- tibble(Entry=tmp$Entry, GO=strsplit(tmp$GO, split = ';'),
                  geneID = str_extract_all(tmp$geneID, '[RYQ][^;]+')) %>%
  unnest(cols= GO, keep_empty = T) %>% unnest(cols = geneID, keep_empty = T) %>%
  mutate(GO = gsub(' $','',gsub('^ ','',GO))) %>% distinct() %>%
  mutate(GO = ifelse(is.na(GO),'',GO)) %>%
  group_by(Entry, geneID) %>%
  summarise(GO = paste(GO, collapse = ';')) %>% ungroup() %>%
  right_join(tmp %>% select(Entry, Entry.name, Length, Protein.names) %>% distinct())

map.SC <- assign.group(anno.SC %>% select(-Entry.name) %>% left_join(gpr.SC))
write.table(map.SC, 'map_group_SC.tsv', sep = '\t', row.names = F, quote = F)

##### I. orientalis #####
## improve annotation using orthologs
orthologs <- read.delim('../../fastas/ortholog_IO_SC.tsv', sep = '\t') %>%
  group_by(pid.IO) %>% filter((e.fwd + e.bwd) == min(e.fwd + e.bwd)) %>%
  filter((identity.fwd + identity.bwd) == min(identity.fwd + identity.bwd)) %>% ungroup()

anno.subsystem.IO <- read.delim('../../abs_proteomics/subsystem_annotation_IO.tsv', sep = '\t')
gpr.IO <- read.delim('../../../code/integration/extfiles/gene_subsystem.tsv',
                     sep = '\t', stringsAsFactors = F) %>% rename(geneID = gene) %>%
  left_join(anno.subsystem.IO) %>%
  mutate(category = ifelse(is.na(category), 'other', as.character(category))) %>%
  select(-rxnname) %>% distinct() %>%
  group_by(geneID) %>% summarise(pathway = paste(category, collapse = ';')) %>% ungroup()

tmp <- read.delim('../../abs_proteomics/annotation_IO.tab', sep = '\t', stringsAsFactors = F) %>% select(-Gene.ontology.IDs) %>%
  pivot_longer(cols = starts_with('Gene.ontology'), names_to = 'GO.category',values_to = 'GO') %>%
  mutate(Entry = as.character(Entry),GO = as.character(GO)) %>% 
  select(Entry, Length, GO, Gene.names,Protein.names)

anno.IO <- tibble(Entry=tmp$Entry, Length = tmp$Length, GO=strsplit(tmp$GO, split = ';')) %>%
  unnest(cols= GO) %>%
  mutate(GO = gsub(' $','',gsub('^ ','',GO)))%>%
  distinct() %>%
  group_by(Entry,Length) %>%
  mutate(GO = ifelse(is.na(GO),'',GO)) %>%
  summarise(GO = paste(GO, collapse = ';')) %>% ungroup() %>%
  right_join(tmp %>% select(Entry, Length, Gene.names, Protein.names) %>% distinct())

# assign group by ortholog mapping to SC, for genes without ortholog mapping, assign based on Io annotation
tmp <- anno.IO %>% mutate(geneID = str_extract(Gene.names, 'JL09\\_g[0-9]+')) %>%
  select(-Gene.names) %>% left_join(gpr.IO) 
# ortholog exist
tmp.orth <- tmp %>%
  inner_join(orthologs %>% select(pid.SC, pid.IO,pname.IO,pname.SC)%>% rename(Entry = pid.IO)) %>%
  mutate(Protein.names = ifelse(is.na(Protein.names),pname.IO,Protein.names)) %>%
  mutate(Protein.names = ifelse(grepl('[uU]ncharacterize',Protein.names),paste0('ortholog of S.cere ',pname.SC),Protein.names)) %>%
  left_join(map.SC %>% rename(pid.SC = Entry,GO.SC = GO) %>% select(pid.SC, group, GO.SC)) %>%
  mutate(GO = ifelse(is.na(GO), GO.SC, GO)) %>%
  select(-pname.IO,-pname.SC,-pid.SC,-GO.SC) %>%
  group_by(Entry, Length, GO, geneID, pathway, group) %>%
  summarise(Protein.names = paste0(Protein.names, collapse = ' or ')) %>% ungroup()

tmp.uniq <- assign.group(tmp %>% anti_join(orthologs %>% rename(Entry = pid.IO) %>% select(Entry)))
map.IO <- rbind(tmp.orth, tmp.uniq)
write.table(map.IO, 'map_group_IO.tsv', sep = '\t', row.names = F, quote = F)

order.group <- c('unannotated','other','transcription','translation/ribosome','other enzymes',
                 'amino acid','nucleic acid','lipid','one carbon','PPP','mitochondrial enzymes','mitochondria','TCA','ox phos','glycolysis')
# color.group <- c(RColorBrewer::brewer.pal(12, 'Set3')[c(2,9,4,1,7,5,6,3,8)],'#C6DBEF','#4292C6','#6BAED6','#FFED6F')
color.group <- c(RColorBrewer::brewer.pal(12, 'Set3')[c(2,9,11,1,8,8,8,8,8,8,8)],'#C6DBEF','#6BAED6','#4292C6','#FFED6F')
names(color.group) <- order.group

save(list=c('map.IO','map.SC','order.group','color.group'), file = 'assign_group_files.RData')
