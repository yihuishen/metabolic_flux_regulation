
require(dplyr)

assign.group <- function(df) { # based on GO, Protein.names, pathway
  df.assigned <- df %>%
    mutate(group = case_when(grepl('TCA',pathway) ~ 'TCA',
                             grepl('ox phos',pathway) ~ 'ox phos',
                             (grepl('one-carbon metabolic process',GO) | grepl('lipid',GO)|(pathway == 'lipid') |
                              grepl('nucleic acid',pathway) | grepl('amino acid',GO)|grepl('amino acid',pathway) | !is.na(pathway)) &
                               (grepl('[mM]itochondria',Protein.names) | grepl('mitochondrial translation', GO))~ 'mitochondrial enzymes',
                             grepl('[mM]itochondria',Protein.names) | grepl('mitochondrial translation', GO)~ 'mitochondria',
                             (grepl('glycolytic',GO) & (grepl('glycolysis',pathway)|is.na(pathway)))|grepl('fermentation',GO)  ~ 'glycolysis',
                             # (grepl('gluconeogenesis',GO) & (grepl('gluconeogenesis',pathway)|is.na(pathway)))  ~ 'gluconeogenesis',
                             grepl('translation',GO)|grepl('ribosom',GO)|grepl('ribosom',Protein.names) ~ 'translation/ribosome',
                             grepl('transcription',GO)|grepl('transcription',Protein.names)|grepl('chromatin',GO)|grepl('tRNA',GO) ~ 'transcription',
                             grepl('mitochondria',GO) ~ 'mitochondria',
                             grepl('pentose-phosphate',GO) | (pathway == 'PPP') ~ 'PPP',
                             grepl('one-carbon metabolic process',GO) ~ 'one carbon',
                             grepl('lipid',GO)|(pathway == 'lipid') ~ 'lipid',
                             grepl('nucleic acid',pathway) ~ 'nucleic acid',
                             grepl('amino acid',GO)|grepl('amino acid',pathway) ~ 'amino acid',
                             !is.na(pathway) ~ 'other enzymes',
                             ((GO == '') | is.na(GO)) & (is.na(Protein.names)|grepl('[uU]ncharacterized',Protein.names)) ~ 'unannotated',
                             TRUE ~ 'other'))
  return(df.assigned )
}

