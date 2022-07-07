## identify statistically significant metabolic regulator and evaluate metabolic leverage

library(dplyr)
library(stringr)
library(egg)
library(rmatio)
library(tidyr)

### 
options(stringsAsFactors = F)
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if(file.exist('report_SIMMER.RData')) {  # load existing data for plot
  load('report_SIMMER.RData')
} else {
  load('report_rxnForms.RData')
  rxnForms <- rxnForms.report
  rm(rxnForms.report)
  rxn.report <- read.delim('report_rxns.tsv', sep = '\t') %>% rename(rxnid.GEM = rxnid, rxnid = rxnidx)
  
} 


#### create summary table from NNLS-MCMC result ####
n <- length(rxnForms)
sum_tb <- data.frame(rxnid = character(length = n),
                     rxnname = character(length = n),
                     metidx  = character(length = n),
                     reg = character(length = n),
                     rMech = character(length = n),
                     mode = character(length = n),
                     maxlik = numeric(length = n),
                     EC = numeric(length = n),
                     r = numeric(length = n),
                     r.enzyme = numeric(length = n),
                     npar = numeric(length = n)) # summarize all rMech
for (k in 1:n) {
  rxnSummary <- rxnForms[[k]]
  if (length(rxnSummary$fit) == 0) {
    next
  }
  fit <- rxnSummary$fit
  sum_tb$rxnid[k] <- rxnSummary$rxnID
  sum_tb$rxnname[k] <- rxnSummary$reaction
  sum_tb$EC[k] <- rxnSummary$EC
  reg.idx <- rxnSummary$rxnFormData$SubstrateID[!rxnSummary$rxnFormData$Subtype %in% c('substrate','product')]
  sum_tb$metidx[k] <- ifelse(length(reg.idx) == 0, NA, unname(reg.idx))
  sum_tb$reg[k] <- ifelse(length(reg.idx) == 0, NA, unname(rxnSummary$metNames[reg.idx]))
  sum_tb$rMech[k] <- rxnSummary$listEntry
  reg.mode <- rxnSummary$rxnFormData$Subtype[!rxnSummary$rxnFormData$Subtype %in% c('substrate','product')]
  sum_tb$mode[k] <- ifelse(length(reg.idx) == 0, NA, reg.mode)
  
  sum_tb$maxlik[k] <- max(fit$lik_track)
  sum_tb$r[k] <- fit$r[which(fit$lik_track == max(fit$lik_track))[1]]
  sum_tb$npar[k] <- sum(rxnSummary$kineticPars$measured, na.rm = T) + 2
  enz <- rxnSummary$omic_data$enzyme_abund
  sum_tb$r.enzyme[k] <- cor(2^rowMeans(enz[,!is.infinite(enz[1,]),drop = F]),rowMeans(rxnSummary$omic_data$flux[,2:3]))
}

#### identification of meaningful regulation ####

# calculate p value from likelihood ratio test between regulatory mechanism with and without activator/inhibitor
# obtain q value with B-H FDR correction

sum_tb_p <- sum_tb %>% filter(!is.na(mode)) %>%
  left_join(sum_tb %>% filter(grepl('-rm$',rMech)) %>% select(rxnid, maxlik) %>% rename(maxlik.woReg = maxlik)) %>%
  mutate(p = 1-pchisq(2*(maxlik-maxlik.woReg), df=1)) %>%
  mutate(q = p.adjust(p, 'fdr'))

# meaningful regulation at FDR = 0.1

reg.FDR <- sum_tb_p %>% filter(q < 0.1)%>%
  mutate(modtype = ifelse(grepl('act',rMech),'act','inh'))%>%
  select(-rMech,-mode) %>%
  group_by(rxnid, rxnname,reg, modtype) %>%
  filter(maxlik == max(maxlik)) %>% ungroup() %>% distinct()

message(paste('# total rMech supported by data, ',nrow(reg.FDR)))
message(paste('# reactions predicted to have regulator, ',length(unique(reg.FDR$rxnid))))

# report, with reference to occurance in Brenda

load('report_all_affinities_brenda.RData')
affn.brenda <- all_affinities %>%
  filter(speciesType == 'regulator',
         EC %in% reg.FDR$EC) %>%
  select(ligandname, EC, modtype, isYeast, nQuant, nQual, metidx) %>%
  mutate(n.report = nQuant + nQual) %>% select(-nQuant,-nQual)

reg.total <- affn.brenda %>%
  # group_by(EC, modtype, isYeast) %>%
  group_by(EC, modtype) %>%
  summarise(n.total = sum(n.report)) %>% ungroup()

reg.predict <- reg.FDR  %>%
  left_join(affn.brenda) %>%
  # group_by(EC, metidx, modtype, isYeast) %>%
  group_by(EC, metidx, modtype) %>%
  summarise(n.predict = sum(n.report)) %>% ungroup()

regulator.report <- reg.FDR %>%
  left_join(reg.predict) %>%
  left_join(reg.total)

regulator.report.table <- regulator.report %>% 
  mutate(reg = paste0(reg,' (',n.predict,'/',n.total,', p = ',
                      format(q, digits = 2, scientific = T),')')) %>%
  group_by(EC, modtype,rxnname) %>%
  summarise(reg = paste0(reg, collapse = '; ')) %>% ungroup()

write.xlsx(regulator.report.table,'report_regulator_table.xlsx')

#### best supported reaction kinetics ####
rMech.best <- rbind(sum_tb %>% filter(!grepl('woE',rMech), is.na(reg)) ,
                    sum_tb_p %>% filter(q<0.1) %>% select(-maxlik.woReg, -p,-q)) %>% group_by(rxnid) %>%
  filter(r == max(r)) %>% ungroup() %>% select(-r.enzyme)

#### Goodness of fit accounting for different reaction components ####
# Pearson's R
r <- rbind(sum_tb %>% filter(grepl('woE',rMech))%>% select(-r.enzyme) %>%mutate(enzyme = F, reactant = T, regulator = F) ,
           sum_tb %>% filter(!grepl('woE',rMech), is.na(reg)) %>% select(-r) %>% rename(r = r.enzyme) %>% mutate(enzyme = T, reactant = F, regulator = F),
           sum_tb %>% filter(!grepl('woE',rMech), is.na(reg)) %>% select(-r.enzyme) %>% mutate(enzyme = T, reactant = T, regulator = F),
           rMech.best %>% mutate(enzyme = T, reactant = T, regulator = T))%>%
  mutate(model = paste0(ifelse(enzyme, 'enzyme ',''),ifelse(reactant, 'reactant ',''),ifelse(regulator,'regulator',''))) %>%
  mutate(model = gsub(' $','',model))

r.report.table <- r %>% select(rxnid,rxnname,EC,model,r) %>% pivot_wider(names_from = 'model', values_from = 'r') %>%
  left_join(rxn.report)

write.xlsx(r.report.table,'report_pearsons_r.xlsx')

model.order <- c('enzyme','reactant','enzyme reactant','enzyme reactant regulator')
color.model <-  c('grey60',"#6BAED6","#3182BD", "#08519C")
names(color.model) <- model.order
r.plot <- r %>% mutate(model = factor(model, levels= model.order))
pdf('r_improvement_violin_plot_filtered.pdf', w = 5, h = 1.8)
ggplot(r.plot,
       aes(x = model, y = r, fill = model)) + 
  geom_violin(color = NA) + ylab("Pearson's R") +
  annotate('text', label = paste0('n = ',length(unique(r.plot$rxnid)),' reactions'), x = 2.5, y = 1.2) +
  scale_fill_manual(values = color.model) + scale_y_continuous(breaks = c(-1,0,1))+
  theme_classic() + theme(axis.text.x = element_blank(),axis.title.x = element_blank())
dev.off()

#### Metabolic Leverage ####

metLev <- list()
regAnal <- list()

for (rMech in rMech.best$rMech) {
  rxnSummary <- rxnForms[[rMech]]
  pos <- which(rxnSummary$fit$lik_track == max(rxnSummary$fit$lik_track))
  par.opt.K <- rxnSummary$fit$markov_par_vals[pos[1],,drop = F]
  colnames(par.opt.K) <- rxnSummary$kineticPars$formulaName
  par.opt.kcat <- rxnSummary$fit$nnls_par_vals[pos[1],,drop = F]
  colnames(par.opt.kcat) <- paste('E',rxnSummary$rxnID,colnames(par.opt.kcat), sep = '_')

  rxnEquations <- rxnSummary$rxnEquations
  omic_data <- rxnSummary$omic_data
  
  # to calculate variance contributed by enzyme and metabolite, var(logE) and var(logM), and regulation coefficient
  {
    flx_fit <- rxnSummary$fit$markov_fit[pos[1],]
    vals <- data.frame(omic_data$met_abund,
                       as.data.frame(omic_data$enzyme_abund),
                       2^par.opt.K, par.opt.kcat)
    occupancy <- as.numeric(sapply(rxnEquations[['l_occupancyExpression']],
                                   function(x) {eval(x,vals)}))
    enzyme <- as.numeric(flx_fit / occupancy)
    dt.reg <- data.frame(log2flx = as.numeric(log2(abs(omic_data$flux$standardQP))),
                         log2E = log2(enzyme),
                         log2M = log2(abs(occupancy)))%>%
      mutate(condition = colnames(flx_fit)) 
    
    # if calculating rho_e and rho_m for each nutrient condition
    # mutate(rho_e = (log2E-log2E[condition == 'B'])/(log2flx-log2flx[condition == 'B']),
    #        rho_m = (log2M-log2M[condition == 'B'])/(log2flx-log2flx[condition == 'B']))
    
    var.log2flx <- var(dt.reg$log2flx)
    var.log2E <- var(dt.reg$log2E)
    var.log2M <- var(dt.reg$log2M)
    {
      # linear regression
      dt.reg <- dt.reg %>% filter(!is.infinite(log2flx))
      rho_e <- lm(log2E~log2flx, dt.reg)[['coefficients']][['log2flx']]
      rho_m <- lm(log2M~log2flx, dt.reg)[['coefficients']][['log2flx']]
    }
  }
  
  vals <- data.frame(omic_data$met_abund %>% summarise_all(mean),
                     as.data.frame(omic_data$enzyme_abund) %>% summarise_all(mean),
                     2^par.opt.K, par.opt.kcat)
  sens <- sapply(rxnEquations[['kinetic_form_partials']],
                 function(x) {eval(x,vals)})
  
  stdev <- data.frame(omic_data$met_abund %>% summarise_all(sd),
                      as.data.frame(omic_data$enzyme_abund) %>% summarise_all(sd))
  
  metLev[[rMech]] <- data.frame(rMech = rMech,
                                metidx = colnames(stdev),
                                variance = t(stdev)^2) %>%
    mutate(metname = rxnSummary$metNames[metidx]) %>%
    mutate(variance = ifelse(is.nan(variance), 0, variance)) %>%
    mutate(sens = sens[metidx]) %>%
    mutate(sens = ifelse(is.na(sens), 0, sens)) %>%
    mutate(rxnStoi = rxnSummary$rxnStoi[metidx]) %>%
    mutate(type = case_when(is.na(rxnStoi) ~ 'enzyme',
                            (rxnStoi == 0) ~ 'regulator',
                            (rxnStoi > 0) ~ 'product',
                            (rxnStoi < 0) ~ 'substrate')) %>%
    mutate(lvr = sens^2*variance) %>%
    mutate(lvr = lvr/sum(lvr))
  regAnal[[rMech]] <- data.frame(var.log2flx, var.log2E, var.log2M, rho_e,rho_m)
}

metLev.df <- do.call(rbind, metLev) %>%
  left_join(rMech.best %>% select(rMech, rxnid, rxnname, EC)) %>%
  left_join(rxn.report)

regAnal.df <- do.call(rbind, regAnal) %>% tibble::rownames_to_column(var = 'rMech') %>%
  mutate(rMech = gsub('\\.[0-9]+$','',rMech)) %>%
  left_join(rMech.best %>% select(rMech, rxnid, rxnname, EC))%>%
  left_join(rxn.report)

# report metabolic leverage
write.xlsx(metLev.df, 'report_metabolic_leverage.xlsx')

# plot
metLev.df.plot <- metLev.df %>%
  mutate(type = factor(type, levels = c('enzyme','substrate','product','regulator'))) %>%
  arrange(type, desc(lvr)) %>%
  mutate(rxnid.GEM = factor(rxnid.GEM, levels = unique(rxnid.GEM)))

color.manual <- c('grey60',"#6BAED6","#6BAED6", "#3182BD")
names(color.manual) <- levels(metLev.df$type)
figure <- ggplot(metLev.df.plot, aes(x = rxnid.GEM, y = lvr, fill = type)) +
  geom_col(position = position_stack(), size= 0.5) +
  coord_cartesian(expand = F) +
  scale_fill_manual(values = color.manual) +
  theme_classic() + ylab('metabolic leverage') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(fill = NA, color = 'black',size = 1))
figure
pdf('regulation_analysis.pdf', w = 9, h = 2.5)
pdf('regulation_analysis_small.pdf', w = 4, h = 2.5)
print(figure)
dev.off()

pdf('regulation_coefficient.pdf', w = 2, h = 1.9)
ggplot(regAnal.df, aes(x = rho_e, y = rho_m)) +
  geom_abline(slope = 1, color = 'grey') +  geom_abline(intercept = 1, slope = -1, color = 'grey', linetype = 'dashed') + 
  # ggrepel::geom_text_repel(aes(label = rxnname))+
  # geom_hline(yintercept = 0, color = 'black') + geom_vline(xintercept= 0, color = 'black') +
  geom_point(pch = 21, color = 'skyblue3') + scale_y_continuous(limits = c(-1.2,2), expand = c(0,0))+ scale_x_continuous(limits = c(-1.2,2), breaks = c(-1,0,1,2), expand = c(0,0)) +
  theme_classic() + theme(panel.background = element_rect(fill = NA))
dev.off()

if (!file.exists('report_SIMMER.RData')) {
  save.image('report_SIMMER.RData')
}


