check.balance <- function(met = 'atp\\_c', mfa, ci =T){
  if (ci) {
    flux <- mfa%>% select(id,reaction,flx,lb,ub) %>% distinct()  %>%
      filter(grepl(paste(paste0(c('^',' '),met), collapse = '|'),reaction)) %>% filter(!flx == 0) %>%
      mutate(dir = ifelse(grepl(paste(paste0(c('^',' '),met), collapse = '|'),unlist(lapply(strsplit(reaction, split = '>'),'[',1))), -1,1)) %>%
      mutate(stoi = str_extract(reaction,paste(paste0(c('^','[0-9. ]+'),met), collapse = '|') )) %>%
      mutate(stoi = as.numeric(gsub(met,'',stoi)), stoi = ifelse(is.na(stoi),1,stoi)*dir) %>%
      mutate(flx = stoi*flx, lb.tmp = stoi* ifelse(stoi > 0,lb, ub), ub.tmp = stoi*ifelse(stoi>0,ub, lb)) %>%
      select(-lb,-ub) %>% rename(lb = lb.tmp, ub = ub.tmp)
  } else {
    flux <- mfa %>% select( id,reaction,flx) %>% distinct()  %>%
      filter(grepl(paste(paste0(c('^',' '),met), collapse = '|'),reaction)) %>% filter(!flx == 0) %>%
      mutate(dir = ifelse(grepl(paste(paste0(c('^',' '),met), collapse = '|'),unlist(lapply(strsplit(reaction, split = '>'),'[',1))), -1,1)) %>%
      mutate(stoi = str_extract(reaction,paste(paste0(c('^','[0-9. ]+'),met), collapse = '|') )) %>%
      mutate(stoi = as.numeric(gsub(met,'',stoi)), stoi = ifelse(is.na(stoi),1,stoi)*dir) %>%
      mutate(flx = stoi*flx)
  }
  return(flux)
}




  

