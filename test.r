rm(list = ls())
library(pgSupertest)
library(tidyverse)

al = read.delim("86402 Array Layout.txt")

getIDTable <- function(df){

    df %>%
    dplyr::select(ID) %>%
    mutate(.ri = 1:n(), .ci = 1)
}


db = function(){
  pgSupertest::upstreamDb()
}
selectMappings = function(df){
  tabs1 = TableS1()
  db() %>%
    ungroup() %>%
    filter(PepProtein_SeqSimilarity >= seqhom) %>%
    filter(PepProtein_Residue == "Y") %>%
    filter(Kinase_Name %in% tabs1$KR_Name) %>%
    filter(ID %in% df$ID)
}

scores = function(mappings){
  dbi = mappings %>%
    scoreIviv(piviv, pOutOfGroup)
  dbp = mappings %>%
    scorePNet(ranks = c(pnet_bp1, pnet_bp2), scores = c(p_pnet0, p_pnet1, p_pnet2), pOutOfGroup)



  result =  dbi %>%
    bind_rows(dbp) %>%
    distinct(ID, Kinase_Name, Database, .keep_all = TRUE) # remove multiple lines from same db type
    result = combinedScores(result, p0)
    result = normalizeScores(result)
    return(result)
}

seqhom = 0.9
piviv = 0.9
pnet_bp1 = 4
pnet_bp2 = 12
p_pnet0 = 0.9
p_pnet1 = 0.7
p_pnet2 = 0.3
p0 = 0.25
pOutOfGroup = 0.5

csdata = al %>%
  getIDTable() %>%
  selectMappings() %>%
  scores() %>%
  arrange(ID, Kinase_Name)


kinase.s1 = TableS1()

kinase.table = csdata %>%
  ungroup() %>%
  distinct(Kinase_Name) %>%
  select(Kinase_Name) %>%
  mutate(.ci = (1:n()) -1) %>%
  left_join(kinase.s1, by = c(Kinase_Name = "KR_Name")) %>%
  select(.ci, Kinase_Name, Kinase_UniprotID = KR_UniprotID, Kinase_GroupName = GroupName, Kinase_Family = Family) %>%
  arrange(Kinase_Family)



id.table = al %>%
  getIDTable() %>%
  select(.ri, ID) %>%
  distinct(ID, .keep_all = TRUE)

result = csdata %>%
  left_join(kinase.table, by = "Kinase_Name") %>%
  left_join(id.table, by = "ID")

result = result %>%
  select(.ri, .ci, ID, s, sc.nor)

db = upstreamDb()
