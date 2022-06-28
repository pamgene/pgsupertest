library(tidyverse)
library(data.table)
library(globaltest)

#' @import tidyverse data.table globaltest
#' @export
upstreamDb = function(){
  load(system.file("extdata", "220512-86402-86412-87102_UpstreamDb.Rdata" , package = "pgSupertest"))
  UpstreamDatabase %>%
    ungroup()
}

#' @export
TableS1 = function(){
  df = read.delim(system.file("extdata", "TableS1 Kinases in Kinome Render.txt", package = "pgSupertest"))
}

#' @export
flagOutOfGroup = function(db){
    db %>%
    left_join(TableS1(), by = c(Kinase_Name = "KR_Name")) %>%
    mutate(outOfGroup = case_when(PepProtein_Residue == "Y"  & GroupName != "TK" ~ TRUE,
                                  PepProtein_Residue != "Y"  & GroupName == "TK" ~ TRUE,
                                  TRUE ~ FALSE))
}

#' @export
scoreIviv = function(db, score, outOfGroup_score = NULL){
  db %>%
    flagOutOfGroup() %>%
    filter(Database == "iviv") %>%
    mutate(s = score) %>%
    mutateOutOfGroup(outOfGroup_score)
}

#' @export
scorePNet = function(db, ranks, scores, outOfGroup_score = NULL){
  db %>%
    flagOutOfGroup() %>%
    filter(Database == "PhosphoNET") %>%
    mutate(
      s  = case_when(
        Kinase_Rank > ranks[2] ~ scores[3],
        Kinase_Rank > ranks[1] ~ scores[2],
        TRUE ~ scores[1])
    ) %>%
    mutateOutOfGroup(outOfGroup_score)
}

#' @export
scoreNoMatch = function(db, score, outOfGroup_score = NULL){
  db = db %>%
    flagOutOfGroup() %>%
    mutate(s = case_when(is.nan(s) ~ score,
                         TRUE ~ s)) %>%
    mutateOutOfGroup(outOfGroup_score)
}

#' @import reshape2
#' @export
combinedScores = function(db, minscore, outOfGroupScore = 0){
  dbc = db %>%
    mutate(Kinase_Name = Kinase_Name %>% as.factor,
           ID = ID %>% as.factor) %>%
    group_by(Kinase_Name, ID) %>%
    summarise(s = 1-prod(1-s), outOfGroup = all(outOfGroup))

  Koog  = dbc %>%
    reshape2::dcast(ID ~ Kinase_Name, value.var = "outOfGroup") %>%
    select(-ID) %>%
    apply(2, any, na.rm = TRUE) %>%
    as.data.frame() %>%
    mutate(Kinase_Name = rownames(.))

  Koog = Koog %>%
    mutate(Kinase_Name = rownames(Koog)) %>%
    select(Kinase_Name, oog = ".")

  db.all = expand.grid(levels(droplevels(dbc$Kinase_Name)), levels(droplevels(dbc$ID)) )
  colnames(db.all) = c("Kinase_Name", "ID")

  db.all %>%
    left_join(dbc, by = c("Kinase_Name", "ID")) %>%
    left_join(Koog, by = "Kinase_Name") %>%
    mutate(s = case_when(
      s < minscore ~ minscore,
      s >= minscore  ~ s,
      is.na(s) & !oog ~ minscore,
      is.na(s) & oog ~ outOfGroupScore)) %>%
    select(Kinase_Name, ID, s, oog)

}

#' @export
normalizeScores = function(db){
  sumdf = db %>% group_by(ID) %>% summarise(sumsc = sum(s))
  db %>%
    left_join(sumdf, by = "ID") %>%
    group_by(ID) %>%
    mutate(sc.nor = s/sumsc)
}

mutateOutOfGroup = function(db, oogs = NULL){
  if(!is.null(oogs)){
    db = db %>%
      mutate(s = case_when(outOfGroup ~ oogs,
                           TRUE ~ s))
  }
  db
}

#' @import data.table
#' @export
gtest = function(kinase.set, X, grp){
  bset = colnames(X) %in% kinase.set$ID
  X = X[, bset,drop = FALSE]
  if(dim(X)[2]>0){
    globtest = gt(grp ~ X, directional = TRUE, standardize = TRUE)
    p = attr(globtest, "result")[1,1]
    delta = mean(gdelta(X, grp))
  } else {
    p = NaN
    delta = NaN
  }
  result = data.table(p = p, N = dim(X)[2], delta = delta, gt = list(globtest))
}

#' @export
gdelta = function(Xin, grp){
  X1 = Xin[grp == levels(grp)[1],,drop = FALSE]
  X2 = Xin[grp == levels(grp)[2],,drop= FALSE]
  M1 =apply(X1, 2, mean)
  M2 =apply(X2, 2, mean)
  delta = M2 - M1
}
