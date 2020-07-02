#' Create listing of AEs for Mayo DSMB reports
#' 
#' This function is design to replicate irbtoxrep.sas macro AE listing output for Mayo DSMB reports.
#' 
#' @param cytox The toxicity dataset containing multiple observations per patient with one observation per AE > grade 0
#' @param id The patient ID variable used in the dataset
#' @param byvar Variable used to separate summary by arm. If not present, "All Patients" printed for arm
#' @param grade.term Character. The variable that contains the AE grade
#' @param relation.term Character. The variable that containst the AE attribution
#' @param tox.term Character. The variable that contains the toxicity (MEDRA) code
#' @return This function creates a listing of the AE events for the Mayo DSMB report
#' @examples 
#' 
#' #With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' irb.listing <- irbtox_listing(cytox = cytox, id = "dcntr_id", byvar = "arm", tox.term = "toxicity")
#' irb.listing$table
#' 
#' @export
#' @author Sawyer Jacobson
#' @author Adam Pettinger
#' 
# Last Edited: 1/20/2020
#' 
#' 

irbtox_listing <-function(cytox, id = "dcntr_id", byvar="arm", grade.term = "grade", 
                          relation.term = "rel_smed", tox.term = "toxicity") {
  
  cytox <- as.data.frame(cytox)
  
  cytox <- cytox %>%
    rename(dcntr_id = id)
  
  cytox <- cytox %>%
    rename(grade = grade.term)
  
  cytox <- cytox %>%
    rename(rel_smed = relation.term)
  
  cytox <- cytox %>%
    rename(toxicity = tox.term)
  
  ######  Make sure the arm exists and to use the correct arm otherwise use All Patients#######
  if (any(colnames(cytox) == byvar)) { 
    cytox <- cytox %>%
      rename(arm = byvar) %>%
      mutate(arm = as.character(arm))
  } else {
    cytox$arm <- rep("All Patients", length(cytox$dcntr_id)) 
  }
  
  ####Will need to get the number eval by each group
  ###take out all the baseline toxicity if there are any
  if (any(colnames(cytox) == "cycle")){
    cytox = cytox[which(cytox$cycle>0),] 
  }
  
  unique.list = cytox[!duplicated(cytox[, "dcntr_id"]),]
  num.eval <- base::table(unique.list$arm) %>% 
    data.frame() %>% 
    mutate(arm = as.character(Var1)) %>% 
    dplyr::select(-Var1)  
  #######Now, lose all the grade 0 and cycle 0#####
  ###IF Cycle exists, which it won't in my theradex data
  if (any(colnames(cytox) == "cycle")){
    cytox = cytox %>% filter(cycle>0 & grade >0)
  }
  
  #######Now, lose all the grade 0#####
  cytox = cytox %>% filter(grade>0)
  
  
  ######Create CRTOX, the worst grade per patient, per toxicity
  crtox = cytox %>%
    group_by(dcntr_id, toxicity) %>% 
    slice(which.max(grade)) %>%
    select(dcntr_id, arm, toxicity, grade) 
  
  ##Count up each tox/arm/grades...
  table2 = crtox %>%
    group_by(toxicity, arm, grade) %>%
    summarise(count=n())
  
  ###seperate by grade to mutate their own column..  may be a better way, but this is good for me :)
  grade1 = table2 %>%
    filter(grade==1) %>%
    mutate(`1` = count) %>%
    select(toxicity, arm, `1`)
  grade2 = table2 %>%
    filter(grade==2) %>%
    mutate(`2` = count) %>%
    select(toxicity, arm, `2`)
  grade3 = table2 %>%
    filter(grade==3) %>%
    mutate(`3` = count) %>%
    select(toxicity, arm, `3`)
  grade4 = table2 %>%
    filter(grade==4) %>%
    mutate(`4` = count) %>%
    select(toxicity, arm, `4`)
  grade5 = table2 %>%
    filter(grade==5) %>%
    mutate(`5` = count) %>%
    select(toxicity, arm, `5`)
  
  ####Create a blank template to fill in####
  #####Here is where you define the cross table to merge with the grades, right now only 1 variable called arm
  ### One observatio for each toxicity*arm pairing
  template.table = data.frame(tidyr::crossing(toxicity=table2$toxicity, arm=table2$arm))
  
  ######Add toxterm and SOC##########
  formats = haven::read_sas("/people/ccs4/ccsicprd/cc-sas-all/sasdata/mart/xstudy/toxcodes.sas7bdat") %>%
    select("v5Meddra_Term", "v5Meddra_Code", "BODYSYS", "soc") %>%
    rename(v5meddra.term=v5Meddra_Term, v5meddra.code=v5Meddra_Code, bodysys = BODYSYS)
  
  formats = formats[!is.na(formats$v5meddra.code) & !duplicated(formats$v5meddra.code),]
  ###formats = formats[!duplicated(formats$v5meddra.code),]
  formats$toxicity=formats$v5meddra.code
  formats$toxterm=formats$v5meddra.term 
  
  ####merge and assign correct term if needed 
  template.table2= left_join(template.table, formats,  by="toxicity")
  
  table3 = crtox %>%
    group_by(toxicity) %>%
    summarise(total=n())
  
  try1 = left_join(template.table2, grade1, by=c("toxicity", "arm")) %>%
    left_join(., grade2, by=c("toxicity", "arm")) %>%
    left_join(., grade3, by=c("toxicity", "arm")) %>%
    left_join(., grade4, by=c("toxicity", "arm")) %>%
    left_join(., grade5, by=c("toxicity", "arm")) %>%
    left_join(num.eval, by = "arm") %>%
    left_join(table3, by = "toxicity")
  
  ##Add 0 to all toxicity arm grade combinations that don't exist.
  try1[is.na(try1)] = 0
  
  #######Add the percentages based on the num.eval for each arm
  
  try3 <- try1 %>%
    mutate(percent1 = as.character(ifelse(`1` == 0, NA, paste0(`1`, " (", round((`1`/Freq)*100, digits = 0), "%)"))),
           percent2 = as.character(ifelse(`2` == 0, NA, paste0(`2`, " (", round((`2`/Freq)*100, digits = 0), "%)"))),
           percent3 = as.character(ifelse(`3` == 0, NA, paste0(`3`, " (", round((`3`/Freq)*100, digits = 0), "%)"))),
           percent4 = as.character(ifelse(`4` == 0, NA, paste0(`4`, " (", round((`4`/Freq)*100, digits = 0), "%)"))),
           percent5 = as.character(ifelse(`5` == 0, NA, paste0(`5`, " (", round((`5`/Freq)*100, digits = 0), "%)")))) 
  
  
  final = try3 %>%
    arrange(soc, toxterm, arm) %>%
    mutate(soc = if_else(duplicated(soc), "", soc)) %>%
    arrange(desc(total), toxicity) %>%
    select(toxterm, arm, percent1, percent2, percent3, percent4, percent5) 
  
  all.listing.table <-  final %>% 
    flextable() %>%
    delete_part(part = "header") %>%
    add_header_row(values = c("Type", "Arm", "N (%)", "N (%)", "N (%)", "N (%)", "N (%)"), colwidths = c(1,1,1,1,1,1,1)) %>%
    add_header_row(values = c("Adverse Event", "1", "2", "3", "4", "5"), colwidths = c(2,1,1,1,1,1)) %>%
    add_header_row(values = c("Adverse Event", "Grade"), colwidths = c(2,5)) %>%
    theme_box() %>% 
    align(align = "center", part = "header") %>%
    merge_v(j = 1, part = "body") %>%
    merge_v(j = 1, part = "header") %>%
    bold(j = c(1,2), bold = TRUE, part = "body" ) %>%
    width(j = c(1,2,3,4,5,6,7), width = c(2.65, .5, .7, .7, .7, .7, .7)) %>%
    align(j = c(1,2), align = "left", part = "body") %>%
    align(j = c(3:7), align = "center", part = "body") %>%
    font(fontname = "Thorndale AMT", part = "all") %>%
    fontsize(j = c(3:7), size = 10, part = "body") %>%
    fontsize(size = 11, part = "header")
  
  out = list(table = all.listing.table, data = final)
  class(out) <- "irbtox_listing"
  return(out)
}


#' @export
print.irbtox_listing <- function(x, ...)
{
  print(x$table)
  invisible(x)
}
