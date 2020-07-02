#' Create summary of AEs for Mayo DSMB reports
#' 
#' This function is design to replicate irbtoxrep.sas macro AE summary output for Mayo DSMB reports.
#' 
#' @param cytox The toxicity dataset containing multiple observations per patient with one observation per AE > grade 0
#' @param id The patient ID variable used in the dataset
#' @param byvar Variable used to separate summary by arm. If not present, "All Patients" printed for arm
#' @param grade.term Character. The variable that contains the AE grade
#' @param relation.term Character. The variable that containst the AE attribution
#' @param tox.term Character. The variable that contains the toxicity (MEDRA) code
#' @param prot_num The study protocol number. Found in protref$dc_num
#' @param s_desc Short description of the study. Found in protref$s_desc.
#' @return This function returns a list of 1: Protocol Number and evaluable patients per Arm as a string, 2: Summary AE table, and 3: Total number of evaluable patients.
#' @examples 
#' 
#' #With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' irb.summary <- irbtox_summary(cytox = cytox, id = "dcntr_id", byvar = "arm", tox.term = "toxicity",
#'  prot_num = protref$dc_num, s_desc = protref$shrt_ttl_m)
#' irb.summary$table
#' 
#' @export
#' @author Sawyer Jacobson
#' @author Adam Pettinger
#' 
# Last Edited: 1/20/2020
#' 
#' 

irbtox_summary <-function(cytox, id = "dcntr_id", byvar = "arm", grade.term = "grade", 
                          relation.term = "rel_smed", tox.term = "toxicity", prot_num = '', s_desc = '') {
  
  prot_num = as.character(prot_num)
  s_desc = as.character(s_desc)
  
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
    cytox$arm = rep("All Patients", length(cytox$dcntr_id)) 
  }
  
  ####Will need to get the number eval by each group
  ###take out all the baseline toxicity if there are any
  if (any(colnames(cytox) == "cycle")){
    cytox = cytox[which(cytox$cycle>0),] 
  }
  
  unique.list = cytox[!duplicated(cytox[, "dcntr_id"]),]
  num.eval = base::table(unique.list$arm) %>% 
    data.frame() %>% 
    mutate(arm = as.character(Var1)) %>%
    dplyr::select(-Var1)
  
  #######Now, lose all the grade 0 and cycle 0#####
  ###IF Cycle exists, which it won't in my theradex data
  if (any(colnames(cytox) == "cycle")){
    cytox = cytox %>% 
      mutate(rel_smed = ifelse(is.na(rel_smed), 1, rel_smed)) %>%
      filter(cycle > 0 & grade >= 3)
  }
  
  if(any(colnames(cytox) == "bodysys")){
    cytox <- cytox %>% 
      dplyr::select(-c(bodysys))
  }
  if(any(colnames(cytox) == "soc")){
    cytox <- cytox %>% 
      dplyr::select(-c(soc))
  }
  
  ######Add Bodysys from toxcodes based on the meddra toxicity code if not there already
  if (!any(colnames(cytox) == "bodysys")){
    formats = haven::read_sas("/people/ccs4/ccsicprd/cc-sas-all/sasdata/mart/xstudy/toxcodes.sas7bdat") %>%
      dplyr::select("v5Meddra_Term", "v5Meddra_Code", "BODYSYS", "soc") %>%
      rename(v5meddra.term=v5Meddra_Term, v5meddra.code=v5Meddra_Code, bodysys = BODYSYS)
    
    formats = formats[!is.na(formats$v5meddra.code) & !duplicated(formats$v5meddra.code),]
    ###formats = formats[!duplicated(formats$v5meddra.code),]
    formats$toxicity = formats$v5meddra.code
    formats$toxterm = formats$v5meddra.term 
    ####merge and assign correct term if needed 
    cytox <- left_join(cytox, formats,  by = "toxicity")
  }
  
  ########Find each of the categories ########
  worst.grade <- cytox %>%
    group_by(dcntr_id) %>% 
    slice(which.max(grade)) %>%
    mutate(category = "Total") %>%
    dplyr::select(dcntr_id, category, grade, arm)
  
  worst.heme <- cytox %>%
    filter (bodysys == "1") %>% 
    group_by(dcntr_id) %>% 
    slice(which.max(grade)) %>%
    mutate(category = "Heme") %>%
    dplyr::select(dcntr_id, category, grade, arm)
  
  worst.nonheme <- cytox %>% 
    filter (bodysys != "1")%>% 
    group_by(dcntr_id) %>% 
    slice(which.max(grade)) %>%
    mutate(category = "Non-Heme") %>%
    dplyr::select(dcntr_id, category, grade, arm)
  
  table1.raw <- full_join(worst.grade, worst.heme, by = c("dcntr_id", "category", "grade", "arm")) %>%
    full_join(worst.nonheme, by = c("dcntr_id", "category", "grade", "arm"))
  
  template <- data.frame(tidyr::crossing(category = c("Total", "Heme", "Non-Heme"), grade = 3:5, arm = as.character(unique(num.eval$arm))))
  
  
  grade3plus_total <- table1.raw %>%
    group_by(category, arm) %>%
    filter(grade >= 3) %>%
    dplyr::count(name = 'count') %>% 
    right_join(template %>% filter(grade == 3), by = c("category", "arm")) %>%
    bind_cols(data.frame(event = c(rep("Grade 3+ Heme Adverse Event", length(unique(num.eval$arm))), 
                                   rep("Grade 3+ Non-Hem Adverse Event", length(unique(num.eval$arm))),
                                   rep("Grade 3+ Adverse Event", length(unique(num.eval$arm)))), 
                         stringsAsFactors = FALSE)) %>%
    ungroup() %>%
    arrange(event) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    select(category, event, arm, count)
  grade4plus_total <- table1.raw %>%
    group_by(category, arm) %>%
    filter(grade >= 4) %>%
    dplyr::count(name = 'count') %>% 
    right_join(template %>% filter(grade == 4), by = c("category", "arm")) %>%
    bind_cols(data.frame(event = c(rep("Grade 4+ Heme Adverse Event", length(unique(num.eval$arm))), 
                                   rep("Grade 4+ Non-Hem Adverse Event", length(unique(num.eval$arm))),
                                   rep("Grade 4+ Adverse Event", length(unique(num.eval$arm)))), 
                         stringsAsFactors = FALSE)) %>%
    ungroup() %>%
    arrange(event) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    select(category, event, arm, count)
  grade5plus_total <- table1.raw %>%
    group_by(category, arm) %>%
    filter(grade >= 5) %>%
    dplyr::count(name = 'count') %>% 
    right_join(template %>% filter(grade == 5), by = c("category", "arm")) %>%
    bind_cols(data.frame(event = c(rep("Grade 5 Heme Adverse Event", length(unique(num.eval$arm))), 
                                   rep("Grade 5 Non-Hem Adverse Event", length(unique(num.eval$arm))),
                                   rep("Grade 5 Adverse Event", length(unique(num.eval$arm)))),
                         stringsAsFactors = FALSE)) %>%
    ungroup() %>%
    arrange(event) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    select(category, event, arm, count)
  
  newtable <- bind_rows(grade3plus_total, grade4plus_total, grade5plus_total) %>%
    left_join(num.eval, by = 'arm') %>%
    mutate(percent = round((count/Freq)*100, digits = 1),
           count_char = as.character(count),
           percent_char = ifelse(count == 0, "0.0", paste0('', percent, '')),
           ord_by = ifelse(event == "Grade 3+ Adverse Event", 1 ,
                           ifelse(event == "Grade 4+ Adverse Event", 2, 
                                  ifelse(event == "Grade 5 Adverse Event", 3, 
                                         ifelse(event == "Grade 3+ Heme Adverse Event", 4, 
                                                ifelse(event == "Grade 4+ Heme Adverse Event", 5, 
                                                       ifelse(event == "Grade 5 Heme Adverse Event", 6, 
                                                              ifelse(event == "Grade 3+ Non-Hem Adverse Event", 7, 
                                                                     ifelse(event == "Grade 4+ Non-Hem Adverse Event", 8, 
                                                                            ifelse(event == "Grade 5 Non-Hem Adverse Event", 9, NA)))))))))) %>%
    arrange(ord_by, arm) %>%
    select(-c(category, ord_by, percent, Freq, count))

  
  
  
  # evaluable <- cytox %>% 
  #   filter(arm != '') %>% 
  #   distinct(arm) %>% 
  #   nrow()
  evaluable <- nrow(num.eval)
  
  eval <- character(evaluable)
  for(i in 1:evaluable){
    eval[i] <- paste0("Arm ", num.eval[i, 2], " Evaluable Patients: ", num.eval[i, 1])
  }
  arms <- paste(eval, sep = '', collapse = " \n ")
  if(nchar(prot_num) > 0 & nchar(s_desc) > 0){
    header <- paste0(prot_num, " - ", s_desc, " \n ", "Adverse Events (Regardless of Attribution) Reported as of ", 
                     format(Sys.time(), format = "%B %d, %Y"), " \n ", arms)
    }
  else{
    header <- paste0("Adverse Events (Regardless of Attribution) Reported as of ", 
                     format(Sys.time(), format = "%B %d, %Y"), " \n ", arms)
  }
  
  new_flex_table <- newtable
  names(new_flex_table) <- c("Patients with at least one:", "Arm", "N", "%")
  all.summary.table <-  new_flex_table %>% 
    flextable() %>%
    delete_part(part = "header") %>%
    add_header_row(values = c("Patients with at least one:", "Arm", '', ''), colwidths = c(1,1,1,1)) %>%
    add_header_row(values = c('', 'N', '%'), colwidths = c(2,1,1)) %>%
    merge_v(j = 1, part = "body") %>%
    bold(j = c(1,2), bold = TRUE, part = "body" ) %>%
    theme_box() %>%
    width(j = c(1,2,3,4), width = c(2.5, .5, .35, .45)) %>%
    align(align = "left", part = "header") %>%
    align(j = c(1,2), align = "left", part = "body") %>%
    align(j = c(4), align = "center", part = "header") %>%
    font(fontname = "Thorndale AMT", part = "all") %>%
    fontsize(j = c(3,4), size = 10, part = "body") %>%
    fontsize(size = 11, part = "header") 
  out <- list(header = header, table = all.summary.table, evaluable = sum(num.eval[, 1]), data = newtable)
  class(out) <- "irbtox_summary"
  return(out)
}


#' @export
print.irbtox_summary <- function(x, ...)
{
  print(x$table)
  invisible(x)
}

