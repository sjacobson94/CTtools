#' Create listing of AEs for Alliance DSMB reports
#' 
#' This function is design to replicate aerpt_alliance_macro.sas AE listing output for Alliance DSMB reports.
#' 
#' @param cytox The toxicity dataset containing multiple observations per patient with one observation per AE > grade 0
#' @param id The patient ID variable used in the dataset
#' @param grade.cutoff Grade cutoff for minimum AE grade. Inclusive
#' @param rel.cutoff Attribution cutoff for minimum AE attribution. Inclusive
#' @param rel.include Value of 0 or 1. If 0, "Regardless of Attribution" printed in header. If 1, "At least possibly related" printed in header
#' @param byvar Variable used to separate summary by arm. If not present, "All Patients" printed for arm. Set to NULL if calculations are desired for all patients together.
#' @param arm.include Include "Arm" in table (1) or not (0). Default is 1.
#' @param arm.labels Option to include additional Arm labels in the table header. Add as vector of characters.
#' @param toxvar Variable containing numeric toxicity codes.
#' @param relvar Variable containing numeric relation values 1-5.
#' @param gradevar Grade variable, numeric.
#' @param bodysysvar Body system variable, character.
#' @return This function creates a listing of the AE events for the Alliance DSMB report
#' @examples 
#' 
#' # With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' all.listing <- aerpt.listing(cytox, id="dcntr_id", byvar="arm", 
#' grade.cutoff = 3, rel.cutoff = 0, rel.include = 0)
#' all.listing$table
#' 
#' @export
#' @author Sawyer Jacobson
#' @author Adam Pettinger
#' 
# Last Edited: 2/29/2020

aerpt.listing <- function(cytox, id, grade.cutoff = 1, rel.cutoff = 0, rel.include = 1, byvar="arm", arm.include = 1, arm.labels = NULL, 
                          toxvar = "toxicity", relvar = "rel_smed", gradevar = "grade", bodysysvar = "bodysys") 
  {
  
  cytox <- as.data.frame(cytox)
  
  arm.labels <- arm.labels
  
  if(is.null(byvar)) arm.include = 0
  
  # Setting up attribution variable
  if(rel.include == 1){
    rel = "At Least Possibly Related\n"
  } else {
    rel = "Regardless of Attribution\n"
  }
  cytox <- cytox %>%
    rename(dcntr_id = tidyselect::all_of(id), 
           toxicity = tidyselect::all_of(toxvar), 
           grade = tidyselect::all_of(gradevar),
           rel_smed = tidyselect::all_of(relvar))
  
  ######  Make sure the arm exists and to use the correct arm otherwise use All Patients#######
  if (any(colnames(cytox) == byvar)) { 
    cytox <- cytox %>%
      rename(arm = tidyselect::all_of(byvar)) %>%
      mutate(arm = as.character(arm))
  } else {
    cytox$arm <- rep("A", length(cytox$dcntr_id)) 
  }
  
  
  ####Will need to get the number eval by each group
  ###take out all the baseline toxicity if there are any
  if (any(colnames(cytox) == "cycle")){
    cytox <- cytox[which(cytox$cycle>0),] 
  }
  
  unique.list <- cytox[!duplicated(cytox[, "dcntr_id"]),]
  
  num.eval <- base::table(unique.list$arm) %>% 
    data.frame() %>% 
    mutate(arm = as.character(Var1)) %>% 
    dplyr::select(-Var1)
  
  if(arm.include == 0) TRUE
  else if(!is.null(arm.labels) & length(arm.labels) != nrow(num.eval)) stop('Number of arm levels not equal to number of arm labels')
  
  #######Now, lose all the grade 0 and cycle 0#####
  ###IF Cycle exists, which it won't in my theradex data
  if (any(colnames(cytox) == "cycle")){
    cytox <- cytox %>% 
      mutate(rel_smed = ifelse(is.na(rel_smed), 1, rel_smed)) %>%
      filter(cycle > 0 & grade >= grade.cutoff & rel_smed >= rel.cutoff)
  }
  
  ######Add toxterm and SOC##########
  data("tox_codes")
  
  #formats <- haven::read_sas("/people/ccs4/ccsicprd/cc-sas-all/sasdata/mart/xstudy/toxcodes.sas7bdat") %>%
  tox_codes <- tox_codes %>%
    select("v5Meddra_Term", "v5Meddra_Code", "BODYSYS", "soc") %>%
    rename(v5meddra.term = v5Meddra_Term, v5meddra.code = v5Meddra_Code)#, bodysys = BODYSYS)
  
  tox_codes <- tox_codes[!is.na(tox_codes$v5meddra.code) & !duplicated(tox_codes$v5meddra.code),]
  tox_codes$toxicity <- tox_codes$v5meddra.code
  tox_codes$toxterm <- tox_codes$v5meddra.term 
  
  if(any(names(cytox) == bodysysvar))
  {
    cytox <- cytox %>%
      rename(bodysys = all_of(bodysysvar))
  } else
  {
    cytox <- cytox %>%
      left_join(tox_codes %>% select(toxicity, bodysys = BODYSYS) %>% distinct(toxicity, .keep_all = TRUE), by = "toxicity")
  } 
  
  cytox1 <- cytox %>% 
    #rename(soc_cytox = soc) %>%
    dplyr::select(-c(bodysys)) 
  
  ######Create CRTOX, the worst grade per patient, per toxicity
  crtox <- cytox1 %>%
    group_by(dcntr_id, toxicity) %>% 
    slice(which.max(grade)) %>%
    dplyr::select(dcntr_id, arm, toxicity, grade) 
  
  bodysys <- cytox %>% 
    dplyr::select(bodysys, toxicity) %>% 
    distinct(toxicity, bodysys)
  ##Count up each tox/arm/grades...
  table2 <- crtox %>%
    group_by(toxicity, arm, grade) %>%
    summarise(count = n()) %>% 
    left_join(bodysys, by = "toxicity")
  
  
  ###seperate by grade to mutate their own column..  may be a better way, but this is good for me :)
  grade1 <- table2 %>%
    filter(grade == 1) %>%
    mutate(`1` = count) %>%
    dplyr::select(toxicity, arm, `1`)
  grade2 <- table2 %>%
    filter(grade == 2) %>%
    mutate(`2` = count) %>%
    dplyr::select(toxicity, arm, `2`)
  grade3 <- table2 %>%
    filter(grade == 3) %>%
    mutate(`3` = count) %>%
    dplyr::select(toxicity, arm, `3`)
  grade4 <- table2 %>%
    filter(grade == 4) %>%
    mutate(`4` = count) %>%
    dplyr::select(toxicity, arm, `4`)
  grade5 <- table2 %>%
    filter(grade == 5) %>%
    mutate(`5` = count) %>%
    dplyr::select(toxicity, arm, `5`)
  
  
  ####Create a blank template to fill in####
  #####Here is where you define the cross table to merge with the grades, right now only 1 variable called arm
  ### One observation for each toxicity*arm pairing
  template.table <- data.frame(tidyr::crossing(toxicity = table2$toxicity, arm = table2$arm))
  
  
  
  ####merge and assign correct term if needed 
  template.table2 <- left_join(template.table, tox_codes %>% distinct(toxicity, .keep_all = TRUE),  by = "toxicity")
  
  try1 <- template.table2 %>%
    left_join(grade1, by = c("toxicity", "arm")) %>%
    left_join(grade2, by = c("toxicity", "arm")) %>%
    left_join(grade3, by = c("toxicity", "arm")) %>%
    left_join(grade4, by = c("toxicity", "arm")) %>%
    left_join(grade5, by = c("toxicity", "arm")) %>%
    select(-BODYSYS) %>%
    left_join(bodysys, by = "toxicity")
  ##Add 0 to all toxicity arm grade combinations that don't exist.
  try1[is.na(try1)] <- 0
  
  # all_formats <- haven::read_sas("/projects/ccsstat/across/R/formats.sas7bdat") %>%
  #   sdms::cn_tolower() %>% 
  #   filter(fmtname == "BODYSYS") %>%
  #   dplyr::select(start, label)
  # #all_formats$fmtname <- tolower(all_formats$fmtname)
  # all_formats1 <- all_formats %>% 
  #   bind_rows(c(start = "**", label = "**"))
  #bodysys_formats <- read.csv("/projects/ccsstat/across/R/bodysys_formats.csv", stringsAsFactors = FALSE)
  data("bodysys_formats")
  
  try2 <- try1 %>%
    left_join(bodysys_formats, by = c("bodysys" = "start"))
  
  try3 <- try2 %>%
    mutate(category = ifelse(toxicity %in% c(10002272,10007839,10025256,10029366,10035528,10048580,10049182,
                                                    10005329,10019150,10019491,10024378,10025258,10028533,10041633,
                                                    10055599,10065973), "Hematologic Adverse Events", "Non-Hematologic Adverse Events"),
                  eventcat = ifelse(toxicity %in% c(10002272,10007839,10025256,10029366,10035528,10048580,10049182,
                                                    10005329,10019150,10019491,10024378,10025258,10028533,10041633,
                                                    10055599,10065973), "Blood/Bone Marrow", 
                                    ifelse(toxicity == 10016288, "Blood and lymphatic sys disord", soc)),
                  label = ifelse(eventcat == "Blood/Bone Marrow", "Blood/Bone Marrow", label)) %>%
    left_join(num.eval, by = "arm")
  
  grade_list <- c("1-Mild\n  n (%)", "2-Mod\n  n (%)", "3-Severe\n  n (%)", "4-LifeThr\n  n (%)", "5-Lethal\n  n (%)")
  
  #######Add the percentages based on the num.eval for each arm
  try3 <- try3 %>%
    mutate(percent1 = paste0(`1`, " (", round((`1`/Freq)*100, digits = 0), "%)"),
                  percent2 = paste0(`2`, " (", round((`2`/Freq)*100, digits = 0), "%)"),
                  percent3 = paste0(`3`, " (", round((`3`/Freq)*100, digits = 0), "%)"),
                  percent4 = paste0(`4`, " (", round((`4`/Freq)*100, digits = 0), "%)"),
                  percent5 = paste0(`5`, " (", round((`5`/Freq)*100, digits = 0), "%)"),
                  category = factor(category, levels = c("Total", "Hematologic Adverse Events", "Non-Hematologic Adverse Events"), 
                                    labels = c("Total", "Hematologic Adverse Events", "Non-Hematologic Adverse Events"))) %>% 
    arrange(category, label, toxterm, arm) %>%
    dplyr::select(category, label, toxterm, arm, paste0("percent", 1:5)) %>%
    rename('Arm\n ' = arm, 'AE Term' = toxterm)
  
  grades <- try3 %>% 
    dplyr::select(paste0("percent", 1:5))
  try4 <- try3 %>%
    dplyr::select(-c(paste0("percent", 1:5)))
  names(grades) <- grade_list
  grade <- grades[as.integer(stringr::str_sub(names(grades),1,1)) >= grade.cutoff]
  final <- try4 %>%
    bind_cols(grade)
  
  evaluable <- cytox %>% 
    filter(arm != '') %>% 
    distinct(arm) %>% nrow()
  
  eval <- character(evaluable)
  # for(i in 1:evaluable){
  #   eval[i] <- paste0("Arm ", num.eval[i, 2], "=", num.eval[i, 1])
  # }
  # arms <- paste(eval, sep = '', collapse = "  ")
  # 
  # summary_title <- paste0("Listing of Grade ", grade.cutoff, "+ Adverse Events\nMax Grade Per Patient Per Event\n",
  #                         rel, 
  #                         "Number of Evaluable Patients:\n",
  #                         arms)
  # 
  for(i in 1:evaluable){
    if(!is.null(arm.labels) & arm.include == 1){
      eval[i] <- paste0(arm.labels[i]," (Arm ", num.eval[i, 2], "=", num.eval[i, 1], ")")
    }else if (arm.include == 1 ){
      eval[i] <- paste0("Arm ", num.eval[i, 2], "=", num.eval[i, 1])
    } else{
      eval[i] <- as.integer(num.eval[i, 1])
    }
  }
  if(!is.null(arm.labels) & arm.include == 1){
    arms <- paste(eval, sep = '', collapse = "\n")
    
    summary_title <- paste0("Listing of Grade ", grade.cutoff, "+ Adverse Events\nMax Grade Per Patient Per Event\n",
                            rel,
                            "Number of Evaluable Patients:\n",
                            arms)
  }else if(is.null(arm.labels) & arm.include == 1){
    arms <- paste(eval, sep = '', collapse = "\t")
    
    summary_title <- paste0("Listing of Grade ", grade.cutoff, "+ Adverse Events\nMax Grade Per Patient Per Event\n",
                            rel,
                            "Number of Evaluable Patients:\n",
                            arms)
  }else if(arm.include == 0){
    arms <- sum(as.numeric(eval))
    summary_title <- paste0("Listing of Grade ", grade.cutoff, "+ Adverse Events\nMax Grade Per Patient Per Event\n",
                            rel,
                            "Number of Evaluable Patients: ",
                            arms)
  }
  
  if(arm.include == 1){
    col_names <- names(final)[-c(1,2,3)]
    newtable <- final %>% 
      as_grouped_data(groups = c('category', 'label')) %>%
      as_flextable(hide_grouplabel = TRUE) %>%
      delete_part(part = "header") %>%
      add_header_row(values = c("", col_names)) %>%
      bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
      bold(j = 1, i = ~ !is.na(label), bold = TRUE, part = "body" ) %>%
      fontsize(size = 10, part = "all") %>% 
      font(fontname = "Arial", part = "all") %>% 
      add_header_row(values = c('', "Grade of Adverse Events"), colwidths = c(2, (6 - grade.cutoff))) %>%
      add_header_lines(values = summary_title) %>%
      bold(bold = TRUE, part = "header") %>% 
      merge_v(j = 3, part = "body") %>%
      border_remove() %>% # theme_box() %>% # if boxes are wanted for better readability
      width(j = c(1:(grade.cutoff + 2)), width = c(2, .5, rep(1, grade.cutoff))) %>%
      align(align = "center", part = "header") %>%
      align(j = 1, align = "left", part = 'body') %>%
      # align(j = 1, i = ~ !is.na(category), align = "left", part = "body") %>%
      # align(j = 1, i = ~ !is.na(label), align = "left", part = "body") %>%
      align(j = c(3:(grade.cutoff + 2)), align = "right", part = "body") %>% 
      align(i = 2, j = 1, align = "left", part = "header") %>%
      valign(j = 1:3, part = "body", valign = "top") %>%
      padding(i = ~ (is.na(category) & is.na(label)), padding.left = 22) %>%
      border_outer(part = "all", border = officer::fp_border(color = "black", style = "solid", width = 1)) %>%
      fix_border_issues()
  }else{
    col_names <- names(final)[-c(1,2,3,4)]
    newtable <- final %>% 
      select(-`Arm\n `) %>%
      as_grouped_data(groups = c('category', 'label')) %>%
      as_flextable(hide_grouplabel = TRUE) %>%
      delete_part(part = "header") %>%
      add_header_row(values = c("", col_names)) %>%
      bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
      bold(j = 1, i = ~ !is.na(label), bold = TRUE, part = "body" ) %>%
      fontsize(size = 10, part = "all") %>% 
      font(fontname = "Arial", part = "all") %>% 
      add_header_row(values = c('', "Grade of Adverse Events"), colwidths = c(1, (6 - grade.cutoff))) %>%
      add_header_lines(values = summary_title) %>%
      bold(bold = TRUE, part = "header") %>% 
      merge_v(j = 3, part = "body") %>%
      border_remove() %>% # theme_box() %>% # if boxes are wanted for better readability
      width(j = c(1:(grade.cutoff + 1)), width = c(2, rep(1, grade.cutoff))) %>%
      align(align = "center", part = "header") %>%
      align(j = 1, align = "left", part = 'body') %>%
      # align(j = 1, i = ~ !is.na(category), align = "left", part = "body") %>%
      # align(j = 1, i = ~ !is.na(label), align = "left", part = "body") %>%
      align(j = c(3:(grade.cutoff + 1)), align = "right", part = "body") %>% 
      align(i = 2, j = 1, align = "left", part = "header") %>%
      valign(j = 1:3, part = "body", valign = "top") %>%
      padding(i = ~ (is.na(category) & is.na(label)), padding.left = 22) %>%
      border_outer(part = "all", border = officer::fp_border(color = "black", style = "solid", width = 1)) %>%
      fix_border_issues()
  }
  out <- list(table = newtable, data = final)
  class(out) <- "aerpt.listing"  
  out
  
}

#' @export
print.aerpt.listing <- function(x, ...)
{
  print(x$table)
  invisible(x)
}

