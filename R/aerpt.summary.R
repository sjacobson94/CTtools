#' Create summary of AEs for Alliance DSMB reports
#' 
#' This function is design to replicate aerpt_alliance_macro.sas AE summary output for Alliance DSMB reports.
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
#' @return This function creates a summary of the AE events for the Alliance DSMB report
#' @examples 
#' 
#' #With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' all.summary <- aerpt.summary(cytox, id="dcntr_id", byvar="arm", 
#' grade.cutoff = 3, rel.cutoff = 0, rel.include = 0)
#' all.summary$table
#' 
#' @export
#' @author Sawyer Jacobson
#' @author Adam Pettinger
#' 
# Last Edited: 2/29/2020

aerpt.summary <- function(cytox, id, grade.cutoff = 1, rel.cutoff = 0, rel.include = 1, byvar = "arm", arm.include = 1, arm.labels = NULL,
                          toxvar = "toxicity", relvar = "rel_smed", gradevar = "grade", bodysysvar = "bodysys") 
  {
  cytox <- as.data.frame(cytox)
  
  arm.labels <- arm.labels
  
  if(is.null(byvar)) arm.include = 0
  
  ######  Make sure to use the correct ID#######
  if(rel.include == 1){
    rel = "At least possibly related\n"
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
    cytox$arm = rep("A", length(cytox$dcntr_id)) 
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
  
  if(arm.include == 0) TRUE
  else if(!is.null(arm.labels) & length(arm.labels) != nrow(num.eval)) stop('Number of arm levels not equal to number of arm labels')
  
  #######Now, lose all the grade 0 and cycle 0#####
  ###IF Cycle exists, which it won't in my theradex data
  if (any(colnames(cytox) == "cycle")){
    cytox = cytox %>% 
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
  
  table1.raw <- full_join(worst.heme, worst.nonheme, by = c("dcntr_id", "category", "grade", "arm")) %>%
    full_join(worst.grade, by = c("dcntr_id", "category", "grade", "arm"))
  
  template <- data.frame(tidyr::crossing(category = c("Heme", "Non-Heme", "Total"), grade = grade.cutoff:5, arm = as.character(unique(cytox$arm))))
  
  
  newtable <- table1.raw %>%
    group_by(category, arm, grade) %>%
    summarise(count = n()) %>% 
    full_join(template, by = c("category", "grade", "arm")) %>% 
    mutate(count = ifelse(is.na(count), as.integer(0), count)) %>% 
    left_join(num.eval, by = 'arm') %>%
    mutate(percent = round((count/Freq)*100, digits = 1),
           grade_char = paste0("Grade ", grade, " Event"),
           percent_char = ifelse(count == 0, "(0.0%)", paste0('(', percent, '%)'))) %>%
    dplyr::select(category, grade_char, arm, count, percent_char) %>% 
    ungroup()
  
  newtable <- newtable %>% 
     mutate(category = factor(category, levels = c("Total", "Heme", "Non-Heme"), 
                                                    labels = c("Total", "Hematologic Adverse Events", "Non-Hematologic Adverse Events"))) %>% 
    arrange(category, grade_char, arm) 
  
  categories <- newtable %>% 
    count(category)
  
  evaluable <- cytox %>% 
    filter(arm != '') %>% 
    distinct(arm) %>% 
    nrow()

  eval <- character(evaluable)
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
  
  summary_title <- paste0("Summary of Grade ", grade.cutoff, "+ Adverse Events\n",
                          rel,
                          "Number of Evaluable Patients:\n",
                          arms)
  }else if(is.null(arm.labels) & arm.include == 1){
    arms <- paste(eval, sep = '', collapse = "\t")
    
    summary_title <- paste0("Summary of Grade ", grade.cutoff, "+ Adverse Events\n",
                            rel,
                            "Number of Evaluable Patients:\n",
                            arms)
  } else if(arm.include == 0){
    arms <- sum(as.numeric(eval))
    summary_title <- paste0("Summary of Grade ", grade.cutoff, "+ Adverse Events\n",
                            rel,
                            "Number of Evaluable Patients: ",
                            arms)
  }
  myHeader <- c(summary_title = 4)
  names(myHeader) <- c(summary_title)
  new_flex_table <- newtable
  names(new_flex_table) <- c("Category", "Patients with a maximum:", "Arm", "n", "(%)")
  if(arm.include == 1){
    all.summary.table <-  new_flex_table %>% 
      as_grouped_data(groups = c('Category')) %>% 
      as_flextable(hide_grouplabel = TRUE) %>% 
      bold(j = 1, i = ~ !is.na(Category), bold = TRUE, part = "body" ) %>%
      bold(bold = TRUE, part = "footer") %>%
      italic(j = 4, italic = TRUE, part = "body") %>%
      fontsize(size = 10, part = "all") %>% 
      font(fontname = "Arial", part = "all") %>% 
      merge_v(j = 2, part = "body") %>%
      add_header_lines(values = summary_title) %>%
      add_footer_lines(values = "Note: Summaries are based on available patient data") %>%
      theme_box() %>% 
      width(j = c(1:4), width = c(3.5, .5, .4, .74)) %>%
      align(align = "center", part = "header") %>%
      align(j = 1, i = ~ is.na(Category), align = "center", part = "body") %>%
      align(j = c(2, 3, 4), align = "right", part = "body") %>% 
      align(i = 2, j = 1, align = "left", part = "header")
  } else if(arm.include == 0){
    all.summary.table <-  new_flex_table %>% 
      select(-Arm) %>%
      as_grouped_data(groups = c('Category')) %>% 
      as_flextable(hide_grouplabel = TRUE) %>% 
      bold(j = 1, i = ~ !is.na(Category), bold = TRUE, part = "body" ) %>%
      bold(bold = TRUE, part = "footer") %>%
      italic(j = 3, italic = TRUE, part = "body") %>%
      fontsize(size = 10, part = "all") %>% 
      font(fontname = "Arial", part = "all") %>% 
      merge_v(j = 2, part = "body") %>%
      add_header_lines(values = summary_title) %>%
      add_footer_lines(values = "Note: Summaries are based on available patient data") %>%
      theme_box() %>% 
      width(j = c(1:3), width = c(3.5, .4, .74)) %>%
      align(align = "center", part = "header") %>%
      align(j = 1, i = ~ is.na(Category), align = "center", part = "body") %>%
      align(j = c(2, 3), align = "right", part = "body") %>% 
      align(i = 2, j = 1, align = "left", part = "header")
  }
  out <- list(table = all.summary.table, data = newtable)
  class(out) <- "aerpt.summary"
  return(out)
}

#' @export
print.aerpt.summary <- function(x, ...)
{
  print(x$table)
  invisible(x)
}
