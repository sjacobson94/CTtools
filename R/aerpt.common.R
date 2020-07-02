#' Create listing of commomly occurring AEs for Alliance DSMB reports
#' 
#' This function is used to create commonly occurring AE listing output for Alliance DSMB reports.
#' 
#' @param cytox The toxicity dataset containing multiple observations per patient with one observation per AE > grade 0
#' @param id The patient ID variable used in the dataset
#' @param threshold Numeric threshold level for commonly occurring Grade 3+ AEs. Default to 10
#' @param grade.cutoff Grade cutoff for minimum AE grade. Inclusive
#' @param rel.cutoff Attribution cutoff for minimum AE attribution. Inclusive
#' @param rel.include Value of 0 or 1. If 0, "Regardless of Attribution" printed in header. If 1, "At least possibly related" printed in header
#' @param byvar Variable used to separate summary by arm. If not present, "All Patients" printed for arm. Set to NULL if calculations are desired for all patients together.
#' @param arm.labels Option to include additional Arm labels in the table header. Add as vector of characters.
#' @param toxvar Variable containing numeric toxicity codes.
#' @param relvar Variable containing numeric relation values 1-5.
#' @param gradevar Grade variable, numeric.
#' @return This function creates a listing of the AE events for the Alliance DSMB report
#' @examples 
#' 
#' # With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' all.common <- aerpt.common(cytox, id = "dcntr_id", byvar = "arm", 
#' grade.cutoff = 3, rel.cutoff = 0, rel.include = 0)
#' all.common$grade3
#' all.common$grade45
#' 
#' @export
#' @author Sawyer Jacobson
#' 
# Last Edited: 2/29/2020

aerpt.common <- function(cytox, id, threshold = 10, grade.cutoff = 3, rel.cutoff = 0, rel.include = 1, 
                         byvar="arm", arm.labels = NULL, toxvar = "toxicity", relvar = "rel_smed", gradevar = "grade") 
  {
  cytox <- as.data.frame(cytox)
  
  arm.labels <- arm.labels
  
  #if(is.null(byvar)) arm.include = 0
  
  # Setting up attribution variable
  if(rel.include == 1)
    {
    rel = "At least possibly related\n"
  } else 
    {
    rel = "Regardless of Attribution\n"
  }
  cytox <- cytox %>%
    rename(dcntr_id = tidyselect::all_of(id), 
           toxicity = tidyselect::all_of(toxvar), 
           grade = tidyselect::all_of(gradevar),
           rel_smed = tidyselect::all_of(relvar))
  
  ######  Make sure the arm exists and to use the correct arm otherwise use All Patients#######
  if (any(colnames(cytox) == byvar)) 
    { 
    cytox <- cytox %>%
      rename(arm = tidyselect::all_of(byvar)) %>%
      mutate(arm = as.character(arm))
  } else 
    {
    cytox$arm <- rep("A", length(cytox$dcntr_id)) 
  }
  
  
  ####Will need to get the number eval by each group
  ###take out all the baseline toxicity if there are any
  if (any(colnames(cytox) == "cycle"))
    {
    cytox <- cytox[which(cytox$cycle>0),] 
  }
  
  unique.list <- cytox[!duplicated(cytox[, "dcntr_id"]),]
  
  num.eval <- base::table(unique.list$arm) %>% 
    data.frame() %>% 
    mutate(arm = as.character(Var1)) %>% 
    dplyr::select(-Var1)
  
  #if(arm.include == 0) TRUE
  if(!is.null(arm.labels) & length(arm.labels) != nrow(num.eval)) stop('Number of arm levels not equal to number of arm labels')
  
  #######Now, lose all the grade 0 and cycle 0#####
  ###IF Cycle exists, which it won't in my theradex data
  if (any(colnames(cytox) == "cycle"))
    {
    cytox <- cytox %>% 
      mutate(rel_smed = ifelse(is.na(rel_smed), 1, rel_smed)) %>%
      filter(cycle > 0 & grade >= grade.cutoff & rel_smed >= rel.cutoff)
  }
  cytox1 <- cytox #%>%
    #rename(soc_cytox = soc) %>%
    #dplyr::select(-c(bodysys))
  
  ######Create CRTOX, the worst grade per patient, per toxicity
  crtox <- cytox1 %>%
    group_by(dcntr_id, toxicity) %>% 
    slice(which.max(grade)) %>%
    dplyr::select(dcntr_id, arm, toxicity, grade) 
  
  # bodysys <- cytox %>% 
  #   dplyr::select(bodysys, toxicity) %>% 
  #   distinct(toxicity, bodysys)
  ##Count up each tox/arm/grades...
  table2 <- crtox %>%
    group_by(toxicity, arm, grade) %>%
    summarise(count = n()) #%>% 
    #left_join(bodysys, by = "toxicity")
  
  
  ###seperate by grade to mutate their own column..  may be a better way, but this is good for me :)
  
  grade3plus <- crtox %>%
    filter(grade >= 3) %>%
    group_by(toxicity, arm) %>%
    summarise(`3` = n()) %>% 
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
  
  ######Add toxterm and SOC##########
  data("tox_codes")
  
  #formats <- haven::read_sas("/people/ccs4/ccsicprd/cc-sas-all/sasdata/mart/xstudy/toxcodes.sas7bdat") %>%
  tox_codes <- tox_codes %>%
    select("v5Meddra_Term", "v5Meddra_Code", "BODYSYS", "soc") %>%
    rename(v5meddra.term = v5Meddra_Term, v5meddra.code = v5Meddra_Code)#, bodysys = BODYSYS)
  
  tox_codes <- tox_codes[!is.na(tox_codes$v5meddra.code) & !duplicated(tox_codes$v5meddra.code),]
  tox_codes$toxicity <- tox_codes$v5meddra.code
  tox_codes$toxterm <- tox_codes$v5meddra.term 
  
  ####merge and assign correct term if needed 
  template.table2 <- left_join(template.table, tox_codes,  by = "toxicity")
  
  total_3 <- grade3plus %>%
    rename(total = `3`) %>%
    dplyr::select(toxicity, arm, total) %>%
    group_by(toxicity) %>%
    summarise(total = sum(total)) %>%
    left_join(tox_codes, by = "toxicity") %>%
    select(toxterm, total) %>%
    mutate(total_display = paste0(total, " (", round(total/sum(num.eval[,1])*100, 1), "%)"),
           total_perc = round(total/sum(num.eval[,1])*100, 1)) %>%
    select(-total)
  
  total_4 <- grade4 %>%
    rename(total = `4`) %>%
    dplyr::select(toxicity, arm, total) %>%
    group_by(toxicity) %>%
    summarise(total = sum(total)) %>%
    left_join(tox_codes, by = "toxicity") %>%
    select(toxterm, total) %>%
    mutate(total_display = paste0(total, " (", round(total/sum(num.eval[,1])*100, 1), "%)"),
           total_perc = round(total/sum(num.eval[,1])*100, 1)) %>%
    select(-total)
  
  total_5 <- grade5 %>%
    rename(total = `5`) %>%
    dplyr::select(toxicity, arm, total) %>%
    group_by(toxicity) %>%
    summarise(total = sum(total)) %>%
    left_join(tox_codes, by = "toxicity") %>%
    select(toxterm, total) %>%
    mutate(total_display = paste0(total, " (", round(total/sum(num.eval[,1])*100, 1), "%)"),
           total_perc = round(total/sum(num.eval[,1])*100, 1)) %>%
    select(-total)
  
  try1 <- template.table2 %>%
    left_join(grade3plus, by = c("toxicity", "arm")) %>%
    left_join(grade4, by = c("toxicity", "arm")) %>%
    left_join(grade5, by = c("toxicity", "arm")) #%>%
    #select(-BODYSYS) %>%
    #left_join(bodysys, by = "toxicity")
  ##Add 0 to all toxicity arm grade combinations that don't exist.
  try1[is.na(try1)] <- 0
  
  #data("bodysys_formats")
  
  try2 <- try1 #%>%
    #left_join(bodysys_formats, by = c("bodysys" = "start"))
  
  try3 <- try2 %>%
    # mutate(category = ifelse(toxicity %in% c(10002272,10007839,10025256,10029366,10035528,10048580,10049182,
    #                                          10005329,10019150,10019491,10024378,10025258,10028533,10041633,
    #                                          10055599,10065973), "Hematologic Adverse Events", "Non-Hematologic Adverse Events"),
    #        eventcat = ifelse(toxicity %in% c(10002272,10007839,10025256,10029366,10035528,10048580,10049182,
    #                                          10005329,10019150,10019491,10024378,10025258,10028533,10041633,
    #                                          10055599,10065973), "Blood/Bone Marrow", 
    #                          ifelse(toxicity == 10016288, "Blood and lymphatic sys disord", soc)),
    #        label = ifelse(eventcat == "Blood/Bone Marrow", "Blood/Bone Marrow", label)) %>%
    left_join(num.eval, by = "arm")
  
  grade_list <- c("3-Severe\n  n (%)", "4-LifeThr\n  n (%)", "5-Lethal\n  n (%)")
  
  #######Add the percentages based on the num.eval for each arm
  grade_3 <- try3 %>%
    arrange(toxterm, arm)  %>%
    select(toxterm, arm, `3`) %>% 
    group_by(toxterm, arm) %>%
    tidyr::spread(arm, `3`) %>%
    as.data.frame()
  
  grade_3_try <- grade_3 
  for(i in 1:nrow(num.eval))
    {
    grade_3[i + 1] <- paste0(grade_3_try[,i+1], " (", round(grade_3_try[,i+1]/num.eval[i,1]*100, 1), "%)") 
    grade_3[i + nrow(num.eval) + 1] <- round(grade_3_try[,i+1]/num.eval[i,1]*100, 1)
  }
  
  grade_3_index <- apply(data.frame(grade_3[, (2+nrow(num.eval)):ncol(grade_3)]), 1, function(x) any(x >= threshold))
  grade_3_thresh <- grade_3[grade_3_index, 1:(1+nrow(num.eval))] %>%
    left_join(total_3, by = "toxterm")  %>% 
    arrange(desc(total_perc)) %>%
    select(-total_perc) %>%
    rename(AE = toxterm, Total = total_display) %>%
    mutate(category = "Commonly Occurring Grade 3+ AE") %>%
    select(category, everything()) 
  
  grade_4 <- try3 %>%
    arrange(toxterm, arm)  %>%
    select(toxterm, arm, `4`) %>% 
    group_by(toxterm, arm) %>%
    tidyr::spread(arm, `4`) %>%
    as.data.frame()
  
  grade_4_try <- grade_4 
  for(i in 1:nrow(num.eval))
    {
    grade_4[i + 1] <- paste0(grade_4_try[,i+1], " (", round(grade_4_try[,i+1]/num.eval[i,1]*100, 1), "%)") 
    grade_4[i + nrow(num.eval) + 1] <- round(grade_4_try[,i+1]/num.eval[i,1]*100, 1)
  }
  
  #grade_4_index <- apply(grade_4[, (2+nrow(num.eval)):ncol(grade_4)], 1, function(x) any(x >= 0))
  grade_4_thresh <- grade_4[, 1:(1+nrow(num.eval))] %>%
    left_join(total_4, by = "toxterm")  %>% 
    filter(!is.na(total_perc)) %>%
    arrange(desc(total_perc)) %>%
    select(-total_perc) %>%
    rename(AE = toxterm, Total = total_display) %>%
    mutate(category = "Grade 4") %>%
    select(category, everything())
  
  grade_5 <- try3 %>%
    arrange(toxterm, arm)  %>%
    select(toxterm, arm, `5`) %>% 
    group_by(toxterm, arm) %>%
    tidyr::spread(arm, `5`) %>%
    as.data.frame()
  
  grade_5_try <- grade_5 
  for(i in 1:nrow(num.eval))
    {
    grade_5[i + 1] <- paste0(grade_5_try[,i+1], " (", round(grade_5_try[,i+1]/num.eval[i,1]*100, 1), "%)") 
    grade_5[i + nrow(num.eval) + 1] <- round(grade_5_try[,i+1]/num.eval[i,1]*100, 1)
  }
  
  #grade_4_index <- apply(grade_4[, (2+nrow(num.eval)):ncol(grade_4)], 1, function(x) any(x >= 0))
  grade_5_thresh <- grade_5[, 1:(1+nrow(num.eval))] %>%
    left_join(total_5, by = "toxterm")  %>% 
    filter(!is.na(total_perc)) %>%
    arrange(desc(total_perc)) %>%
    select(-total_perc) %>%
    rename(AE = toxterm, Total = total_display) %>%
    mutate(category = "Grade 5/Deaths on Treatment") %>%
    select(category, everything())
  
  final45 <- bind_rows(grade_4_thresh, grade_5_thresh) #%>%
    #mutate(AE = paste0("  ", AE))
  
  evaluable <- cytox %>% 
    filter(arm != '') %>% 
    distinct(arm) %>% nrow()
  
  eval <- character(evaluable)

  for(i in 1:evaluable)
    {
    if(!is.null(arm.labels))
      {
      eval[i] <- paste0(arm.labels[i]," (Arm ", num.eval[i, 2], "=", num.eval[i, 1], ")")
    } else if (is.null(arm.labels) )
      {
      eval[i] <- paste0("Arm ", num.eval[i, 2], "=", num.eval[i, 1])
    } else
      {
      eval[i] <- as.integer(num.eval[i, 1])
    }
  }
  if(!is.null(arm.labels))
    {
    arms <- paste(eval, sep = '', collapse = "\n")
    
    grade_3_title <- paste0("Listing of Grade 3+ Adverse Events\nMax Grade Per Patient Per Event\nOccurrences Across Any Level > ",
                            threshold, "%\n",
                            rel,
                            "Number of Evaluable Patients:\n",
                            arms)
    grade_45_title <-  paste0("Listing of all Grade 4 or 5 Adverse Events\nMax Grade Per Patient Per Event\n",
                              rel,
                              "Number of Evaluable Patients:\n",
                              arms)
  } else if(is.null(arm.labels)) 
    {
    arms <- paste(eval, sep = '', collapse = "\t")
    grade_3_title <- paste0("Listing of Grade 3+ Adverse Events\nMax Grade Per Patient Per Event\nOccurrences Across Any Level > ",
                            threshold, "%\n",
                            rel,
                            "Number of Evaluable Patients:\n",
                            arms)
    grade_45_title <-  paste0("Listing of all Grade 4 or 5 Adverse Events\nMax Grade Per Patient Per Event\n",
                              rel,
                              "Number of Evaluable Patients:\n",
                              arms)
  }
  # else if(is.null(arm.labels)){
  #   arms <- sum(as.numeric(eval))
  #   grade_3_title <- paste0("Listing of Grade ", grade.cutoff, "+ Adverse Events\nMax Grade Per Patient Per Event\nOccurences Across Any Level ≥ ",
  #                           threshold, "%\n",
  #                           rel,
  #                           "Number of Evaluable Patients: ",
  #                           arms)
  # }

  if(nrow(num.eval) == 1)
    {
    if(nrow(grade_3_thresh) == 0) grade3_table <- grade_3_thresh
    else
    {
      col_names <- names(grade_3_thresh)[-c(1, 2, ncol(grade_3_thresh))]
      grade3_table <- grade_3_thresh %>% 
        select(-Total) %>%
        as_grouped_data(groups = c('category')) %>%
        as_flextable(hide_grouplabel = TRUE) %>%
        delete_part(part = "header") %>%
        add_header_row(values = c("", col_names), colwidths = rep(1, ncol(grade_3_thresh) - 2)) %>%
        bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
        fontsize(size = 10, part = "all") %>% 
        font(fontname = "Arial", part = "all") %>% 
        add_header_lines(values = grade_3_title) %>%
        bold(bold = TRUE, part = "header") %>% 
        theme_box() %>%
        width(j = c(1:(nrow(num.eval) + 1)), width = c(2.63, rep(.95, nrow(num.eval)))) %>%
        align(align = "center", part = "header") %>%
        align(align = "left", part = "body") %>%
        padding(j = 1, i = ~ is.na(category), padding.left = 14, part = 'body')
    }
    
    if(nrow(final45) == 0) grade45_table <- final45
    else
    {
      col_names <- names(final45)[-c(1, 2, ncol(final45))]
      grade45_table <- final45 %>% 
        select(-Total) %>%
        as_grouped_data(groups = c('category')) %>%
        as_flextable(hide_grouplabel = TRUE) %>%
        delete_part(part = "header") %>%
        add_header_row(values = c("", col_names), colwidths = rep(1, ncol(final45) - 2)) %>%
        bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
        fontsize(size = 10, part = "all") %>% 
        font(fontname = "Arial", part = "all") %>% 
        add_header_lines(values = grade_45_title) %>%
        bold(bold = TRUE, part = "header") %>% 
        theme_box() %>%
        width(j = c(1:(nrow(num.eval) + 1)), width = c(2.63, rep(.95, nrow(num.eval)))) %>%
        align(align = "center", part = "header") %>%
        align(align = "left", part = "body")%>%
        padding(j = 1, i = ~ is.na(category), padding.left = 14, part = 'body')
    }
  } else
  {
    if(nrow(grade_3_thresh) == 0) grade3_table <- grade_3_thresh
    else
    {
      col_names <- names(grade_3_thresh)[-c(1, 2)]
      grade3_table <- grade_3_thresh %>% 
        as_grouped_data(groups = c('category')) %>%
        as_flextable(hide_grouplabel = TRUE) %>%
        delete_part(part = "header") %>%
        add_header_row(values = c("", col_names), colwidths = rep(1, ncol(grade_3_thresh) - 1)) %>%
        bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
        fontsize(size = 10, part = "all") %>% 
        font(fontname = "Arial", part = "all") %>% 
        add_header_lines(values = grade_3_title) %>%
        bold(bold = TRUE, part = "header") %>% 
        theme_box() %>%
        width(j = c(1:(nrow(num.eval) + 2)), width = c(2.63, rep(.95, nrow(num.eval) + 1))) %>%
        align(align = "center", part = "header") %>%
        align(align = "left", part = "body")%>%
        padding(j = 1, i = ~ is.na(category), padding.left = 14, part = 'body')
    }
    if(nrow(final45) == 0) grade45_table <- final45
    else
    {
      col_names <- names(final45)[-c(1, 2)]
      grade45_table <- final45 %>% 
        as_grouped_data(groups = c('category')) %>%
        as_flextable(hide_grouplabel = TRUE) %>%
        delete_part(part = "header") %>%
        add_header_row(values = c("", col_names), colwidths = rep(1, ncol(final45) - 1)) %>%
        bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body" ) %>%
        fontsize(size = 10, part = "all") %>% 
        font(fontname = "Arial", part = "all") %>% 
        add_header_lines(values = grade_45_title) %>%
        bold(bold = TRUE, part = "header") %>% 
        theme_box() %>%
        width(j = c(1:(nrow(num.eval) + 2)), width = c(2.63, rep(.95, nrow(num.eval) + 1))) %>%
        align(align = "center", part = "header") %>%
        align(align = "left", part = "body")%>%
        padding(j = 1, i = ~ is.na(category), padding.left = 14, part = 'body')
    }
  }
  
  out <- list(grade3plus = grade3_table, grade45 = grade45_table)
  class(out) <- "aerpt.common"
  out
}
  