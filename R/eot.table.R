#' @title
#' End of Treatment Table
#' @description 
#' This function creates the end of treatment table for Alliance DSMB reports
#' 
#' @details 
#' 
#' If endatrsn column is a factor, the given values will be used. However, if the column is 
#' numeric, the values will be mapped to pre-specified descriptions.
#' 
#' Value map for numeric input:
#' 
#' NA = On Treatment
#' 
#' 1 = Treatment Completed Per Protocol
#' 
#' 2 = Patient Withdrawal/Refusal After Beginning Protocol Therapy
#' 
#' 3 = Adverse Events/Side Effects/Complications
#' 
#' 4 = Disease Progression, Relapse During Active Treatment
#' 
#' 5 = Alternative Therapy
#' 
#' 6 = Patient Off-Treatment For Other Complicating Disease
#' 
#' 7 = Death On Study
#' 
#' 8 = Other
#' 
#' 10 = Disease Progression Before Active Treatment
#' 
#' 24 = Patient Withdrawal/Refusal Prior To Beginning Protocol Therapy
#'
#' @param data A dataset with 1 row for each patient accrued. 
#' @param id Patient ID variable
#' @param endatrsn End of treatment reason variable. May be numeric or factored.
#' @param arm Arm variable.
#' @return This function returns an end of treatment table formatted for the Alliance DSMB report.
#' @examples
#' 
#' crtlibn(d = "A041501")
#' 
#' x = eot.table(casecrse, arm = "curr_arm")
#' x$table
#'
#' @export
#' @author Sawyer Jacobson
#'
# Created 3/10/2020

eot.table <- function(data, id = "dcntr_id", endatrsn = "endatrsn", arm = "arm")
{
  endat <- as.data.frame(data)
  
  endat <- endat %>%
    select(dcntr_id = tidyselect::all_of(id), 
           endatrsn = tidyselect::all_of(endatrsn), 
           arm = tidyselect::all_of(arm)) %>%
    droplevels()
  endat$arm <- as.character(endat$arm)
  num.eval <- base::table(endat$arm)
  
  if(!is.factor(endat$endatrsn))
  {
    endat$endatrsn <- factor(endat$endatrsn, levels = c(1,2,3,4,5,6,7,8,10,24), labels = c("Treatment Completed Per Protocol", 
                                                                                           "Patient Withdrawal/Refusal After Beginning Protocol Therapy",
                                                                                           "Adverse Events/Side Effects/Complications", 
                                                                                           "Disease Progression, Relapse During Active Treatment", 
                                                                                           "Alternative Therapy", 
                                                                                           "Patient Off-Treatment For Other Complicating Disease", 
                                                                                           "Death On Study", 
                                                                                           "Other", 
                                                                                           "Disease Progression Before Active Treatment", 
                                                                                           "Patient Withdrawal/Refusal Prior To Beginning Protocol Therapy"))
  }
  
  endat$endatrsn <- forcats::fct_explicit_na(endat$endatrsn, na_level = "On Treatment")
  endat$endatrsn <- forcats::fct_relevel(endat$endatrsn, "On Treatment")
  
  endat <- endat %>% 
    group_by(arm, endatrsn) %>% 
    summarise(count = n()) %>% 
    tidyr::spread(arm, count) 
  endat[is.na(endat)] <- 0
  
  endat$total <- apply(endat[2:ncol(endat)], 1, sum, na.rm = TRUE)
  
  tots <- apply(endat[2:ncol(endat)], 2, sum)
  
  off_tx_count <- apply(endat[2:nrow(endat), 2:ncol(endat)], 2, sum)
  
  for(i in 1:length(tots))
  {
    endat[[i + 1]][1] = paste0(endat[[i + 1]][1], " (", round(100*endat[[i + 1]][1]/as.numeric(tots[i]), 1), "%)")
    endat[[i + 1]][2:nrow(endat)] = paste0(as.numeric(endat[[i + 1]][2:nrow(endat)]), " (", round(100*as.numeric(endat[[i + 1]][2:nrow(endat)])/as.numeric(off_tx_count[i]), 1), "%)")
  }
  
  off_tx <- c(endatrsn = "Off Treatment", off_tx_count)

  off_tx[2:length(off_tx)] <- paste0(off_tx_count, " (", round(100*off_tx_count/as.numeric(tots), 1), "%)")
  
  endat <- endat %>%
    mutate(endatrsn = as.character(endatrsn))
  
  endat <- bind_rows(endat[1, ], off_tx, endat[2:nrow(endat), ]) %>%
    rename(Total = total)
  
  header <- c("", paste0(colnames(endat[2:ncol(endat)]), "\n", "(N=", as.numeric(tots), ")"))
  
  if(length(num.eval) > 1)
  {
    final <- endat %>% 
      flextable()  %>%
      delete_part(part = "header") %>%
      add_header_row(values = header, colwidths = rep(1, 1 + length(tots))) %>%
      theme_box() %>%
      width(width = c(3.06, rep(1.08, length(num.eval)), 1.18)) %>%
      align(i = 1:2, align = "left") %>%
      align(align = "left") %>%
      align(align = "center", part = "header") %>%
      font(fontname = "Arial") %>%
      fontsize(size = 10) %>%
      padding(i = 3:nrow(endat), padding.left = 22)
  }
  else
  {
    final <- endat %>% 
      select(-Total) %>%
      flextable()  %>%
      delete_part(part = "header") %>%
      add_header_row(values = header[-length(header)], colwidths = rep(1, length(tots))) %>%
      theme_box() %>%
      width(width = c(3.06, rep(1.08, length(num.eval)))) %>%
      align(i = 1:2, align = "left") %>%
      align(align = "left") %>%
      align(align = "center", part = "header") %>%
      font(fontname = "Arial") %>%
      fontsize(size = 10) %>%
      padding(i = 3:nrow(endat), padding.left = 22)
  }
    
  out <- list(table = final, data = endat)
  class(out) <- "eot.table"
  out
  
}

#' @export
print.eot.table <- function(x, ...)
{
  print(x$table)
  invisible(x)
}