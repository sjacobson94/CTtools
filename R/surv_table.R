#' Create survfit table
#' 
#' Use this to create a table for survival analysis. Creates a data.frame containing sample size, number of events, median survival time, and 95\% CI between by groups.
#'
#' @param data data frame that contains time to event, event, and by variables required to create survfit object 
#' @param time Time variable
#' @param event Event variable
#' @param by By variable
#'
#' @return This function returns a data.frame containing a column for each level of the by variable and the pvalue with rows "N", "Events" (# events),
#'  "Median" (median survival time), "95\% CI" (for median survival time).
#' @export
#'
# Created 02/12/2020

surv_table <- function(data, time, event, by){
  # saving the needed variables in a data.frame()
  data1 <- data.frame(time = data[[time]], event = data[[event]], by_var = data[[by]])
  # creating the survfit object
  surv <- data1 %>% 
    survival::survfit(survival::Surv(time, event) ~ by_var, data = .)
  # getting the logrank pvalue between byvars 
  logrank <- survminer::surv_pvalue(surv, data=data1, method="survdiff")
  # creating the data.frame() that contains N, median, and the 95% CIs
  sum.table <- as.data.frame(summary(surv)$table[,c("records","events","median","0.95LCL","0.95UCL")]) %>%
    mutate(median = round(median, 2), 
           CI = paste0("(", round(`0.95LCL`, 2), ", ", round(`0.95UCL`, 2), ")")) %>% 
    select(-c(`0.95LCL`, `0.95UCL`)) %>%
    t() %>%
    as.data.frame()
  # adding column names based on the by_var
  colnames(sum.table) <- unique(data1$by_var)
  # adding the rownames
  rownames(sum.table) <-  c("N", "Events", "Median", "95% CI")
  # adding the pvalue column
  sum.table$`p value` <- c(as.character(round(logrank$pval, 3)), rep("", nrow(sum.table)-1))
  out <- list(table = sum.table, surv = surv)
  class(out) <- "surv_table"
  out
}
