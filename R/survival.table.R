#' Survival table
#'
#' This function creates a table/dataset with number of events/N, 95\% median (months) survival CI, 95\% survival rate CI (eg. 12-month PFS rate), and 95\% HR CI.
#'
#' @param data A dataset containing columns necessary for survival analysis.
#' @param stat Survival event variable
#' @param time Survival time variable, currently only supports months
#' @param label Label for survival variable
#' @param arm Arm variable
#' @param ref Option to specify the reference arm. Must be entered how it is in the data.
#' @param rate_span Vector of time levels desired for survival rate. Used in conjunction with time_level.
#' @param time_level Time level desired. Currently supports "months" and "years".
#' @param rate_label Label for survival rate in the table (eg. rate_label = "PFS" gives "12-month PFS rate).
#' @param one.sided Logical. If TRUE returns one-sided p-values. Default is FALSE
#'
#' @return This function returns a nicely formatted table of survival analysis statistics.
#' @export
#'
survival.table <- function(data, stat, time, label="Rate", arm="arm", ref=as.character(levels(data2$arm)[1]), rate_span = c(12), 
                  time_level = "month", rate_label = "PFS", one.sided = FALSE)
{
  ##Create dataset to work with
  data2 <- data.frame(time=data[[time]], stat=data[[stat]], arm=factor(data[[arm]]))
  data2$arm <- relevel(data2$arm, ref=ref)
  ##Create survfit object
  km.fit <- survival::survfit(survival::Surv(time, as.numeric(stat))~arm, data=data2)
  c.fit <- survival::coxph(survival::Surv(time, as.numeric(stat))~arm, data=data2)
  ##Output p-value and hazard ratio to global variables
  logrank <- survminer::surv_pvalue(km.fit, data=data2, method="survdiff")
  hazardratio <- exp(c.fit$coefficients)
  beta <- coef(c.fit)
  se   <- sqrt(diag(c.fit$var))
  p    <- round(1 - pchisq((beta/se)^2, 1), 5)
  rate <- summary(km.fit, times = c(if(tolower(time_level) == "year") rate_span*12 else rate_span))
  w <- length(rate_span) 
  
  list.df <- list(data.frame(arm = rep(c(levels(data2$arm)), c(rep(w, w)))), 
                  data.frame(time = rep(if(tolower(time_level) == "year") rate$time * 12 else rate$time, if(length(rate$time) < w*length(levels(data2$arm))) length(levels(data2$arm)) else 1)), 
                  data.frame(survival = rate$surv), data.frame(lower = rate$lower), data.frame(upper = rate$upper))
  
  max.rows <- max(unlist(lapply(list.df, nrow), use.names = F))
  
  list.df <- lapply(list.df, function(x) {
    na.count <- max.rows - nrow(x)
    if (na.count > 0L) {
      na.dm <- matrix(NA, na.count, ncol(x))
      colnames(na.dm) <- colnames(x)
      rbind(x, na.dm)
    } else {
      x
    }
  })
  
  rateCI <- do.call(cbind, list.df)
  
  rateCI <- rateCI %>%
    mutate(ci = paste0(round(survival, 2), "\n(", round(lower, 3), ", ", round(upper, 3), ")")) %>%
    select(arm, time, ci) %>% 
    tidyr::spread(arm, ci) %>% 
    mutate(pval = "") %>%
    as.matrix()
  
  ##Create summary table
  sum.table<-as.table(summary(km.fit)$table[,c("records","events","median","0.95LCL","0.95UCL")])
  sum.table<-as.table(rbind(cbind(paste0(sum.table[,2],"/",sum.table[,1]), 
                                  paste0(sprintf('%.2f', sum.table[,3]), "\n(", sprintf('%.2f', sum.table[,4]), ", ", 
                                         sprintf('%.2f', sum.table[,5]), ")"), 
                                  rbind("Ref", matrix(paste0(round(hazardratio,2), "\n(", paste0(round(summary(c.fit)$conf.int[,3],3)), 
                                                             ", ", paste0(round(summary(c.fit)$conf.int[,4],3)), ")")))), 
                            cbind("", matrix(c(round(logrank$pval / if(one.sided) 2 else 1, 3), round(unname(p), 3)), ncol = 2))))
  
  out <-  t(sum.table)
  colnames(out) <-  c(levels(data2$arm), "P")
  
  out <- rbind(rep("", ncol(out)), out[1:2, ], rateCI[, 2:ncol(rateCI)], out[3, ])
  
  rownames(out) <-  c(label, "Events/N", "Median (months)\n(95% CI)", paste0(rate_span, "-", time_level, " ", rate_label, " rate\n (95% CI)"), "HR\n(95% CI)")
  
  out <- out %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column(var = " ") %>%
    mutate_all(as.character)
  out1 <- out %>% tidy_tableby()
  final <- list(table = out1, data = out)
  class(final) <- "survival.table"
  final
}

#' @export
print.survival.table <- function(x, ...)
{
  print(x$table)
  invisible(x)
}


# tidy_tableby <- function(x) 
# {
#   header <- colnames(x)
#   bold_ind <- which(x[, 2] == '')
#   padd_ind <- which(x[, 2] != '')
#   colnames(x)[1] <- "var"
#   x$var <- trimws(str_remove_all(x$var, "[\\[\\]]"))
#   out <- x %>%
#     flextable() %>%
#     bold(j = 1, i = bold_ind, part = "body") %>%
#     delete_part(part = "header") %>%
#     autofit() %>%
#     add_header_row(values = header, colwidths = rep(1, ncol(x))) %>%
#     theme_box() %>%
#     padding(i =padd_ind, j = 1, padding.left = 20) %>%
#     align(j = 2:ncol(x), align = "center", part = "all") %>%
#     align(j = 1, align = "left")
#   out
# }
