#' Forest Plot
#'
#' This function creates the data for and plots a forest plot comparing hazard rates for categorical variables by a 2 level by variable.
#'
#' @param data A dataset containing the necessary survival variables (event, time, other covariables).
#' @param stat The event survival variable, numeric, with 1 coded as an event and 0 a censor.
#' @param time The time survival variable, numeric.
#' @param byvar The by variable with 2 distinct levels.
#' @param covars A vector of categorical covariates that are desired for analysis.
#' @param ref Reference by variable level (optional). Use this if specific by variable level is desired as the reference.
#' @param title Title for the plot.
#' @param labels An optional vector, soon to be list, of labels for the covariates in the 'covars' variable. Currently must be the same length as the 'covars' vector.
#' @param graph.pos Column position for the graph.
#' @param xaxis.breaks Distance between the x-axis breaks of the forest plot.
#' @param x.lab Label for x-axis of plot.
#' @param plot.include Logical to have plot return from function call.
#' @param x.max Max value for x-axis. CI bands for intervals that pass this value will have an arrow on the right indicating a continuing band.
#' @param ...  Additional arguments for the 'forestplot' function.
#'
#' @return This function returns a list of the forestplot produced and tables necessary to reproduce that plot.
#' @export
#'
#' @examples
#' 
#' library(arsenal)
#' data(mockstudy)
#' mockstudy %>%
#'  mutate(fu.stat = fu.stat - 1) %>%
#'  filter(arm != "A: IFL" & age.ord != "10-19" & race != "Other") %>%
#'  forest_plot(data = ., stat = "fu.stat", time = "fu.time", byvar = "arm", 
#'              covars = c("sex", "age.ord"), title = "Overall Survival", 
#'              labels = c("Gender", "Age Groups"), 
#'              xaxis.breaks = .5)
#' 
forest_plot <- function(data, stat, time, byvar, covars, ref = as.character(levels(data2$arm)[1]), labels = NULL, graph.pos = 3,
                        xaxis.breaks = .25, title = "Hazard Rates", x.lab = "Hazard Ratios", plot.include = TRUE, x.max = max(values$upper, na.rm = TRUE),
                        ...)
{
  data2 <- data.frame(time=data[[time]], stat=data[[stat]], arm=factor(data[[byvar]])) %>% 
    bind_cols(data[covars] %>%
                mutate_all(as.character))
  data2$arm <- relevel(data2$arm, ref=ref)
  
  nevent <- c(); hr <- c(); var <- c(); ci <- c(); upper <- c(); lower <- c(); pval <- c(); hrCI <- c()
  i <- 0
  for(variable in covars)
  {
    i <- i + 1
    nevent = c(nevent, "")
    hr = c(hr, "")
    if(!is.null(labels) & length(labels) == length(covars)) 
    {
      var = c(var, labels[i])
    }
    else{
      var = c(var, variable)
    }
    ci = c(ci, "")
    upper = c(upper, "")
    lower = c(lower, "")
    pval = c(pval, "")
    hrCI = c(hrCI, "")
    for(level in unique(data2[[variable]]))
    {
      if(!is.na(level))
      {
        data3 <- data2[data2[[variable]] == level, ]
        fit = survival::coxph(survival::Surv(time, stat) ~ arm, data = data3)
        tidy_fit = broom::tidy(fit, exponentiate = TRUE)
        nevent = c(nevent, paste0(fit$nevent, "/", fit$n))
        hr = c(hr, round(tidy_fit$estimate, 2))
        ci = c(ci, paste0("(", round(tidy_fit$conf.low, 3), ", ", round(tidy_fit$conf.high, 3), ")"))
        hrCI = c(hrCI, paste0(round(tidy_fit$estimate, 2), " ", "(", round(tidy_fit$conf.low, 3), ", ", round(tidy_fit$conf.high, 3), ")"))
        upper = c(upper, round(tidy_fit$conf.high, 3))
        lower = c(lower, round(tidy_fit$conf.low, 3))
        var = c(var, paste0("    ", level))
        pval = c(pval, round(tidy_fit$p.value, 3))
      }
    }
  }
  forest_data <- data.frame(
    variables = c(paste0(levels(data2$arm)[2], " vs ", levels(data2$arm)[1], " (ref)"), "Subgroup", var),
    nevent = c("", "Events/N", nevent),
    hazard = c("", "Hardard Rate", hr),
    ci = c("", "95% CI", ci),
    lower = c("", "Lower bound", lower),
    upper = c("", "Upper bound", upper), 
    pval = c("", "P-value", pval),
    hrCI = c("", "HR (95% CI)", hrCI)
  ) %>%
    filter(upper != "Inf")
  
  values <- forest_data[3:nrow(forest_data), c("hazard", "lower", "upper")] %>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric)
  
  table_text <- forest_data[, c("variables", 'nevent', "hrCI", "pval")]

  if(plot.include == TRUE)
  {
  forestplot::forestplot(labeltext = table_text, graph.pos = graph.pos, 
             mean=c(NA, NA, values$hazard), 
             lower=c(NA, NA,values$lower), upper=c(NA, NA, values$upper),
             title = title,
             hrzl_lines = TRUE,
             zero = 1,
             clip = c(0, x.max),
             xticks = seq(plyr::round_any(min(values$lower, na.rm = TRUE), xaxis.breaks), 
                          plyr::round_any(max(pmin(x.max, values$upper, na.rm = TRUE)), xaxis.breaks), by = xaxis.breaks),
             xticks.digits = 2,
             txt_gp = forestplot::fpTxtGp(label=grid::gpar(cex=.8),
                            ticks=grid::gpar(cex=.8),
                            xlab=grid::gpar(cex = 1.2),
                            title=grid::gpar(cex = 1.2)),
             col=forestplot::fpColors(box="black", lines="black", zero = "gray50"),
             cex=0.9, 
             lineheight = "auto",
             boxsize=0.15, 
             colgap=unit(1,"mm"),
             lwd.ci=2, 
             ci.vertices=TRUE, 
             ci.vertices.height = 0.1,
             xlab = x.lab,
             ...
  )
  }
  out <- list(table = table_text, values = values)
  class(out) <- "forest_plot"
  invisible(out)
}


# @export
# print.forest_plot <- function(x, ...)
# {
#   print(x$table)
#   invisible(x)
# }
