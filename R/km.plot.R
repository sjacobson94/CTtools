
###Create Function to consistently plot KM curves
#' KM plotting function
#' 
#' @description This is a function to consistently plot aesthetically pleasing KM plots using survminer::ggsurvplot() that include risk table below and survival table within.
#' 
#' @param data The data.
#' @param stat The survival status variable used in survfit()
#' @param time The survival time variable used in survfit()
#' @param arm By variable used in survfit()
#' @param x.lab Label for x-axis of KM plot
#' @param x.time Specific unit for the survival time variable in x-axis label, character (eg. "Days", "Months", "Years").
#' @param y.lab Label for y-axis of KM plot
#' @param ref Specify reference arm for HR values in survival table.
#' @param x.max Max value for the x-axis
#' @param x.by By value for tick marks on the x-axis based on the scaled time variable.
#' @param surv.legend.labs KM plot legend arm labels
#' @param table.legend.labs Survival table arm labels
#' @param color Color for KM plot lines
#' @param lty Type of lines for KM plot
#' @param table.font.size Font size for the survival table
#' @param y.scale Scale for y-axis. Options are "default" (decimals, 0.20) and "percent" (percent, 20\%)
#' @param one.sided TRUE for one-sided p-value, FALSE for two-sided p-value
#' @param risk.table Include risk table (TRUE/FALSE)
#' @param risk.table.height Title for the at risk table below plot
#' @param table.include Include survival table (TRUE/FALSE)
#' @param pval.include Include p-value. if TRUE and pval = NULL, logrank p-value from survfit() is included
#' @param pval User specified p-value. This can be used to specify p-value used if something different than logrank is desired. Can be character or numeric.
#' @param pval.loc Vector of length 2 to specify the position for the p-value if necessary.
#' @param size Overall size of the plot?
#' @param risk.table.height Percent of plot size taken up by the risk table.
#' @param surv.plot.height Percent of plot size taken up by the KM plot.
#' @param table.x.min Minimum x-axis value for the survival table if table.loc = "topright".
#' @param table.x.max Maximum x-axis value for the survival table if table.loc = "bottomleft".
#' @param x.scale Numeric variable to allow for scaling of the time variable. Will divide the time variable by the value given.
#' @param surv.order Character vector. Allows for the display order of the variables in the legend, plot, survival table, and risk table to be specified. Must include all levels of the "arm" variable.
#' @param surv.legend.title Legend title to appear at the top of the plot (i.e. "Arm: ", "Strata: ').
#' @param legend.font Font size for the legend.
#' @param theme Desired theme for the KM plot and risk table.
#' @param table.loc Character. Specify where the survival table should go inside the plot. Currently support "topright" and "bottomleft".
#' @param ... Additional arguments passed onto ggsurvplot()
#'
#' @return This function returns a publication quality KM curve.
#' @export
#'
#' @examples
#' 
#' library(arsenal)
#' data(mockstudy)
#' mockstudy %>% 
#'     km.plot(data = ., stat = "fu.stat", time = "fu.time", arm = "arm", 
#'             table.loc = "topright", x.scale = 30.44)
#' 
km.plot <- function(data, stat, time, arm = arm, x.lab = "Time (Months)", x.time = NULL, y.lab = "% Event-Free", ref = as.character(levels(data2$arm)[1]), 
                    x.max = max(data2$time), x.by = 12, surv.legend.labs = levels(data2$arm2), table.legend.labs = levels(data2$arm2), 
                    color = RColorBrewer::brewer.pal(8,"Set1"), lty = 1, table.font.size = 9, y.scale = "percent", 
                    one.sided = FALSE, risk.table = TRUE, risk.table.title = "Patients at Risk", table.include = TRUE, 
                    pval.include = TRUE, pval = NULL, pval.loc = NULL, size = .5, risk.table.height = .175, 
                    surv.plot.height = .825, table.x.min = x.max/5, table.x.max = x.max/1.5, x.scale = 1,
                    surv.order = NULL, surv.legend.title = "Arm: ", legend.font = 9, theme = survminer::theme_survminer(),
                    table.loc = "topright", ...){
  ##Create dataset to work with
  time <- rlang::enquo(time)
  stat <- rlang::enquo(stat)
  arm <- rlang::enquo(arm)
  
  data2 <- data %>%
    dplyr::select(time = !!time, stat = !!stat, arm = !!arm) %>% 
    dplyr::mutate_at(vars(arm), as.factor)
  # data2 <- data.frame(time=data[[time]], stat=data[[stat]], arm=factor(data[[arm]]))
  if(is.numeric(x.scale)) 
  {
    data2$time <- data2$time/x.scale
  } 
  
  data2$arm2 <- relevel(data2$arm, ref=ref)
  ##Create survfit object
  km.fit <- survival::survfit(survival::Surv(time, as.numeric(stat))~arm2, data=data2)
  c.fit <- survival::coxph(survival::Surv(time, as.numeric(stat))~arm2, data=data2)
  ##Output p-value and hazard ratio to global variables
  logrank <- round(survminer::surv_pvalue(km.fit, data=data2, method="survdiff")$pval / if(one.sided) 2 else 1, 3)
  hazardratio <- exp(c.fit$coefficients)
  ##Create summary table
  sum.table <- as.table(summary(km.fit)$table[, c("records","events","median","0.95LCL","0.95UCL")])
  sum.table <- as.table(cbind(paste0(sum.table[,2],"/",sum.table[,1]), 
                              paste0(sprintf('%.2f', sum.table[,3])," (", sprintf('%.3f', sum.table[,4]), ",", sprintf('%.3f', sum.table[,5]), ")"),
                              rbind("ref", matrix(paste0(round(hazardratio,2), " (", paste0(round(summary(c.fit)$conf.int[,3],3)), ",",
                                                         paste0(round(summary(c.fit)$conf.int[,4],3)), ")"))))) %>%
    as.data.frame.matrix() 
  rownames(sum.table) <- levels(data2$arm2)
  if(!is.null(surv.order) & length(surv.order) == length(unique(data2$arm)))
  {
    sum.table = sum.table[surv.order, ]
    data2$arm2 <- factor(data2$arm, levels = surv.order)
    km.fit <- survival::survfit(survival::Surv(time, as.numeric(stat))~arm2, data=data2)
    
  }

  if(!is.null(x.time)){
    x.lab = paste0("Time (", x.time, ")")
  }
  
  pval.coord <- if(table.loc == "topright" & is.null(pval.loc)) # change this so it works better, probably by leaving this as the default in the function call and letting this be a vector
  {
    c(x.max*(2/3), 1 - (.175 * length(levels(data2$arm))))
  } 
  else if(table.loc == "bottomleft" & is.null(pval.loc)){
    c(x.max*(1/3), (.175 * length(levels(data2$arm))))
  } else pval.loc
  
  ##Create KM Curve
  km_plot <- survminer::ggsurvplot(km.fit, 
                                   data=data2,  
                                   censor.shape="|",
                                   censor.size=4, 
                                   xscale=1, 
                                   break.x.by=x.by, 
                                   xlim=c(0, x.max), 
                                   xlab=x.lab, 
                                   ylab=y.lab, 
                                   palette=color, 
                                   linetype=lty,
                                   size = size, 
                                   surv.scale = y.scale,
                                   #pval.method=TRUE,
                                   pval = if (!is.null(pval)) pval else if(pval.include) paste0("Log-rank\n p = ", logrank) else FALSE, 
                                   pval.coord = pval.coord, 
                                   #pval.method.coord=c(x.max*(2/3),.9-(.2*length(levels(data2$arm)))), 
                                   pval.size=2.75, 
                                   #pval.method.size=2.75, 
                                   risk.table=risk.table, 
                                   risk.table.col="strata", 
                                   #legend="none", 
                                   legend.title = surv.legend.title,
                                   legend.labs = surv.legend.labs, 
                                   risk.table.title = risk.table.title, 
                                   risk.table.height = risk.table.height, 
                                   surv.plot.height = surv.plot.height,
                                   fontsize = 4,
                                   font.x = 10, 
                                   font.y = 10,
                                   font.tickslab = 10,
                                   font.legend = legend.font,
                                   ggtheme = theme,
                                   tables.theme = survminer::theme_cleantable(), 
                                   ...)
  ####Add in stats table
  if(table.include){
    km_plot$plot <- km_plot$plot + 
      ggplot2::annotation_custom(
        gridExtra::tableGrob(
          sum.table, 
          rows=table.legend.labs, 
          cols=c("Events/N", "Median (95% CI)", "HR (95% CI)"), 
          theme=gridExtra::ttheme_minimal(
            base_size=table.font.size, 
            rowhead=list(
              fg_params=list(col=c("black",color), fontface="bold")))), 
        xmin = if(table.loc == "topright") table.x.min else 0, 
        xmax = if(table.loc == "topright") x.max else table.x.max, 
        ymin = if(table.loc == "topright") .65 else 0,
        ymax = if(table.loc == "topright") 1 else .35
      )
  }
  ####edit PAR table
  if(risk.table)
  {
    km_plot$table <- km_plot$table +
      theme(
        axis.ticks.y=element_blank(), 
        axis.text.y=element_text(size = 9, face="bold"), 
        plot.title = element_text(size = 10, face = "bold"))
  }
  ##Print KM Curve
  return(km_plot)
}

