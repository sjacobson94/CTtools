#' Create monthly accrual graph
#' 
#' This is a function to create accrual graphs for trials. Similar output to the \%accgph() macro.
#' 
#' @param protnum The study number as a character.
#' @param data The dataset with the patients you want to use. Must contain patient IDs and patient "on date"
#' @param id The patient ID variable in quotes used in the dataset 
#' @param date.on The variable that contains the patient "on date". Entered as character.
#' @param open.dt The study open date. Alternatively, the lower limit for date you want to use for the start of the graph. Entered in "mm/dd/yyyy" format.
#' @param froz.dt The upper limit for the patient "on dates". Entered in "mm/dd/yyyy" format.
#' @param expmth The number of expected patients accrued per month. Available in protocol.
#' @param y.breaks The distance between each y-axis break. Default is y.lim/10 rounded to the nearest 10.
#' @param y.lim Limit for y-axis. Default is the number of rows in the given dataset.
#' @param exp.pts The number of patient expected to be accrued for the study
#' @return This function creates a month by month accrual graph.
#' @examples 
#' 
#' #With data loaded for A041501
#' 
#' crtlibn(d = 'A041501')
#' 
#' accrual <- accgph(protnum = "A041501", 
#'     data = crse, id = "dcntr_id", 
#'     date.on = "date_on",
#'     open.dt = "06/15/2018",
#'     froz.dt = format(Sys.Date(), "%m/%d/%Y"), 
#'     expmth = 5, exp.pts = 324) 
#' accrual
#' 
#' @export
#' @author Sawyer Jacobson
#' @author Adam Pettinger
#' 
# Last Edited: 1/22/2020
#               Fixed issue with February not appearing if the trial opened after the 28th of the month
# Edited: 5/27/2020
#           Updated the y-axis feature to allow for the limit to be specified and the distance between breaks to be specified as well.
#' 
#' 


# Creating the function
accgph <- function(data, protnum, id = "dcntr_id", date.on, open.dt, froz.dt, expmth = 0, y.breaks = round(y.lim/10, -1), exp.pts, y.lim = nrow(one)){
  # Bringing in the dataset and setting the on study date
  one <- as.data.frame(data)
  one$date.on <- one[, date.on]
  if(class(date.on) != "Date"){
    one$date.on <- as.Date(one$date.on)
  }
  
  temp <- one %>%
    filter(lubridate::mdy(open.dt) <= date.on & date.on <=lubridate::mdy(froz.dt)) %>%
    mutate(
      prot = protnum,
      accrual = 0,
      num = 1,
      date_open = lubridate::mdy(open.dt),
      date_close = lubridate::mdy(froz.dt),
      year_mon = zoo::as.yearmon(date.on)
    ) %>% 
    arrange(date.on) %>%
    group_by(year_mon) %>%
    dplyr::count()
  
  year_mon <- data.frame(date = lubridate::floor_date(seq(lubridate::floor_date(lubridate::mdy(open.dt), unit = "month"), lubridate::mdy(froz.dt), 
                                   by = "1 month"), unit = "month")) %>%
    mutate(year_mon = zoo::as.yearmon(date))
  min_date <- year_mon[1, "date"]
  max_date <- year_mon[nrow(year_mon), "date"]
  final <- year_mon %>%
    left_join(temp, by = "year_mon") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  final$exp_totals <- seq(0, nrow(final)*expmth - expmth, by = expmth)
  
  final$totals <- purrr::accumulate(final$n, `+`)
  year_mon_row <- nrow(year_mon)
  
  if(year_mon_row < 24){
    n_month = "1"
  }else if(dplyr::between(year_mon_row, 24, 36)){
    n_month = "2"
  }else if(year_mon_row >36){
    n_month = "3"
  }
  
  p <- ggplot(data = final) + 
    geom_point(aes(date, totals, color = "actual"), size = 2) + 
    geom_text(aes(date, totals, label=n), hjust=0, vjust=1.35) + 
    geom_line(aes(date, totals, color = "actual")) + 
    geom_line(aes(date, exp_totals, color = "expected"), lty = 2, size = 1) + 
    scale_x_date(date_labels = "%b %y", date_breaks = paste0(n_month," month"), limits = c(min_date, max_date)) + 
    scale_y_continuous(n.breaks = y.breaks, 
                       breaks = seq(0, y.lim, by = y.breaks),
                       limits = c(0, y.lim)) +
    # scale_y_continuous(breaks=pretty(seq(0, nrow(one)+5, by=1), n = y.breaks),
    #                    labels=pretty(seq(0, nrow(one)+5, by=1), n = y.breaks), limits=c(0, nrow(one)+.05*nrow(one))) +
    scale_linetype_manual(values = c("actual" = 1, "expected" = 2)) +
    scale_colour_manual(name = NULL, 
                                 values =c('actual'='black','expected'='black'), 
                                 labels = c('Actual Accrual', paste0('Expected Accrual ', exp.pts, ' Patients'))) +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.line.x.top = element_line(colour = 'black', size = 1), 
          axis.line.y.right = element_line(colour = 'black', size = 1),
          axis.line.x.bottom = element_line(colour = 'black', size = 1), 
          axis.line.y.left = element_line(colour = 'black', size = 1),
          axis.text.x=element_text(angle=45, hjust=1, vjust = 1, face = "bold", size = 8),
          axis.text.y=element_text(face = "bold"),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = c(.23, .95),
          axis.ticks.length = grid::unit(.2, "cm"),
          legend.key.width = grid::unit(1, "cm")) +
    guides(colour = guide_legend(override.aes = list(lty = c(1, 2), size = c(.75, .75), shape = c(19, NA)))) + 
    ylab("Total Number of Patients") +
    xlab("Month") + 
    ggtitle(paste0(protnum, " Cumulative Accrual: Actual Versus Expected\nNumbers on the Accrual Line are Monthly Accrual"))
 out <- list(graph = p, data = final)
 class(out) <- "accgph"
 out
}

#' @export
print.accgph <- function(x, ...)
{
  print(x$graph)
  invisible(x)
}

