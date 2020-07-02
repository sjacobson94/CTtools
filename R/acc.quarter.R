#' Quarterly Accrual
#'
#' This function creates the quarterly accrual table for Alliance DSMB reports.
#'
#' @param data The data containing one row for each patient with on study date.
#' @param id Patient identifier
#' @param date.on Variable for on study date
#' @param open.dt Study open date
#' @param froz.dt Frozen date
#' @param expmth Expected patients per month for an amount of quarters set in expmth.quarters. Can be an integer 
#' value or vector of integers. Vector must be same length as expmth.quarters.
#' @param expmth.quarters Number of quarters for each expmth value. Can be integer or vector of values. 
#' Each vector represents the number of quarters for each expmth value. Eg. expmth = c(5,10,15) and 
#' expmth.quarters = c(3,3,99) would expect 5 patients accrued per month for the first 3 quarters, 
#' 10 patients accrued per month for quarters 4-6, and 15 patients per month thereafter.
#' @param reg.stat Registration status for accrual table. "reg" for registered, "pre" for pre-registered, 
#' and "rand" for randomized.
#' @param stars Option to add asterisks to the first 2 quarter start times and the last quarter stop time for 
#' comments or footnotes. Set to FALSE if no asterisks are desired.
#' @param date.fmt Format for dates in the table. Default "\%b \%d, \%Y".
#' @return This function creates a quarterly accrual table.
#' @examples
#' 
#' crtlibn(d = "A041501")
#' 
#' x <- acc.quarter(casecrse, open.dt = "6/15/2018", froz.dt = "3/20/2020", expmth = 8)
#' x
#' 
#' y <- acc.quarter(casecrse, open.dt = "6/15/2018", froz.dt = "3/20/2020", expmth = c(5,10,15),
#'  expmth.quarters = c(2,2,99))
#' y
#'
#' @export
#' @author Sawyer Jacobson

acc.quarter <- function(data, id = "dcntr_id", date.on = "date_on", open.dt, froz.dt, expmth = c(0), 
                        expmth.quarters = c(99), reg.stat = "reg", stars = TRUE, date.fmt = "%b %d, %Y")
{
  one <- as.data.frame(data)
  reg.stat <- tolower(reg.stat)
  if(reg.stat == "reg")
  {
    reg <- "Registered"
  } else if(reg.stat == "pre")
  {
    reg <- "Pre-Registered"
  } else if(reg.stat == "rand")
  {
    reg <- "Randomized" 
  }
  
  if(length(expmth) != length(expmth.quarters)) stop("Lengths of expected accrual and accrual quarters must be the same.")
  
  one$date.on <- one[, date.on]
  if(class(date.on) != "Date"){
    one$date.on <- as.Date(one$date.on)
  }
  
  open_dt <- lubridate::mdy(open.dt)
  froz.dt <- lubridate::mdy(froz.dt)
  
  temp <- one %>%
    filter(lubridate::mdy(open.dt) <= date.on & date.on <= froz.dt) %>%
    select(dcntr_id, date.on) %>%
    arrange(date.on) 
  
  quarters.begin <- numeric(0)
  class(quarters.begin) <- "Date"
  quarters.end <- numeric(0)
  class(quarters.end) <- "Date"
  i <- 2
  #start.q <- lubridate::floor_date(open_dt, unit = "month") 
  quarters.begin[1] <- open_dt
  while ((open_dt %m+% months(3)) < froz.dt){
    quarters.begin[i] <- quarters.begin[i - 1] %m+% months(3)
    i <- i + 1
    open_dt <- open_dt %m+% months(3)
  }
  quarters.end <- quarters.begin %m+% months(3) %m-% lubridate::days(1)
  #quarters.begin[1] <- lubridate::mdy(open.dt)
  quarters.end[i - 1] <- froz.dt
  quarters <- data.frame(quarter = 1:length(quarters.begin), start = quarters.begin, end = quarters.end, year = lubridate::year(quarters.begin))
  
  counting <- function(x, start, end)
  {
    ret <- numeric(length(start))
    for(i in 1:length(start)){
      x.start <- start[i]
      x.end <- end[i]
      ret[i] <- sum(as.numeric(between(x, x.start, x.end)))
    }
    ret
  }
  
  quarters$total <- counting(temp$date.on, quarters$start, quarters$end)
  expmth.acc <- rep(expmth, expmth.quarters)
  quarter_acc <- data.frame(quarter = 1:length(expmth.acc), expmth = expmth.acc)
  
  quarters <- quarters %>%
    left_join(quarter_acc, by = "quarter")

  quarters <- quarters %>%
    mutate(
      exp = 3*expmth, 
      start_f = format(start, date.fmt), 
      end_f = format(end, date.fmt))
  quarters[nrow(quarters), "exp"] <- ceiling(expmth[length(expmth)] * as.numeric(quarters[nrow(quarters), "end"] - quarters[nrow(quarters), "start"])/30)
  
  if(stars)
  {
    quarters[1, "exp"] <- floor(expmth[1] * as.numeric(quarters[1, "end"] - quarters[1, "start"])/30)
    quarters[1, "start_f"] <- paste0(quarters[1, "start_f"], "*")
    quarters[2, "start_f"] <- paste0(quarters[2, "start_f"], "**")
    quarters[nrow(quarters), "end_f"] <- paste0(quarters[nrow(quarters), "end_f"], "***") 
  }
  quarters <- quarters %>%
    rowwise() %>%
    mutate(
      percent = paste0(round(total/max(exp, .1) * 100, digits = 1), "%"),
      qrange = paste0(start_f, " - ", end_f)
      )  %>%
    select(quarter, qrange, total, exp, percent)
  final <- quarters %>%
    flextable() %>%
    delete_part(part = "header") %>%
    add_header_row(values = c("Quarter", "Time Period", paste0('Actual Number ', reg), paste0('Expected Number ', reg), 
                              paste0("% Expected Number of ", reg)), colwidths = c(1,1,1,1,1)) %>%
    theme_box() %>%
    width(j = 1:5, width = c(.66, 2.38, .88, .91, 1)) %>%
    align(align = "center", part = "all") %>%
    fontsize(size = 10, part = "all") %>% 
    font(fontname = "Arial", part = "all")
  
  out = list(table = final, data = quarters)
  class(out) <- "acc.quarter"
  out
}

#' @export
print.acc.quarter <- function(x, ...)
{
  print(x$table)
  invisible(x)
}

