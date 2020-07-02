

#' Tidy Word Tableby
#'
#' @param x Dataframe of a summary() of tableby object 
#'
#' @return This function outputs a tidy version of a dataframe summary(tableby) object to be output in word.
#' @export
#'
tidy_tableby <- function(x) 
{
  x <- as.data.frame(x)
  header <- str_replace(names(x), "\\(N=", replacement = "\n\\(N=")
  bold_ind <- which(x[, 2] == '')
  padd_ind <- which(x[, 2] != '')
  colnames(x)[1] <- "var"
  x$var <- trimws(str_remove_all(x$var, "[\\[\\]]"))
  out <- x %>%
    flextable() %>%
    bold(j = 1, i = bold_ind, part = "body") %>%
    delete_part(part = "header") %>%
    autofit() %>%
    add_header_row(values = header, colwidths = rep(1, ncol(x))) %>%
    theme_box() %>%
    padding(i =padd_ind, j = 1, padding.left = 20) %>%
    align() %>%
    align(align = "center", part = "header") %>%
    fontsize(size = 10, part = "all") 
  out
}

# tidy_tableby_pdf <- function(x)
# {
#   header <- names(x) #str_replace(names(x), "\\(N=", replacement = "\newline\\(N=")
#   bold_ind <- which(x[, 2] == '')
#   padd_ind <- which(x[, 2] != '')
#   colnames(x)[1] <- "var"
#   x$var <- trimws(str_remove_all(x$var, "[\\[\\]]"))
#   x[bold_ind ,1] <- cell_spec(x[bold_ind, 1], bold = TRUE)
#   #x[padd_ind ,1] <- paste0("\\quad ", x[padd_ind ,1])
#   out <- x %>%
#     kable(format = "latex", row.names = FALSE, align = c("lc"), col.names = header, escape = FALSE) %>%
#     kable_styling(position = "center", full_width = FALSE, protect_latex = TRUE) %>%
#     column_spec(column = 1, border_left = TRUE) %>%
#     column_spec(column = ncol(x), border_right = TRUE) %>%
#     add_indent(padd_ind) 
#   out
# }
