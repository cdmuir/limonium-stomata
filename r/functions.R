# Calculations of gsmax follow Sack and Buckley 2016
biophysical_constant <- function(D_wv, v) D_wv / v

morphological_constant <- function(c, h, j) {
  (pi * c ^ 2) / (j ^ 0.5 * (4 * h * j + pi * c))
} 

# Custom theme
theme_cdm <- function(base_size = 10, base_family = "") {
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.2)),
      legend.text = element_text(size = rel(1)),
      legend.title = element_text(size = rel(1)),
      panel.grid.minor = element_blank(), 
      strip.text = element_text(size = rel(1)),
      aspect.ratio = 1,
      complete = TRUE
    )
}

# Get p-value from posterior
get_p <- function(.post) {
  
  m <- median(.post)
  n <- length(.post)
  p <- if (m > 0) {
    1 - length(which(.post > 0)) / n
  } else {
    1 - length(which(.post < 0)) / n
  }
  
  2 * p
    
}

# Functions for manuscript

add_zeroes <- function(s, digits) {
  
  s %<>% ifelse(str_detect(., "\\."), ., str_c(., "."))
  
  n <- s %>%
    str_remove("\\.") %>%
    str_remove("^0+") %>%
    str_remove("0+$") %>%
    nchar()
  
  map2_chr(n, s, ~ str_c(.y, str_c(rep("0", digits - .x), collapse = "")))
  
}  

summarize_trait <- function(x, digits) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  n <- length(na.omit(x))
  se <- sigma / sqrt(n)
  str_c(formatC(mu, digits = digits, format = "fg"), " $\\pm$ ", 
        formatC(se, digits = digits, format = "fg"), " (", n, ")")
}
