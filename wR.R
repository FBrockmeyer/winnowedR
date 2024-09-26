# --------------------------- #
####### Sometimes useful ###### 
# --------------------------- #

#' @export
`fb_-+` = \(a, b) c(a - b, a + b)

#' @export 
fb_identifyOS = \() {
  switch(.Platform$OS.type, 
         unix = { print("This OS is either 'Linux' or 'Mac'.") }, 
         windows  = { print("This OS is 'Windows'.") }
  )
} 

#' @export
fb_DTput = \(x) {
  attr(x, which = ".internal.selfref") = NULL 
  dput(x)
}

#' @export
fb_library = \(ls) {
  if(!is.list(ls)) stop("Argument 'ls' is not a list.")
  invisible(lapply(ls, library, character.only = TRUE))
}

#' @export
fb_dfBn_Scc = function(r, rho, n) {
  # implementation of exact density function f(r) 
  # of sample correlation coefficient r, 
  # for data that follow a bivariate normal distribution: 
  nr <- (n - 2L) * gamma(n - 1L) * (1L - rho ^ 2L) ^ ((n - 1L) / 2L) * 
    (1L - r ^ 2L) ^ ((n - 4L) / 2L)
  dr <- sqrt(2L * pi) * gamma(n - .5) * (1L - rho * r) ^ (n - 1.5)
  ghf <- hypergeo::hypergeo(A = .5, B = .5, C = .5 * (2L * n - 1L), 
                            z = .5 * (rho * r + 1L)) 
  ghf * nr / dr
}

#' @export
fb_DataEllipseFdist <- function(x, y, levels, col, legend = TRUE, ...) {
  stopifnot(is.vector(x), is.vector(y), length(x) == length(y))
  if(missing(levels)) levels <- c(.2, .4, .6, .8)
  if(missing(col)) col <-  colorRampPalette(colors = c("green", "yellow"))(length(levels)) 
  message("Notice, function expects two vectors of equal length.\n Assumes data tio be bivariate Normal.")
  stopifnot("Lengths of colour and level(s) vectors differ." =
              length(col) == length(levels))
  # car:::ellipse #
  nseg <- 100L # := number of line segments; determination?
  cm <- cov(cbind(x, y))
  centre <- c(mean(x), mean(y))
  if (max(abs(cm - t(cm))) / max(abs(cm)) > 1e-10) 
    stop("The covariance matrix of `cbind(x,y)` must be a symmetric.")
  angles <- (0L:nseg) * 2L * pi / nseg
  uc <- cbind(cos(angles), sin(angles))
  Q <- chol(cm, pivot = TRUE)
  ord <- order(attr(Q, "pivot"))
  # car:::ellipse #
  lapply(seq_along(levels), \(i) {
    lines(t(centre + 
              sqrt(2L * qf(levels[[i]], 2L, length(x) - 1L)) * 
              t(uc %*% Q[, ord])), 
          type = "l", col = col[[i]], 
          lwd = 2L, lty = "dashed")
  })
  if(legend) {
    legend("bottomright",
           legend = lapply(seq_along(levels), \(i) levels[[i]]),
           col = col, lwd = 2L, lty = "dashed", bty = "n", ...)
  }
  x <- lapply(levels, \(i) {
    r <- sqrt(2L * qf(i, 2L, length(x) - 1L))
    t(centre + r * t(uc %*% Q[, ord])) })
  invisible(x)
}

#' @export
#' returns positions in vec
fb_ordinaryOutlierD <- \(vec, max_diff) {
  if(missing(vec))
    stop("Data vector is missing.")
  if(missing(max_diff))
    stop("Please set a value of maximum difference between consecutive data points.")
  if(length(na.omit(vec)) < 2L) return(numeric())
  for(i in seq_along(vec)[-1]) {
    last_valid <- tail(!is.na(head(vec, i - 1L)), 1L)
    if(abs(vec[[i]] - last_valid) > max_diff) {
      vec[[i]] = NA
    }
  }
  which(is.na(vec))
}

#' @export
fb_rollsumr = \(x, k, ...) {
  stopifnot(is.numeric(x), is.integer(k))
  res = rep(NA, length(x))
  for (i in seq_along(x)) # this should be right-aligned
    if (i > k) res[i] = sum( x[(i-k+1L):i], ...)
  res
}

#' @export
fb_split_df = \(x, f, ...) {
  # preserves order  
  f = factor(f, levels = unique(f))
  split(x, f, ...)
}

#' export 
fb_cbind_list = \(l, fill=NA, names=NULL, transpose=FALSE) {
  # Note. If names are set and transpose=TRUE, names are not used yet. 
  stopifnot("Input is not in list format." = is.list(l))
  ll = lengths(l);  ml = max(ll)
  d = if(length(unique(ll))==1L) list2DF(l) else 
    lapply(l, \(j) c(j, rep(fill,  ml-length(j)))) |> list2DF()
  d = if(!(missing(names))) 
    if(length(names) == ncol(d)) setNames(d, names) else { 
      warning("Length of names does not match number of columns."); d  } else d 
  if(transpose) as.data.frame(t.data.frame(d), row.names=FALSE) else d 
}

#' export 
fb_NApadToSameLength = \(l) {
  stopifnot(is.list(l))
  lapply(l, `length<-`, max(lengths(l))) |> data.frame()
}

#' export 
fb_rowVars = \(x, ...) rowSums((x - rowMeans(x, ...)) ^ 2L, ...) / (nrow(x) - 1L) 

#' export 
fb_rowMins = \(D, k=1L) {
  if(!(is.data.frame(D) | is.matrix(D))) stop("Object is neither data.frame nor matrix.")
  if(!all(sapply(D, is.numeric))) stop("Not all columns are numeric.")
  if(is.data.frame(D)) D = as.matrix(D)
  matrix(D[order(row(D), D)], ncol = ncol(D), byrow = TRUE)[, k]
}

#' export 
quickdf = \(l) { 
  # https://adv-r.hadley.nz/perf-improve.html#as.data.frame
  class(l) = "data.frame"; attr(l, "row.names") = .set_row_names(length(l[[1L]])); l 
}

#' export 
#' multiple numbers from string
fb_numsFstr = \(y) Map(\(x) x[nchar(x) > 0L], strsplit(y, "\\D+"))

#' export
#' Vectorized seq()
fb_Vseq = Vectorize(\(from, to) seq(from, to, by = 1L))

#' export 
#' replace values with first value by group 
#' https://stackoverflow.com/questions/66376333/replace-all-values-with-first-observation-by-group
fb_replaceAllButFirstByGroup = \(df, x, g) {
  df[(i <- !duplicated(df[g])), x][cumsum(i)]
}

#' export 
fb_Mode = \(x) {
  if(length(x) <= 2L) return(x[[1L]])
  if(anyNA(x)) x = x[!is.na(x)]
  ux = unique(x); tab = tabulate(match(x, ux))
  ux[tab == max(tab)]
}

#' export 
fb_consecutive_id = \(x) with(rle(x), rep(seq_along(values), lengths)) 

#' export 
fb_lead = \(x, k, fill = NA) c(x[-seq(k)], rep(fill, k))

#' export 
fb_lag = \(x, k, fill = NA) c(rep(fill, k), head(x, -k))

#' export 
#' Freedman-Diaconis rule
#' https://stats.stackexchange.com/a/862/400850
fb_FD = \(x, f=2L) f * IQR(x) / length(x)^(1L/3L)

#' export 
fb_QuarterDaysOfDate = \(x) {
  ## imports data.table() functions 
  stopifnot(inherits(x, "Date"))
  y = data.table::year(x)
  y = y%%4L == 0L & y%%100L != 0L | y%%400L == 0L 
  q = data.table::quarter(x)
  d = c(90L:92L, 92L) 
  r = d[q]
  r[q==1L & y] = 91L
  r
}

#' export 
fb_rename = \(df, new, old, returndf = TRUE) {
  stopifnot(is.data.frame(df), length(old) == length(new), is.character(new))
  # this check is porbably bad practice:
  if(! inherits(old, c("numeric", "character", "integer"))) 
    stop("Specify 'old' either as character or integer/numeric.") 
  if(is.character(old)) 
    if(all(sapply(old, \(x) x %in% names(df))==FALSE)) 
      stop("Names in 'old' do not match with names of df.") else 
        names(df)[match(old, names(df))] = new 
      else names(df)[old] = new 
      if(returndf) df
}
