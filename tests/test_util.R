## Test utils: to be sourced from other test files

if (length(strsplit(packageDescription("tramME")$Version, "\\.")[[1]]) > 3) {
  Sys.setenv("NOT_CRAN" = "true")
}


chk <- function(x, y) stopifnot(isTRUE(all.equal(x, y)))
chkerr <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
chkwarn <- function(expr, wm) stopifnot(tryCatch(expr, warning = function(w) grepl(wm, w)))
