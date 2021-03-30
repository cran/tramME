## Test utils: to be sourced from other test files

if (length(strsplit(packageDescription("tramME")$Version, "\\.")[[1]]) > 3) {
  Sys.setenv("NOT_CRAN" = "true")
}


chkeq <- function(x, y, ...) stopifnot(isTRUE(all.equal(x, y, ...)))
chkid <- function(x, y, ...) stopifnot(isTRUE(identical(x, y, ...)))

chkerr <- function(expr, em = NULL) {
  stopifnot(
    tryCatch(expr,
      error = function(e) {
        if (!is.null(em)) return(grepl(em, e))
        else return(TRUE)
      }
      )
  )
}

chkwarn <- function(expr, wm = NULL) {
  stopifnot(
    tryCatch(expr,
      warning = function(w) {
        if (!is.null(wm)) return(grepl(wm, w))
        else return(TRUE)
      }
      )
  )
}
