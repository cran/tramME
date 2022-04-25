## Test utils: to be sourced from other test files

## Decide based on the package version number whether this is a CRAN version
if (length(strsplit(packageDescription("tramME")$Version, "\\.")[[1]]) > 3) {
  Sys.setenv("NOT_CRAN" = "true")
}

..sumfail <- 0

chkeq <- function(x, y, ..., chkdiff = FALSE) {
  fail <- !isTRUE(all.equal(x, y, ...))
  if (chkdiff) {
    fail <- !fail
    msg <- "The arguments are not different."
  } else msg <- "The arguments are not equal."
  fail_action(fail, match.call(), msg = msg)
}

chkid <- function(x, y, ..., chkdiff = FALSE) {
  fail <- !isTRUE(identical(x, y, ...))
  if (chkdiff) {
    fail <- !fail
    msg <- "The arguments are not different."
  } else msg <- "The arguments are not identical."
  fail_action(fail, match.call(), msg = msg)
}

chkerr <- function(expr, em = NULL) {
  fail <- tryCatch({expr; 1L},
      error = function(e) {
        if (!is.null(em) && !grepl(em, e)) return(2L)
        else return(0L)
      }
      )
  msg <- if (fail < 2L) "No error was raised."
         else "Error message doesn't match."
  fail_action(fail > 0L, match.call(), msg = msg)
}

chkwarn <- function(expr, wm = NULL) {
  fail <- tryCatch({expr; 1L},
      warning = function(w) {
        if (!is.null(wm) && !grepl(wm, w)) return(2L)
        else return(0L)
      })
  msg <- if (fail < 2L) "No warning was raised."
         else "Warning message doesn't match."
  fail_action(fail > 0L, match.call(), msg = msg)
}

fail_action <- function(fail, call,
                        raise_error = identical(Sys.getenv("NOT_CRAN"), "true"),
                        msg = NULL) {
  if (fail) {
    message("\n==== TEST FAILED: ========\n",
            "\t", deparse(call), "\n",
            if (length(msg)) paste0(msg, "\n") ,
            "==========================\n")
    if (raise_error) {
      msg <- if (!is.null(msg)) paste0(": ", msg) else "!"
      stop(paste0("Test failed", msg))
    }
    if (exists("..sumfail")) ..sumfail <<- ..sumfail + 1
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

summarize_tests <- function() {
  if (exists("..sumfail"))
    message("==========================\n",
            "Number of failed tests: ", ..sumfail,"\n",
            "==========================")
}
