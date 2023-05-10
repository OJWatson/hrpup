#------------------------------------------------
# length(x) equals n
#' @noRd
assert_length <- function(x, n, message = "%s must be of length %s", name = deparse(substitute(x))) {
  assert_pos_int(n)
  if (length(x) != n[1]) {
    stop(sprintf(message, name, n), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is positive integer (with or without zero allowed)
#' @noRd
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}

#------------------------------------------------
# x is integer
#' @noRd
assert_int <- function(x, message = "%s must be integer valued", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (!isTRUE(all.equal(x, as.integer(x), check.attributes = FALSE))) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is numeric
#' @noRd
assert_numeric <- function(x, message = "%s must be numeric", name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is positive (with or without zero allowed)
#' @noRd
assert_pos <- function(x, zero_allowed = TRUE, message1 = "%s must be greater than or equal to zero", message2 = "%s must be greater than zero", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (zero_allowed) {
    if (!all(x >= 0)) {
      stop(sprintf(message1, name), call. = FALSE)
    }
  } else {
    if (!all(x > 0)) {
      stop(sprintf(message2, name), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# x is single positive (with or without zero allowed)
#' @noRd
assert_single_pos <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_length(x, n = 1, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}

#------------------------------------------------
# x and nms have the same length and the set of names of x is same as names
#' @noRd
assert_same_names <- function(x, nms,
                              message = "%s must have the same names as %s",
                              name_x = deparse(substitute(x)),
                              name_y = deparse(substitute(y))) {

  assert_length(x, length(nms))
  if(!identical(sort(names(x)), sort(nms))) {
    stop(sprintf(message, name_x, name_y), call. = FALSE)
  }

}

#------------------------------------------------
# checks length and orders by nms
#' @noRd
order_names <- function(x, nms,
                       name_x = deparse(substitute(x)),
                       name_y = deparse(substitute(nms))) {
  assert_same_names(x, nms, name_x = name_x, name_y = name_y)
  x[match(nms, names(x))]
}
