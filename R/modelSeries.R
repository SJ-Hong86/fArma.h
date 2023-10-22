#' @title .modelSeries
#'
#' @description Models a time Series object to use formulas
#' @param formula some formula
#' @param data a timeSeries, a data.frame or a numeric vector.
#' @param fake default is `FALSE`
#' @param lhs default is `FALSE`
#' @return no return value.
#' @details
#' Time Series Modelling
#' Regression Modelling
#' Data Management
#' @export



.modelSeries <-
  function(formula, data, fake = FALSE, lhs = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   formula -
    #   data - a timeSeries, a data.frame or a numeric vector
    #   fake -
    #   lhs -

    # Details:
    #   Time Series Modelling
    #   Regression Modelling
    #   Data Management

    # FUNCTION:

    # If no respnonse is pecified:
    if (length(formula) == 2) {
      formula = as.formula(paste("x", formula[1], formula[2], collapse = ""))
      stopifnot(!missing(data))
    }

    # If data is missing, take the first data set from the search path:
    if (missing(data)) {
      data = eval(parse(text = search()[2]), parent.frame())
    }

    if (is.numeric(data)) {
      data = data.frame(data)
      colnames(data) = all.vars(formula)[1]
      lhs = TRUE
    }

    # If we consider a faked formula:
    if (fake) {
      response = as.character(formula)[2]
      Call = as.character(match.call()[[2]][[3]])
      method = Call[1]
      predictors = Call[2]
      formula = as.formula(paste(response, "~", predictors))
    }

    # If only left hand side is required:
    if (lhs) {
      response = as.character(formula)[2]
      formula = as.formula(paste(response, "~", 1))
    }

    # Create Model Data:
    x = model.frame(formula, data)

    # Convert:
    if (class(data) == "timeSeries") x = timeSeries(x)
    if (fake) attr(x, "control") <- method

    # Return value:
    x

  }


################################################################################
