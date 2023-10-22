#' @title Integrated ARMA Time Series Modelling
#'
#' @description
#' A collection and description of simple to use functions to model univariate
#' autoregressive moving average time series processes, including time series simulation,
#' parameter estimation, diagnostic analysis of the fit, and predictions of future values.
#'
#' The Functions are:
#'
#' `armaSim`	Simulates an artificial ARMA time series process.
#' `armaFit`	Fits the parameters of an ARMA time series process.
#' `print`	  Print Method.
#' `plot`   	Plot Method.
#' `summary`	Summary Method.
#' `predict`	Forecasts and optionally plots an ARMA process.
#' `fitted` 	Method, returns fitted values.
#' `coef`   	Method, returns coefficients.
#' `coefficients` 	Method, returns coefficients.
#' `residuals`Method, returns residuals.`
#'
#' @usage
#' armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
#'         innov = NULL, n.start = 100, start.innov = NULL,
#'         rand.gen = rnorm, rseed = NULL, addControl = FALSE, ...)
#' armaFit(formula, data, method = c("mle", "ols"), include.mean = TRUE,
#'         fixed = NULL, title = NULL, description = NULL, ...)
#'
#' # S4 method for signature 'fARMA'
#' show(object)
#'
#' # S3 method for class 'fARMA'
#' fitted(object, ...)
#'
#' # S3 method for class 'fARMA'
#' coef(object, ...)
#'
#' # S3 method for class 'fARMA'
#' residuals(object, ...)
#'
#' @param addControl [armaSim] -
#'   a logical value. Should control parameters added to the returned series
#'   as a control attribute?
#' @param data an optional timeSeries or data frame object containing
#'   the variables in the model. If not found in data, the variables are taken from
#'   `environment(formula)`, typically the environment from which `armaFit` is
#'   called. If data is an univariate series, then the series is converted into a
#'   numeric vector and the name of the response in the formula will be neglected.
#' @param description a character string which allows for a brief description.
#' @param doplot [armaRoots] -
#'   a logical. Should a plot be displayed?
#'   [predict] -
#'   is used by the predict and `summary` methods. By default, this value is set
#'   to `TRUE` and thus the function calls generate beside written also
#'   graphical printout. Additional arguments required by underlying
#'   functions have to be passed through the `dots` argument.
#' @param fixed [armaFit] -
#'   is an optional numeric vector of the same length as the total number of
#'   parameters. If supplied, only `NA` entries in `fixed` will be varied. In this
#'   way subset ARMA processes can be modeled. ARIMA modelling
#'   supports this option. Thus for estimating parameters of subset ARMA
#'   and AR models the most easiest way is to specify them by the formulas
#'   `x~ARIMA(p, 0, q)` and `x~ARIMA(p, 0, 0)`, respectively.
#' @param formula [armaFit] -
#'   a formula specifying the general structure of the ARMA form. Can have
#'   one of the forms `x ~ ar(q)`, `x ~ ma(p)`, `x ~ arma(p, q)`,
#'   `x ~ arima(p, d, q)`, or `x ~ arfima(p, q)`.
#'   `x` is the response variable optionally to appear in the formula expression.
#'   In the first case R's function `ar` from the ts package will be used
#'   to estimate the parameters, in the second case R's function `arma` from
#'   the `tseries` package will be used, in the third case R's function
#'   `arima` from the ts package will be used, and in the last case R's function
#'   `fracdiff` from the `fracdiff` package will be used. The state space modelling
#'   based `arima` function allows also to fit ARMA models using `arima(p, d=0, q)`,
#'   and AR models using `arima(q, d=0, q=0)`, or pure MA models using `arima(q=0, d=0, p)`
#'   (Exogenous variables are also allowed and can be passed through the `...` argument.).
#' @param gof.lag [print][plot][summary][predict] -
#'   the maximum number of lags for a goodness-of-fit test.
#' @param include.mean [armaFit] -
#'   Should the ARIMA model include a mean term? The default is `TRUE`,
#'   note that for differenced series a mean would not affect the fit nor predictions.
#' @param innov [armaSim] -
#'   is a univariate time series or vector of innovations to produce the series.
#'   If not provided, `innov` will be generated using the random number
#'   generator specified by `rand.gen`. Missing values are not allowed.
#'   By default the normal random number generator will be used.
#' @param method [armaFit] -
#'   a character string denoting the method used to fit the model. The
#'   default method for all models is the log-likelihood parameter estimation
#'   approach, `method="mle"`. In the case of an AR model the parameter
#'   estimation can also be done by ordinary least square estimation, `"ols"`.
#' @param model [armaSim] -
#'   a list with one (AR), two (ARMA) or three (ARIMA, FRACDIFF) elements.
#'   `ar` is a numeric vector giving the AR coefficients, `d` is an integer value
#'   giving the degree of differencing, and `ma` is a numeric vector giving the
#'   MA coefficients. Thus the order of the time series process is
#'   (F)ARIMA(p, d, q) with `p=length(ar)` and `q=length(ma)`.
#'   `d` is a positive integer for ARIMA models and a numeric value
#'   for FRACDIFF models. By default an ARIMA(2, 0, 1) model with coefficients
#'   `ar=c(0.5, -0.5)` and `ma=0.1` will be generated.
#' @param n [armaSim] -
#'   an integer value setting the length of the series to be simulated
#'   (optional if `innov` is provided). The default value is 100.
#' @param n.ahead [print][plot][summary][predict] -
#'   are presetted arguments for the `predict` method. `n.ahead` determines
#'   how far ahead forecasts should be evaluated together with errors on the
#'   confidence intervals given by the argument `conf`. If a forecast plot is
#'   desired, which is the default and expressed by `doplot=TRUE`, then
#'   `n.back` sets the number of time steps back displayed in the graph.
#' @param n.back [print][plot][summary][predict] -
#'   are presetted arguments for the `predict` method. `n.ahead` determines
#'   how far ahead forecasts should be evaluated together with errors on the
#'   confidence intervals given by the argument `conf`. If a forecast plot is
#'   desired, which is the default and expressed by `doplot=TRUE`, then
#'   `n.back` sets the number of time steps back displayed in the graph.
#' @param conf [print][plot][summary][predict] -
#'   are presetted arguments for the `predict` method. `n.ahead` determines
#'   how far ahead forecasts should be evaluated together with errors on the
#'   confidence intervals given by the argument `conf`. If a forecast plot is
#'   desired, which is the default and expressed by `doplot=TRUE`, then
#'   `n.back` sets the number of time steps back displayed in the graph.
#' @param n.start [armaSim] -
#'   gives the number of start-up values discarded when simulating non-stationary
#'   models. The start-up innovations will be generated by `rand.gen`
#'   if `start.innov` is not provided.
#' @param object [summary][predict] -
#'   is an object of class `fARMA` returned by the fitting function `armaFit` and
#'   serves as input for the `summary`, and predict methods. Some methods
#'   allow for additional arguments.
#' @param rand.gen [armaSim] -
#'   is the function which is called to generate the innovations. Usually,
#'   `rand.gen` will be a random number generator. Additional arguments
#'   required by the random number generator `rand.gen`, usually the
#'   location, scale and/or shape parameter of the underlying distribution
#'   function, have to be passed through the `dots` argument.
#' @param rseed [armaSim] -
#'   the random number seed, by default `NULL`. If this argument is set to an
#'   integer value, then the function `set.seed(rseed)` will be called.
#' @param start.innov [armaSim] -
#'   is a univariate time series or vector of innovations to be used as start up
#'   values. Missing values are not allowed.
#' @param title a character string which allows for a project title.
#' @param which [plot][summary] -
#'   if `which` is set to `"ask"` the function will interactively ask which plot
#'   should be displayed. This is the default value for the `plot` method. If
#'   `which="all"` is specified all plots will be displayed. This is the default
#'   setting for the `summary` method. On the other hand, if a vector of
#'   logicals is specified, then those plots will be displayed for which the
#'   elements of the vector are set to `TRUE`.
#' @param x [print][plot] -
#'   is an object of class `fARMA` returned by the fitting function `armaFit` and
#'   serves as input for the `predict`, `print`, `print.summary`, and `plot` methods.
#'   Some methods allow for additional arguments.
#' @param ... additional arguments to be passed to the output timeSeries.
#'
#' @details
#' AR - Auto-Regressive Modelling:
#'
#' The argument x~ar(p) calls the underlying functions ar.mle or ar.ols depending on
#' the method's choice. For definiteness, the AR models are defined through
#'
#' $$x_t-\mu=a_1(x_{t-1}-\mu)+\cdots+a_p(x_{t-p}-\mu)+e_t$$
#'
#' Order selection can be achieved through the comparison of AIC values for different
#' model specifications. However this may be problematic, as of the methods here only
#' `ar.mle` performs true maximum likelihood estimation. The AIC is computed as if the
#' variance estimate were the MLE, omitting the determinant term from the likelihood.
#' Note that this is not the same as the Gaussian likelihood evaluated at the estimated
#' parameter values. With `method="yw"` the variance matrix of the innovations is computed
#' from the fitted coefficients and the autocovariance of `x`. Burg's method allows for two
#' alternatives `method="burg1"` or `method="burg2"` to estimate the innovations variance and
#' hence AIC. Method 1 is to use the update given by the Levinson-Durbin recursion
#' (Brockwell and Davis, 1991), and follows S-PLUS. Method 2 is the mean of the sum of
#' squares of the forward and backward prediction errors (as in Brockwell and Davis, 1996).
#' Percival and Walden (1998) discuss both.
#'
#' MA - Moving-Average Modelling:
#'
#' The argument `x~ma(q)` maps the call to the argument `x ~ arima(0, 0, q)`.
#'
#' ARMA - Auto-Regressive Moving-Average Modelling:
#'
#' The argument `x~arma(p,q)` maps the call to the argument `x~arima(p, 0, q)`.
#'
#' ARIMA - Integrated ARMA Modelling:
#'
#' The argument `x~arima()` calls the underlying function `arima` from R's `ts` package.
#' For definiteness, the AR models are defined through
#'
#' $$x_t=a_1x_{t-1}+\cdots+a_px_{t-p}+e_t+b_1e_{t-1}+\cdots+b_qe_{t-q}$$
#'
#' and so the MA coefficients differ in sign from those of S-PLUS. Further,
#' if `include.mean` is `TRUE`, this formula applies to $x-m$ rather than $x$. For
#' ARIMA models with differencing, the differenced series follows a zero-mean ARMA model.
#' The variance matrix of the estimates is found from the Hessian of the log-likelihood,
#' and so may only be a rough guide.
#' Optimization is done by `optim`. It will work best if the columns in `xreg` are roughly
#' scaled to zero mean and unit variance, but does attempt to estimate suitable scalings.
#' The exact likelihood is computed via a state-space representation of the ARIMA
#' process, and the innovations and their variance found by a Kalman filter.
#' The initialization of the differenced ARMA process uses stationarity. For a differenced
#' process the non-stationary components are given a diffuse prior (controlled by `kappa`).
#' Observations which are still controlled by the diffuse prior (determined by having
#' a Kalman gain of at least 1e4) are excluded from the likelihood calculations.
#' (This gives comparable results to `arima0` in the absence of missing values, when
#' the observations excluded are precisely those dropped by the differencing.)
#' Missing values are allowed, and are handled exactly in method `"ML"`.
#' If `transform.pars` is true, the optimization is done using an alternative parametrization
#' which is a variation on that suggested by Jones (1980) and ensures that the model is
#' stationary. For an AR(p) model the parametrization is via the inverse tanh of the partial
#' autocorrelations: the same procedure is applied (separately) to the AR and seasonal AR
#' terms. The MA terms are not constrained to be invertible during optimization, but they
#' will be converted to invertible form after optimization if `transform.pars` is true.
#' Conditional sum-of-squares is provided mainly for expositional purposes. This computes
#' the sum of squares of the fitted innovations from observation `n.cond` on, (where `n.cond`
#' is at least the maximum lag of an AR term), treating all earlier innovations to be zero.
#' Argument `n.cond` can be used to allow comparability between different fits. The “part
#' log-likelihood” is the first term, half the log of the estimated mean square. Missing
#' values are allowed, but will cause many of the innovations to be missing.
#' When regressors are specified, they are orthogonalized prior to fitting unless any
#' of the coefficients is fixed. It can be helpful to roughly scale the regressors
#' to zero mean and unit variance.
#' Note from `arima`: The functions parse their arguments to the original time series
#' functions available in R's time series library `ts`.
#' The results are likely to be different from S-PLUS's `arima.mle`, which computes a
#' conditional likelihood and does not include a mean in the model. Further, the
#' convention used by `arima.mle` reverses the signs of the MA coefficients.
#'
#' ARFIMA/FRACDIFF Modelling:
#'
#' The argument `x~arfima()` calls the underlying functions from R's `fracdiff` package.
#' The estimator calculates the maximum likelihood estimators of the parameters of a
#' fractionally-differenced ARIMA (p,d,q) model, together (if possible) with their estimated
#' covariance and correlation matrices and standard errors, as well as the value of the
#' maximized likelihood. The likelihood is approximated using the fast and accurate
#' method of Haslett and Raftery (1989). Note, the number of AR and MA coefficients
#' should not be too large (say < 10) to avoid degeneracy in the model.
#' The optimization is carried out in two levels: an outer univariate unimodal optimization
#' in d over the interval [0,.5], and an inner nonlinear least-squares optimization in the AR
#' and MA parameters to minimize white noise variance.
#'
#' @return `armaFit` returns an S4 object of class `"fARMA"`,
#'   with the following slots:
#'
#'   `call`	the matched function call.
#'   `data` the input data in form of a data.frame.
#'   `description` allows for a brief project description.
#'   `fit` the results as a list returned from the underlying time series model function.
#'   `method` the selected time series model naming the applied method.
#'   `formula` the formula expression describing the model.
#'   `parameters` named parameters or coefficients of the fitted model.
#'   `title` a title string.
#'
#' @export
#' @examples
#'
#' # armaSim - Simulation
#' x = armaSim(model = list(ar = c(0.5, -0.5), ma = 0.1), n = 1000)
#'
#' # armaFit - Estimate the Parameters
#' fit = armaFit(~ arma(2, 1), data = x)
#' print(fit)
#'
#' # summary - Diagnostic Analysis
#' par(mfrow = c(2, 2), cex = 0.7)
#' summary(fit, which =  "all")
#'
#' # predict - Forecast 5 Steps Ahead
#' par(mfrow = c(1, 1))
#' predict(fit, 5)
#'
#' # armaFit - Alternative Calls
#' TS = MSFT
#' armaFit(formula = diff(log(Close)) ~ ar(5), data = TS)
#' armaFit(Close ~ ar(5), data = returns(TS, digits = 12))
#' TS.RET = returns(TS, digits = 12)
#' armaFit(Close ~ ar(5), TS.RET)
#' armaFit(Close ~ ar(5), as.data.frame(TS.RET))
#' armaFit(~ ar(5), as.vector(TS.RET[, "Close"]))
#' armaFit(~ ar(5), as.ts(TS.RET)[, "Close"])
#' attach(TS.RET)
#' armaFit(Close ~ ar(5))
#' detach(TS.RET)


setClass("fARMA",
         representation(
           call = "call",
           formula = "formula",
           method = "character",
           parameter = "list",
           data = "list",
           fit = "list",
           residuals = "list",
           fitted = "list",
           title = "character",
           description = "character"
         )
)


# ------------------------------------------------------------------------------
#' @export

setMethod("show", "fARMA",
          function(object)
          {   # A function implemented by Diethelm Wuertz

            # Description:
            #   Old S3 print method for a fitted ARMA timeSeries object

            # FUNCTION:

            # Title:
            cat("\nTitle:\n ")
            cat(object@title, "\n")

            # Call:
            cat("\nCall:\n ")
            cat(paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n", sep = "")

            # Model:
            cat("\nModel:\n ", object@fit$tstitle, "\n", sep = "")

            # Coefficients:
            cat("\nCoefficient(s):\n")
            digits = max(4, getOption("digits") - 4)
            print.default(format(object@fit$coef, digits = digits), print.gap = 2,
                          quote = FALSE)

            # Description:
            cat("\nDescription:\n ")
            cat(object@description, "\n\n")

            # Return Value:
            invisible()
          })


# ------------------------------------------------------------------------------
#' @export

summary.fARMA <-
  function (object, doplot = TRUE, which = "all", ...)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Analyzes a Fitted ARMA timeSeries Object

    # FUNCTION:

    # Initialize:
    if (object@fit$tsmodel == "arfima" & doplot) {
      warning(" Plot Method for arfima Models not yet implemented")
      doplot = FALSE
    }
    ans = NULL

    # Fit Call and Model:
    x <- object
    object = x@fit
    ans$call = object$call
    ans$tsmodel = object$tstitle

    # Calculate Residuals and Variance:
    # ans$residuals = na.remove(object$residuals)
    ans$residuals = as.vector(na.omit(object$residuals))
    if (length(ans$residuals) == 0) {
      ans$var = 0 }
    if (length(ans$residuals) > 0) {
      ans$var = var(ans$residuals) }
    ans$sigma2 = object$sigma2

    # Generate Coef Matrix:
    tval = object$coef/object$se.coef
    prob = 2 * (1 - pnorm(abs(tval)))
    ans$coefmat = cbind(object$coef, object$se.coef, tval, prob)
    dimnames(ans$coefmat) = list(names(object$coef),
                                 c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))

    # More Parameters: aic, etc ...
    if (object$tsmodel == "ar") {
      ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used *
                   log(ans$var) + 2 * length(object$coef)) }
    if (object$tsmodel == "arma") {
      ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used *
                   log(ans$var) + 2 * length(object$coef))
      ans$css = object$css }
    if (object$tsmodel == "arima") {
      ans$aic = object$aic
      ans$loglik = object$loglik }
    if (object$tsmodel == "fracdiff") {
      doplot = FALSE }

    # Print Part:

    # Title:
    cat("\nTitle:\n ")
    cat(x@title, "\n")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    # Model:
    cat("\nModel:\n ", object$tstitle, "\n", sep = "")

    # Coefficients:
    cat("\nCoefficient(s):\n")
    digits = max(4, getOption("digits") - 4)
    print.default(format(object$coef, digits = digits), print.gap = 2,
                  quote = FALSE)

    # Residuals:
    digits = max(4, getOption("digits") - 4)
    if (length(object$residuals) > 2) {
      cat("\nResiduals:\n")
      rq = structure(quantile(ans$residuals),
                     names = c("Min", "1Q", "Median", "3Q", "Max"))
      print(rq, digits = digits)
      # Moments:
      cat("\nMoments: \n")
      skewness = sum((ans$residuals - mean(ans$residuals))^3 /
                       sqrt(var(ans$residuals))^3)/length(ans$residuals)
      kurtosis = sum((ans$residuals - mean(ans$residuals))^4 /
                       var(ans$residuals)^2)/length(ans$residuals) - 3
      stats = structure(c(skewness, kurtosis),
                        names = c("Skewness", "Kurtosis"))
      print(stats, digits = digits) }

    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    signif.stars = getOption("show.signif.stars")
    printCoefmat(ans$coefmat, digits = digits,
                 signif.stars = signif.stars, ...)

    # Fit:
    cat("\n")
    if (x@fit$tsmodel == "ar") {
      cat("sigma^2 estimated as:       ",
          format(object$var, digits = digits), "\n")
      cat("AIC Criterion:              ",
          format(round(object$aic, 2)), "\n") }
    if (x@fit$tsmodel == "arma") {
      cat("sigma^2 estimated as:       ",
          format(object$sigma2, digits = digits), "\n")
      cat("Conditional Sum-of-Squares: ",
          format(round(object$css, digits=2)), "\n")
      ## cat("AIC Criterion:              ",
      ##    format(round(object$aic, digits=2)), "\n")
    }
    if (x@fit$tsmodel == "arima") {
      cm = object$call$method
      if (is.null(cm) || cm != "CSS")
        cat(
          "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
          "\nlog likelihood:       ", format(round(object$loglik, 2)),
          "\nAIC Criterion:        ", format(round(object$aic, 2)),
          "\n", sep = "")
      else
        cat(
          "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
          "\npart log likelihood:  ", format(round(object$loglik,2)),
          "\n", sep = "") }

    # Doplot:
    if (doplot) plot.fARMA(x, which = which, ...)

    # Description:
    cat("\nDescription:\n ")
    cat(x@description, "\n\n")

    # Return Value:
    invisible()
  }


# ------------------------------------------------------------------------------
#' @export

plot.fARMA <-
  function(x, which = "ask", gof.lag = 10, ...)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fARMA'

    # FUNCTION:

    # Check:
    if (x@fit$tsmodel == "arfima") {
      warning(" Plot method for ARFIMA models not yet implemented")
      return()
    }

    # Store Lag:
    x@fit$gof.lag = gof.lag

    # Plot:
    .interactiveArmaPlot(
      x,
      choices = c(
        "Standardized Residuals",
        "ACF of Residuals",
        "QQ Plot of Residuals",
        "Ljung-Box p Values"),
      plotFUN = c(
        ".plot.arma.1",  ".plot.arma.2",  ".plot.arma.3", ".plot.arma.4"),
      which = which)

    # Return Value:
    invisible(x)
  }


# ------------------------------------------------------------------------------
#' @export

.plot.arma.1 <-
  function(x, ...)
  {
    # 1. Standardized Residuals Plot:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    stdres = rs/sqrt(object$sigma2)
    plot(stdres, type = "h",
         main = "Standardized Residuals",
         ylab = "Residuals", col = "steelblue", ...)
    grid()
    abline(h = 0, col = "grey")
  }


# ------------------------------------------------------------------------------
#' @export

.plot.arma.2 <-
  function(x, ...)
  {
    # 2. ACF of Residuals:
    object = x@fit
    acf(object$residuals, plot = TRUE, main = "ACF of Residuals",
        na.action = na.pass, ...)
    grid()
  }


# ------------------------------------------------------------------------------
#' @export

.plot.arma.3 <-
  function(x, ...)
  {
    # 3. QQ Plot of Residuals:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    stdres = rs/sqrt(object$sigma2)
    qqnorm(stdres,
           xlab = "Normal Quantiles",
           ylab = "Residual Quantiles",
           main = "QQ-Plot of Residuals",
           pch = 19, col = "steelblue", ...)
    qqline(stdres, col = "grey")
    grid()
  }


# ------------------------------------------------------------------------------
#' @export

.plot.arma.4 <-
  function(x, ...)
  {
    # 4. Ljung-Box p Values:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    nlag = x@fit$gof.lag
    pval = numeric(nlag)
    for (i in 1:nlag)
      pval[i] = Box.test(rs, i, type = "Ljung-Box")$p.value
    plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1),
         pch = 19, col = "steelblue", main = "Ljung-Box p-values", ...)
    abline(h = 0.05, lty = 2, col = "grey")
    grid()
  }


# ------------------------------------------------------------------------------
#' @export

.interactiveArmaPlot =
  function(x, choices = paste("Plot", 1:19),
           plotFUN = paste("plot.", 1:19, sep = ""), which = "all", ...)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class "template".

    # Arguments:
    #   x - an object to be plotted
    #   choices - the character string for the choice menu
    #   plotFUN - the names of the plot functions
    #   which - plot selection, which graph should be
    #     displayed. If a character string named "ask" the
    #     user is interactively asked which to plot, if
    #     a logical vector of length N, those plots which
    #     are set "TRUE" are displayed, if a character string
    #     named "all" all plots are displayed.

    # Note:
    #   At maximum 19 plots are supported.

    # FUNCTION:

    # Some cecks:
    if (length(choices) != length(plotFUN))
      stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices))
      stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN))
      stop("Arguments which has incorrect length.")
    if (length(choices) > 19)
      stop("Sorry, only 19 plots at max are supported.")

    # Plot:
    if (is.numeric(which)) {
      Which = rep(FALSE, times = length(choices))
      Which[which] = TRUE
      which = Which
    }
    if (which[1] == "all") {
      which = rep(TRUE, times = length(choices))
    }
    if (which[1] == "ask") {
      .multArmaPlot(x, choices, ...)
    } else {
      for ( i in 1:length(which) ) {
        FUN = match.fun(plotFUN[i])
        if (which[i]) FUN(x)
      }
    }

    # Return Value:
    invisible(x)
  }


# ------------------------------------------------------------------------------
#' @export

.multArmaPlot <-
  function (x, choices, ...)
  {
    # Internal "askPlot" Function:

    # Match Functions, up to nine ...
    if (length(plotFUN) < 19) plotFUN =
        c(plotFUN, rep(plotFUN[1], times = 19 - length(plotFUN)))
    plot.1  = match.fun(plotFUN[1]);  plot.2  = match.fun(plotFUN[2])
    plot.3  = match.fun(plotFUN[3]);  plot.4  = match.fun(plotFUN[4])
    plot.5  = match.fun(plotFUN[5]);  plot.6  = match.fun(plotFUN[6])
    plot.7  = match.fun(plotFUN[7]);  plot.8  = match.fun(plotFUN[8])
    plot.9  = match.fun(plotFUN[9]);  plot.10 = match.fun(plotFUN[10])
    plot.11 = match.fun(plotFUN[11]); plot.12 = match.fun(plotFUN[12])
    plot.13 = match.fun(plotFUN[13]); plot.14 = match.fun(plotFUN[14])
    plot.15 = match.fun(plotFUN[15]); plot.16 = match.fun(plotFUN[16])
    plot.17 = match.fun(plotFUN[17]); plot.18 = match.fun(plotFUN[18])
    plot.19 = match.fun(plotFUN[19])
    pick = 1
    while (pick > 0) { pick = menu (
      ### choices = paste("plot:", choices),
      choices = paste(" ", choices),
      title = "\nMake a plot selection (or 0 to exit):")
    # up to 19 plot functions ...
    switch (pick,
            plot.1(x),  plot.2(x),  plot.3(x),  plot.4(x),  plot.5(x),
            plot.6(x),  plot.7(x),  plot.8(x),  plot.9(x),  plot.10(x),
            plot.11(x), plot.12(x), plot.13(x), plot.14(x), plot.15(x),
            plot.16(x), plot.17(x), plot.18(x), plot.19(x))
    }
  }


################################################################################




