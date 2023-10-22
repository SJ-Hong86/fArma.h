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



test.fARMA =
  function()
  {
    # Slot Names:
    Names = getSlots("fARMA")
    target = names(Names)
    current = c(
      "call",
      "formula",
      "method",
      "parameter",
      "data",
      "fit",
      "residuals",
      "fitted",
      "title",
      "description")
    checkIdentical(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.armaSim =
  function()
  {
    # armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
    #   positions = NULL, innov = NULL, n.start = 100, start.innov = NULL,
    #   rand.gen = rnorm, rseed = NULL, addControl = FALSE, ...)

    # ts: Simulate ARMA(2,1):
    # Note, if "positions=NULL",
    #    then a 'ts' object will be returned ...
    ts = armaSim(n = 25)
    class(ts)
    print(ts)
    ts = armaSim(n = 25, addControl = TRUE)
    print(ts)

    # timeSeries: Simulate ARMA(2,1):
    # Note, if "positions" is a timeDate object,
    #   then a 'timeSeries' object will be returned ...
    tS = armaSim(n = 12)
    class(tS)
    print(tS)
    tS = armaSim(n = 12, addControl = TRUE)
    print(tS)

    # ts: Simulate t4-ARMA(2,1):
    ts = armaSim(n = 25, rand.gen = rt, df = 4, rseed = 4711)
    print(ts)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.ar2Fit =
  function()
  {
    # Simulate AR(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5)), n = 1000)
    head(x)

    # method = c("mle", "ols")
    args(ar)

    # AR(2) - OLS Fit:
    # R's internal function ar() will be used ...
    object = armaFit(formula = ~ ar(2), data = x, method = "ols")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # AR(2) - MLE Fit:
    # R's internal function ar() will be used ...
    object = armaFit(formula = ~ ar(2), data = x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # For the expert ...
    # Note, also other methods can be used supported by ar():

    # Yule-Walker Fit:
    object = armaFit(formula = ~ ar(2), data = x, method = "yw")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Burg 1 Fit:
    object = armaFit(formula = ~ ar(2), data = x, method = "burg1")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Burg 2 Fit:
    object = armaFit(formula = x ~ ar(2), data = x, method = "burg2")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Note, also arma() or arima() formulas can be applied:
    # In these cases R's internal function arima() will be used ...
    args(arima)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data =  x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Alternatively you can use also method=CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data =  x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # or method=CSS Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # or method=ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.ar2Report =
  function()
  {
    # Simulate AR(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5)), n = 1000)
    head(x)

    # Fit:
    object = armaFit(formula = ~ ar(2), data =  x, method = "mle")

    # Short Print Report:
    print(object)

    # Plot: Standardized Residuals, ACF, QQ-Plot, Ljung-Box p-Values
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Try:
    par(mfrow = c(1, 1))
    plot(object, which = 1)
    plot(object, which = 2)
    plot(object, which = 3)
    plot(object, which = 4)

    # Interactive Plot:
    # plot(object)

    # Summary Method - No Plot:
    summary(object, doplot = FALSE)

    # Summary Method - Including All Plots:
    par(mfrow = c(2, 2), cex = 0.7)
    summary(object, doplot = TRUE, which = "all")

    # Extractor Functions - Get Values:
    coefficients(object)
    coef(object)
    fitted = fitted(object)
    class(fitted)
    head(fitted)
    residuals = residuals(object)
    class(residuals)
    head(residuals)

    # Should we have ?
    # getCoef
    # getResiduals
    # getFitted

    # Predict Method:
    # predict.fARMA is now in namespace as an hidden method
    # args(predict.fARMA)
    predict(object)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.ma2Fit =
  function()
  {
    # Simulate MA(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0, ma = c(0.5, -0.5)), n = 5000)

    # To Fit a MA Model use ma(q), arma(0,q) or arima(0, 0, q):
    object = armaFit(formula = ~ ma(2), data = x)
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Note, also arma() or arima() formulas can be applied:

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # ML fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.ma2Report =
  function()
  {
    # Simulate MA(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0, ma = c(0.5, -0.5)), n = 5000)

    # To Fit a MA Model use ma(q), arma(0,q) or arima(0, 0, q):
    object = armaFit(formula = ~ ma(2), data = x)

    # Report:
    print(object)

    # Plot: Standardized Residuals, ACF, QQ-Plot, Ljung-Box p-Values
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary:
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arma21Fit =
  function()
  {
    # Simulate ARMA(2,1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), ma = 0.1), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arma(2, 1), data =  x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # Note, also arima() formulas can be applied:

    # "mle" == "ML" Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x)
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arma21Report =
  function()
  {
    # Simulate ARMA(2, 1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), ma = 0.1), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arma(2, 1), data = x)

    # Report:
    print(object)

    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary:
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arima211Fit =
  function()
  {
    # Simulate ARIMA(2, 1, 1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 1, ma = 0.1), n = 1000)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # ML Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arima211Report =
  function()
  {
    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 1, ma = 0.1), n = 1000)

    # mle Integrated ARMA Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x)

    # Report:
    print(object)

    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arfima00Fit =
  function()
  {
    # Simulate ARFIMA(0, 0):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0.3), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arfima(0, 0), data = x)
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = 0.3
    checkEqualsNumeric(target, current)

    # Parameter:
    target = unlist(object@parameter)
    print(target)
    current = c(include.mean = 1, M = 100, h = -1)
    checkIdentical(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.arfima00Report =
  function()
  {
    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0.3), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arfima(0, 0), data = x, M = 50, h = -1)

    # Report:
    print(object)

    # Plot:
    # plot(object, which = "all")         # CHECK not yet implemented

    # Summary:
    summary(object, doplot = FALSE)       # CHECK use always doplot=FALSE

    # Get Values:
    coefficients(object)
    coef(object)
    fitted = fitted(object)
    class(fitted)
    tail(fitted)                          # CHECK head yields NA's
    residuals = residuals(object)
    class(residuals)
    tail(residuals)

    # Predict:
    # predict(object)                     # CHECK not yet implemented

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.armaFormula =
  function()
  {
    # This example shows variants of alternative usage
    #   of the argument formula ...

    # Load From Ecofin Package:
    TS = MSFT
    head(TS)
    class(TS)
    colnames(TS)

    # Fit:
    armaFit(formula = diff(log(Close)) ~ ar(5), data = TS)

    # Fit:
    # Note, data may be a timeSeries object ...
    armaFit(Close ~ ar(5), data = returns(TS, digits = 12))

    # Fit:
    TS.RET = returns(TS, digits = 12)
    # Note, data may be a timeSeries object ...
    armaFit(Close ~ ar(5), TS.RET)

    # Fit:
    # Note, data may be a 'data.frame' ...
    armaFit(Close ~ ar(5), as.data.frame(TS.RET))

    # Fit:
    # Note, data may be a 'numeric' vector ...
    armaFit(~ ar(5), as.vector(TS.RET[, "Close"]))

    # Fit:
    # Note, data may be an object of class 'ts' ...
    armaFit(~ ar(5), as.ts(TS.RET)[, "Close"])

    # Fit:
    TS.RET = returns(TS)
    colnames(TS.RET)
    # attach(TS.RET)                    # CHECK doesn't work under RUnit ...
    attach(TS.RET)
    head(Close)
    armaFit(Close ~ ar(5))

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.armaArguments =
  function()
  {
    # armaFit(
    #   formula, data, method = c("mle", "ols"), include.mean = TRUE,
    #   fixed = NULL, title = NULL, description = NULL, ...)

    # arima(
    #   x, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA),
    #   xreg = NULL, include.mean = TRUE, transform.pars = TRUE,
    #   fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"),
    #   n.cond, optim.control = list(), kappa = 1e+06)

    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 1000)

    # Include Mean - by default:
    armaFit(~ arma(2, 1), data = x)

    # Don't include Mean:
    armaFit(~ arma(2, 1), data = x, include.mean = FALSE)

    # Full Model:
    armaFit(~ arma(2, 1), data = x)

    # Fixed - AR(2[2]) Subset Model:
    #                                         ar1 ar2 ma1 intercept
    armaFit(~ arma(2, 1), data = x, fixed = c(0.5, NA, NA, NA))

    # Return Value:
    return()
  }


################################################################################





