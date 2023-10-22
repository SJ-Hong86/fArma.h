#' @title Long Range Dependence Modelling
#'
#' @description
#' A collection and description of functions to investigate
#' the long range dependence or long memory behavior of an univariate
#' time series process. Included are functions to simulate factional
#' Gaussian noise and fractional ARMA processes, and functions to estimate
#' the Hurst exponent by several different methods.
#'
#' The Functions and methods are:
#'
#' Functions to simulate long memory time series processes:
#'
#' `fnmSim`	Simulates fractional Brownian motion:
#' * `mvn`	from the numerical approximation of the stochastic integral.
#' * `chol`	from the Choleski's decomposition of the covariance matrix.
#' * `lev`	using the method of Levinson.
#' * `circ`	using the method of Wood and Chan.
#' * `wave`	using the wavelet synthesis.
#'
#' `fgnSim`	Simulates fractional Gaussian noise:
#' * `beran`	using the method of Beran.
#' * `durbin`	using the method Durbin and Levinson.
#' * `paxson`	using the method of Paxson.
#'
#' `farimaSim` simulates FARIMA time series processes.
#'
#' Functions to estimate the Hurst exponent:
#'
#' `aggvarFit`	Aggregated variance method.
#' `diffvarFit`	Differenced aggregated variance method.
#' `absvalFit`	aggregated absolute value (moment) method.
#' `higuchiFit`	Higuchi's or fractal dimension method.
#' `pengFit`	Peng's or variance of residuals method.
#' `rsFit`	    R/S Rescaled Range Statistic method.
#' `perFit`	    periodogram method.
#' `boxperFit`	boxed (modified) periodogram method.
#' `hurstSlider`Interactive Display of Hurst Estimates.
#'
#' Function for the wavelet estimator:
#'
#' * `waveletFit`wavelet estimator.
#'
#' @usage
#' fbmSim(n = 100, H = 0.7, method = c("mvn", "chol", "lev", "circ", "wave"),
#'        waveJ = 7, doplot = TRUE, fgn = FALSE)
#' fgnSim(n = 1000, H = 0.7, method = c("beran", "durbin", "paxson"))
#' farimaSim(n = 1000, model = list(ar = c(0.5, -0.5), d = 0.3, ma = 0.1),
#'           method = c("freq", "time"), ...)
#'
#' aggvarFit(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
#'           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' diffvarFit(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
#'            doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' absvalFit(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5), moment = 1,
#'           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' higuchiFit(x, levels = 50, minnpts = 2, cut.off = 10^c(0.7, 2.5),
#'            doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' pengFit(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
#'         method = c("mean", "median"),
#'         doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' rsFit(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
#'       doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#' perFit(x, cut.off = 0.1, method = c("per", "cumper"),
#'        doplot = FALSE, title = NULL, description = NULL)
#' boxperFit(x, nbox = 100, cut.off = 0.10,
#'           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
#'
#' hurstSlider(x = fgnSim())
#'
#' waveletFit(x, length = NULL, order = 2, octave = c(2, 8),
#'            doplot = FALSE, title = NULL, description = NULL)
#'
#' ## S4 method for signature 'fHURST'
#' show(object)
#'
#' @param cut.off [*Fit] -
#'   a numeric vector with the lower and upper cut-off points for the estimation.
#'   They should be chosen to define a linear range. The default values are c(0.7, 2.5),
#'   i.e. 10\^0.7 and 10\^2.5, respectively.
#' @param description [*Fit] -
#'   a character string which allows for a brief description.
#' @param doplot [*Fit] -
#'   a logical flag, by default FALSE. Should a plot be displayed?
#' @param fgn [fbmSim] -
#'   a logical flag, if FALSE, the functions returns a FBM series otherwise a FGN series.
#' @param H [fgnSim] -
#'   the Hurst exponent, a numeric value between 0.5 and 1, by default 0.7.
#' @param length [waveletFit] -
#'   the length of data to be used, must be power of 2. If set to NULL,
#'   the previous power will be used.
#' @param levels [*Fit] -
#'   the number of aggregation levels or number of blocks from which the variances
#'   or moments are computed.
#' @param method [fbmSim] -
#'   the method how to generate the time series sequence, one of the following
#'   five character strings: "mvn", "chol", "lev", "circ", or "wave".
#'   [fgnSim] -
#'   the method how to generate the FGN time series sequence, one of the following
#'   three character strings: "beran", "durbin", or "paxson".
#'   [farimaSim] -
#'   the method how to generate the time series sequence, one of the following tow
#'   character strings: "freq", or "time".
#'   [pengFit] -
#'   a string naming the method how to do the averaging, either calculating the "mean"
#'   or the "median".
#'   [perFit] -
#'   a string naming the method how to fit the data, either using the peridogram itself
#'   "per", or using the cumulated periodogram "cumper".
#' @param minnpts [*Fit] -
#'   the minimum number of points or blocksize to be used to estimate the variance
#'   or moments at any aggregation level.
#' @param model a list with model parameters `ar`, `ma` and `d`. `ar` is a numeric vector
#'   giving the AR coefficients, `d` is an integer value giving the degree of differencing,
#'   and `ma` is a numeric vector giving the MA coefficients. Thus the order
#'   of the time series process is FARMA(p, d, q) with `p=length(ar)` and `q=length(ma)`.
#'   `d` is a fractional value for FARMA models. By default an FARMA(2, d, 1) model
#'   with coefficients `ar=c(0.5, -0.5)`, `ma=0.1`, and `d=0.3` will be generated.
#' @param moment [absvalHurst] -
#'   an integer value, by default 1 which denotes absolute values. For values larger than
#'   one this argument determines what absolute moment should be calculated.
#' @param n [fgnSim][farimaSim] -
#'   number of data points to be simulated, a numeric value, by default 1000.
#' @param nbox [boxperFit] -
#'   is the number of boxes to divide the data into. A numeric value, by default 100.
#' @param object an object of class `fHurst`.
#' @param octave [waveletFit] -
#'   beginning and ending octave for estimation. An integer vector with two elements.
#'   By default c(2, 8). If the upper value is too large, it will be replaced by
#'   the maximum allowed value.
#' @param order [waveletFit] -
#'   the order of the wavelet. An integer value, by default 2.
#' @param title a character string which allows for a project title.
#' @param trace a logical value, by defaul FALSE. Should the estimation process be traced?
#' @param waveJ [fbmSim] -
#'   an integer parameter for the simulation of FBM using the wavelet method.
#' @param X [*Fit] -
#'    the numeric vector of data, an object of class `timeSeries`, or any other object
#'    which can be transofrmed into a numeric vector by the function `as.vector`.
#' @param ... arguments to be passed.
#' @keywords LrdModelling {fArma}
#' @return `fgnSim` and `farimaSim` return a numeric vector of length `n`,
#'   the FGN or FARIMA series.
#'
#'   `Fit` returns an S4 object of class `fHURST` with the following slots:
#'   * `@call`      the function call.
#'   * `@method`    a character string with the selected method string.
#'   * `@hurst`     a list with at least one element, the Hurst exponent named `H`.
#'                  Optional values may be the value of the fitted slope `beta`,
#'                  or information from the fit.
#'   * `@parameters`a list with a varying number of elements describing the input
#'                  from the argument list.
#'   * `@data`      a list holding the input data.
#'   * `@fit`       a list object with all information of the fit.
#'   * `@plot`      a list object which holds information to create a plot of the fit.
#'   * `@title`     a character string with the name of the test.
#'   * `@description`a character string with a brief description of the test.
#'
#'   * waveletFit
#'
#' @details
#' Functions to Simulate Long Memory Processes:
#'
#' Fractional Gaussian Noise:
#' The function `fgnSim` simulates a series of fractional Gaussian noise, FGN. FGN provides a
#' parsimonious model for stationary increments of a self-similar process parameterized by
#' the Hurst exponent H and variance. Fractional Gaussian noise with H < 0.5
#' demonstrates negatively autocorrelated or anti-persistent behaviour, and FGN with
#' H > 0.5 demonstrates 1/f , long memory or persistent behaviour, and the special case.
#' The case H = 0.5 corresponds to the classical Gaussian white noise.
#' One can select from three different methods. The first generator named "beran" uses the
#' fast Fourier transform to generate the series based on SPLUS code written
#' originally by J. Beran [1994]. The second generator named "durbin" produces a FGN series
#' by using the Durbin-Levinson coefficients.
#' The algorithm was reimplemented in pure S based on the C source code written by V. Teverovsky [199x].
#' The third generator named "paxson" was proposed by V. Paxson [199x],
#' this approaximate method is a very fast and requires low storage.
#' However, the algorithm reveals some weakness in the method which was
#' discussed by D.A. Rolls [2001].
#'
#' Fractional ARIMA Processes:
#' The function `farimaSim` is a generator for fractional ARIMA time series processes.
#' A Gaussian FARIMA(0,d,0) series can be created, where d is related to the the Hurst
#' exponent H through d=H-0.5. This is a particular case of the more general Gaussian
#' FARIMA(p,d,q) process which follows the same asymptotic relations for their
#' autocovariance and the spectral density as do the Gaussian FARIMA(0,d,0) processes.
#' Two different generators are implement in S. The first named "freq" works in the
#' frequence domain and generates the series from the fast Fourier transform based on
#' SPLUS code written originally by J. Beran [1994]. The second method creates the series
#' in the time domain, therefore named "time". The algorithm was reimplemented in pure
#' S based on the Fortran source code from the R's `fracdiff` package originally written by
#' C. Fraley [1991]. Details for the algorithm are given in J Haslett and A.E. Raftery [1989].
#'
#' Functions to Estimate the Hurst Exponent:
#'
#' These are 9 functions as described by M.S. Taqqu, V. Teverovsky, and W. Willinger [1995]
#' to estimate the self similarity parameter and/or the intensity of long-range dependence
#' in a time series.
#'
#' Aggregated Variance Method:
#' The function `aggvarFit` computes the Hurst exponent from the variance of an
#' aggregated FGN or FARIMA time series process. The original time series is divided into
#' blocks of size `m`. Then the sample variance within each block is computed. The slope
#' `beta=2*H-2` from the least square fit of the logarithm of the sample variances versus the
#' logarithm of the block sizes provides an estimate for the Hurst exponent `H`.
#'
#' Differenced Aggregated Variance Method:
#' To distinguish jumps and slowly decaying trends which are two types of non-stationary,
#' from long-range dependence, the function `diffvarFit` differences the sample variances
#' of successive blocks. The slope `beta=2*H-2` from the least square fit of the logarithm of
#' the differenced sample variances versus the logarithm of the block sizes provides an
#' estimate for the Hurst exponent `H`.
#'
#' Aggregated Absolute Value/Moment Method:
#' The function `absvalFit` computes the Hurst exponent from the moments `moment=M` of
#' absolute values of an aggregated FGN or FARIMA time series process. The first moment
#' `M=1` coincides with the absolute value method, and the second moment `M=2` with the
#' aggregated variance method. Again, the slope `beta=M*(H-1)` of the regression line of the
#' logarithm of the statistic versus the logarithm of the block sizes provides an estimate
#' for the Hurst exponent `H`.
#'
#' Higuchi or Fractal Dimension Method:
#' The function `higuchiFit` implements a technique which is very similar to the absolute
#' value method. Instead of blocks a sliding window is used to compute the aggregated
#' series. The function involves the calculation the calculation of the length of a path and,
#' in principle, finding its fractal Dimension `D`. The slope `D=2-H` from the least square fit of
#' the logarithm of the expected path lengths versus the logarithm of the block (window)
#' sizes provides an estimate for the Hurst exponent `H`.
#'
#' Peng or Variance of Residuals Method:
#' The function `pengFit` uses the method described by Peng. In Peng's variance of
#' residuals method the series is also divided into blocks of size `m`. Within each block the
#' cumulated sums are computed up to t and a least-squares line `a+b*t` is fitted to the
#' cumulated sums. Then the sample variance of the residuals is computed which is
#' proportional to `m^(2*H)`. The "mean" or "median" are computed over the blocks. The
#' slope `beta=2*H` from the least square provides an estimate for the Hurst exponent H.
#'
#' The R/S Method:
#' The function `rsFit` implements the algorithm named rescaled range analysis which is
#' dicussed for example in detail by B. Mandelbrot and Wallis [199x], B. Mandelbrot [199x]
#' and B. Mandelbrot and M.S. Taqqu [199x].
#'
#' The Periodogram Method:
#' The function `perFit` estimates the Hurst exponent from the periodogram. In the finite
#' variance case, the periodogram is an estimator of the spectral density of the time series.
#' A series with long range dependence will show a spectral density with a lower law
#' behavior in the frequency. Thus, we expect that a log-log plot of the periodogram
#' versus frequency will display a straight line, and the slopw can be computed as 1-2H.
#' In practice one uses only the lowest 10% of the frequencies, since the power law behavior
#' holds only for frequencies close to zero. Varying this cut off may provide additional
#' information. Plotting H versus the cut off, one should select that cut off where the curve
#' flattens out to estimate `H`. This approach can be selected by the argument `method="per"`.
#' Alternatively we can select `method="cumper"`. In this case, instead of using the
#' periodgram itself, the cmulative periodgram will be investigated. The slope of the
#' double logarithmic fit is given by 2-2H. More details can be found in the work of J.
#' Geweke and S. Porter-Hudak [1983] and in Taqqu [?].
#'
#' The Boxed or Modified Periodogram Method:
#' The function `boxperFit` is a modification of the periodogram method. The algorithm
#' devides the frequency axis into logarithmically equally spaced boxes and averages the
#' periodogram values corresponding to the frequencies inside the box.
#'
#' The original functions were written by V. Teverovsky and W. Willinger for SPLUS calling
#' internal functions written in C. The software can be found on M. Taqqu's home page:
#' http://math.bu.edu/people/murad/
#' In addition the Whittle estimator uses SPlus functions written by J. Beran. They can be
#' found in the appendix of his book or on the StatLib server:
#' http://lib.stat.cmu.edu/S/
#' Note, all nine R functions and internal utility functions are reimplemented entirely in S.
#'
#' Functions to perform a Wavelet Analysis:
#'
#' The function `waveletFit` computes the Discrete Wavelet Transform, averages the squares
#' of the coefficients of the transform, and then performs a linear regression on the
#' logarithm of the average, versus the log of the scale parameter of the transform. The
#' result should be directly proportional to H providing an estimate for the Hurst exponent.
#'
#' @examples
#'
#' # fgnSim
#' par(mfrow = c(3, 1), cex = 0.75)
#'
#' # Beran's Method
#' plot(fgnSim(n = 200, H = 0.75), type = "l",
#'      ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Beran")
#'
#' # Durbin's Method
#' plot(fgnSim(n = 200, H = 0.75, method = "durbin"),
#'      type = "l", ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Durbin")
#'
#' # Paxson's Method
#' plot(fgnSim(n = 200, H = 0.75, method = "paxson"), type = "l",
#'      ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Paxson")



################################################################################
# FUNCTIONS:          HURST EXPONENT:
#  'fHURST'            S4 Class Representation
#  show.fHURST         S3 Print Method
#  plot.fHURST         S3 Plot Method
#  aggvarFit           3.1 Aggregated variance method
#  diffvarFit          3.2 Differenced aggregated variance method
#  absvalFit           3.3 Absolute values (moments) method
#  higuchiFit          3.4 Higuchi's method
#  pengFit             3.5 Peng's or Residuals of Regression method
#  rsFit               3.6 R/S method
#  perFit              3.7 Periodogram and cumulated periodogram method
#  boxperFit           3.8 Boxed (modified) peridogram method
#  whittleFit          3.9 Whittle estimator -> PART II
#  hurstSlider         Hurst Slider
################################################################################


################################################################################
# Reimplemented functions from
#   Taqqu M.S, Teverovsky V, Willinger W.
#   Estimators for Long-Range Dependence: An Empirical Study
#   Fractals, Vol 3, No. 4, 785-788, 1995

# ------------------------------------------------------------------------------
#' @export

setClass("fHURST",
         representation(
           call = "call",
           method = "character",
           hurst = "list",
           parameter = "list",
           data = "list",
           fit = "list",
           plot = "list",
           title = "character",
           description = "character")
)


# ------------------------------------------------------------------------------
#' @export

setMethod("show", "fHURST",
          function(object)
          {   # A function implemented by Diethelm Wuertz

            # Description:
            #   Prints a fHURST Object

            # FUNCTION:

            # Setting:
            x = object
            doplot = TRUE

            # Title:
            cat("\nTitle:\n ", x@title, "\n", sep = "")

            # Call:
            cat("\nCall:\n ")
            cat(paste(deparse(x@call), sep = "\n", collapse = "\n"), "\n", sep = "")

            # Method:
            cat("\nMethod:\n ", x@method, "\n", sep = "")

            # Hurst Exponent:
            cat("\nHurst Exponent:\n")
            H = as.numeric(unlist(x@hurst)[1:2])
            names(H) = names(unlist(x@hurst)[1:2])
            output = capture.output(H)
            cat(paste(" ", output), sep = "\n")

            # Hurst Exponent Diagnostic:
            if (!is.null(x@hurst$diag)) {
              cat("\nHurst Exponent Diagnostic:\n ")
              print(x@hurst$diag[2, ])
            }

            # Parameter Settings:
            cat("\nParameter Settings:\n")
            parms = unlist(x@parameter)
            integer.parms = as.integer(parms)
            names(integer.parms) = names(parms)
            # output = capture.output(integer.parms)
            # cat(paste(" ", output), sep = "\n")
            print(integer.parms)

            # Description:
            cat("\nDescription:\n ", x@description, sep = "")
            cat("\n\n")

            # Plot:
            fit = object
            if (x@plot$doplot) {
              labels = TRUE
              if (labels) {
                xlab = fit@plot$xlab
                ylab = fit@plot$ylab
                H = as.character(round(fit@hurst$H, digits = 4))
                main = paste(fit@method, "\n H =", H)
                M = fit@plot$m[fit@plot$weights == 1]
                min.M = as.character(min(M))
                max.M = as.character(max(M))
                M = as.character(length(M))
                gridlines = TRUE
              } else {
                xlab = ""
                ylab = ""
                main = ""
                gridlines = FALSE
              }
              # Plot:
              x = c(0, log10(fit@plot$m))
              y = c(0, log10(fit@plot$data))
              wt = fit@plot$weights
              plot(x, y, type = "n", xlab = xlab, ylab = ylab)
              title(main = main)
              if (gridlines) grid()
              x = x[-1]
              y = y[-1]
              points(x[wt == 1], y[wt == 1], pch = 19, cex = fit@plot$cex)
              points(x[wt == 0], y[wt == 0], pch = 3, cex = fit@plot$cex)
              # y = mx + c
              m = fit@fit$coeff[[2]]
              a = fit@fit$coeff[[1]]
              x.fit = x[wt == 1]
              y.fit = m * x.fit + a
              lines(x.fit, y.fit, lwd = 2)
              if (is.numeric(fit@plot$abline))
                abline(fit@plot$abline[1], fit@plot$abline[2], lty = 3)
            }

            # Return Value:
            invisible()
          })


# ------------------------------------------------------------------------------
#' @export

aggvarFit <-
  function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Aggregated Variance - [Taqqu 3.1]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, aggregated variances
    #       'AGGVAR', and weights for the fit ,'wt', the numeric
    #       values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M > 0]

    # Create Data:
    AGGVAR = NULL
    for (m in M) {
      nCols = n %/% m
      X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
      STATS = var(colMeans(X))
      AGGVAR = c( AGGVAR, STATS )
      if (trace) cat("\n\tm = \t", m, "\tAggVar = \t", STATS  )
    }
    if(trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(AGGVAR), wt)
    fitH = lsfit(log10(M), 0.5*log10(AGGVAR*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (beta + 2) / 2

    # Return Value:
    plot = list(m = M, data = AGGVAR, weights = wt,
                abline = c(0, -1), cex = 0.7, doplot = doplot,
                xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Aggregated Variances"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Aggregated Variance Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

diffvarFit <-
  function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   # Differenced Aggregated Variance - [Taqqu 3.2]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    n = length(as.vector(x))
    data = list(x = x)
    x = as.vector(x)

    # Compute Aggregated Variances:
    ans = aggvarFit(x, levels, minnpts, cut.off)

    # Create Differenced Data:
    DIFFVAR = -diff(ans@plot$data)

    # What M's to use?
    # M = ( ans@plot$data[-levels, 1] + ans@plot$data[-1] ) / 2
    # M = sqrt ( ans@plot$data[-levels] * ans@plot$data[-1] )
    # M = ans@plot$data[-1]
    M = ans@plot$m[-levels]

    # Remove negative and zero values:
    M = M[DIFFVAR > 0]
    wt = (ans@plot$weights[-levels])[DIFFVAR > 0]
    DIFFVAR = DIFFVAR[DIFFVAR > 0]

    if (trace) {
      for ( i in 1:length(M) )
        cat("\n\tm = \t", M[i], "\tDiffVar = \t", DIFFVAR[i]  )
      cat("\n")
    }

    # Fit:
    fit = lsfit(log10(M), log10(DIFFVAR), wt)
    fitH = lsfit(log10(M), 0.5*log10(DIFFVAR*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (beta + 2) / 2

    # Return Value:
    plot = list(m = M, data = DIFFVAR, weights = wt,
                abline = FALSE, cex = 0.7, doplot = doplot,
                xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Differenced Aggregated Variances"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Differenced Aggregated Variance",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

absvalFit <-
  function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5), moment = 1,
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Absolute Value/Moments Method - [Taqqu 3.3]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))

    # Compute Absolute Moments:
    ABSVAL = NULL
    for (m in M) {
      nCols = n %/% m
      X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
      Y = colMeans(X)
      MEAN = mean(Y)
      STATS = sum( (abs(Y-MEAN))^moment ) / (length(Y) - 1)
      ABSVAL = c( ABSVAL, STATS )
      if (trace) cat("\n\tm = \t", m, "\tAbsVal = \t", STATS  )
    }
    if(trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(ABSVAL), wt)
    fitH = lsfit(log10(M), log10(ABSVAL*M^moment)/moment, wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta/moment + 1

    # Return Value:
    plot = list(m = M, data = ABSVAL, weights = wt,
                abline = c(0, -0.5), cex = 0.7, doplot = doplot,
                xlab = "log10(m)", ylab = "log10(variances)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Absolute Values"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = paste("Absolute Moment - No.", as.character(moment)),
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off, moment = moment),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

higuchiFit =
  function(x, levels = 50, minnpts = 2, cut.off = 10^c(0.7, 2.5),
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Higuchi Method / Fratal Dimension Method - [Taqqu 3.4]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    y = cumsum(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))

    # Higuchi Method:
    if (trace) cat("\nHiguchi Iteration Path:")
    HIGUCHI = NULL
    for ( m in M ) {
      k.max = max(floor((n-(1:m))/m) )
      X = matrix(rep(0, length = m*k.max), byrow = FALSE, ncol = k.max)
      for ( i in 1:m ) {
        for ( k in 1:(floor((n-i)/m)) ) {
          X[i, k] = abs(y[1+k*m] - y[i+(k-1)*m]) / floor((n-i)/m)
        }
      }
      Y = sum(X) * (n-1) / m^3
      STATS = Y / n
      HIGUCHI = c( HIGUCHI, STATS )
      if (trace) cat("\n\tm = \t", m, "\tHiguchi = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(HIGUCHI), wt)
    fitH = lsfit(log10(M), log10(HIGUCHI*M*M), wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta + 2

    # Return Value:
    plot = list(m = M, data = HIGUCHI, weights = wt,
                abline = c(0, -0.5), cex = 0.7, doplot = doplot,
                xlab = "log10(m)", ylab = "log10(curve length)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Higuchi Method"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Higuchi Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

pengFit <-
  function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
           method = c("mean", "median"),
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Peng's Method / Variance of Residuals - [Taqqu 3.5]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M>2]

    # Averaging Function:
    stats = match.fun(method[1])

    # Peng's Method:
    PENG = NULL
    for (m in M) {
      nCols = n %/% m
      X = matrix(x[1:(m*nCols)], byrow = FALSE, ncol = nCols)
      Y = colCumsums(X)
      V = NULL
      t = cbind(1, as.matrix(1:m))
      for (i in 1:nCols ) {
        y = Y[, i]
        nrx = nry = NROW(t)
        ncx = NCOL(t)
        ncy = NCOL(y)
        # 17-09-2012 (YC): reverted external .Fortran call to R
        # function lm.fit to comply with new CRAN policy
        res <- lm.fit(t, y)$residuals
        V = c(V, var(res))
      }
      STATS = stats(V)
      PENG = c(PENG, STATS)
      if (trace) cat("\n\tm = \t", m, "\tPENG = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(PENG), wt)
    fitH = lsfit(log10(M), log10(PENG)/2, wt)
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta/2

    # Return Value:
    plot = list(m = M, data = PENG, weights = wt,
                abline = FALSE, cex = 0.7, doplot = doplot,
                xlab = "log10(m)", ylab = "log10(var[residuals])")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Peng Method"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Peng Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

rsFit <-
  function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5),
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   R/S Statistic Method - [Taqqu 3.6]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    increment = (log10(n/minnpts))/levels
    M = floor(10^((1:levels)*increment))
    M = M[M > 1]

    # R/S Method:
    Y = cumsum(x)
    Y2 = cumsum(x*x)
    RS = NULL
    for (m in M) {
      S = sqrt(Y2[m]/m - (Y[m]/m)^2)
      Z = Y[1:m]-(1:m)*Y[m]/m
      STATS = (max(Z) - min(Z))/S
      RS = c(RS, STATS)
      if (trace) cat("\n\tm = \t", m, "\tR/S = \t", STATS  )
    }
    if (trace) cat("\n")

    # Fit:
    wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
    fit = lsfit(log10(M), log10(RS), wt)
    fitH = fit
    fitH$wt = NULL
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = beta

    # Plot Values:
    plot = list(m = M, data = RS, weights = wt,
                abline = FALSE, cex = 0.7, doplot = doplot,
                xlab = "log10(d)", ylab = "log10(r/s)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from R/S Method"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "R/S Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, levels = levels, minnpts = minnpts,
                         cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

perFit <-
  function(x, cut.off = 0.10,
           method = c("per", "cumper"), doplot = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Periodogram Method - [Taqqu 3.7]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.
    #   levels - the number of aggregation levels
    #   minnpts - the minimum block size
    #   cut.off - the lower and upper cut off for the fit

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    n = length(x)
    FFT = Mod(fft(x))^2/(2*pi*n)
    pgram = FFT[1:(n %/% 2+1)]
    N = length(pgram)

    # Periodogram Method:
    if (method[1] == "per") {
      Method = "Periodogram Method"
      X = (pi/n)*c(2:((n*cut.off)))
      Y = pgram[2:((n*cut.off))]
      fit = lsfit(x = log10(X), y = log10(Y))
      fitH = lsfit(log10(X), log10(X/Y)/2)
      diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
      beta = fit$coef[[2]]
      H = (1-beta)/2
      U = (pi/n)*(1:n)
      V = FFT
    }

    # Cumulated Periodogram Method:
    if (method[1] == "cumper") {
      Method = "Cumulated Periodogram Method"
      PGRAM = cumsum(pgram[2:n])
      U = (pi/n)*c(1:(n-1))
      V = PGRAM[1:(n-1)]
      X = (pi/n)*c(1:(((n-1)*cut.off)))
      Y = PGRAM[1:(((n-1)*cut.off))]
      fit = lsfit(x = log10(X), y = log10(Y))
      fitH = lsfit(log10(X), log10(X*X/Y)/2)
      diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
      beta = fit$coef[[2]]
      H = (2-beta)/2
      U = (pi/n)*(1:n)
      V = cumsum(FFT)
    }

    # Plot Values:
    plot = list(m = U, data = V, weights = rep(1, times = n),
                abline = FALSE, cex = 0.25, doplot = doplot,
                xlab = "log10(frequency)", ylab = "log10(periodogram)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Periodgram Method"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = Method,
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, cut.off = 100*cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------
#' @export

boxperFit <-
  function(x, nbox = 100, cut.off = 0.10,
           doplot = FALSE, trace = FALSE, title = NULL, description = NULL)
  {
    # A functions implemented by Diethelm Wuertz

    # Description:
    #   Boxed (Modified) Periodogram Method - [Taqqu 3.8]

    # Arguments:
    #   x - a numeric vector, a 'timeSeries' object, or any other
    #       object which can be transformed into a vector by the
    #       function 'as.vector'.

    # Value:
    #   Returns a list with the folllowing elements:
    #   data - a data frame with blocksizes M, differenced aggregated
    #       variances 'DIFFVAR', and weights for the fit ,'wt', the
    #       numeric values of 'beta' and the Hurst exponent 'H'.

    # FUNCTION:

    # Settings:
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    len = length(x)
    pgram = (Mod(fft(x))^2/(2*pi*len)) [1:(len %/% 2+1)]
    n = length(pgram)

    # Calculate fractions from percentage:
    per1 = cut.off
    per2 = 1 - per1
    m = log10(per2 * n) / nbox

    # Do the boxes (except for beginning few points):
    padj = z = NULL
    for (i in 1:nbox) {
      m1 = floor(10^(m * i - m) + per1 * n)
      m2 = floor(10^(m * i) + per1 * n)
      padj[i] = sum(pgram[m1:m2])/(m2 - m1 + 1)
      z[i] = log10((pi * (m2 + m1))/(2 * n))
    }

    # x|y points:
    X = c( 0, log10((pi/n) * (2:floor(per1 * n))) )
    Y = c( 0, log10(pgram[2:floor(per1 * n)]) )
    i = (floor(per1 * n) + 1):(floor(per1 * n) + nbox)
    X = c(X, z[i - floor(per1 * n)] )
    Y = c(Y, log10(padj[i - floor(per1 * n)]) )

    # Fit:
    XN = 10^X
    YN = 10^Y
    fit = lsfit(log10(XN), log10(YN))
    fitH = lsfit(log10(XN), log10(XN/YN)/2)
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (1-beta)/2

    # Plot Values:
    plot = list(m = XN, data = YN, weights = rep(1, times = length(XN)),
                abline = FALSE, cex = 0.5, doplot = doplot,
                xlab = "log10(frequency)", ylab = "log10(periodogram)")

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Boxed Periodgram Method"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Boxed Periodogram",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(n = n, nbox = nbox, cut.off = cut.off),
        data = data,
        fit = fit,
        plot = plot,
        title = title,
        description = description
    )
  }


################################################################################









