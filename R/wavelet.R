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
#' @export
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

waveletFit <-
  function(x, length = NULL, order = 2, octave = c(2, 8),
           doplot = FALSE, title = NULL, description = NULL)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Function to do the Wavelet estimator of H.

    # Arguments:
    #   x - Data set.
    #   length - Length of data to be used (must be power of 2)
    #       if NULL, the previous power will be used
    #   octave - Beginning and ending octave for estimation.

    # Details:
    #   This method computes the Discrete Wavelet Transform, averages the
    #   squares of the coefficients of the transform, and then performs a
    #   linear regression on the logarithm of the average, versus the log
    #   of j, the scale parameter of the transform. The result should be
    #   directly proportional to H.
    #   There are several options available for using this method: method.
    #   1.  The length of the data must be entered (power of 2).
    #   2.  c(j1, j2) are the beginning and ending octaves for the estimation.
    #   3.  'order' is the order of the wavelet. (2 default)
    #   5.  Calls functions from R's Wavelet package. ( wd, accessD ).
    #   6.  Inside function, a bound.effect is used in the estimation to
    #       avoid boundary effects on the coefficients.

    # Authors:
    #   Based on work by Ardry and Flandrin.
    #   Originally written by Vadim Teverovsky 1997.

    # Notes:
    #   Calls functions from R's package 'wavethresh'

    # FUNCTION:

    # Settings:
    N = order
    call = match.call()
    data = list(x = x)
    x = as.vector(x)
    j1 = octave[1]
    j2 = octave[2]
    if(is.null(length)) length = 2^floor(log(length(x))/log(2))
    noctave = log(length, base = 2) - 1
    bound.effect = ceiling(log(2*N, base = 2))

    # Calculate:
    transform = .wd(x[1:(length)], filter.number = N)
    statistic = rep(0, noctave)
    if (j2 > noctave - bound.effect) {
      # cat("Upper bound too high, resetting to ", noctave-bound.effect, "\n")
      j2 = noctave - bound.effect
      octave[2] = j2
    }
    for (j in 1:(noctave - bound.effect)) {
      statistic[j] = log(mean((.accessD(transform,
                                        level = (noctave+1-j))[N:(2^(noctave+1-j)-N)])^2), base = 2)
    }

    # Fit:
    X = 10^c(j1:j2)
    Y = 10^statistic[j1:j2]
    fit = lsfit(log10(X), log10(Y))
    fitH = lsfit(log10(X), log10(Y*X)/2)
    diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
    beta = fit$coef[[2]]
    H = (beta+1)/2

    # Plot Values:
    plot = list(
      m = 10^c(1:(noctave - bound.effect)),
      data = 10^statistic[1:(noctave - bound.effect)],
      weights = rep(1, times = length(X)),
      abline = FALSE, cex = 0.5, doplot = doplot,
      xlab = "Octave", ylab = "Statistic")

    # Add Slope:
    # abline(fit$coef[[1]], fit$coef[[2]])}

    # Add:
    if (is.null(title))
      title = "Hurst Exponent from Wavelet Estimator"
    if (is.null(description))
      description = description()

    # Return Value:
    new("fHURST",
        call = call,
        method = "Wavelet Method",
        hurst = list(H = H, beta = beta, diag = diag),
        parameter = list(length = length, order = order, octave = octave),
        data = data,
        plot = plot,
        fit = fit,
        title = title,
        description = description
    )
  }


# ------------------------------------------------------------------------------


# Following auxilliary code copied from:


# Package: wavethresh
# Version: 2.2-8
# Note:    ----- Version also in wvrelease() in ./R/release.R
# Date: 2004-03-08
# Author: Guy Nason <G.P.Nason@Bristol.ac.uk>
#   of R-port: Arne Kovac (1997) and Martin Maechler (1999)
# Maintainer: Martin Maechler <maechler@stat.math.ethz.ch>
# Title: Software to perform wavelet statistics and transforms.
# Description: Software to perform 1-d and 2-d wavelet statistics and transforms
# Depends: R (>= 1.4)
# License: GPL version 2 or later

# ------------------------------------------------------------------------------
#' @export

.accessD <-
  function(wd.obj, level, boundary = FALSE)
  {   # A function copied from R-package wavethresh

    # FUNCTION:

    # Settings:
    ctmp = class(wd.obj)
    if (is.null(ctmp) || all(ctmp != "wd"))
      stop("argument `wd.obj' is not of class \"wd\"")
    if (level < 0)
      stop("Must have a positive level")
    else if(level > wd.obj$nlevels - 1)
      stop("`level' must be less than resolution (= nlevels)")
    level = level + 1
    first.last.d = wd.obj$fl.dbase$first.last.d
    first.level = first.last.d[level, 1]
    last.level  = first.last.d[level, 2]
    off.l = first.last.d[level, 3]

    if (boundary) {
      n = last.level - first.level + 1
      wd.obj$D[(off.l + 1):(off.l + n)]
    } else {
      n = 2^(level - 1)
      wd.obj$D[(off.l + 1 - first.level):(off.l + n - first.level)]
    }
  }


# ------------------------------------------------------------------------------
#' @export

.wd <-
  function(data, filter.number = 2, family = c("DaubExPhase", "DaubLeAsymm"),
           bc = c("periodic", "symmetric"), verbose = getOption("verbose"))
  {
    # A function copied from R-package wavethresh

    # FUNCTION:

    # Settings:
    if (verbose) cat("Argument checking...")
    family = match.arg(family)
    bc = match.arg(bc)
    DataLength = length(data)

    # Check that we have a power of 2 data elements
    nlevels = log(DataLength)/log(2)
    if(round(nlevels) != nlevels)
      stop("The length of data is not a power of 2")

    # Select the appropriate filter
    if (verbose) cat("...done\nFilter...")
    filter = .filter.select(filter.number = filter.number, family = family)

    # Build the first/last database
    if (verbose) cat("...selected\nFirst/last database...")
    fl.dbase = .first.last(LengthH = length(filter$H), DataLength =
                             DataLength, bc = bc)

    # Put in the data
    C = rep(0, fl.dbase$ntotal)
    C[1:DataLength] = data
    error = if (verbose) 1 else 0
    if (verbose) cat("built\nObject code...")

    # Compute the decomposition
    if(verbose) cat("Decomposing...\n")
    nbc = switch(bc, periodic = 1, symmetric = 2)
    wavelet.decomposition = .C("wavedecomp",
                               C = as.double(C),
                               LengthC = as.integer(fl.dbase$ntotal),
                               D = double(fl.dbase$ntotal.d),
                               LengthD = as.integer(fl.dbase$ntotal.d),
                               H = as.double(filter$H),
                               LengthH = as.integer(length(filter$H)),
                               nlevels = as.integer(nlevels),
                               firstC = as.integer(fl.dbase$first.last.c[, 1]),
                               lastC = as.integer(fl.dbase$first.last.c[, 2]),
                               offsetC = as.integer(fl.dbase$first.last.c[, 3]),
                               firstD = as.integer(fl.dbase$first.last.d[, 1]),
                               lastD = as.integer(fl.dbase$first.last.d[, 2]),
                               offsetD = as.integer(fl.dbase$first.last.d[, 3]),
                               nbc = as.integer(nbc),
                               error = as.integer(error),
                               PACKAGE = "fArma")

    if (verbose) cat("done\n")
    error = wavelet.decomposition$error
    if (error != 0) stop(paste("Error", error, " occured in wavedecomp"))

    # Result:
    l = list(C = wavelet.decomposition$C, D = wavelet.decomposition$D,
             nlevels = wavelet.decomposition$nlevels, fl.dbase = fl.dbase,
             filter = filter, bc = bc)
    class(l) = "wd"

    # Return Value:
    l
  }


# ------------------------------------------------------------------------------
#' @export

.filter.select <-
  function(filter.number, family = c("DaubExPhase", "DaubLeAsymm"),
           constant = 1)
  {
    # A function copied from R-package wavethresh

    # FUNCTION:

    # Settings:
    family = match.arg(family)# one of the two, maybe abbrev. in call
    if ((filter.number = as.integer(filter.number)) <= 0)
      stop("invalid `filter.number'")
    if (family == "DaubExPhase") {

      # The following wavelet coefficients are taken from
      #   Daubechies, I (1988)
      #   Orthonormal Bases of Wavelets
      #   Communications on Pure and Applied Mathematics. Page 980
      # or
      #   Ten Lectures on Wavelets, Daubechies, I, 1992
      #   CBMS-NSF Regional Conference Series, page 195, Table 6.1
      #   Comment from that table reads:
      #   "The filter coefficients for the compactly supported wavelets
      #   with extremal phase and highest number of vanishing moments
      #   compatible with their support width".
      filter.name =
        switch(filter.number,
               {
                 # filter.number  -- 1 --
                 # This is for the Haar basis. (not in Daubechies).
                 H = rep(1/sqrt(2), 2)
                 "Haar wavelet"
               },
               { # filter.number  -- 2 --
                 H = c(0.48296291314500001,   0.83651630373800001,
                       0.22414386804200001,
                       -0.12940952255099999)
                 "Daub cmpct on ext. phase N=2"
               },
               { # filter.number  -- 3 --
                 H = c(0.33267055294999998,   0.80689150931099995,
                       0.45987750211799999,  -0.13501102001000001,
                       -0.085441273882,        0.035226291882000001)
                 "Daub cmpct on ext. phase N=3"
               },
               { # filter.number  -- 4 --
                 H = c(0.23037781330900001,    0.71484657055300005,
                       0.63088076793000003,   -0.027983769417,
                       -0.18703481171899999,    0.030841381836,
                       0.032883011667000001,  -0.010597401785)
                 "Daub cmpct on ext. phase N=4"
               },
               { # filter.number  -- 5 --
                 H = c(0.160102397974,         0.60382926979700002,
                       0.72430852843799998,    0.138428145901,
                       -0.242294887066,        -0.032244869585000002,
                       0.07757149384,         -0.006241490213,
                       -0.012580751999,         0.003335725285)
                 "Daub cmpct on ext. phase N=5"
               },
               { # filter.number  -- 6 --
                 H = c(0.11154074335,          0.49462389039799998,
                       0.751133908021,         0.31525035170900001,
                       -0.22626469396500001,   -0.12976686756700001,
                       0.097501605586999995,   0.02752286553,
                       -0.031582039318000001,   0.000553842201,
                       0.004777257511,        -0.001077301085)
                 "Daub cmpct on ext. phase N=6"
               },
               {  # filter.number -- 7 --
                 H = c(0.077852054084999997,   0.396539319482,
                       0.72913209084599995,    0.469782287405,
                       -0.14390600392899999,   -0.22403618499399999,
                       0.071309219267,         0.080612609151000006,
                       -0.038029936935000001,  -0.016574541631,
                       0.012550998556,         0.000429577973,
                       -0.001801640704,         0.0003537138)
                 "Daub cmpct on ext. phase N=7"
               },
               { # filter.number  -- 8 --
                 H = c(0.054415842243000001,   0.31287159091400002,
                       0.67563073629699999,    0.58535468365400001,
                       -0.015829105256,        -0.28401554296199999,
                       0.000472484574,         0.12874742662,
                       -0.017369301002000001,  -0.044088253931,
                       0.013981027917,         0.008746094047,
                       -0.004870352993,        -0.000391740373,
                       0.000675449406,        -0.000117476784)
                 "Daub cmpct on ext. phase N=8"
               },
               { # filter.number  -- 9 --
                 H = c(0.038077947363999998,   0.24383467461300001,
                       0.60482312369000002,    0.65728807805099998,
                       0.13319738582499999,   -0.293273783279,
                       -0.096840783222999993,   0.14854074933799999,
                       0.030725681479000001,  -0.067632829061000002,
                       0.000250947115,         0.022361662124,
                       -0.004723204758,        -0.004281503682,
                       0.001847646883,         0.000230385764,
                       -0.000251963189,         3.934732e-05)
                 "Daub cmpct on ext. phase N=9"
               },
               { # filter.number -- 10 --
                 H = c(0.026670057901000001,   0.188176800078,
                       0.52720118893199996,    0.688459039454,
                       0.28117234366100002,   -0.24984642432699999,
                       -0.19594627437699999,    0.127369340336,
                       0.093057364604000006,  -0.071394147165999997,
                       -0.029457536822,         0.033212674058999997,
                       0.003606553567,        -0.010733175483,
                       0.001395351747,         0.001992405295,
                       -0.000685856695,        -0.000116466855,
                       9.358867e-05,          -1.3264203e-05)
                 "Daub cmpct on ext. phase N=10"
               }
        ) #  switch ( filter.number )
      if (is.null(filter.name))
        stop(paste("Unknown filter number (not in {1:10}) for", family,
                   "(Daubechies wavelets with extremal phase...)"))

    } else {
      # if(family == "DaubLeAsymm")
      # The following wavelet coefficients are taken from
      # Ten Lectures on Wavelets, Daubechies, I, 1992
      # CBMS-NSF Regional Conference Series, page 198, Table 6.3
      # Comment from that table reads:
      #   "The low pass filter coefficients for the "least-asymmetric"
      #   compactly supported wavelets with maximum number of
      #   vanishing moments, for N = 4 to 10
      filter.name =
        switch(filter.number,
               NULL, NULL, NULL,
               { # filter.number  -- 4 --
                 H = c(-0.107148901418, -0.041910965125,   0.703739068656,
                       1.136658243408,   0.421234534204,  -0.140317624179,
                       -0.017824701442,   0.045570345896)
                 "Daub cmpct on least asymm N=4"
               },
               { # filter.number  -- 5 --
                 H = c(0.038654795955,   0.041746864422,  -0.055344186117,
                       0.281990696854,   1.023052966894,   0.89658164838,
                       0.023478923136,  -0.247951362613,  -0.029842499869,
                       0.027632152958)
                 "Daub cmpct on least asymm N=5"
               },
               { # filter.number  -- 6 --
                 H = c(0.021784700327,   0.004936612372,  -0.166863215412,
                       -0.068323121587,   0.694457972958,   1.113892783926,
                       0.477904371333,  -0.102724969862,  -0.029783751299,
                       0.06325056266,    0.002499922093,  -0.011031867509)
                 "Daub cmpct on least asymm N=6"
               },
               { # filter.number  -- 7 --
                 H = c(0.003792658534,  -0.001481225915,  -0.017870431651,
                       0.043155452582,   0.096014767936,  -0.070078291222,
                       0.024665659489,   0.758162601964,   1.085782709814,
                       0.408183939725,  -0.198056706807,  -0.152463871896,
                       0.005671342686,   0.014521394762)
                 "Daub cmpct on least asymm N=7"
               },
               { # filter.number  -- 8 --
                 H = c(0.002672793393,  -0.0004283943,    -0.021145686528,
                       0.005386388754,   0.069490465911,  -0.038493521263,
                       -0.073462508761,   0.515398670374,   1.099106630537,
                       0.68074534719,   -0.086653615406,  -0.202648655286,
                       0.010758611751,   0.044823623042,  -0.000766690896,
                       -0.004783458512)
                 "Daub cmpct on least asymm N=8"
               },
               { # filter.number  -- 9 --
                 H = c(0.001512487309,  -0.000669141509,  -0.014515578553,
                       0.012528896242,   0.087791251554,  -0.02578644593,
                       -0.270893783503,   0.049882830959,   0.873048407349,
                       1.015259790832,   0.337658923602,  -0.077172161097,
                       0.000825140929,   0.042744433602,  -0.016303351226,
                       -0.018769396836,   0.000876502539,   0.001981193736)
                 "Daub cmpct on least asymm N=9"
               },
               { # filter.number  -- 10 --
                 H = c(0.001089170447,   0.00013524502,   -0.01222064263,
                       -0.002072363923,   0.064950924579,   0.016418869426,
                       -0.225558972234,  -0.100240215031,   0.667071338154,
                       1.0882515305,     0.542813011213,  -0.050256540092,
                       -0.045240772218,   0.07070356755,    0.008152816799,
                       -0.028786231926,  -0.001137535314,   0.006495728375,
                       8.0661204e-05,   -0.000649589896)
                 "Daub cmpct on least asymm N=10"
               }
        ) # switch ( filter.number )

      if(is.null(filter.name))
        stop(paste("Unknown filter number (not in {4:10}) for", family,
                   "\n (Daubechies wavelets with least asymmetry...)"))
      H = H/sqrt(2)
    } # ""DaubLeAsymm" family

    H = H/constant

    # Return Value:
    list(H = H, name = filter.name, family = family,
         filter.number = filter.number)
  }


# ------------------------------------------------------------------------------
#' @export

.first.last <-
  function(LengthH, DataLength, bc = c("periodic", "symmetric"))
  {

    # A function copied from R-package wavethresh

    # FUNCTION:

    # Settings:
    bc = match.arg(bc)
    levels = log(DataLength)/log(2)
    first.last.c = matrix(0, nrow = levels + 1, ncol = 3,
                          dimnames = list(NULL, c("First", "Last", "Offset")))
    first.last.d = matrix(0, nrow = levels, ncol = 3,
                          dimnames = list(NULL, c("First", "Last", "Offset")))
    if (bc == "periodic") {
      # Periodic boundary correction
      first.last.c[, 1] = rep(0, levels + 1)
      first.last.c[, 2] = 2^(0:levels) - 1
      first.last.c[, 3] = rev(c(0, cumsum(rev(1 + first.last.c[, 2])
      )[1:levels]))
      first.last.d[, 1] = rep(0, levels)
      first.last.d[, 2] = 2^(0:(levels - 1)) - 1
      first.last.d[, 3] = rev(c(0, cumsum(rev(1 + first.last.d[, 2])
      )[1:(levels - 1)]))
      ntotal = 2 * DataLength - 1
      ntotal.d = DataLength - 1
    } else {
      # (bc == "symmetric")
      # Symmetric boundary reflection
      first.last.c[levels + 1, 1] = 0
      first.last.c[levels + 1, 2] = DataLength - 1
      first.last.c[levels + 1, 3] = 0
      ntotal = first.last.c[levels + 1, 2] - first.last.c[levels + 1, 1] + 1
      ntotal.d = 0
      for (i in levels:1) {
        first.last.c[i, 1] = trunc(0.5 * (1 - LengthH +
                                            first.last.c[i + 1, 1]))
        first.last.c[i, 2] = trunc(0.5 * first.last.c[i + 1, 2])
        first.last.c[i, 3] = first.last.c[i + 1, 3] +
          first.last.c[i + 1, 2] - first.last.c[i + 1, 1] + 1
        first.last.d[i, 1] = trunc(0.5 * (first.last.c[i + 1, 1] - 1))
        first.last.d[i, 2] = trunc(0.5 * (first.last.c[i + 1, 2] +
                                            LengthH - 2))
        if(i != levels) {
          first.last.d[i, 3] = first.last.d[i + 1, 3] +
            first.last.d[i + 1, 2] - first.last.d[i + 1, 1] + 1
        }
        ntotal = ntotal + first.last.c[i, 2] - first.last.c[i, 1] + 1
        ntotal.d = ntotal.d + first.last.d[i, 2] - first.last.d[i, 1] + 1
      }
    }
    names(ntotal) = NULL
    names(ntotal.d) = NULL

    # Return Value:
    list(first.last.c = first.last.c, ntotal = ntotal,
         first.last.d = first.last.d, ntotal.d = ntotal.d)
  }


################################################################################
# Statistical Tests:
#' @export

.beranTest <-
  function()
  {
    # A function implemented by Diethelm Wuertz

    # Result:
    ans = NA

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.rsTest <-
  function(x, q)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the statistic of the modified R/S test

    # Arguments:
    #   x - time series
    #   q - number of lags included for calculation of covariances

    # Details:
    #   significance level: 0.050,  0.10
    #   critical value:     1.747,  1.62

    # Notes:
    #   This functions uses partly code from Christoph Helwig
    #   presented in the R-help list, 2004.

    # References:
    #   Lo (1991), Long-term Memory in Stock Market Prices,
    #   Econometrica 59, 1279--1313

    xbar = mean(x)
    N = length(x)
    r = max(cumsum(x-xbar)) - min(cumsum(x-xbar))

    covar = NULL
    for (i in 1:q) {
      covar = c(covar, sum((x[1:(N-i)]-xbar)*(x[(1+i):N]-xbar)))
    }

    if (q > 0) {
      s = sum((x-xbar)^2)/N + sum((1-(1:q)/(q+1))*covar)*2/N
    } else {
      s = sum((x-xbar)^2)/N
    }

    rs = r/(sqrt(s)*sqrt(N))
    method = "R/S Test for Long Memory"
    names(rs) = "R/S Statistic"
    names(q) = "Bandwidth q"

    # Result:
    ans = structure(list(statistic = rs, parameter = q, method = method,
                         data.name = deparse(substitute(x))), class = "htest")

    # Return Value:
    ans
  }



# ------------------------------------------------------------------------------
#' @export

.vsTest <-
  function(x, q)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the statistic of the modified V/S test

    # Arguments:
    #   x - time series
    #   q - number of lags included for calculation of covariances

    # Details:
    #   significance level: 0.01,   0.05,   0.10
    #   critical value:     0.2685, 0.1869, 0.1518

    # Notes:
    #   This functions uses partly code from Christoph Helwig
    #   presented in the R-help list, 2004.

    # References:
    #   Giraitis, Kokoszka und Leipus (2000), Rescaled variance
    #   and related tests for long memory in volatility and levels

    xbar = mean(x)
    N = length(x)
    v = sum((cumsum(x-xbar))^2) - (sum(cumsum(x-xbar)))^2/N
    covar = NULL

    for (i in 1:q) {
      covar = c(covar,
                sum((x[1:(N-i)]-xbar)*(x[(1+i):N]-xbar)))
    }

    if (q > 0) {
      s = sum((x-xbar)^2)/N + sum((1-(1:q)/(q+1))*covar)*2/N
    } else {
      s = sum((x-xbar)^2)/N
    }

    vs = v/(s*N^2)
    method = "V/S Test for Long Memory"
    names(vs) = "V/S Statistic"
    names(q) = "Bandwidth q"

    # Result:
    ans = structure(list(statistic = vs, parameter = q, method = method,
                         data.name = deparse(substitute(x))), class = "htest")

    # Return Value:
    ans
  }


################################################################################
# Slider:


## .xHurst = NA

# ------------------------------------------------------------------------------
#' @export

hurstSlider <-
  function(x = fgnSim())
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively Hurst exponent estimates

    # Arguments:
    #   x - a numerical vector or any other object which can
    #       transformed into a numeric vector.

    # FUNCTION:

    # Transform and Save Series:
    .xHurst <- as.vector(x)

    # Internal Function:
    refresh.code = function(...)
    {
      # Sliders:
      method = .sliderMenu(no = 1)
      levels = .sliderMenu(no = 2)
      minnpts = .sliderMenu(no = 3)
      lower  = .sliderMenu(no = 4)
      range = .sliderMenu(no = 5)

      # Graph Frame:
      par(mfrow = c(1, 1))

      # Plot:
      description = paste("Method", method, description())
      if (method == 1) ans = aggvarFit(x = .xHurst, levels = levels,
                                       minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                       doplot = TRUE, description = description)

      if (method == 2) ans = diffvarFit(x = .xHurst, levels = levels,
                                        minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                        doplot = TRUE, description = description)

      if (method == 3) ans = absvalFit(x = .xHurst, levels = levels,
                                       minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                       doplot = TRUE, description = description)

      if (method == 4) ans = higuchiFit(x = .xHurst, levels = levels,
                                        minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                        doplot = TRUE, description = description)

      if (method == 5) ans = pengFit(x = .xHurst, levels = levels,
                                     minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                     method = "mean", doplot = TRUE, description = description)

      if (method == 6) ans = rsFit(x = .xHurst, levels = levels,
                                   minnpts = minnpts, cut.off = 10^c(lower, lower+range),
                                   doplot = TRUE, description = description)

      if (method >= 7) ans = perFit(x = .xHurst, cut.off = lower,
                                    doplot = TRUE, description = description)

      # Add Legend:
      show(ans)
      if (method == 7) {
        mtext(text = paste(
          "cut.off = ", lower,
          sep = ""), line = -1.5, side = 3, cex = 0.8)
      } else {
        mtext(text = paste(
          "levels = ", levels, " | ",
          "minnpts = ", minnpts, " | ",
          "cut.off = 10^[", lower, ", ", lower+range, "]",
          sep = ""), line = -1.5, side = 3, cex = 0.8)
      }
      what = paste("Method:",
                   "1:aggvar | 2:diffvar | 3:absval | 4:higuchi | 5:peng | 6:rs | 7:per")
      mtext(what, side = 4, cex = 0.55, line = 0.9, adj = 0,
            col = "steelblue")

      # Reset Frame:
      par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
                names       = c("method", "levels", "minnpts", "lower", "range"),
                minima      = c(       1,       10,         1,     0.1,     0.1),
                maxima      = c(       7,      200,        10,     1.5,     3.0),
                resolutions = c(       1,        5,         1,     0.1,     0.1),
                starts      = c(       1,       50,         3,     0.7,     1.8))
  }


################################################################################


