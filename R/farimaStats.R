#' @title Statistics of the True LDR Process
#'
#' @description
#' A collection and description of functions to investigate the true statistics of the long
#' range dependence or long memory behavior of an univariate time series process.
#'
#' The Functions are:
#'
#' Functions to model the true autocorrelation function and spectrum:
#'
#' `fgnTrueacf`	Returns true FGN covariances.
#' `fgnTruefft`	returns true FGN fast Fourier transform.
#' `fgnStatsSlider`	returns a plot of true FGN Statistics.
#' `farimaTrueacf`	returns true FARIMA covariances.
#' `farimaTruefft`	returns true FARIMA fast Fourier transform.
#' `farimaStatsSlider`	returns a plot of true FARIMA Statistics.
#'
#' @usage
#' fgnTrueacf(n = 100, H = 0.7)
#' fgnTruefft(n = 100, H = 0.7)
#' fgnStatsSlider()
#'
#' farimaTrueacf(n = 100, H = 0.7)
#' farimaTruefft(n = 100, H = 0.7)farimaTruefft(n = 100, H = 0.7)
#' farimaStatsSlider()
#'
#' @param H the Hurst exponent, a numeric value between 0.5 and 1, by default 0.7.
#' @param n number of data points to be simulated, a numeric value, by default 100.
#'
#' @details
#' Functions to Model True Correlations and Spectrum:
#'
#' The functions `fgnTrueacf` and `farimaTrueacf` return the true covariances of
#' an FGN and Gaussian FARIMA(0,d,0) time series process.
#'
#' The functions `fgnTruefft` and `farimaTruefft` return the true fast Fourier
#' transform of an FGN and Gaussian FARIMA(0,d,0) time series process.
#'
#' @return `fgnTrueacf` and `farimaTrueacf` return the true covariance of
#'   an FGN or FARIMA time series process.
#'
#'   `fgnTruefft` and `farimaTruefft` return the true spectrum of an FGN
#'   or FARIMA time series process.
#'
#' @export
#' @examples
#'
#' # fgnTrueacf
#' fgnTrueacf(n = 20, H = 0.8)
#'
#' # fgnTruefft
#' fgnTruefft(n = 20, H = 0.8)
#'
#' # farimaTrueacf
#' farimaTrueacf(n = 20, H = 0.8)
#'
#' # farimaTruefft
#' farimaTruefft(n = 20, H = 0.8)



farimaTrueacf <-
  function(n = 100, H = 0.7)
  {
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # ACF:
    ans = .ckFARIMA0(n = n, H = H)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

farimaTruefft =
  function(n = 100, H = 0.7)
  {   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # FFT:
    ans <- .gkFARIMA0(n = n, H = H)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.ckFARIMA0 <-
  function(n, H)
  {
    # Description:
    #   Computes the covariances of a fractional ARIMA(0,d,0) process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Covariances up to lag n-1

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Covariances:
    # result = (0:(n-1))
    # k = 1:(n-1)
    # d = H - 0.5
    # result[1] = gamma(1-2*d)/gamma(1-d)**2
    # result[k+1] = result[1]*(k**(2*H-2))*gamma(1-d)/gamma(d)
    # result[k+1] = result[1]*gamma(k+d)* gamma(1-d)/(gamma(k-d+1)*gamma(d))
    # ans = drop(result)

    # Covariances:
    res = numeric(n)
    d = H - 0.5
    g1d = gamma(1 - d)
    gd = pi/(sin(pi * d) * g1d)
    res[1] = gamma(1 - 2 * d)/g1d^2
    k = 1:min(50, n - 1)
    res[k + 1] = res[1] * gamma(k + d) * g1d/(gamma(k - d + 1) * gd)
    if (n > 51) {
      k <- 51:(n - 1)
      res[k + 1] <- res[1] * g1d/gd * k^(2 * H - 2)
    }

    # Return Value:
    res
  }


# ------------------------------------------------------------------------------
#' @export

farimaStatsSlider <-
  function()
  {
    # A function implemented by Diethelm Wuertz

    # Description
    #   Displays farima true Statistics: ACF and FFT

    # Example:
    #   farimaStatsSlider()

    # FUNCTION:

    # Internal Function:
    refresh.code <- function(...)
    {
      # Sliders:
      n <- .sliderMenu(no = 1)
      H <- .sliderMenu(no = 2)

      # Frame:
      par(mfrow = c(2, 1), cex = 0.7)

      # FGN ACF:
      ans <- farimaTrueacf(n = n, H = H)
      plot(ans, type = "h", col = "steelblue")
      title(main = "FARIMA True ACF")
      grid()
      abline(h = 0, col = "grey")

      # FGN FFT:
      ans <- farimaTruefft(n = n, H = H)
      plot(Re(ans), type = "h", col = "steelblue")
      title(main = "FARIMA True FFT")
      grid()
      abline(h=0, col = "grey")

      # Reset Frame:
      par(mfrow = c(1, 1), cex = 0.7)
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
                names =       c(  "n",    "H"),
                minima =      c(   10,   0.01),
                maxima =      c(  200,   0.99),
                resolutions = c(   10,   0.01),
                starts =      c(  100,   0.70))
  }


# ------------------------------------------------------------------------------
#' @export

.gkFARIMA0 <-
  function(n, H)
  {
    # Description:
    #   Calculates  gk=fft of V=(r(0),...,r(n-2),r(n-1),r(n-2),...,r(1)),
    #   where r = the autocovariances of a fractional ARIMA with innovation
    #   variance 0

    # Arguments:
    #   n <- length of time series
    #   H <- self-similarity parameter

    # Value:
    #   gk = Fourier transform of V at the Fourier frequencies

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # FFT:
    gammak = .ckFARIMA0(n,H)
    ind = c(0:(n - 2), (n - 1), (n - 2):1)
    gk = gammak[ind+1]
    ans = drop(fft(c(gk), inverse = TRUE))

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.simFARIMA0 <-
  function(n, H)
  {
    # Description:
    #   Simulates a series X(1),...,X(n) of a fractional ARIMA(0,d,0)
    #   process (d=H-1/2)

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Simulated series X(1),...,X(n)

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Simulate:
    z = rnorm(2*n)
    zr = z[c(1:n)]
    zi = z[c((n+1):(2*n))]
    zic = -zi
    zi[1] = 0
    zr[1] = zr[1]*sqrt(2)
    zi[n] = 0
    zr[n] = zr[n]*sqrt(2)
    zr = c(zr[c(1:n)], zr[c((n-1):2)])
    zi = c(zi[c(1:n)], zic[c((n-1):2)])
    z = complex(real = zr, imaginary = zi)
    cat("n = ", n, "h = ", H)
    gksqrt = Re(.gkFARIMA0(n, H))
    if (all(gksqrt > 0)) {
      gksqrt = sqrt(gksqrt)
      z = z*gksqrt
      z = fft(z, inverse = TRUE)
      z = 0.5*(n-1)**(-0.5)*z
      ans = drop(Re(z[c(1:n)]))
    } else {
      gksqrt = 0*gksqrt
      stop("Re(gk)-vector not positive")
    }

    # Return Value:
    ans
  }


################################################################################







