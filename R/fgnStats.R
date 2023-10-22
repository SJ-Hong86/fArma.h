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



fgnTrueacf =
  function(n =100, H = 0.7)
  {
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # True ACF:
    ans <- .ckFGN0(n = n, H = H)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

fgnTruefft <-
  function(n = 100, H = 0.7)
  {
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # True FFT:
    ans <- .gkFGN0(n = n, H = H)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

fgnStatsSlider <-
  function()
  {
    # A function implemented by Diethelm Wuertz

    # Description
    #   Displays FGN true statistics: ACF and FFT

    # Example:
    #   fgnStatsSlider()

    # FUNCTION:

    # Internal Function:
    refresh.code <- function(...)
    {
      # Sliders:
      n = .sliderMenu(no = 1)
      H = .sliderMenu(no = 2)

      # Frame:
      par(mfrow = c(2, 1), cex = 0.7)

      # FGN ACF:
      ans = fgnTrueacf(n = n, H = H)
      plot(ans, type = "h", col = "steelblue")
      title(main = "FGN True ACF")
      grid()
      abline(h = 0, col = "grey")

      # FGN FFT:
      ans = fgnTruefft(n = n, H = H)
      plot(Re(ans), type = "h", col = "steelblue")
      title(main = "FGN True FFT")
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

.ckFGN0 <-
  function(n, H)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes covariances of a fractional Gaussian  process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Covariances upto lag n-1

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:
    k = 0:(n-1)
    H2 = 2 * H
    ans = drop((abs(k-1)**H2-2*abs(k)**H2+abs(k+1)**H2)/2)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.gkFGN0 <-
  function(n, H)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates gk = fft of V=(r(0),...,r(n-2),
    #   r(n-1), r(n-2),...,r(1), where r=the autocovariances
    #   of a fractional Gaussian process with variance 1

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   gk = Fourier transform of V at Fourier frequencies

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # FFT:
    gammak = .ckFGN0(n, H)
    ind = c(0:(n - 2), (n - 1), (n - 2):1)
    gk = gammak[ind+1]
    ans = drop(fft(c(gk), inverse = TRUE))

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.simFGN0 =
  function(n, H)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a series X(1),...,X(n) of a fractional Gaussian process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   simulated series X(1),...,X(n)

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Simulation:
    z = rnorm(2*n)
    zr = z[c(1:n)]
    zi = z[c((n+1):(2*n))]
    zic = -zi
    zi[1] = 0
    zr[1] = zr[1]*sqrt(2)
    zi[n] = 0
    zr[n] = zr[n]*sqrt(2)
    zr = c(zr[c(1:n)],zr[c((n-1):2)])
    zi = c(zi[c(1:n)],zic[c((n-1):2)])
    z = complex(real = zr,imaginary = zi)
    gksqrt = Re(.gkFGN0(n, H))
    if (all(gksqrt > 0)) {
      gksqrt = sqrt(gksqrt)
      z = z*gksqrt
      z = fft(z, inverse = TRUE)
      z = 0.5*(n-1)**(-0.5)*z
      ans = drop(Re(z[c(1:n)]))
    } else {
      gksqrt = 0*gksqrt
      cat("Re(gk)-vector not positive") }

    # Return Value:
    ans
  }

# ------------------------------------------------------------------------------
#' @export

test.fgnTrueacf =
  function()
  {
    # Try:
    #   farimaStatsSlider()
    #   fgnStatsSlider()

    # FGN True ACF:
    ans = fgnTrueacf(n = 1000, H = 0.7)[1:5]
    print(ans)
    target = round(sum(ans), 2)
    print(target)
    current = 1.78
    checkEquals(target, current)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.fgnTruefft =
  function()
  {
    # FGN True FFT:
    ans = fgnTruefft(n = 1000, H = 0.7)[1:3]
    print(ans)
    target = Re(round(sum(ans), 2))
    print(target)
    current = 40.19
    checkEquals(target, current)

    # Return Value:
    return()
  }

################################################################################


