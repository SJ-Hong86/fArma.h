#' @title Statistics of the True ARMA Process
#'
#' @description
#' A collection and description of functions to compute statistics of
#' a true ARMA time series process.
#'
#' The Functions are:
#'
#' `armaRoots`	Roots of the characteristic ARMA polynomial.
#' `armaTrueacf`	True autocorrelation function of an ARMA process.
#'
#' @usage
#' armaRoots(coefficients, n.plot = 400, digits = 4, ...)
#' armaTrueacf(model, lag.max = 20, type = c("correlation", "partial", "both"),
#'             doplot = TRUE)
#'
#' @param coefficients [armaRoots] -
#'   a numeric vector with the coefficients of the characterisitic polynomial.
#' @param digits [armaRoots] -
#'   output precision, an integer value.
#' @param doplot [armaRoots] -
#'   a logical. Should a plot be displayed?
#' @param lag.max [armaTrueacf] -
#'   maximum number of lags at which to calculate the acf or pacf, an
#'   integer value by default 20.
#' @param model [armaTrueacf] -
#'   a specification of the ARMA model with two elements: `model$ar` is the
#'   vector of the AR coefficients, and `model$ma` is the vector of the MA
#'   coefficients.
#' @param n [armaSim] -
#'   an integer value setting the length of the series to be simulated
#'   (optional if `innov` is provided). The default value is 100.
#' @param n.plot [armaRoots] -
#'   the number of data points to plot the unit circle; an integer value.
#' @param type [armaTrueacf] -
#'   a character string, "correlation" to compute the true autocorrelation
#'   function, "partial" to compute the true partial autocorrelation function,
#'   or "both" if both functions are desired. The start of one of the strings
#'   will suffice.
#' @param ... additional arguments to be passed.
#'
#' @return `armaRoots` and returns a three column data frame with the real,
#'   the imaginary part and the radius of the roots.
#'   The number of rows corresponds to the coefficients.
#'
#'   `armaTrueacf` returns a two column data frame with the lag and
#'   the correlation function.
#'
#' @export
#' @examples
#'
#' # armaRoots - Calculate and plot the roots of an ARMA process:
#' par(mfrow = c(2, 2), cex = 0.7)
#' coefficients = c(-0.5, 0.9, -0.1, -0.5)
#' armaRoots(coefficients)
#'
#' # armaTrueacf
#' model = list(ar = c(0.3, +0.3), ma = 0.1)
#' armaTrueacf(model)
#' model = list(ar = c(0.3, -0.3), ma = 0.1)
#' armaTrueacf(model)


armaTrueacf <-
  function(model, lag.max = 20, type = c("correlation", "partial", "both"),
           doplot = TRUE)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   A synonyme to ARMAacf

    # Notes:
    #   A synonyme for arma.tacf under R. See R's .First.lib.
    #   Implemented from ARMAacf

    # FUNCTION:

    # Settings:
    lag <- 0:lag.max
    result <- NA

    # Partial:
    if (type == "partial" || type == "both") {
      main = ylab = "True PACF"
      lag = 1:lag.max
      pacf = ARMAacf(model$ar, model$ma, lag.max = lag.max, pacf = TRUE)
      result = data.frame(cbind(lag, pacf))
      if (doplot) {
        plot(x = lag, y = pacf, type = "n", xlab = "Lag",
             ylab = ylab, main = main,
             ylim = c(min(c(pacf, 0)), 1) )
        lines(x = lag, y = pacf, type = "h")
        abline(h = 0)
      }
    }

    # Correlation:
    if (type == "correlation" || type == "both") {
      main = ylab = "True ACF"
      lag = 0:lag.max
      acf = ARMAacf(model$ar, model$ma, lag.max = lag.max, pacf = FALSE)
      result = data.frame(cbind(lag, acf))
      if (doplot) {
        plot(x=lag, y = acf, type = "n", xlab = "Lag",
             ylab = ylab, main = main,
             ylim = c(min(c(acf, 0)), 1) )
        lines(x = lag, y = acf, type = "h")
        abline(h = 0)
      }
    }

    # Return Value:
    result
  }


# ------------------------------------------------------------------------------
#' @export

armaRoots =
  function(coefficients, n.plot = 400, digits = 4, ...)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the roots of a characteristc polynomial

    # FUNCTION:

    # Algorithm:
    root = polyroot(c(1, -coefficients))
    real.root = Re(root)
    im.root = Im(root)
    xrange = range(real.root)
    xrange = c(xrange[1] - 1.2*abs(xrange[1]),
               xrange[2]+1.2 * abs(xrange[2]))
    xplot = seq(xrange[1], xrange[2], length = n.plot)
    fpoly = 1
    for(i in 1:length(coefficients)) {
      fpoly = fpoly - xplot^i * coefficients[i]
    }
    plot(xplot, fpoly, type = "l", xlab = "B", ylab = "Function",
         col = "steelblue", pch = 19, ...)
    title(main = "Polynomial Function vs. B")
    abline(h = 0)
    distance = sqrt(real.root^2 + im.root^2)
    root.mat = cbind(round(real.root, digits = digits),
                     round(im.root, digits = digits),
                     round(distance, digits = digits))
    dimnames(root.mat) = list(1:nrow(root.mat), c("re", "im", "dist"))
    size.limit = max(abs(real.root), 1.5, abs(im.root))
    plot(root, xlim = c( - size.limit, size.limit),
         ylim = c( - size.limit, size.limit), xlab = "", ylab = "",
         col = "steelblue", pch = 19, ...)
    x = (2*pi/360)*(0:360)
    # symbols(0, 0, circles = 1, add = TRUE, inches = FALSE, col = 6)
    lines(sin(x), cos(x))
    abline(h = 0)
    abline(v = 0)
    title("Roots and Unit Circle",
          xlab = "Real Part", ylab = "Imaginary Part")
    result = root.mat

    # Return Value:
    data.frame(result)
  }


# ------------------------------------------------------------------------------
#' @export

.armaToeplitz =
  function(x)
  {

    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x - a vector of autocovariances. The returned Toeplitz matrix is
    #       the corresponding covariance matrix of the observatons.

    # FUNCTION:

    # Wraps:
    .toeplitzARMA(x)
  }


# ------------------------------------------------------------------------------
#' @export

.armaFischer =
  function(model = list(ar = c(0.5, -0.5), ma = 0.1))
  {
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Wraps:
    .iARMA(phi = model$ar, theta = model$ma)
  }



# -----------------------------------------------------------------------------
#' @export

.schurTest <-
  function(phi)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Tests for invertibility

    # Details:
    #   Tests if all roots of
    #   1 - phi[1] B - ... - phi[length(phi)] B^length(phi) = 0
    #   are outside the unit circle.

    # References:
    #   Pagano M., (1973);
    #       When is an autoregressive process stationary?
    #       Commun. in Statist. 1, pp. 533-544.
    #   McLeod I., (1974);
    #       Derivation of the Theoretical Autocovariance Function of
    #       Autoregressive-Moving Average Time Series
    #       Appl. Statist. 24, pp. 255--256.

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Case 1 - Return Value:
    if (length(phi) == 0) {
      return(TRUE)
    }

    # Case 2 - Return Value:
    p = length(phi)
    phi = c(-1, phi)
    A = matrix(numeric(p^2), nrow = p, ncol = p)
    for (j in 1:p) {
      for (i in 1:p) {
        if (j > i) {
          A[i, j] = A[j, i]
        } else {
          k = 1:min(i, j)
          A[i, j] = sum(phi[1 + i - k] * phi[1 + j - k] -
                          phi[1 + p + k - i] * phi[1 + p + k - j])
        }
      }
    }

    if (dim(A)[1] == attr(chol(A, pivot = TRUE), "rank")) {
      return(TRUE)
    }

    # Case 3 - Return Value:
    return(FALSE)
  }


# -----------------------------------------------------------------------------
#' @export

.toeplitzARMA =
  function(x)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Toeplitz Matrix

    # Details:
    #   Given a vector x, the Toeplitz matrix is a square matrix of
    #   order length(x) and with [i,j] entries given by x[abs(i-j)].
    #   If x is a vector of autocovariances, the Toeplitz matrix is
    #   the corresponding covariance matrix of the observatons.

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Compute:
    ans = matrix(x[1 + abs(outer(seq(along = x), seq(along = x),
                                 FUN = "-"))], byrow = TRUE, ncol = length(x))

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.iARMA =
  function(phi = numeric(0), theta = numeric(0))
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Information Matrix of an ARMA Process

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Check:
    if (!(.schurTest(phi) & .schurTest(theta))) {
      cat("Model is non-causal or non-invertible\n")
      return(NULL)
    }

    # Compute:
    p = length(phi)
    q = length(theta)
    unames = vnames = character(0)
    if (p > 0) {
      if (p > 1) {
        vvmatrix = (.tccfAR(phi, phi)[ - (1:(p - 1))])[ - (p + 1)] }
      else if (p == 1) {
        vvmatrix = .tccfAR(phi, phi)[ - (p + 1)] }
      vvmatrix = .toeplitzARMA(vvmatrix)
      imatrix = vvmatrix
      vnames = paste("phi(", 1:p, ")", sep = "")
    }
    if (q > 0) {
      if (q > 1) {
        uumatrix = (.tccfAR(theta, theta)[ - (1:(q - 1))])[ - ( q + 1)]
      } else if (q == 1) {
        uumatrix = .tccfAR(theta, theta)[ - (q + 1)]
      }
      uumatrix = .toeplitzARMA(uumatrix)
      imatrix = uumatrix
      unames = paste("theta(", 1:q, ")", sep = "")
    }
    if (p > 0 && q > 0) {
      uvmatrix = matrix(numeric(1), nrow = p, ncol = q)
      tuv =  -.tccfAR(phi, theta)
      for (i in 1:p) {
        for (j in 1:q) {
          uvmatrix[i, j] = tuv[q + i - j]
        }
      }
      imatrix = cbind(rbind(vvmatrix, t(uvmatrix)), rbind(uvmatrix,
                                                          uumatrix))
    }
    inames = c(vnames, unames)
    dimnames(imatrix) = list(inames, inames)

    # Return Value:
    imatrix
  }


# ------------------------------------------------------------------------------
#' @export

.iFARMA =
  function(phi = numeric(0), theta = numeric(0), maxlag = 128)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Information Matrix of a Fractional ARMA Process

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Internal Functions:
    .psiwtsAR = function(phi, maxlag){
      p = length(phi)
      x = numeric(maxlag + 1)
      x = 1
      for (i in 1:p) {
        x[i + 1] = crossprod(phi[1:i], (rev(x))[1:i])
      }
      if (maxlag > p) {
        for (i in (p + 1):maxlag) {
          x[i + 1] = crossprod(phi, (rev(x))[1:p])
        }
      }
      x
    }

    .jFARMA = function(theta, maxlag) {
      psis = .psiwtsAR(theta, maxlag = maxlag)
      q = length(theta)
      J = numeric(q)
      for (k in 1:q) {
        J[k] = sum(psis/(k + 0:maxlag))
      }
      J
    }

    # Check
    if (!(.schurTest(phi) & .schurTest(theta))) {
      cat("Model is non-causal or non-invertible\n")
      return(NULL)
    }

    # Compute:
    I22 = (pi^2)/6
    if ((length(phi) == 0) && (length(theta) == 0)) return(I22)
    I11 = .iARMA(phi = phi, theta = theta)
    J11 = numeric(0)
    if (length(phi) > 0) J11 = -.jFARMA(phi, maxlag)
    J12 = numeric(0)
    if (length(theta) > 0) J12 = .jFARMA(theta, maxlag)
    J = c(J11, J12)
    I = rbind(I11, J)
    J = c(J, I22)
    I = cbind(I, J)
    inames = c(dimnames(I11)[[1]], "d")
    dimnames(I) = list(inames, inames)

    # Return Value:
    I
  }


# ------------------------------------------------------------------------------
#' @export

.psiwtsARMA <-
  function(phi = numeric(0), theta = numeric(0), maxlag = 20)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes MA expansion coefficients

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Compute:
    r = max((p = length(phi)), (q = length(theta)))
    phi2 = theta2 = numeric(r)
    maxlagp1 = 1 + maxlag
    if (q > 0) theta2[1:q] = theta
    if (p > 0) phi2[1:p] = phi
    x = numeric(maxlagp1)
    x[1] = 1
    if (r == 0) return(x[1:maxlagp1])
    for (i in 1:r) {
      x[i + 1] = crossprod(phi2[1:i], rev(x[1:i])) - theta2[i]
    }
    if (p == 0) return(x[1:maxlagp1])
    if (maxlag > r) {
      for (i in (r + 1):maxlag) {
        x[i + 1] = crossprod(phi, x[i - (0:(p - 1))])
      }
    }
    ans = x[1:maxlagp1]

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.tacvfARMA <-
  function(phi = numeric(0), theta = numeric(0), maxlag = 20)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes theoretical autocovariance function of an ARMA
    #   Process

    # Reference:
    #   McLeod A.I. (1975);
    #   Derivation of the theoretical autocovariance function of
    #   autoregressive-moving average models,
    #   Applied Statistics 24, pp. 255-256.

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Check:
    if (!(.schurTest(phi) & .schurTest(theta))) {
      cat("Model is non-causal or non-invertible\n")
      return(NULL)
    }

    # Compute:
    p = length(phi)
    q = length(theta)
    maxlagp1 = maxlag + 1
    if (max(p, q) == 0) {
      return(c(1, numeric(maxlagp1)))
    }
    r = max(p, q) + 1
    b = numeric(r)
    C = numeric(q + 1)
    C[1] = 1
    theta2 = c(-1, theta)
    phi2 = numeric(3 * r)
    phi2[r] = -1
    if (p > 0) {
      phi2[r + 1:p] = phi
    }
    if (q > 0) {
      for (k in 1:q) {
        C[k + 1] =  - theta[k]
        if (p > 0) {
          for (i in 1:min(p, k)) {
            C[k + 1] = C[k + 1] + phi[i] * C[k + 1 - i]
          }
        }
      }
    }
    for (k in 0:q) {
      for (i in k:q) {
        b[k + 1] = b[k + 1] - theta2[i + 1] * C[i - k + 1]
      }
    }
    if (p == 0) {
      g = c(b, numeric(maxlagp1))[1:maxlagp1]
      return(g)
    } else if (p > 0) {
      a = matrix(numeric(r^2), ncol = r)
      for (i in 1:r) {
        for (j in 1:r) {
          if (j == 1) {
            a[i, j] = phi2[r + i - 1]
          } else if (j != 1) {
            a[i, j] = phi2[r + i - j] + phi2[r + i + j - 2]
          }
        }
      }
      g = solve(a, -b)
      if (length(g) <= maxlag) {
        g = c(g, numeric(maxlag - r))
        for (i in (r + 1):maxlagp1) {
          g[i] = phi %*% g[i - 1:p]
        }
        return(g)
      } else if (length(g) >= maxlagp1) {
        return(g[1:maxlagp1])
      }
    }
  }


# ------------------------------------------------------------------------------
#' @export

.tccfAR <-
  function(phi, theta)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the theoretical cross-covariance function of two
    #   autoregressions
    #   z[t] - phi[1] z_[t-1] --- phi[p] z[t-p]     = a[t]
    #   z[t] - theta[1] z_[t-1] --- theta[q] z[t-q] = a[t]
    #   where p, q are length(phi), length(theta)

    # Notes:
    #   Auxilary function used with iarma

    # Author:
    #   Original Version from "iarma" R library: A.I. McLeod, July 1998

    # FUNCTION:

    # Compute:
    p = length(phi)
    q = length(theta)
    if (p == 0 || q == 0) {
      return(numeric(0))
    }
    k = p + q
    rhs = c(-1, rep(0, k - 1))
    A = matrix(numeric(k^2), nrow = k, ncol = k)
    for (i in 1:k) {
      for (j in 1:k) {
        imj = i - j
        ijq = i + j - q - 1
        if (i > q) {
          if (i > j && imj <= q) {
            A[i, j] = theta[imj]
          } else if (i > q && imj == 0) A[i, j] = -1
        } else {
          if (ijq > 0 && ijq <= p) {
            A[i, j] = phi[ijq]
          } else if (ijq == 0)
            A[i, j] = -1
        }
      }
    }
    ans <- solve(A, rhs)

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

test.armaTrueacf =
  function()
  {
    # armaTrueacf: Returns True ARMA autocorrelation function

    # armaTrueacf(model, lag.max = 20, type = "correlation", doplot = TRUE)
    model = list(ar = c(0.5, -0.5))
    armaTrueacf(model, lag.max = 10)

    # Return Value:
    return()
  }


# ------------------------------------------------------------------------------
#' @export

test.armaRoots =
  function()
  {
    # armaRoots:   Returns Roots of the ARMA characteristic polynomial

    # armaRoots(coefficients, n.plot = 400, digits = 4, ...)
    coefficients = c(-0.5, 0.9, -0.1, -0.5)
    ans = armaRoots(coefficients)
    target = round(sum(ans), 2)
    checkSum = 4.58
    checkEqualsNumeric(target, checkSum)

    # Return Value:
    return()
  }


################################################################################










