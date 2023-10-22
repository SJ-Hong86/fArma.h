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



fbmSim <-
  function(n = 100, H = 0.7, method = c("mvn", "chol", "lev", "circ", "wave"),
           waveJ = 7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulation of fractional Brownian motion by five different methods

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion

    # Details:
    #   The underlying functions were ported from SPlus code written
    #   by J.F. Couerjolly. They are documented in the reference given
    #   below.

    # Reference:
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian
    #       Motion: A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Initialization:
    method <- match.arg(method)

    # Match Function:
    fun = paste(".fbmSim.", method, sep = "")
    funFBM = match.fun(fun)

    # Simulate:
    if (method == "wave") {
      ans = funFBM(n, H, waveJ, doplot, fgn)
    } else {
      ans = funFBM(n, H, doplot, fgn)
    }

    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSim.mvn <-
  function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n }
    #   by numerical approximation of stochastic integral

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Mandelbrot B. and Van Ness,
    #       Fractional brownian motions,
    #       fractional noises and applications, SIAM Review, 10, n.4, 1968.
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Initialization:
    dB1 = rnorm(n)
    borne = trunc(n^1.5)
    dB2 = rnorm(borne)
    fBm = rep(0, n)
    CH = sqrt(gamma(2 * H + 1) * sin(pi * H))/gamma(H + 1/2)    ##

    # Main Program:
    ind1 = (1:n)^(H - 1/2)
    ind2 = (1:(borne + n))^(H - 1/2)
    for(i in (2:n)) {
      I1 = dB1[(i - 1):1] * ind1[1:(i - 1)]
      I2 = (ind2[(i + 1):(i + borne)] - ind2[1:borne]) * dB2
      fBm[i] = sum(I1) + sum(I2)
    }
    fBm = fBm * n^( - H) * CH
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "mvnFBM Path"
    if (fgn) {
      ans = c(fBm[1], diff(fBm))
      Title = "mvnFGN Path"
    }

    # Plot of fBm
    if (doplot) {
      time = 1:n
      Nchar = as.character(n)
      Nleg = paste("N=", Nchar, sep = "")
      Hchar = as.character(round(H, 3))
      Hleg = paste(", H=", Hchar, sep = "")
      NHleg = paste(c(Nleg, Hleg), collapse = "")
      leg = paste(c(Title, NHleg), collapse = " - ")
      plot(time, ans, type = "l", main = leg, col = "steelblue")
      grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "mvn", H = H)
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSim.wave =
  function(n = 1000, H = 0.7, J = 7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   J : resolution
    #   doplot : = TRUE ----> plot of path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by wavelet synthesis

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Abry P. and Sellan F.,
    #       The wavelet-based synthesis
    #       for fractional Brownian motion, Applied and computational
    #       harmonic analysis, 1996 + Matlab scripts from P. Abry
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Daubechies filter of length 20
    Db20 = c(0.026670057901000001, 0.188176800078)
    Db20 = c(Db20, 0.52720118932000004, 0.688459039454)
    Db20 = c(Db20, 0.28117234366100002, -0.24984642432699999)
    Db20 = c(Db20, -0.19594627437699999, 0.127369340336)
    Db20 = c(Db20, 0.093057364604000006, -0.071394147165999997)
    Db20 = c(Db20, -0.029457536821999999, 0.033212674058999997)
    Db20 = c(Db20, 0.0036065535670000001, -0.010733175483)
    Db20 = c(Db20, 0.001395351747, 0.0019924052950000002)
    Db20 = c(Db20, -0.00068585669500000003, -0.000116466855)
    Db20 = c(Db20, 9.3588670000000005e-05, -1.3264203000000001e-05)
    secu = 2 * length(Db20)

    # Quadrature mirror filters of Db20
    Db20qmf = (-1)^(0:19) * Db20
    Db20qmf = Db20qmf[20:1]
    nqmf = -18

    # Truncated fractional coefficients appearing in fractional integration
    # of the multiresolution analysis
    prec = 0.0060000000000000001
    hmoy = c(1, 1)
    s = H + 1/2
    d = H - 1/2
    if (H == 1/2) {
      ckbeta = c(1, 0)
      ckalpha = c(1, 0)
    } else {
      # Truncature at order prec
      ckalpha = 1
      ckbeta = 1
      ka = 2
      kb = 2
      while(abs(ckalpha[ka - 1]) > prec) {
        g = gamma(1 + d)/gamma(ka)/gamma(d + 2 - ka)
        if (is.na(g))
          g = 0
        ckalpha = c(ckalpha, g)
        ka = ka + 1
      }
      while(abs(ckbeta[kb - 1]) > prec) {
        g = gamma(kb - 1 + d)/gamma(kb)/gamma(d)
        if (is.na(g))
          g = 0
        ckbeta = c(ckbeta, g)
        kb = kb + 1
      }
    }
    lckbeta = length(ckbeta)
    lckalpha = length(ckalpha)  ##

    # Number of starting points
    nbmax = max(length(ckbeta), length(ckalpha))
    nb0 = n/(2^(J)) + 2 * secu  ##

    # Sequence fs1:
    fs1 = .convol(ckalpha, Db20)
    fs1 = .convol(fs1, hmoy)
    fs1 = 2^( - s) * fs1
    fs1 = fs1 * sqrt(2) #

    # Sequence gs1:
    gs12 = .convol(ckbeta, Db20qmf)
    gs1 = cumsum(gs12)
    gs1 = 2^(s) * gs1
    gs1 = gs1 * sqrt(2) ##

    # Initialization:
    nb1 = nb0 + nbmax
    bb = rnorm(nb1)
    b1 = .convol(bb, ckbeta)
    bh = cumsum(b1)
    bh = bh[c(nbmax:(nb0 + nbmax - 1))]
    appro = bh
    tappro = length(appro)  ##

    # Useful function:
    dilatation = function(vect) {
      # dilates one time vector vect
      ldil = 2 * length(vect) - 1
      dil = rep(0, ldil)
      dil[seq(1, ldil, by = 2)] = vect
      drop(dil)
    }

    # Synthese's algorithm:
    for(j in 0:(J - 1)) {
      appro = dilatation(appro)
      appro = .convol(appro, fs1)
      appro = appro[1:(2 * tappro)]
      detail = rnorm(tappro) * 2^(j/2) * 4^( - s) * 2^( - j * s)
      detail = dilatation(detail)
      detail = .convol(detail, gs1)
      detail = detail[( - nqmf + 1):( - nqmf + 2 * tappro)]
      appro = appro + detail
      tappro = length(appro)
    }
    debut = (tappro - n)/2
    fBm = appro[c((debut + 1):(debut + n))]
    fBm = fBm - fBm[1]
    fGn = c(fBm[1], diff(fBm))
    fGn = fGn * 2^(J * H) * n^( - H)    # path on [0,1]
    fBm = cumsum(fGn)
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "waveFBM Path"
    if (fgn) {
      ans = c(fBm[1], diff(fBm))
      Title = "waveFGN Path"
    }

    # Plot of fBM/FGN:
    if (doplot) {
      time = 1:n
      Nchar = as.character(n)
      Nleg = paste("N=", Nchar, sep = "")
      Hchar = as.character(round(H, 3))
      Hleg = paste(", H=", Hchar, sep = "")
      NHleg = paste(c(Nleg, Hleg), collapse = "")
      leg = paste(c(Title, NHleg), collapse = " - ")
      plot(time, ans, type = "l", main = leg, col = "steelblue")
      grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "wave", H = H)
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.convol <-
  function(x, y)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x,y : vectors

    # Value:
    #   convolution of vectors x and y

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # FUNCTION:

    # Convolution:
    if (missing(x) | missing(y)) {
      break
    } else {
      a = c(x, rep(0, (length(y) - 1)))
      b = c(y, rep(0, (length(x) - 1)))
      a = fft(a, inverse = FALSE)
      b = fft(b, inverse = FALSE)
      conv = a * b
      conv = Re(fft(conv, inverse = TRUE))
      conv = conv/length(conv)
      drop(conv)
    }
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSim.chol =
  function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n }
    #   by Choleki's decomposition of the covariance matrix of the fBm

    # Author:
    #   Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Construction of covariance matrix of fBm
    H2 = 2 * H
    matcov = matrix(0, n - 1, n - 1)
    for(i in (1:(n - 1))) {
      j = i:(n - 1)
      r = 0.5 * (abs(i)^H2 + abs(j)^H2 - abs(j - i)^H2)
      r = r/n^H2
      matcov[i, j] = r
      matcov[j, i] = matcov[i, j]
    }
    L = chol(matcov)
    Z = rnorm(n - 1)
    fBm = t(L) %*% Z
    fBm = c(0, fBm)

    # Result:
    ans = drop(fBm)
    Title = "cholFBM Path"
    if (fgn) {
      ans = c(fBm[1], diff(fBm))
      Title = "cholFGN Path"
    }

    # Plot of fBm:
    if (doplot) {
      time = 1:n
      Nchar = as.character(n)
      Nleg = paste("N=", Nchar, sep = "")
      Hchar = as.character(round(H, 3))
      Hleg = paste(", H=", Hchar, sep = "")
      NHleg = paste(c(Nleg, Hleg), collapse = "")
      leg = paste(c(Title, NHleg), collapse = " - ")
      plot(time, ans, type = "l", main = leg, col = "steelblue")
      grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "chol", H = H)
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSim.lev <-
  function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   plotfBm :  =1 ---> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by Levinson's method

    # Author:
    #   Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Peltier R.F.,
    #       Processus stochastiques fractals avec
    #       applications en finance, these de doctorat, p.42, 28.12.1997
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Covariances of fGn:
    k = 0:(n - 1)
    H2 = 2 * H
    r = (abs((k - 1)/n)^H2 - 2 * (k/n)^H2 + ((k + 1)/n)^H2)/2

    # Initialization of algorithm:
    y = rnorm(n)
    fGn = rep(0, n)
    v1 = r
    v2 = c(0, r[c(2:n)], 0)
    k =  - v2[2]
    aa = sqrt(r[1]) #

    # Levinson's algorithm:
    for(j in (2:n)) {
      aa = aa * sqrt(1 - k * k)
      v = k * v2[c(j:n)] + v1[c((j - 1):(n - 1))]
      v2[c(j:n)] = v2[c(j:n)] + k * v1[c((j - 1):(n - 1))]
      v1[c(j:n)] = v
      bb = y[j]/aa
      fGn[c(j:n)] = fGn[c(j:n)] + bb * v1[c(j:n)]
      k =  - v2[j + 1]/(aa * aa)
    }
    fBm = cumsum(fGn)
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "levFBM Path"
    if (fgn) {
      ans = c(fBm[1], diff(fBm))
      Title = "levFGN Path"
    }

    # Plot of fBm:
    if (doplot) {
      time = 1:n
      Nchar = as.character(n)
      Nleg = paste("N=", Nchar, sep = "")
      Hchar = as.character(round(H, 3))
      Hleg = paste(", H=", Hchar, sep = "")
      NHleg = paste(c(Nleg, Hleg), collapse = "")
      leg = paste(c(Title, NHleg), collapse = " - ")
      plot(time, ans, type = "l", main = leg, col = "steelblue")
      grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "lev", H = H)
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSim.circ <-
  function(n = 100, H = 0.7, doplot = TRUE, fgn = FALSE)
  {
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ---> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by Wood-Chan's method

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Wood A. and Chan G.,
    #       Simulation of stationnary Gaussian processes,
    #       Journal of computational and grahical statistics, Vol.3, 1994.
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # First line of the circulant matrix, C, built via covariances of fGn
    lineC =
      function(n, H, m)
      {
        k = 0:(m - 1)
        H2 = 2 * H
        v = (abs((k - 1)/n)^H2 - 2 * (k/n)^H2 + ((k + 1)/n)^H2)/2
        ind = c(0:(m/2 - 1), m/2, (m/2 - 1):1)
        v = v[ind + 1]
        drop(v)
      }

    # Next power of two > n:
    m = 2
    repeat {
      m = 2 * m
      if (m >= (n - 1)) break
    }
    stockm = m  ##

    # Research of the power of two (<2^18) such that C is definite positive:
    repeat {
      m = 2 * m
      eigenvalC = lineC(n, H, m)
      eigenvalC = fft(c(eigenvalC), inverse = FALSE)
      ### DW: That doesn't work on complex vectors !
      ### if ((all(eigenvalC > 0)) | (m > 2^17)) break
      ### We use:
      if ((all(Re(eigenvalC) > 0)) | (m > 2^17)) break
    }
    if (m > 2^17) {
      cat("----> exact method, impossible!!", fill = TRUE)
      cat("----> can't find m such that C is definite positive", fill = TRUE)
      break
    } else {
      # Simulation of W=(Q)^t Z, where Z leads N(0,I_m)
      # and   (Q)_{jk} = m^(-1/2) exp(-2i pi jk/m):
      ar = rnorm(m/2 + 1)
      ai = rnorm(m/2 + 1)
      ar[1] = sqrt(2) * ar[1]
      ar[(m/2 + 1)] = sqrt(2) * ar[(m/2 + 1)]
      ai[1] = 0
      ai[(m/2 + 1)] = 0
      ar = c(ar[c(1:(m/2 + 1))], ar[c((m/2):2)])
      aic =  - ai
      ai = c(ai[c(1:(m/2 + 1))], aic[c((m/2):2)])
      W = complex(real = ar, imaginary = ai)  ##

      # Reconstruction of the fGn:
      W = (sqrt(eigenvalC)) * W
      fGn = fft(W, inverse = FALSE)
      fGn = (1/(sqrt(2 * m))) * fGn
      fGn = Re(fGn[c(1:n)])
      fBm = cumsum(fGn)
      fBm[1] = 0

      # Result:
      ans = drop(fBm)
      Title = "circFBM Path"
      if (fgn) {
        ans = c(fBm[1], diff(fBm))
        Title = "circFGN Path"
      }

      # Plot of fBm:
      if (doplot) {
        time = 1:n
        Nchar = as.character(n)
        Nleg = paste("N=", Nchar, sep = "")
        Hchar = as.character(round(H, 3))
        Hleg = paste(", H=", Hchar, sep = "")
        NHleg = paste(c(Nleg, Hleg), collapse = "")
        leg = paste(c(Title, NHleg), collapse = " - ")
        plot(time, ans, type = "l", main = leg, col = "steelblue")
        grid()
      }
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "circ", H = H)
    ans
  }


# ------------------------------------------------------------------------------
#' @export

.fbmSlider =
  function()
  {
    # A function implemented by Diethelm Wuertz

    # Description
    #   Displays fbm simulated time Series

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
      # Sliders:
      n      <- .sliderMenu(no = 1)
      H      <- .sliderMenu(no = 2)
      method <- .sliderMenu(no = 3)

      # Select Method:
      Method = c("mvn", "chol", "lev", "circ", "wave")
      Method = Method[method]

      # Frame:
      par(mfrow = c(2, 1))

      # FBM TimeSeries:
      fbmSim(n = n, H = H, method = Method, doplot = TRUE)

      # FGN TimeSeries:
      fbmSim(n = n, H = H, method = Method, doplot = TRUE, fgn = TRUE)

      # Reset Frame:
      par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
                names =       c(  "n",    "H", "method"),
                minima =      c(   10,   0.01,       1),
                maxima =      c(  200,   0.99,       5),
                resolutions = c(   10,   0.01,       1),
                starts =      c(  100,   0.70,       1))
  }


################################################################################



