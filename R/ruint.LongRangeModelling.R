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


test.fbmSim =
    function()
    {
        # Simulation:
        
        # Try:
        #   .fbmSimSlider()
        
        # fbmSim(
        #   n = 100, H = 0.7,
        #   method = c("mvn", "chol", "lev", "circ", "wave"),
        #   waveJ = 7,
        #   doplot = TRUE, fgn = FALSE)
        
        # Graphics Frame:
        par(mfrow = c(3, 2), cex = 0.7)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fbmSim(n = 50, method = "mvn")
        print(x)
        target = round(mean(x), 2)
        print(target)
        current = +0.05
        checkEquals(target, current)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fbmSim(n = 50, method = "chol")
        print(x)
        target = round(mean(x), 2)
        print(target)
        current = +0.45
        checkEquals(target, current)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fbmSim(n = 50, method = "lev")
        print(x)
        target = round(mean(x), 2)
        print(target)
        current = 0.66
        checkEquals(target, current)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fbmSim(n = 50, method = "circ")
        print(x)
        target = round(mean(x), 2)
        print(target)
        current = -0.45
        checkEquals(target, current)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fbmSim(n = 50, method = "wave")
        print(x)
        target = round(mean(x), 2)
        print(target)
        current = 0.09
        checkEquals(target, current)
        
        # Return Value:
        return()
    }


# ------------------------------------------------------------------------------
#' @export

test.fgnSim =
    function()
    {
        # Simulation:
        
        # Try:
        #   .fgnSimSlider()
        
        # fbmSim(
        #   n = 100, H = 0.7,
        #   method = c("mvn", "chol", "lev", "circ", "wave"),
        #   waveJ = 7, doplot = TRUE, fgn = FALSE)
        
        # Graphics Frame:
        par(mfrow = c(3, 2), cex = 0.7)
        
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        
        
        # Return Value:
        return()
    }


# ------------------------------------------------------------------------------
#' @export

test.farimaSim =
    function()
    {
        # Simulation:
        
        # Try:
        #   .farimaSimSlider()
        
        # farimaSim(
        #   n = 1000,
        #   model = list(ar = c(0.5, -0.5), d = 0.3, ma = 0.1),
        #   method = c("freq", "time"),
        #   ...)
        
        # Frequency Method:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = farimaSim(n = 50, method = "freq")
        print(x)
        target = round(var(x), 2)
        print(target)
        current = 1.6
        print(current)
        checkEquals(target, current)
        
        # Time Methhod:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = farimaSim(n = 50, method = "time")
        print(x)
        target = round(var(x), 2)
        print(target)
        current = 1.6
        print(current)
        checkEquals(target, current)
        
        # Return Value:
        return()
    }


# ------------------------------------------------------------------------------
#' @export

test.whittleFit =
    function()
    {
        # whittleFit(
        #   x, order = c(1, 1), subseries = 1,
        #   method = c("fgn", "farma"),
        #   trace = FALSE, spec = FALSE, title = NULL, description = NULL)


        # Beran - Simulate:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fgnSim(n = 1000, H = 0.7)

        # Fit:
        Hurst = whittleFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.701
        checkEquals(target, current)

        # Durbin - Simulate:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fgnSim(n = 1000, H = 0.7, method = "durbin")

        # Fit:
        Hurst = whittleFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.679
        checkEquals(target, current)

        # Paxson - Simulate:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fgnSim(n = 1000, H = 0.7, method = "paxson")

        # Fit:
        Hurst = whittleFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.728
        checkEquals(target, current)

        # Return Value:
        return()
    }


# ------------------------------------------------------------------------------
#' @export

test.hurstFit =
    function()
    {
        # Try:
        #   hurstSlider(x = fbmSim(n = 1000, H = 0.7, fgn = TRUE))
        
        # Methods:
        #   method = 1: aggvarFit
        #   method = 2: diffvarFit
        #   method = 3: absvalFit
        #   method = 4: higuchiFit
        #   method = 5: pengFit
        #   method = 6: rsFit
        #   method = 7: perFit
        
        # Simulate:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fgnSim(n = 1000, H = 0.7)
        
        Hurst = aggvarFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.731
        checkEqualsNumeric(target, current)
        
        Hurst = diffvarFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.869
        checkEqualsNumeric(target, current)
        
        Hurst = absvalFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.787
        checkEqualsNumeric(target, current)
        
        Hurst = higuchiFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.714
        checkEqualsNumeric(target, current)
        
        Hurst = pengFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.620
        checkEqualsNumeric(target, current)
        
        Hurst = rsFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.594
        checkEqualsNumeric(target, current)
        
        Hurst = perFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.602
        checkEqualsNumeric(target, current)
        
        # More Estimators:
        # ...
        
        # Return Value:
        return()
    }


# ------------------------------------------------------------------------------
#' @export

test.waveletFit =
    function()
    {
        # waveletFit(
        #   x, length = NULL, order = 2, octave = c(2, 8),
        #   doplot = FALSE, title = NULL, description = NULL)
        
        # Simulate:
        RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
        set.seed(4711, kind = "Marsaglia-Multicarry")
        x = fgnSim(n = 1000, H = 0.7)
        
        # Fit:
        Hurst = waveletFit(x)@hurst$H
        print(Hurst)
        target = round(Hurst, 3)
        print(target)
        current = 0.696
        checkEquals(target, current)
        
        # Return Value:
        return()
    }


################################################################################
