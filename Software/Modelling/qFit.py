import scipy
import scipy.integrate
import scipy.optimize
import numpy as np
import pylab
import matplotlib.pyplot as plt
import copy
import os
import sys
import scipy.stats as stats


class Param:
    '''Parameter for model fitting.

    Grabbed Code from http://www.scipy.org/Cookbook/FittingData
    '''
    def __init__(self, value, name='unknown'):
        self.value = value
        self.uncertainty = None
        self.covindex = None  # the corresponding index in the cov matrix
        self.initial_value = value  # this will be unchanged after the fit
        self.name = name

    def set(self, value):
        self.value = value

    def __call__(self):
        # When paramter is called, return itself as the value
        return self.value


def fminfit(function, parameters, y, yerr, x=None, algorithm=None):
    '''
    Note: this code is not well fleshed out compared to fit() in terms of the
    object-orienty nature or returning estimates of uncertainties. Right now
    all it will solve/return is the best fit values.

    Use fmin or a variant thereof to optimize/find the minimum of a the
    chi-squared statistic
    ((y-function(x))**2/(yerr**2)).sum()

    where yerr is the known variance of the observation, y is the observed data
    and function(x) is the theoretical data.[1] This definition is only useful
    when one has estimates for the error on the measurements, but it leads to a
    situation where a chi-squared distribution can be used to test goodness of
    fit, provided that the errors can be assumed to have a normal distribution.

    if algorithm == bfgs then use the bfgs fmin fit,
        if algorithm == None use the default simplex

    from fmin guide:
        This algorithm has a long history of successful use in applications.
        But it will usually be slower than an algorithm that uses first or
        second derivative information. In practice it can have poor performance
        in high-dimensional problems and is not robust to minimizing
        complicated functions. Additionally, there currently is no complete
        theory describing when the algorithm will successfully converge to the
        minimum, or how fast it will if it does.
    '''
    def errfunc(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return ((y-function(x))**2/(yerr**2)).sum()  # chi sq

    paraminit = [param() for param in parameters]
    print "Initial Parameters: ", paraminit
    if not algorithm:
        p = scipy.optimize.fmin(errfunc, paraminit, maxfun=None,
                                maxiter=None, full_output=1)
    elif algorithm == 'bfgs':
        p = scipy.optimize.fmin_bfgs(errfunc, paraminit)
    solved_values = p[0]
    func_value = p[1]  # final function value, chi-sq in this case
    niter = p[2]  # number of iterations
    nfunc = p[3]  # number of function calls

    print "Solved Values: ", solved_values
    # don't really need to return this as the values are going to be
    # stored in the Param objects
    return solved_values


def fit(function, parameters, y, yerr, x=None, return_covar=False):
    '''Fit performs a simple least-squares fit on a function.  To use:

    Give initial paramaters:
        mu = Param(7)
        sigma = Param(3)
        height = Param(5)

    Define your function:
        def f(x): return height() * exp(-((x-mu())/sigma())**2)

    Make simulated data:
        xvals = np.arange(100)
        gaussian = lambda x: 3*exp(-(30-x)**2/20.)
        ydata = gaussian(xvals)
        ydata = scipy.randn(100)*.1+ydata #adding noise
        yerr = np.zeros(100)+.1 #array of uncertainties

    Fit the function (provided 'data' is an array with the data to fit):
        fit(f, [mu, sigma, height], ydata, yerr, xvals)

    Plot the fitted model over the data if desired
        simxvals = np.arange(10000)/100. # 10000 points from 0-100
        plot(simxvals, f(simxvals))
        scatter(xvals, ydata)

    Your input parameters will now have the fitted values assigned to them.
    '''
    # This errfunc assumes a uniform error for each datapoint!
    def errfunc(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return (y - function(x)) / yerr

    # errfunc = lambda x, y, err: (y-function(x))/err

    pr = """Fitting {0} free parameters on {1} datapoints...
         """.format(len(parameters), len(y))
    print pr

    # If no x axis given, set x array as integral steps starting with
    # zero for the length of the y array.  For instance, if
    # y = array([1.32, 2.15, 3.01, 3.92]),
    # then x will be x = array([0, 1, 2, 3])
    if x is None:
        x = np.arange(y.shape[0])  # not sure if this works
    paraminit = [param() for param in parameters]
    print 'Initial Parameters:', paraminit
    # print "Parameter list:", paraminit
    # print "error function:", errfunc
    fitout = scipy.optimize.leastsq(errfunc, paraminit, full_output=1)
    # print fitout

    # paramfinal = [param() for param in parameters]
    paramfinal = fitout[0]
    covarmatrix = fitout[1]
    info = fitout[2]
    mesg = fitout[3]
    errint = fitout[4]
    # errint of 1-4 means a solution was found
    if errint not in np.arange(1, 5):
        raise Exception(mesg)
    print info['nfev'], ' function calls required to find solution'
    # print errint
    # print mesg
    print 'Final Parameters:', paramfinal
    print 'Covariance Matrix', covarmatrix
    print ''
    # If paramfinal is not an array, make it one to avoid scripting errors
    if not isinstance(paramfinal, np.ndarray):
        paramfinal = np.array([paramfinal])

    # Calculate the chi-squared value, and the reduced chi-squared
    # http://mail.scipy.org/pipermail/scipy-user/2005-June/004632.html
    chi2 = sum(np.power(errfunc(paramfinal), 2))
    degrees_of_freedom = y.shape[0] - len(paramfinal)
    chi2r = chi2/degrees_of_freedom

    print "chi^2 / dof = %.2f / %i" % (chi2, degrees_of_freedom)
    print "reduced chi^2 = %.3f  \n" % (chi2r)

    retdict = {'parameters': parameters, 'covarmatrix': covarmatrix,
               'chi2': chi2, 'dof': degrees_of_freedom
               }

    count = 0
    fitstrlist = []
    for param in retdict['parameters']:
        param.covindex = count  # corresponding index in the covmatrix
        uncertainty = np.sqrt(retdict['covarmatrix'].diagonal()[count])
        # update the uncertainty in the param object
        param.uncertainty = uncertainty
        fitstr = '%s: %.3f +/- %.3f' % (param.name, param.value, uncertainty)
        print fitstr
        fitstrlist.append(fitstr)
        count += 1
    retdict.update({'strings': fitstrlist})

    if return_covar:
        return covarmatrix
    else:
        return retdict  # return the dictonary of outputs


def plot_marg_from_fitdict(fitdict, paramnames):
    '''Given a fit dictionary from fit(), and a tuple of parameter names (from
    param.name), get the covariance matrix and plot the marginalization e.g.
    paramnames = ('Av_1','beta_1')
    '''
    allvalues = np.zeros(len(fitdict['parameters']))
    indices = [-1, -1]
    values = [0, 0]
    names = ['a1', 'a2']
    covmat = fitdict['covarmatrix']
    count = 0
    for param in fitdict['parameters']:
        allvalues[count] = param.value
        count += 1
        if paramnames[0] == param.name:
            indices[0] = param.covindex
            names[0] = param.name
            values[0] = param.value
        if paramnames[1] == param.name:
            indices[1] = param.covindex
            names[1] = param.name
            values[1] = param.value

    # print indices
    # return allvalues
    ret = plot_marginalization(covmat=covmat,
                               indices=indices,
                               names=names,
                               values=values)
    return ret


def plot_marginalization(covmat=None, indices=None,
                         names=None, values=None,
                         plot_delta_values=False,
                         invert_yticks=True,
                         storepath='./'):
    '''Suppose we dont care much about some parameters and want to explore
    the uncertainties involved in just two. We can marginalize over the
    other parameters by first extracting the relevant values from the
    covariance matrix above to form a new 2x2 covariance matrix:
    '''
    if covmat is None:  # default just for illustration
        covmat = np.matrix([
            [5.29626719, 0.57454987, -0.73125854],
            [0.57454987, 1.16079146, -0.28095744],
            [-0.73125854, -0.28095744, 0.23075755]])

    if indices is None:  # default for illustration
        ind1 = 1
        ind2 = 2
    else:
        ind1 = indices[0]
        ind2 = indices[1]
    if names is None:  # default for names
        names = ['a1', 'a2']

    unc_1 = np.sqrt(covmat[ind1, ind1])
    unc_2 = np.sqrt(covmat[ind2, ind2])

    # slice the covariance matrix
    cov_i = np.matrix([[covmat[ind1, ind1], covmat[ind1, ind2]],
                       [covmat[ind2, ind1], covmat[ind2, ind2]]])

    # And we invert this to get a new curvature matrix:
    curv_i = cov_i.getI()

    # Now from this curvature matrix we can write (c.f. Eq.
    # 9.3 of Aficionados):
    # \[
    # \Delta \chi^2_{\mathbf{a_i}} = \delta \mathbf{a_i^T} \cdot
    #               [\alpha_\chi]_\mathbf{i} \cdot \delta \mathbf{a_i},
    # \]
    # where $\delta \mathbf{a_i^T} = [\delta a_1 \; \delta a_2]$.
    #
    # To visualize these $\chi^2$ values, we create a 256 by 256 grid of
    #  $\delta a_1, \delta a_2$ values and calculate the $\chi^2$ for each.

    scale = max((unc_1, unc_2))*3.5

    dx = np.linspace(-1*scale, scale, 256)
    dy = np.linspace(-1*scale, scale, 256)
    delta_chi_sq = np.zeros((256, 256))

    # calculate grid of delta_chi_squared:
    x_ind = 0
    for dx_i in dx:
        y_ind = 0
        for dy_j in dy:
            delta_a = np.matrix([[dx_i], [dy_j]])
            delta_chi_sq[x_ind, y_ind] = delta_a.getT()*curv_i*delta_a
            y_ind += 1
        x_ind += 1

    # Now with this data grid, we can plot the contours corresponding to
    # $\Delta \chi^2_{\mathbf{a_i}} = 1.0$ and
    # $\Delta \chi^2_{\mathbf{a_i}} = 2.3.$
    # The latter value is where 68.3\% of measured values should lie inside for
    # two degrees of freedom.
    # $\Delta \chi^2_{\mathbf{a_i}} = 6.17.$ is where 95.4% of measured values
    # should lie for 2 degrees of freedom (corresponding to 4sigma for 1dof)
    # 9.21 is 99% confidence
    levels = np.array([1.0, 2.3, 9, 11.8])

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    CS = ax.contour(dx, dy, delta_chi_sq, levels)
    ax.clabel(CS, inline=1, fontsize=10)
    # label these axes only if desired; could add confusion
    if plot_delta_values:
        yname = '$\delta %s$' % names[0]
        xname = '$\delta %s$' % names[1]
        ax.set_ylabel(yname)
        ax.set_xlabel(xname)
    else:  # get rid of the tickmarks
        ax.set_xticks([])
        ax.set_yticks([])

    # marginalize over two parameters and plot the corresponding lines;
    # these values are equivalent to the uncertainties from the diagonals
    # of the original covariance matrix above.  It can be seen that they bound
    # the error ellipse corresponding to$\Delta \chi^2 = 1.0,$ as they should
    # (hence the term 'marginalization' - this contour, projected into the
    # margins, gives the uncertainty for a single parameter of interest).

    ax.axvline(x=0, linestyle='dotted', color='grey')
    ax.axhline(y=0, linestyle='dotted', color='grey')

    # 1 sigma
    ax.axvline(x=unc_2, linestyle='dashed')
    ax.axvline(x=-1*unc_2, linestyle='dashed')

    ax.axhline(y=unc_1, linestyle='dashed')
    ax.axhline(y=-1*unc_1, linestyle='dashed')
    # 3 sigma
    ax.axvline(x=3*unc_2, linestyle='dashed', color='orange')
    ax.axvline(x=-3*unc_2, linestyle='dashed', color='orange')

    ax.axhline(y=3*unc_1, linestyle='dashed', color='orange')
    ax.axhline(y=-3*unc_1, linestyle='dashed', color='orange')

    if not plot_delta_values:  # get rid of the tickmarks
        ax.set_xticks([])
        ax.set_yticks([])

    if values is not None:
        xlim = ax.get_xlim() + values[1]
        ylim = ax.get_ylim() + values[0]
        ax2 = ax.twiny()
        ax3 = ax.twinx()
        ax2.set_xlim(xlim)
        ax3.set_ylim(ylim)
        # HACK
        if not "Av_1" and "beta_1" in names:
            yname = '$%s$' % names[0]
            xname = '$%s$' % names[1]
        else:     # HACK changing Av_1 to \Delta A_V
            yname = '$\Delta A_V$'
            xname = '$\Delta \\beta$'
        # END HACK
        ax3.set_ylabel(yname)
        ax2.set_xlabel(xname)

    # ### HACK ### Av needs to be inverted
    if invert_yticks:
        yticks = ax3.get_yticks()
        ax3.set_yticklabels(yticks * -1)
    # ### END HACK ###
    if not plot_delta_values:  # get rid of the tickmarks
        ax.set_xticks([])
        ax.set_yticks([])

    path = storepath + 'marginalization.png'
    fig.savefig(path)
    path = storepath + 'marginalization.pdf'
    fig.savefig(path)

    return (dx, dy, delta_chi_sq, levels)


def test_fit():
    '''
    Test the leastsq algorithm on a toy model

    It seems that the standard leastsq fit is more sensitive to this choice
    of initial parameters (for mu in particular) than the fmin fit below.
    When I had mu = 20, it sometimes converged to noise with a terrible chisq,
    '''
    import matplotlib.pyplot as plt

    # Make simulated data:
    xvals = np.arange(100)
    zeros = np.zeros(100)
    zipxvals = zip(xvals, zeros)

    gaussian = lambda x: 3*np.exp(-(30-x)**2/20.)
    # true values: mu = 30, height = 3, sigma = sqrt(20) = 4.472
    ydata = gaussian(xvals)
    ydata = scipy.randn(100)*.05+ydata  # adding noise
    yerr = np.zeros(100)+.05  # array of uncertainties

    # Give initial paramaters:

    mu = Param(34, name='mu')
    sigma = Param(4, name='sigma')
    height = Param(5, name='height')

    # Define your function:
    def f(x):
        ''' here we are using proof of concept of a multi-input function
        e.g. y = f(x, t)
        'x' is a list of tuples, zipped with the zip function (see below)
        we unzip them and then evaluate.
        x needs to be a single parameter because of the way the errfunc
        in fit() is defined
        e.g., x = [(0, 0.0), (1, 0.0), (2, 0.0),
                   (3, 0.0), (4, 0.0), (5, 0.0)]
        '''
        xval, zero = zip(*x)  # unzip the x feature vector
        return height() * np.exp(-((xval-mu())/sigma())**2) + zero

    # Fit the function (provided 'data' is an array with the data to fit):
    retdict = fit(f, [mu, sigma, height], ydata, yerr, zipxvals)

    # Plot the fitted model over the data if desired
    simxvals = np.arange(10000)/100.  # 10000 points from 0-100
    simzeros = np.zeros(len(simxvals))
    zipsimxvals = zip(simxvals, simzeros)

    fig2 = plt.figure()
    ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(simxvals, f(zipsimxvals))

    ax.scatter(xvals, ydata)
    fig2.show()
    return retdict


def test_linear(fix_intercept=False):
    '''
    Test the leastsq algorithm on a toy linear model
    '''
    import matplotlib.pyplot as plt

    model1color = 'green'
    model2color = 'red'
    truthcolor = 'black'

    # Make simulated data:
    xvals = np.arange(20)*4
    # gaussian = lambda x: 3*np.exp(-(30-x)**2/20.)
    # true values: mu = 30, height = 3, sigma = sqrt(20) = 4.472
    truedict = {'slope': 0.50, 'intercept': 0.0}

    linear = lambda x: truedict['slope'] * x + truedict['intercept']
    ydata_true = linear(xvals)
    est_y_err = 5000*(xvals+10.0)**-1.8  # making some heteroskedastic noise
    ydata = scipy.randn(20)*est_y_err+ydata_true  # adding noise to the data
    yerr = np.zeros(20)+est_y_err  # array of uncertainties

    # set up your dictonary of values to fit
    fitdict = {}
    names = ['slope', 'intercept']
    for name in names:
        fitdict.update({name: {'init': 1, 'fixed': False}})

    if fix_intercept:
        # fit the intercept at 0
        fitdict['intercept']['init'] = 0
        fitdict['intercept']['fixed'] = True

    # BUILD UP PARAMETERS
    fullparamlist = []
    fitparamlist = []
    fixparamlist = []

    # set parameters
    for key, val in fitdict.iteritems():
        param = Param(val['init'], name=key)
        fullparamlist.append(param)
        if val['fixed'] is False:
            fitparamlist.append(param)
        elif val['fixed'] is True:
            fixparamlist.append(param)
        else:
            raise ValueError('Invalid value for Fixed. Must be True/False.')

    # creating parameters to fit
    # beta1 = Param(3, name = 'slope')
    # beta2 = Param(10, name = 'intercept')
    # paramlist = [beta1, beta2]
    #
    # Define your function:
    def myfunc(x, paramlist):
        for param in paramlist:
            if param.name == 'slope':
                slope = param.value
            elif param.name == 'intercept':
                intercept = param.value
            else:
                raise Exception('Invalid parameter')
        return slope*x + intercept

    def f(x):
        return myfunc(x, fullparamlist)

    # Plot the fitted model over the data if desired
    simxvals = np.array([-5, 50, 85])  # 10000 points from 0-100
    simzeros = np.zeros(len(simxvals))

    fig2 = plt.figure()
    ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

    # plot the underlying truth
    ax.plot(simxvals, linear(simxvals), lw=7, alpha=0.4, color=truthcolor)

    # Fit the function (provided 'data' is an array with the data to fit):
    retdict = fit(f, fitparamlist, ydata, yerr, xvals)
    ax.plot(simxvals, f(simxvals), lw=2, color=model1color)

    # Now assuming tiny/no error
    # Fit the function (provided 'data' is an array with the data to fit):
    retdictnoerror = fit(f, fitparamlist, ydata, np.zeros(20)+.001, xvals)
    ax.plot(simxvals, f(simxvals), lw=2, color=model2color)

    # Plot the data with error bars
    ax.errorbar(xvals, ydata, yerr=yerr, fmt='.', color=truthcolor)

    ax.set_xlim((-5, 85))
    ax.set_ylim((-120, 120))

    # bottom left annotations
    xloc = 0.2
    textoffset = 0.16
    textincrement = 0.04
    textsize = 12
    color = model1color
    modeldict = retdict
    title = 'Including measurement error'
    # model1 annotations
    string = 'chi2 / dof = %.2f / %i' % (modeldict['chi2'], modeldict['dof'])
    fig2.text(xloc, textoffset, string)
    textoffset += textincrement
    for string in modeldict['strings']:
        fig2.text(xloc, textoffset, string, color=color)
        textoffset += textincrement
    fig2.text(xloc, textoffset, title, size=textsize)

    # bottom right annotations
    xloc = 0.55
    textoffset = 0.16
    textincrement = 0.04
    textsize = 12
    color = model2color
    modeldict = retdictnoerror
    title = 'Assuming no measurement error :'
    # model2 annotations
    string = 'chi2 / dof = %.2f / %i' % (modeldict['chi2'], modeldict['dof'])
    fig2.text(xloc, textoffset, string)
    textoffset += textincrement
    for string in modeldict['strings']:
        fig2.text(xloc, textoffset, string, color=color)
        textoffset += textincrement
    fig2.text(xloc, textoffset, title, size=textsize)

    # True annotations
    textoffset = 0.7
    for key, val in truedict.iteritems():
        string = key + ': ' + str(val)
        fig2.text(0.2, textoffset, string, color=truthcolor,
                  alpha=0.4, size=14)
        textoffset += 0.042
    fig2.text(0.2, textoffset, 'True Inputs:', size=14)

    string = 'Fixed params'
    if fixparamlist:
        fig2.text(0.7, 0.8, string)
        textoffset = 0.76
    for myparam in fixparamlist:
        string = myparam.name + ': ' + str(myparam.value)
        fig2.text(0.7, textoffset, string)
        textoffset -= 0.04

    fig2.show()
    return retdict


def sample_from_multivariate_normal_test(retdict, plot_every_model=True):
    '''Proof of concept of sampling from a multivariate gaussian to generate
    a distribution of model fits.
    docs.scipy.org/doc/numpy/reference/generated/numpy.random.multivariate_normal.html
    The multivariate normal, multinormal or Gaussian distribution is a
    generalization of the one-dimensional normal distribution to higher
    dimensions. Such a distribution is specified by its mean and covariance
    matrix. These parameters are analogous to the mean (average or "center")
    and variance (standard deviation, or "width," squared) of the
    one-dimensional normal distribution.
    '''

    # definition of the underlying functional model
    def f(x, mu, sigma, height):
        return height * np.exp(-((x-mu)/sigma)**2)

    # extracting the mean values and covariance matrix from the fit retdict
    mean = [param.value for param in retdict['parameters']]
    mean = tuple(mean)
    cov = retdict['covarmatrix']

    nxvals = 1000
    nsimulations = 1500

    # nxvals samples from normal distribution
    samples = np.random.multivariate_normal(mean, cov, (nsimulations))

    simxvals = np.linspace(0, 100, nxvals)  # 1000 points from 0-100

    fig2 = plt.figure()
    ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

    # Array of simulated y values
    simymat = np.zeros((nsimulations, nxvals))

    count = 0
    for sample in samples:
        simyvals = f(simxvals, sample[0], sample[1], sample[2])
        simymat[count] = simyvals
        # plot every simulated model with low opacity. This takes a long time
        # but gives a nice illustration of the underlying possible models
        if plot_every_model:
            ax.plot(simxvals, simyvals, lw=2, color='black', alpha=0.01)
        count += 1

    # Two sigma is 95.4%
    simylower02point3 = stats.scoreatpercentile(simymat, 2.3)
    simyupper97point7 = stats.scoreatpercentile(simymat, 97.7)

    if not plot_every_model:
        ax.fill_between(simxvals, simylower02point3,
                        simyupper97point7, color='#CCCCCC')
        simyupper84 = stats.scoreatpercentile(simymat, 84.1)
        simylower16 = stats.scoreatpercentile(simymat, 15.9)
        ax.fill_between(simxvals, simylower16, simyupper84, color='#888888')
        ax.plot(simxvals, f(simxvals, mean[0], mean[1], mean[2]),
                color='black')

    else:
        ax.plot(simxvals, simyupper97point7, color='red', lw=2)
        ax.plot(simxvals, simylower02point3, color='blue', lw=2)
    return simymat


def test_fit_fmin():
    import matplotlib.pyplot as plt

    # Make simulated data:
    xvals = np.arange(100)
    zeros = np.zeros(100)
    zipxvals = zip(xvals, zeros)

    gaussian = lambda x: 3*np.exp(-(30-x)**2/20.)
    # True values: mu = 30, height = 3, sigma = sqrt(20) = 4.472
    ydata = gaussian(xvals)
    ydata = scipy.randn(100)*.05+ydata  # adding noise
    yerr = np.zeros(100)+.05  # array of uncertainties

    # Give initial paramaters:
    mu = Param(34, name='mu')
    sigma = Param(4, name='sigma')
    height = Param(5, name='height')

    # Define your function:
    def f(x):
        ''' here we are using proof of concept of a multi-input function
        e.g. y = f(x, t)
        'x' is a list of tuples, zipped with the zip function (see below)
        we unzip them and then evaluate.
        x needs to be a single parameter because of the way the errfunc
        in fit() is defined
        e.g., x = [(0, 0.0), (1, 0.0), (2, 0.0),
                   (3, 0.0), (4, 0.0), (5, 0.0)]
        '''
        xval, zero = zip(*x)  # unzip the x feature vector
        return height() * np.exp(-((xval-mu())/sigma())**2) + zero

    # Fit the function (provided 'data' is an array with the data to fit):
    fminfit(f, [mu, sigma, height], ydata, yerr, zipxvals)

    # Plot the fitted model over the data if desired
    simxvals = np.arange(10000)/100.  # 10000 points from 0-100
    simzeros = np.zeros(len(simxvals))
    zipsimxvals = zip(simxvals, simzeros)

    fig2 = plt.figure()
    ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(simxvals, f(zipsimxvals))

    ax.scatter(xvals, ydata)
    fig2.show()


def test_multi():
    # Borrowed from SEDfit code
    y0 = np.array([12, 13, 14, 15, 16, 17]*10)
    x0 = np.arange(300).reshape(60, 5)
    names = ['b1', 'b2', 'b3', 'b4', 'b5']
    fitdict = {}
    for name in names:
        fitdict.update({name: {'init': 1, 'fixed': False}})
    # BUILD UP PARAMETERS
    fullparamlist = []
    fitparamlist = []
    fixparamlist = []

    # Set parameters
    for key, val in fitdict.iteritems():
        param = Param(val['init'], name=key)
        fullparamlist.append(param)
        if val['fixed'] is False:
            fitparamlist.append(param)
        elif val['fixed'] is True:
            fixparamlist.append(param)
        else:
            raise ValueError('Invalid value for Fixed. Must be True/False.')

    def myf(x):
        return np.dot(x, np.array([beta.value for beta in fitparamlist]).T)

    retdict = fit(myf, fitparamlist, y0, np.array([0.1]*60), x0)
    return retdict
