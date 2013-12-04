from pylab import *
from scipy import *
from scipy import optimize
import numpy

# I think something is wrong with how this calculates uncertainties in the fit

powerlaw = lambda x, amp, index: amp * (x**index)


xdata=numpy.array([5805.5770000000002, 9605.8690000000006, 11560.356, 15643.641])
ydata=numpy.array([6.802,5.108,4.287,3.908])
yerr=numpy.array([0.833,0.103,0.433,0.075])


logx = log10(xdata)
logy = log10(ydata)
logyerr = yerr / ydata

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit,
                       args=(logx, logy, logyerr), full_output=1)

pfinal = out[0]
covar = out[1]
print pfinal
print covar

index = pfinal[1]
amp = 10.0**pfinal[0]

indexErr = sqrt( covar[0][0] )
ampErr = sqrt( covar[1][1] ) * amp

print indexErr

fake=linspace(5000,20000,100)

clf()
subplot(2, 1, 1)
plot(fake, powerlaw(fake, amp, index))     # Fit
errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
text(10000, 6.5, 'Ampli = %5.3f +/- %5.3f' % (amp, ampErr))
text(10000, 5.5, 'Index = %5.3f +/- %5.3f' % (index, indexErr))
title('Best Fit Power Law')
xlabel('Time Since Burst (s)')
ylabel('UVOT Count rate (cts/s)')
xlim(5000, 16000)


subplot(2, 1, 2)
loglog(fake, powerlaw(fake, amp, index))
errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
xlabel('Time Since Burst (s)')
ylabel('UVOT Count rate (cts/s)')
xlim(5000.0, 16000)