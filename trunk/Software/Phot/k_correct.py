from scipy import integrate
def wavelength_to_freq(wavelength):
    freq = 3*10**8/(wavelength) # in SI
    return freq

def k_correct_integral(z,a, filter):
    #Define spectrum: power law with index a
    def spect(nu, z, a):
        return nu*(1+z)**(1-a)
    #Bandpass
    if filter == 'j':
        r1 = wavelength_to_freq(11500*10^-9) # wavelength in angstrom: 11500
        r2 = wavelength_to_freq(13300*10^-9) # wavelength in angstrom: 13300
    elif filter == 'h':
        r1 = wavelength_to_freq(14900*10^-9) # wavelength in angstrom: 14900
        r2 = wavelength_to_freq(17800*10^-9) # wavelength in angstrom: 17800
    elif filter == 'k':
        r1 = wavelength_to_freq(20300*10^-9) # wavelength in angstrom: 20300
        r2 = wavelength_to_freq(23600*10^-9) # wavelength in angstrom: 23600
        
    #Let's integrate
    args = (z,a)
    results = integrate.quad(spect,r1,r2, args)
    print 'Integral = ', results[0], ' with error = ', results[1]
    return results

def correct(z, a):
    K = 2.5*(a-1)*log10(1+z)
    return K

def time_correct(z):
    T = 1+z
    return T
