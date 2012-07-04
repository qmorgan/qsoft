import numpy as np

def DecayingExponentialAv(timearr,Av_0,Av_1,tau):
    # constraining Av_0 to be negative
    # constraining Av_1 to be negative - cant have dust creation
    Av = -1*abs(Av_0) + -1*abs(Av_1)*np.exp(-1*timearr/tau)
    return Av

def DecayingExponentialbeta(timearr,beta_0,beta_1,tau):
    beta = beta_0 + beta_1*np.exp(-1*timearr/tau)
    return beta

def BrokenPowerLaw(time,Av_0,Av_1,Av_2):
    # Av_2 = tbreak
    try: #try if its an array
        Avlist=[] 
        for tim in time:
            if tim <= Av_2:
                Av = Av_0/(Av_2**Av_1)*tim**Av_1
            if tim > Av_2:
                Av = Av_0 
            Avlist.append(Av)
        Av = Avlist
    except:
        if time <= Av_2:
            Av = Av_0/(Av_2**Av_1)*time**Av_1
        if time > Av_2:
            Av = Av_0 
    return Av
    ###
