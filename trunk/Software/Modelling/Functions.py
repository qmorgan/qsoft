import numpy as np

def DecayingExponential(timearr,Av_0,Av_1,tau):
    Av = Av_0 + Av_1*np.exp(-1*timearr/tau)
    return Av

def BrokenPowerLaw(timearr,Av_0,Av_1,Av_2):
    # Av_2 = tbreak
    Avlist=[]
    for time in timearr:
        if time <= Av_2:
            Av = Av_0/(Av_2**Av_1)*time**Av_1
        if time > Av_2:
            Av = Av_0 
        Avlist.append(Av)
    Avarr = np.array(Avlist)    
    return Avlist
    ###
