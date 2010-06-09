import datetime
import pyfits
from RedshiftMachine import LoadGCN

def t_mid(filepath, GRBid=None):
    '''Given a fitts file and the GRB time, this program determines the t-mid for PAIRITEL'''
    header = pyfits.open(filepath)
    starttime = header[0].header['STRT_CPU']
    stoptime = header[0].header['STOP_CPU']
    
    start = datetime.datetime.strptime(starttime.split('.')[0], "%Y-%m-%d %H:%M:%S")
    stop = datetime.datetime.strptime(stoptime.split('.')[0], "%Y-%m-%d %H:%M:%S")

    
    durdiv2 = (start - stop)/2 + start
    print "durdiv2 is " +  str(durdiv2)

    if not GRBid:
        trg = int(header[0].header['TRGTNAME'][6:])
        dict = LoadGCN.LoadGCN(trg)
    else:
        dict = LoadGCN.LoadGCN(GRBid)
    
    GRBtime = dict.pdict['grb_time_str']
    GRBdate = dict.pdict['grb_date_str']
    
    GRBcomb = GRBdate + ' ' + GRBtime
    
    GRB = datetime.datetime.strptime(GRBcomb.split('.')[0], "%y/%m/%d %H:%M:%S")
    
    t_mid = durdiv2 - GRB  
    t_mid_str = str(t_mid)

    print t_mid_str
    
    if 'days' in t_mid_str:
        t_mid_days = float(t_mid_str.split(' days')[0]) * (24.)
    else:
        t_mid_days = 0.


#edited for photloop-------

    t_mid_time_list = t_mid_str.split(':')

    print t_mid_time_list
    t_mid_time_hour = float(t_mid_time_list[0]) + float(t_mid_time_list[1])*(1/60.) + float(t_mid_time_list[2])*(1/3600.)
    t_mid_hour = t_mid_days + t_mid_time_hour
    

    
#original-------

    #t_mid_time_list = t_mid_str[::-1][0:8][::-1].split(':')

    #print t_mid_time_list
    #t_mid_time_hour = float(t_mid_time_list[0]) + float(t_mid_time_list[1])*(1/60.) + float(t_mid_time_list[2])*(1/3600.)
    #t_mid_hour = t_mid_days + t_mid_time_hour
    


    print t_mid_hour
    return t_mid_hour



