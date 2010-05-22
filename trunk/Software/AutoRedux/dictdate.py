#!/usr/bin/python
# Filename : dictdate.py

def formattime(mydict):
    for i in mydict:
        time1 = mydict[i]['ptel_time']
        time2 = time1[0:4]+time1[5:7]+time1[8:10]
        formated = 'http://skycam.mmto.arizona.edu/skycam/'+time2+'/night_movie.avi'
        mydict[i]['formated_time']=formated





