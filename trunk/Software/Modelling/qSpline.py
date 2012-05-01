def qSpline(xvals,yvals,yerrvals,xgrid):
    
    fitpath = storepath + 'regrSpline_fit.txt'
    fiterrpath = storepath + 'regrSpline_fiterr.txt'
    write_ypath = storepath + 'regrSpline_y_in.txt'
    write_xpath = storepath + 'regrSpline_x_in.txt'
    
    ### WRITE OUT THE DATA TO FILE FOR R TO READ IN
    xout=file(write_xpath,'w')
    for val in xgrid:
        outstr = str(val) + '\n'
        xout.write(outstr)
    xout.close()
    
    yout=file(write_ypath,'w')
    ind = 0
    for xval in xvals:
        yval = yvals[ind]
        yeval = yerrvals[ind]
        outstr = str(xval) + ' ' + str(yval) + ' ' + str(yeval) + '\n'
        yout.write(outstr)
        ind += 1
    yout.close()
    
    ### CALL THE R FUNCTION TO CREATE THE SPLINE OUTPUT
    functionstr='R < %scall_regrSpline.R --no-save' % splinedir
    os.system(functionstr)
    
    
    # read in the new lines 
    f=file(fitpath,'r')
    g=file(fiterrpath,'r')
    
    # assuming line format of
    # '"10" 17.719301369 \n'
    newylist=[]
    for line in f.readlines():
        linesplit = line.split()
        try:
            newylist.append(float(linesplit[1]))
        except:
            print "skipping line %s" % line
    f.close()
    
    
    spline_model_errlist=[]
    for line in g.readlines():
        linesplit = line.split()
        try:
            spline_model_errlist.append(float(linesplit[1]))
        except:
            print "skipping line %s" % line
    g.close()
    
    newyarr = np.array(newylist)
    spline_model_errarr = np.array(spline_model_errlist)
    
    return newyarr, spline_model_errarr