import numpy
import time
    
def add(num=2,max_num=1000,do='rand'):
    dorand = False 
    randsign = 1
    if do.upper() == 'RAND': dorand = True
    elif do.upper() == 'ADD': randsign = 1
    elif do.upper() == 'SUBTRACT': randsign = 0
    else: raise ValueError 
    numarr = numpy.random.random(num)*max_num
    numarr = numarr.round()
    numarr=list(numarr)  # for some reason, index(item) fails for numpy arrays
    if randsign == 0 and num == 2: numarr.sort(reverse=True) # big number first for 2 item subtraction
    pstr ='\n %4d' % (numarr[0])
    for item in numarr[1:]:
        if dorand:
            randsign = numpy.random.randint(2)
        if randsign == 0: 
            sign = '-'
            numarr[numarr.index(item)]*=-1
        else: 
            sign = '+'
        pstr += '\n%s%4d' % (sign,item)
    pstr += '\n-----'
    print pstr
    ans = sum(numarr)
    myans = raw_input(' ')
    if myans.upper().strip(' ') == 'EXIT': return -1 
    try: myans = int(myans)
    except: myans = 0
    if myans == int(ans): 
        print 'Correct!'
        return 1
    else: 
        print 'Wrong.  Answer is %i' % int(ans)
        return 0
    
def add_loop(do='add'):
    print 'Type the answer or EXIT to quit'
    totnumber=-1
    correct=0
    thereturn=0
    t1 = time.time()
    while thereturn != -1:
        correct+=thereturn
        thereturn = add(do=do)
        totnumber+=1 
    percent_correct = float(correct)/float(totnumber)
    tot_time = time.time() - t1
    avg_time_per_correct = tot_time/correct
    print 'You got %i out of %i correct. (%f%%)' % (correct,totnumber,percent_correct) 
    print 'Total time is %i seconds, which is on average %i s per correct answer.' %\
        (int(tot_time), int(avg_time_per_correct))