'''
Given a set of odds, figure out how to place bets to maximize the chance of
return.

'''
import numpy

def EstReturn(odds,stake):
    if type(odds).__name__ != 'tuple':
        print 'Was Epecting Tuple odds'
        return
    elif len(odds) != 2:
        print 'Was expecting tuple of length two for odds'
        return
    return float(odds[0])/float(odds[1]) * stake
    
def test2(odds1,odds2,outstr1="Outcome 1",outstr2="Outcome 2",startbet=10.):
    total = float(startbet)    
    print "Starting with $%.2f" % total
    print "Odds are %i:%i on %s, %i:%i on %s" % (odds1[0],odds1[1],outstr1,
        odds2[0],odds2[1],outstr2)
    for i in numpy.linspace(0,total,101):
        bet1 = i
        bet2 = total - i
        print 
        print "Putting %.2f on %s and %.2f on %s: " % (bet1, outstr1, bet2, outstr2)
        winnings_1 = EstReturn(odds1,bet1) - bet2
        winnings_2 = EstReturn(odds2,bet2) - bet1
        if winnings_1 >= 0: 
            team1_out = 'earn'
        else:
            team1_out = '*lose*'
        if winnings_2 >= 0: 
            team2_out = 'earn'
        else:
            team2_out = '*lose*'
            
        print "If %s, you will %s %.2f" % (outstr1, team1_out, abs(winnings_1))
        print "If %s, you will %s %.2f" % (outstr2, team2_out, abs(winnings_2))
        