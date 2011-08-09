import urllib2
from BeautifulSoup import BeautifulSoup


def parseNatWebTable(html, expected_cols=20):
    table = soup.find('table')
    rows = table.findAll('tr')
    rowcount = 0
    grbname = ''
    grbdict = {}
    for tr in rows:
        subdict={}
        if rowcount == 0:
            headerlist = tr.findAll(text=True)                
        colcount = 0
        # Make sure there are the right number of rows in the table
        if not len(tr) == expected_cols:
            print 'Skipping this burst'
            continue
        cols = tr.findAll('td')
        # loop through every column
        for td in cols:
            currentval = td.findAll(text=True)
            if td.find('a'): # if there is a link in the val
                if currentval[1].strip() and colcount == 0:
                    grbname = currentval[1].strip()
                    colcount += 1
                    continue
                else: 
                    continue    
            if currentval: 
                val = currentval[0].strip()
                # If there are parentheticals, these are the upper and lower bounds.
                # parse them accordingly: '1.8e-06 (1.6e-06,1.9e-06)'
                print val
                if val.find('(') != -1 and val.find(',') != -1:
                    value = val.split('(')[0].strip()
                    print value
                    remainder = val.split('(')[1].strip().split(',')
                    lowlimit = remainder[0].strip()
                    uplimit = remainder[1].strip().strip(')')
                    valname = headerlist[colcount].strip()
                    upname = valname + '_upper'
                    lowname = valname + '_lower'
                    subdict.update({valname:value,upname:uplimit,lowname:lowlimit})
                else:
                    subdict.update({headerlist[colcount].strip():val})
                
                if grbname:
                    grbdict.update({grbname:subdict})
                    
                    
            
                #print subdict
            else:
                pass
            
            colcount += 1
        rowcount += 1
            
    print 'Done'
    print headerlist
    return grbdict

def TestParse():
    sock = urllib2.urlopen('http://astro.berkeley.edu/~nat/swift/bat_spec_table.html')
    html = sock.read()
    sock.close()
    soup=BeautifulSoup(html)
    a = parseNatWebTable(html)
    return a 
