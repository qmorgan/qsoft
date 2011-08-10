import urllib2
from BeautifulSoup import BeautifulSoup


def parseNatWebTable(soup, expected_cols=20,errtype='ul'):
    '''errtype is either:
        ul: '1.8e-06 (1.6e-06,1.9e-06)'
        pm: 41.920 +/-2.017
    '''
    errtypelist = ['ul','pm']
    if errtype not in errtypelist:
        raise ValueError('Error type not acceptable.')
    table = soup.find('table')
    rows = table.findAll('tr')
    rowcount = 0
    grbname = ''
    grbdict = {}
    for tr in rows:
        subdict={}
        if rowcount == 0:
            headerlist = tr.findAll(text=True)
            headerlist = headerlist[0:expected_cols]
            
            # rename any duplicates in the header list by appending a number
            duplicates = [x for x in headerlist if headerlist.count(x)>=2]

            for duplicate in duplicates:
                addd = headerlist.count(duplicate)
                headerlist[headerlist.index(duplicate)] = str(addd) + duplicate
                            
        colcount = 0
        trstrip = [a for a in tr.findAll(text=True) if a != u' ']
        
        # Make sure there are the right number of rows in the table
        # But don't skip the first row no matter what because there is a malformed
        # table where the first row has more 
        if len(trstrip) == expected_cols + 1 and tr.find('a') != None:
            print trstrip.pop(1)
        if not len(trstrip) == expected_cols and rowcount != 0:
            print 'Skipping this burst'
            continue
        cols = tr.findAll('td')
        # loop through every column
        for td in cols:
            currentval = td.findAll(text=True)
            if currentval and colcount == 0:
                if td.find('a') and currentval[1].strip():
                    grbname = currentval[1].strip()
                    colcount += 1
                    continue
                elif currentval[0].split('(')[0].strip():
                    grbname = currentval[0].split('(')[0].strip()
                    colcount += 1
                    continue
                else: 
                    continue    
            elif currentval and colcount != 0:
                val = currentval[0].strip()
                valname = headerlist[colcount].strip()
                # If there are parentheticals, these are the upper and lower bounds.
                # parse them accordingly: '1.8e-06 (1.6e-06,1.9e-06)'
                if errtype == 'ul' and val.find('(') != -1 and val.find(',') != -1:
                    value = val.split('(')[0].strip()
                    print value
                    remainder = val.split('(')[1].strip().split(',')
                    lowlimit = remainder[0].strip()
                    uplimit = remainder[1].strip().strip(')')
                    upname = valname + '_upper'
                    lowname = valname + '_lower'
                    subdict.update({valname:value,upname:uplimit,lowname:lowlimit})
                elif errtype == 'pm' and val.find('+/-') != -1:
                    value = val.split('+/-')[0].strip()
                    error = val.split('+/-')[1].strip()
                    errname = 'D'+valname
                    subdict.update({valname:value,errname:error})
                else:
                    subdict.update({valname:val})
        
                if grbname:
                    grbdict.update({grbname:subdict})
                    
                    
            
                #print subdict
            else:
                pass
            
            colcount += 1
            # Don't try to parse blank columns at the end
            if colcount >= expected_cols:
                break
        rowcount += 1
            
    print 'Done'
    print headerlist
    return grbdict

def ParseBATSpec():
    sock = urllib2.urlopen('http://astro.berkeley.edu/~nat/swift/bat_spec_table.html')
    html = sock.read()
    sock.close()
    soup=BeautifulSoup(html)
    batspec = parseNatWebTable(soup,expected_cols=20,errtype='ul')
    return batspec 
    
def ParseBATTiming():
    sock = urllib2.urlopen('http://astro.berkeley.edu/~nat/swift/bat_time_table.html')
    html = sock.read()
    sock.close()
    soup=BeautifulSoup(html)
    battiming = parseNatWebTable(soup,expected_cols=15,errtype='pm')
    return battiming
    
def ParseXRTSpec():
    sock = urllib2.urlopen('http://astro.berkeley.edu/~nat/swift/xrt_spec_table.html')
    html = sock.read()
    sock.close()
    soup=BeautifulSoup(html)
    xrtspec = parseNatWebTable(soup,expected_cols=17,errtype='ul')
    return xrtspec
