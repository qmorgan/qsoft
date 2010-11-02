'''
Extracts 

Requirements:

urllib2
BeautifulSoup
'''
# Perform initial imports
import urllib2
import BeautifulSoup
import numpy

def obtainlinks(url='http://tanga.com/bottom_feeder/archive?page=1',search_string='/bottom_feeder/'):
    # First, load up the webpage
    myurl = urllib2.urlopen(url)
    soup = BeautifulSoup.BeautifulSoup(myurl)

    # Extract all the links from the webpage
    anchors = soup.findAll('a')
    links = []
    for a in anchors:
        links.append(a['href'])

    # Parse the links to find the desired ones
    parsedlinks = []
    for link in links:
        # check if the search string is in the link before continuing
        if search_string in link:
            
            if len(link.lstrip(search_string)) <= 4:
                parsedlinks.append(link)
                
    return parsedlinks
    
def findTangaInfo(url='http://tanga.com/bottom_feeder/7360'):
    bf1 = urllib2.urlopen(url)
    bf1read = bf1.readlines()
    bf1 = urllib2.urlopen(url)
    soupbf1 = BeautifulSoup.BeautifulSoup(bf1)
    
    numlist = []
    namelist = []
    # obtain the columns
    anchors = soupbf1.findAll('td')
    for anchor in anchors:
        stripanch = str(anchor).lstrip('<td>').rstrip('</td>')
        # if small, likely to be a number
        try:
            numlist.append(int(stripanch))
        except:
            namelist.append(stripanch)
    
    bidline = ''
    # Find the number of bids by straight up parsing
    for line in bf1read:
        if 'Number of bids' in line:
            bidline = line
    
    numberofbids = bidline.strip().split()[-1].rstrip('.')
    
    numlistarr = numpy.array(numlist)
    minval = numlistarr.min()
            
    retdict = {'min':minval,'numofbids':numberofbids,'numlist':numlist,'namelist':namelist}
    return retdict

def loopThruTanga():
    id_list = numpy.arange(25)+1  # list of [1,2,...24,25]
    full_linklist=[]
    
    for id in id_list:
        pageurl = 'http://tanga.com/bottom_feeder/archive?page=%i' % (id)
        linklist = obtainlinks(url=pageurl,search_string='/bottom_feeder/')
        for link in linklist:
            newlink = 'http://tanga.com' + link
            full_linklist.append(newlink)
    print len(full_linklist)
    # Remove duplicates
    full_linklist = list(set(full_linklist))
    print 'looping through ' + str(len(full_linklist)) + ' links'
    dictlist = []
    for link in full_linklist:
        try:
            newdict = findTangaInfo(link)
            print 'For %s, #of bids=%s, lowest bid=%s' % (link, str(newdict['numofbids']), str(newdict['min']))    
            dictlist.append(newdict)
        except:
            print 'Cannot extract information from %s' % (link)
    return dictlist