#! /usr/bin/python

## arXiv script version 0.2

## Copyright 2008 Tom Brown

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## See http://www.stringwiki.org/wiki/ArXiv_script for more usage
## instructions

'''arXiv script
Usage:
python arxiv.py reference [ -htabcjdps ] [ --help ]
"reference" must be a standard arXiv reference, e.g. hep-th/9711200, 0705.0303.
Options:
-h, --help
displays this help message
-t
displays the title
-a
displays the author(s)
-b
displays the aBstract
-c
displays the comments
-j
displays the journal reference
-d
downloads the PDF
-p
downloads the PS
-s
downloads the source file
'''

__version__ = "0.2"
__author__ = "Tom Brown"
__copyright__ = "Copyright 2008 Tom Brown, GNU GPL 3"


import sys, os, getopt, re, urllib,gzip


def findRefType(ref):
    ref = ref.replace('arxiv:','')
    if re.search(r'^[a-zA-Z\-]+/\d{7}$',ref):
        type = 'old-style eprint'
    elif re.search(r'^\d{7}$',ref):
        type = 'old-style eprint'
        ref = 'hep-th/' + ref
    elif re.search('^\d{4}\.\d{4}$',ref):
        type = 'new-style eprint'
    else:
        type = 'not arXiv'

    return type, ref




def downloadPDF(ref,type,downloadPath):
    downloadPath = os.path.expanduser(downloadPath)
    if type == 'old-style eprint':
        urllib.urlretrieve('http://arxiv.org/pdf/' + ref, downloadPath + ref.replace('/','-') + '.pdf')
    elif type == 'new-style eprint':
        urllib.urlretrieve('http://arxiv.org/pdf/' + ref, downloadPath + ref + '.pdf')


def downloadPS(ref,type,downloadPath):
    downloadPath = os.path.expanduser(downloadPath)
    filename = downloadPath + ref.replace('/','-')
    urllib.urlretrieve('http://arxiv.org/ps/' + ref, filename)
    gzipFile =  gzip.GzipFile(filename)
    psFile = open(filename + ".ps","w")
    psFile.write(gzipFile.read())
    psFile.close()
    gzipFile.close()
    os.remove(filename)

def downloadSource(ref,type,downloadPath):
    downloadPath = os.path.expanduser(downloadPath)
    filename = downloadPath + ref.replace('/','-')
    urllib.urlretrieve('http://arxiv.org/e-print/' + ref, filename + ".dum")
    gzipFile =  gzip.GzipFile(filename + ".dum")
    sourceFile = open(filename,"w")
    sourceFile.write(gzipFile.read())
    sourceFile.close()
    gzipFile.close()
    os.remove(filename + ".dum")


def getTitle(html):
    title = html[html.find(">Title:</span>")+15:]
    title = title[:title.find("</h1>")]
    return title


def getAuthors(html):
    authors = html[html.find(">Authors:</span>"):]
    authors = authors[authors.find("\">")+2:]
    authors = authors[:authors.find("</div>")]
    authors = re.sub('<[^>]*>','',authors)
    authors = authors.replace("\n","")
    return authors


def getAbstract(html):
    abstract = html[html.find("Abstract:</span>")+17:]
    abstract = abstract[:abstract.find("</blockquote>")-1]
    return abstract

def getComments(html):
    if html.count("comments") == 0:
        return "no comments"
    else:
        comments = html[html.find("comments\">")+10:]
        comments = comments[:comments.find("</td>")]
        return comments


def getJref(html):
    if html.count("jref") == 0:
        return "no journal reference"
    else:
        jref = html[html.find("jref\">")+6:]
        jref = jref[:jref.find("</td>")]
        return jref




if __name__ == "__main__":

    authorOpt = 0
    titleOpt = 0
    abstractOpt = 0
    commentsOpt = 0
    jrefOpt = 0
    pdfOpt = 0
    psOpt = 0
    sourceOpt = 0

    try:
        options, arguments = getopt.gnu_getopt(sys.argv[1:], 
        'hatbcjdpsv', ['help'])
    except getopt.error:
        print 'error: you tried to use an unknown option or the argument for an option that requires it was missing; try \'arxiv.py -h\' for more information'
        sys.exit(0)

    for o,a in options:
        if o in  ('-h','--help'):
            print __doc__
            sys.exit(0)

        elif o == '-a':
            authorOpt = 1

        elif o == '-t':
            titleOpt = 1

        elif o == '-b':
            abstractOpt = 1

        elif o == '-c':
            commentsOpt = 1

        elif o == '-j':
            jrefOpt = 1

        elif o == '-d':
            pdfOpt = 1

        elif o == '-p':
            psOpt = 1

        elif o == '-s':
            sourceOpt = 1


    if len(options) == 0:
        authorOpt = 1
        titleOpt = 1
        abstractOpt = 1
        commentsOpt = 1
        jrefOpt = 1



    if len(arguments) != 1:
        print 'you didn\'t specify an arXiv reference; try \'arxiv.py -h\' for more information'
        sys.exit(0)
    else:
        ref=arguments[0]





    type, ref = findRefType(ref)

    if type=="not arXiv":
        print "type not of arXiv form"
        sys.exit(0)

    if (authorOpt+titleOpt+abstractOpt+commentsOpt+jrefOpt > 0):
        htmlObject = urllib.urlopen('http://arxiv.org/abs/' + ref)
        html = htmlObject.read()

    if titleOpt:
        title = getTitle(html)
        print title

    if authorOpt:
        authors = getAuthors(html)
        print authors


    if abstractOpt:
        abstract = getAbstract(html)
        print abstract


    if commentsOpt:
        comments = getComments(html)
        print comments


    if jrefOpt:
        jref = getJref(html)
        print jref

    if pdfOpt:
        downloadPDF(ref,type,"")

    if psOpt:
        downloadPS(ref,type,"")

    if sourceOpt:
        downloadSource(ref,type,"")