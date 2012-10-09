#!/bin/env python
"""
Author: Adam Morgan (qmorgan@gmail.com)

You need a twilio account to use this: sign up for a free trial here:
 * https://www.twilio.com/
 * pick a new phone number to use for outgoing twilio calls
 * input your current phone number to verify your new account
 * Note - for free trial calls, you can only call to this verified phone number;
   to call to other numbers, you'll need a paid account.

Running this script requires twilio_python and updated version of httplib2

sudo pip install twilio
sudo pip install --upgrade  httplib2

The code needs your twilio id and token to be used.
(While you can call up your ID and Token directly within TwilioRestClient, 
its better to set them up as environment variables so you do not accidentally
send them to people)

You can get these numbers from https://www.twilio.com/user/account after you 
have signed up.  

Set up your environment variables permanently in e.g. your ~/.bashrc file:
export TWILIO_ACCOUNT_SID='yourIDhere'
export TWILIO_AUTH_TOKEN='yourTokenHere'
"""

from twilio.rest import TwilioRestClient
import os, sys

# Check to make sure the environment variables are set up correctly
if not os.environ.has_key("TWILIO_ACCOUNT_SID"):
    print "You need to set the environment variable TWILIO_ACCOUNT_SID to"
    print "point to your twilio account id"
    sys.exit(1)

if not os.environ.has_key("TWILIO_AUTH_TOKEN"):
    print "You need to set the environment variable TWILIO_AUTH_TOKEN to"
    print "point to your twilio authorization token"
    sys.exit(1)

# set up the client using the environment variables (it searches for these by
# default, but you can define them explicitly in this function if you want)
client=TwilioRestClient()


def _usage():
    print "[USAGE] qCall.py callto callfrom message"
    print "[EXAMPLE] qCall.py 5551234567 5558675309 hello how are you today?"
    return

def call(callto,callfrom,messages):
    '''
    Place a phone call from a twilio number to another number, and have the 
    body of the call be user-defined text-to-speech.
    
    This uses the "message" twimlet web application; see http://twimlets.com
    
    callto: number to place the call to (for free twilio accounts, this must 
        be your verified #)
    callfrom: your twilio number associated with your ID and AUTH_TOKEN
    messages: either a string or a list of strings containing either:
        * a text message to read aloud using text-to-speech
        * a URL to an audio message to play
    
    '''
    # if messages is just a single string
    if isinstance(messages,str): 
        outstring="Message%5B0%5D="+messages.replace(" ","+")
    # if messages is a list of strings
    elif isinstance(messages,list):
        outstring=''
        count=0
        for message in messages:
            message=str(message) #convert to string if not already
            outstring+="Message%5B"
            outstring+=str(count)
            outstring+="%5D="
            outstring+=message.replace(" ","+")+"&"
            count+=1
    else:
        raise ValueError("messages needs to be a list or a string")
        
    url="http://twimlets.com/message?"+outstring
    print url
    call=client.calls.create(to=callto,from_=callfrom,url=url)
    
    ## old version using the echo twimlet
    # textconverted=text.replace(" ","+") 
    # url="http://twimlets.com/echo?Twiml=%3CResponse%3E%3CSay%3E" +\
    #     textconverted + "%3C%2FSay%3E%3C%2FResponse%3E"
    

def sms(callto,callfrom,body):
    '''
    Send a text message from a twilio number to another number.
        
    callto: number to send the sms (for free twilio accounts, this must 
        be your verified #)
    callfrom: your twilio number associated with your ID and AUTH_TOKEN
    body: a string, or something that can be converted to a string. 
        Must be fewer than 160 characters.
    
    '''
    message = client.sms.messages.create(to=callto, from_=callfrom, body=body)

def test_call(callto="5102297683",callfrom="5102503825"):
    """
    Test the calling functionality call()
    Replace the callto and callfrom with your numbers here to test the 
    application"""
    # text = "A+G.R.B.+has+occured+at+coordinates+R.A.+21+hours+32+minutes+23"
    # text += "+seconds+declination+positive+23+degrees+23+minutes+55+seconds."
    text = "Here is a nice song for you to enjoy."
    mp3 = "http://dl.dropbox.com/u/2484822/boil.mp3"
    messages = [text,mp3]
    call(callto,callfrom,messages)

def test_sms(callto="5102297683",callfrom="5102503825"):
    """
    Test the sms sending capability sms()
    Replace the callto and callfrom with your numbers here to test the 
    application"""
    text = "How doth the little crocodile improve his shining tail"
    sms(callto,callfrom,text)


if __name__ == '__main__':
        
    if ((sys.argv[1] == '-h') or (sys.argv[1] == "--h")):
        print _usage()
        sys.exit(0)
    
    # Change if want to be able to send attachments from command line    
    call(sys.argv[1],sys.argv[2],sys.argv[3:])
    sys.exit(0)