'''
q_insert.py

For use on betsy.  
'''

import os
import MySQLdb
import sys
import ephem


def grab_script(filepath,rm=True):
    '''Grab one of the reduction scripts made by PTEL_data.py, loop around it,
    and insert its contents into the reduction queue.  If rm=True, move the
    script to the finished_scripts directory
    '''
    
    print "Inserting into queue: " + os.path.basename(filepath)
    
    finished_scripts_path = '/home/ptelreducer/storage/amorgan/redux/scripts/finished_scripts/'

    if not os.path.exists(filepath):
        print 'File path does not exist. Exiting.'
        return
    
    f = file(filepath,'r')
    for line in f:
        if line[0:5] == 'mkdir':
            os.system(line)
        elif line[0:5] == 'pytho':
            queue_insert(line)
        else:
            print 'Unrecognized line in script!'
            print line
    if rm:
        if not os.path.exists(finished_scripts_path):
            os.system('mv ' + filepath + ' ' + finished_scripts_path + '.')
        else:
            print "**FINISHED SCRIPTS PATH ALREADY EXISTS.  Not moving."
        

def queue_insert(reduction_command):
    '''Insert a Pipeline 3 reduction_command into the reduction queue'''
    que_conn = MySQLdb.connect(host="localhost",user="ptelreducer",
        passwd="ilove2mass",db="reduction_que",connect_timeout=60)
    current_datetime = ephem.date(ephem.now())
    current_datetime = str(current_datetime).replace("/","-")
    insert_str = ('''INSERT INTO que (reduction_cmd, request_entered_time,
        request_started_time, request_completed_time, request_priority) values
        ("''' + reduction_command + '''", "''' + current_datetime + '''",
        "0000-00-00 00:00:00", "0000-00-00 00:00:00", 3)''')
    cursor = que_conn.cursor()
    cursor.execute(insert_str)
    cursor.close()
    que_conn.close()
    
def queue_insert_override(reduction_command):
    '''Insert a Pipeline 3 reduction_command into the reduction queue
    Makes it highest priority for immediate reduction.
    ONLY USE FOR SUPER IMPORTANT TARGETS.
    '''
    que_conn = MySQLdb.connect(host="localhost",user="ptelreducer",
        passwd="ilove2mass",db="reduction_que",connect_timeout=60)
    current_datetime = ephem.date(ephem.now())
    current_datetime = str(current_datetime).replace("/","-")
    insert_str = ('''INSERT INTO que (reduction_cmd, request_entered_time,
        request_started_time, request_completed_time, request_priority) values
        ("''' + reduction_command + '''", "''' + current_datetime + '''",
        "0000-00-00 00:00:00", "0000-00-00 00:00:00", 1)''')
    cursor = que_conn.cursor()
    cursor.execute(insert_str)
    cursor.close()
    que_conn.close()