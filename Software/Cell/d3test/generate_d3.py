import glob
import os
import pandas as pd
import numpy as np

df=pd.read_csv("/Users/amorgan/qsoft/Software/Cell/d3test/Detections.csv")

columnlist = ['id','x','y','r_kron','A','B','theta']  

data_str = '['
for ind in np.arange(int(df.describe().z_index.loc['max'])):
    subdf = df[df.z_index == ind]
    if len(subdf) == 0:
        data_str += '[],\n' # no data for this slice
    else:
        #define formatters to print
        id_fmt  = lambda x: '\t{id:%i,' % x
        a_fmt = lambda x: 'a:%.2f,' % x
        b_fmt = lambda x: 'b:%.2f,' % x
        th_fmt = lambda x: 'th:%.1f},' % x
        x_fmt = lambda x: 'x:%.2f,' % x
        y_fmt = lambda x: 'y:%.2f,' % x
        r_fmt = lambda x: 'r:%.2f,' % x
    
        formatterdict = {'id':id_fmt,'x':x_fmt,'y':y_fmt,'r_kron':r_fmt,
                'A':a_fmt,'B':b_fmt,'theta':th_fmt}
    
    
        out_str = subdf.to_string(index=False,index_names=False,
                header=False,columns=columnlist,formatters=formatterdict)
            
        out_str = '[' + out_str.rstrip(',') + '],\n'# get rid of trailing comma
        data_str += out_str

data_str = data_str.rstrip(',')+']'