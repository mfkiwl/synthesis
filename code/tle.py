"""
Created on Fri Oct  7 11:59:53 2022

@author: 3Q
"""

import pandas as pd
import requests

def tle_json(file_path,url):
    
    tle = requests.get(url)
    data = tle.text
    with open(file_path,'w') as tle:
        tle.write(data)
    
    data = []
    with open(file_path,'r') as tle:
        data = tle.readlines()
        
    tle = []
    i = 0
    for line in data:
        while i < len(data):
            j = i
            tle.append({data[i].rstrip('\n').rstrip(' '):[data[j+2].rstrip('\n'),data[j+4][:len(data[j+4])-1]]})
            i += 6
    
    return tle


