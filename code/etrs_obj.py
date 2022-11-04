# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 00:55:52 2022

@author: TM
"""
import os
import pandas as pd
from pyproj import Transformer


#print(os.getcwd())
os.chdir('../data/')
tri = []
with open('LoD2_32_418_5653_1_NW.obj','r') as obj:
    tri = obj.readlines()

tri = list(filter(lambda i: i.startswith('#') == False, tri))
vertices = list(filter(lambda i: i.startswith('v') == True, tri))
faces = list(filter(lambda i: i.startswith('f') == True, tri))
vertices = list(map(lambda st: st.strip('\n').split(' ')[1:],vertices))
vertices = list(map(lambda st: [float(st[0]),float(st[1]),float(st[2])],vertices))

d = {'xyz': [[1,2]]*len(vertices), 'v': [[1,2]]*len(vertices)}
df = pd.DataFrame(d)
print(df.head())

for i in range(0,len(vertices)):
    df.at[i,'xyz'] = vertices[i]
    
print(df.head())    
print('x {}'.format(df['xyz'][0]))
print('v {}'.format(df['v'][0])) 
transformer = Transformer.from_pipeline('''+proj=pipeline
                                          +step +inv +proj=utm +zone=32 +ellps=GRS80
                                          +step +proj=cart +ellps=GRS80''')
df['v'] = df['xyz'].apply(lambda x: transformer.transform(x[0],x[1],x[2]))
lst = list(df['v'])
with open('testmesh.obj','w') as mesh:
    for ver in lst:
        mesh.write('v {} {} {}\n'.format(ver[0],ver[1],ver[2]))
    for f in faces:
        mesh.write(f)




