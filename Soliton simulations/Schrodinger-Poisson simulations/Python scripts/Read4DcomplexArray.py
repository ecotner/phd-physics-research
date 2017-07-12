# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 18:20:02 2015

@author: 27182_000
"""

from numpy import *

def read_timestepC(g,M):
    "reads in a single timestep from the data stored in a 4D array in a text file. This data is stored in a complex ndarray object"
    s = ","
    nlist = ones((M+1)**3, dtype=complex128)
    for i in range((M+1)**3):
        buff = ""
        while (s == "[" or s == "]" or s == ","):
            s = g.read(1)
        while (s != "[" and s != "]" and s != ","):
            buff = buff + s
            s = g.read(1)
        nlist[i] = complex(buff)
    nlist = nlist.reshape(M+1,M+1,M+1)
    return nlist