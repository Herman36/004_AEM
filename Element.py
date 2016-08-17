# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:46:39 2016

@author: herman
"""

class Element:
    'class for all element information and properties'
    ele_no = 0
    b = 10
    a = 10
    edge_len = {1:a,2:b,3:a,4:b}
    T = 10
    E = 200e+9
    nu = 0
    G = E/(2*(1+nu))
    x = 0           #matrix of x co-ordinates on edge of element where springs will be attached, rows are edges
    y = 0           #matrix of y co-ordinates on edge of element where springs will be attached, rows are edges
    nodes = [0, 0, 0, 0]    #corner nodes of element
    spring_num = 2       #number of spring contact points on element
    r = 0           #rotational angle in radians of element (boundary condition)
    u = 0           #displacements of element (boundary condition)
    
    def __init__(self):
        Element.ele_no += 1
    
    def getPlotMat(self):
        Element.nodes = []
        Element.x = 
        Element.y = 