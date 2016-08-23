# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:47:36 2016

@author: herman
"""

# Basic AEM for two elements
    # element numbering from 1,2...
    # Note  all comments marked as #& indicate changes to be made in AEM_2D.py
    #       all comments marked as #! are concerns to be addressed
    #       all comments marked as ## indicates start and end of sections
    #       all comments marked as #@ working comments

import numpy as np
from math import cos, sin

## Start of input deck
# Material properties
E = 207e+9
nu = 0.3
G = E/(2*(1+nu))
T = 0.15

# Element properties
#& Should be adjusted to load from hex_mesher
grid = [10, 2, 0]
num_ele = grid[0]*grid[1]
el_center = np.zeros(shape=(num_ele, 3))
el_count = 0
for y_c in range(50, 200, 100):
    for x_c in range(100, 2000, 200):
        el_center[el_count, 0] = el_count
        el_center[el_count, 1] = float(x_c)/1000
        el_center[el_count, 2] = float(y_c)/1000
        el_count += 1


a1 = 0.2         # width of element 1
b1 = 0.100         # height of element 1
theta1 = 0.      # radians (for time being not applicable)

#@ Seems redundent unless different size and shape elements
#& Remove in future releases
a2 = 0.2         # width of element 2
b2 = 0.100         # height of element 2
theta2 = 0.      # radians (for time being not applicable)

gap = [0.00, 0.00, 0]         # size of gap between elements

# Spring properties
#@ Assuming same number of connecting springs for horizontal and vetical
num_spring = 10

# Prescribed displacements
    # element number, DOF, value
Up = np.array([0, 0, 0, 0, 0, 0])
#@ Needs to be done to read from input
dof = np.ones(num_ele*3)
dof[0] = 0
dof[1] = 0
dof[2] = 0
dof[30] = 0
dof[31] = 0
dof[32] = 0
#dof[241] = 0
#dof[360] = 0
#dof[361] = 0
#dof[480] = 0
#dof[481] = 0
#dof[600] = 0
#dof[601] = 0
#dof[720] = 0
#dof[721] = 0
#dof[840] = 0
#dof[841] = 0
#dof[960] = 0
#dof[961] = 0
#dof[1080] = 0
#dof[1081] = 0


# Prescribed forces
    # element number, DOF, value
F = np.zeros(shape=(num_ele*3, 1))
#F[118, 0] = -1.0e+3
#F[238, 0] = -0.50e+3
#F[358, 0] = -0.50e+3
#F[118, 0] = -1.0e+3
F[27, 0] = 0.50e+5
F[57, 0] = 0.50e+5

## End of input deck
def AEM_basic(E, nu, G, T, grid, num_ele,
              el_center, a1, a2, b1, b2, theta1, theta2, gap,
              num_spring, Up, dof, F):
    #-----------------------------------------------------------------------------#
    ## Start of calculating stiffness matricies [Prashidha 2014]
    
    #@ Insert incrementation here
    #! Include [if] for element with smallest height
    #! Question 1&2: Are Kn and Ks different if element sizes are different?
        # which element determines the value of d and z if element sizes differ?
        # unless calculating seperate Kn and Ks like for vertical and horizontal?
    #@ Assume all elements of same size currently
    #@ For horizontal element interaction
    dh = float(b1)/(num_spring)
    #! Not sure how gap plays role here
    ah = a1/2 + a2/2
    Knh = (E*dh*T)/ah
    Ksh = (G*dh*T)/ah
    
    #@ For vertical element interaction
    dv = float(a1)/(num_spring)
    #! Not sure how gap plays role here
    av = b1/2 + b2/2
    Knv = (E*dv*T)/av
    Ksv = (G*dv*T)/av
    #-----------------------------------------------------------------------------#
    
    
    def Kte_row(Knh, Ksh, a1, a2, bn1, bn2):
        Kte = np.array([[Knh, 0, -Knh*bn1, -Knh, 0, Knh*bn2],
                       [0, Ksh, Ksh*(a1/2), 0, -Ksh, Ksh*(a2/2)],
                       [-Knh*bn1, Ksh*(a1/2),
                        Knh*bn1**2+Ksh*(a1/2)**2, Knh*bn1,
                        -Ksh*(a1/2), -Knh*bn1*bn2+Ksh*(a1/2)*(a2/2)],
                       [-Knh, 0, Knh*bn1, Knh, 0, -Knh*bn2],
                       [0, -Ksh, -Ksh*(a1/2), 0, Ksh, -Ksh*(a2/2)],
                       [Knh*bn2, Ksh*(a2/2),
                        (Knh*bn1*bn2)+Ksh*(a1/2)*(a2/2),
                        -Knh*bn2, -Ksh*(a2/2),
                        Knh*bn2**2+Ksh*(a2/2)**2]])
        return Kte
    #-----------------------------------------------------------------------------#
    
    
    def Kte_column(Knv, Ksv, an1, an2, b1, b2):
        Kte = np.array([[Ksv, 0, -Ksv*(b1/2), -Ksv, 0, -Ksv*(b1/2)],
                        [0, Knv, -Knv*an1, 0, -Knv, Knv*an2],
                        [-Ksv*(b1/2), -Knv*an1,
                         Knv*an1**2+Ksv*(b1/2)**2, Ksv*b1,
                         Knv*an1, -Knv*an1*an2+Ksv*(b1/2)*(b2/2)],
                        [-Ksv, 0, Ksv*(b1/2), Ksv, 0, Ksv*(b2/2)],
                        [0, -Knv, Knv*an1, 0, Knv, -Knv*an2],
                        [-Ksv*(b2/2), Knv*an2,
                         -Ksv*an1*an2+Ksv*(b1/2)*(b2/2),
                         Ksv*(b2/2), -Knv*an2,
                         Knv*an2**2+Ksv*(b2/2)**2]])
        return Kte
    #-----------------------------------------------------------------------------#
    #! Unsure of the workings of theta1 and theta2
    #! Question 3: Setting up the stiffness matrices, these do not change unless a
                    # spring is no longer active
    #@ Note stiffness matrices are compiled per spring thus per element pair
    ## Start loops over element pairs
    Kg = np.zeros(shape=(num_ele*3, num_ele*3))
    elemn = 0
    # Translation matrix:
    #! Check use of theta2
    L = np.array([[cos(theta1), sin(theta1), 0, 0, 0, 0],
                  [-sin(theta1), cos(theta1), 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cos(theta2), sin(theta2), 0],
                  [0, 0, 0, -sin(theta2), cos(theta2), 0],
                  [0, 0, 0, 0, 0, 1]])
    LT = L.T
    for row in range(0, grid[1]-1):
        for column in range(0, grid[0]-1):
            ## Start of horizontal element pair interaction
            Kel = np.zeros(shape=(6, 6))
            for sign in range(1, 3):
                for n in range(1, (num_spring/2)+1):
                    bn1 = ((-1)**sign)*(dh/2 + dh*(n-1))
                    bn2 = ((-1)**sign)*(dh/2 + dh*(n-1))
                    Kte = Kte_row(Knh, Ksh, a1, a2, bn1, bn2)
                    Kel = Kel + Kte     # adding each spring's contribution
            Ke = LT.dot(Kel.dot(L))
            ## End of horizontal element pair interaction
            ## Adding element pair stiffness to global stiffness
            pg_r = ([3*elemn, 3*elemn+1, 3*elemn+2,
                     3*(elemn+1), 3*(elemn+1)+1, 3*(elemn+1)+2])
            Kg[np.ix_(pg_r, pg_r)] = Kg[np.ix_(pg_r, pg_r)] + Ke
    
            ## Start of vertical element pair interaction
            Kel = np.zeros(shape=(6, 6))
            for sign in range(1, 3):
                for n in range(1, (num_spring/2)+1):
                    an1 = ((-1)**sign)*(dv/2 + dv*(n-1))
                    an2 = ((-1)**sign)*(dv/2 + dv*(n-1))
                    Kte = Kte_column(Knv, Ksv, an1, an2, b1, b2)
                    Kel = Kel + Kte     # adding each spring's contribution
    
            Ke = LT.dot(Kel.dot(L))
            ## End of vertical element pair interaction
            ## Adding element pair stiffness to global stiffness
            pg_c = ([3*elemn, 3*elemn+1, 3*elemn+2,
                     3*(elemn+grid[0]), 3*(elemn+grid[0])+1, 3*(elemn+grid[0])+2])
            Kg[np.ix_(pg_c, pg_c)] = Kg[np.ix_(pg_c, pg_c)] + Ke
            elemn += 1
    
        ## Start of vertical element pair interaction
        Kel = np.zeros(shape=(6, 6))
        for sign in range(1, 3):
            for n in range(1, (num_spring/2)+1):
                an1 = ((-1)**sign)*(dv/2 + dv*(n-1))
                an2 = ((-1)**sign)*(dv/2 + dv*(n-1))
                Kte = Kte_column(Knv, Ksv, an1, an2, b1, b2)
                Kel = Kel + Kte     # adding each spring's contribution
    
        Ke = LT.dot(Kel.dot(L))
        ## End of vertical element pair interaction
        ## Adding element pair stiffness to global stiffness
        pg_c = ([3*elemn, 3*elemn+1, 3*elemn+2,
                 3*(elemn+grid[0]), 3*(elemn+grid[0])+1, 3*(elemn+grid[0])+2])
        Kg[np.ix_(pg_c, pg_c)] = Kg[np.ix_(pg_c, pg_c)] + Ke
        elemn += 1
    for column in range(0, grid[0]-1):
        ## Start of horizontal element pair interaction
        Kel = np.zeros(shape=(6, 6))
        for sign in range(1, 3):
            for n in range(1, (num_spring/2)+1):
                bn1 = ((-1)**sign)*(dh/2 + dh*(n-1))
                bn2 = ((-1)**sign)*(dh/2 + dh*(n-1))
                Kte = Kte_row(Knh, Ksh, a1, a2, bn1, bn2)
                Kel = Kel + Kte     # adding each spring's contribution
        Ke = LT.dot(Kel.dot(L))
        ## End of horizontal element pair interaction
        ## Adding element pair stiffness to global stiffness
        pg_r = ([3*elemn, 3*elemn+1, 3*elemn+2,
                 3*(elemn+1), 3*(elemn+1)+1, 3*(elemn+1)+2])
        Kg[np.ix_(pg_r, pg_r)] = Kg[np.ix_(pg_r, pg_r)] + Ke
        elemn += 1
    ## End of calculating stiffness matrices
    #-----------------------------------------------------------------------------#
    
    ## Start of calculating unknown degrees of freedom and reaction forces
    
    # Next the stiffness matrix will be devided so that the free and
        # prescribed degree of freedoms can be idependantly used to calculate
        # the unknown variables
    
    # Loads are applied at the centre of an element on the degrees of freedom
    
    # Solve displacements
    U = np.zeros(shape=(num_ele*3, 1))
    # To solve for the free displacements the stiffness matrix needs to be devided
        # this is done by clasifying the free and prescribed degrees of freedom
    pdof = np.where(dof == 0)[0]
    fdof = np.nonzero(dof)[0]
    
    Kff = Kg[fdof][:, fdof]
    Kfp = Kg[fdof][:, pdof]
    Kpp = Kg[pdof][:, pdof]
    
    # Now we include the effects of prescribed displacements
    Fp1 = F[fdof, 0] - Kfp.dot(Up)
    # Solve unknown degrees of freedom
    Uf = np.linalg.solve(Kff, Fp1)
    # Finally solve for the support reactions
    KfpT = Kfp.T
    Fp = KfpT.dot(Uf) + Kpp.dot(Up)
    # Placing calculated values back into global F and U
    U[fdof, 0] = Uf
    U[pdof, 0] = Up
    F[pdof, 0] = Fp
    ## End of calculating unknown degrees of freedom and reaction forces
    #-----------------------------------------------------------------------------#
    
    
    #-----------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------#
    # Function to calculate the new co-ordinates of each spring
    def new_coordinates(num_spring, ele, x_co, y_co, U, dof1, dof2):
        Li = np.zeros(shape=(num_spring, 2))
        x_co_n = np.zeros(shape=(num_spring, 2))  # [spring num,element num]
        y_co_n = np.zeros(shape=(num_spring, 2))  # [spring num,element num]
        x_co_n[:, 0] = x_co[:, 0] + U[dof1, 0]
        x_co_n[:, 1] = x_co[:, 1] + U[dof2, 0]
        y_co_n[:, 0] = y_co[:, 0] + U[dof1+1, 0]
        y_co_n[:, 1] = y_co[:, 1] + U[dof2+1, 0]
        for spr in range(0, num_spring):
            x_co_n[spr, 0] = (x_co_n[spr, 0]*cos(U[dof1+2, 0])
                              - y_co_n[spr, 0]*sin(U[dof1+2, 0]))
            x_co_n[spr, 1] = (x_co_n[spr, 1]*cos(U[dof2+2, 0])
                              - y_co_n[spr, 1]*sin(U[dof2+2, 0]))
            y_co_n[spr, 0] = (x_co_n[spr, 0]*sin(U[dof1+2, 0])
                              + y_co_n[spr, 0]*cos(U[dof1+2, 0]))
            y_co_n[spr, 1] = (x_co_n[spr, 1]*sin(U[dof2+2, 0])
                              + y_co_n[spr, 1]*cos(U[dof2+2, 0]))
            Li[spr, 0] = x_co_n[spr, 1] - x_co_n[spr, 0]    # New length of x
            Li[spr, 1] = y_co_n[spr, 1] - y_co_n[spr, 0]    # New length of y
        return Li, x_co_n, y_co_n
    #-----------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------#
    
    
    #-----------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------#
    # Function to add element strains
    def element_strain_column(strain_ele, strain_column, ele, grid):
        strain_ele[0, ele] = (strain_ele[0, ele]
                              + np.average(strain_column[0, :]))
        strain_ele[1, ele] = (strain_ele[1, ele]
                              + np.average(strain_column[1, :]))
        strain_ele[2, ele] = (strain_ele[2, ele]
                              + np.average(strain_column[2, :]))
    
        strain_ele[0, ele+grid[0]] = (strain_ele[0, ele+grid[0]]
                                      + np.average(strain_column[0, :]))
        strain_ele[1, ele+grid[0]] = (strain_ele[1, ele+grid[0]]
                                      + np.average(strain_column[1, :]))
        strain_ele[2, ele+grid[0]] = (strain_ele[2, ele+grid[0]]
                                      + np.average(strain_column[2, :]))
        return strain_ele
    
    
    def element_strain_row(strain_ele, strain_row, ele, grid):
        strain_ele[0, ele] = (strain_ele[0, ele]
                              + np.average(strain_row[0, :]))
        strain_ele[1, ele] = (strain_ele[1, ele]
                              + np.average(strain_row[1, :]))
        strain_ele[2, ele] = (strain_ele[2, ele]
                              + np.average(strain_row[2, :]))
    
        strain_ele[0, ele+1] = (strain_ele[0, ele+1]
                                + np.average(strain_row[0, :]))
        strain_ele[1, ele+1] = (strain_ele[1, ele+1]
                                + np.average(strain_row[1, :]))
        strain_ele[2, ele+1] = (strain_ele[2, ele+1]
                                + np.average(strain_row[2, :]))
        return strain_ele
    
    #-----------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------#
    
    ## Start of calculating stresses and strains
    ele = 0
    strain_ele = np.zeros(shape=(3, grid[0]*grid[1]))   # [strain; ele]
    
    L0_row = np.zeros(shape=(num_spring, 2))
    dL_row = np.zeros(shape=(num_spring, 2))
    strain_row = np.zeros(shape=(3, num_spring))
    
    L0_column = np.zeros(shape=(num_spring, 2))
    dL_column = np.zeros(shape=(num_spring, 2))
    strain_column = np.zeros(shape=(3, num_spring))
    ## Loop over elements
    for row in range(0, grid[1]-1):
        for column in range(0, grid[0]-1):
            # Calculating the original co-ordinates of each spring for horizontal
            x_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
            y_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
            spr = 0
            for sign in range(1, 3):
                for n in range(1, (num_spring/2)+1):
                    bn1 = ((-1)**sign)*(dh/2 + dh*(n-1))
                    bn2 = ((-1)**sign)*(dh/2 + dh*(n-1))
                    x_co[spr, 0] = el_center[ele, 1] + a1/2
                    x_co[spr, 1] = el_center[ele+1, 1] - a2/2
                    y_co[spr, 0] = el_center[ele, 2] + bn1
                    y_co[spr, 1] = el_center[ele+1, 2] + bn2
                    L0_row[spr, 0] = x_co[spr, 1] - x_co[spr, 0]    # Original L N
                    L0_row[spr, 1] = y_co[spr, 1] - y_co[spr, 0]    # Original L S
                    spr += 1
            dof1 = ele*3
            dof2 = (ele+1)*3
            Li_row, x_co_n, y_co_n = new_coordinates(num_spring, ele, x_co, y_co, U, dof1, dof2)
            dL_row[:, 0] = Li_row[:, 0] - L0_row[:, 0]  # x direction
            dL_row[:, 1] = Li_row[:, 1] - L0_row[:, 1]  # y direction
            strain_row[0, :] = dL_row[:, 0]/(L0_row[:, 0] + ah)
            strain_row[2, :] = dL_row[:, 1]/(L0_row[:, 0] + ah)
    
            # Calculating the original co-ordinates of each spring for vertical
            x_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
            y_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
            spr = 0
            for sign in range(1, 3):
                for n in range(1, (num_spring/2)+1):
                    an1 = ((-1)**sign)*(dv/2 + dv*(n-1))
                    an2 = ((-1)**sign)*(dv/2 + dv*(n-1))
                    x_co[spr, 0] = el_center[ele, 1] + an1
                    x_co[spr, 1] = el_center[ele+grid[0], 1] + an2
                    y_co[spr, 0] = el_center[ele, 2] + b1/2
                    y_co[spr, 1] = el_center[ele+grid[0], 2] - b2/2
                    L0_column[spr, 0] = x_co[spr, 1] - x_co[spr, 0]  # Original L S
                    L0_column[spr, 1] = y_co[spr, 1] - y_co[spr, 0]  # Original L N
                    spr += 1
            dof1 = ele*3
            dof2 = (ele+grid[0])*3
            Li_column, x_co_n, y_co_n = new_coordinates(num_spring, ele, x_co, y_co, U, dof1, dof2)
            dL_column[:, 0] = Li_column[:, 0] - L0_column[:, 0]     # x direction
            dL_column[:, 1] = Li_column[:, 1] - L0_column[:, 1]     # y direction
            strain_column[1, :] = dL_column[:, 1]/(L0_column[:, 0] + av)
            strain_column[2, :] = dL_column[:, 0]/(L0_column[:, 0] + av)
    
            strain_ele = element_strain_column(strain_ele,
                                               strain_column, ele, grid)
            strain_ele = element_strain_row(strain_ele, strain_row, ele, grid)
            ele += 1
        # Calculating the original co-ordinates of each spring for vertical
        x_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
        y_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
        spr = 0
        for sign in range(1, 3):
            for n in range(1, (num_spring/2)+1):
                an1 = ((-1)**sign)*(dv/2 + dv*(n-1))
                an2 = ((-1)**sign)*(dv/2 + dv*(n-1))
                x_co[spr, 0] = el_center[ele, 1] + an1
                x_co[spr, 1] = el_center[ele+grid[0], 1] + an2
                y_co[spr, 0] = el_center[ele, 2] + b1/2
                y_co[spr, 1] = el_center[ele+grid[0], 2] - b2/2
                L0_column[spr, 0] = x_co[spr, 1] - x_co[spr, 0]  # Original L S
                L0_column[spr, 1] = y_co[spr, 1] - y_co[spr, 0]  # Original L N
                spr += 1
        dof1 = ele*3
        dof2 = (ele+grid[0])*3
        Li_column, x_co_n, y_co_n = new_coordinates(num_spring, ele, x_co, y_co, U, dof1, dof2)
        dL_column[:, 0] = Li_column[:, 0] - L0_column[:, 0]
        dL_column[:, 1] = Li_column[:, 1] - L0_column[:, 1]
        strain_column[1, :] = dL_column[:, 1]/(L0_column[:, 0] + av)
        strain_column[2, :] = dL_column[:, 0]/(L0_column[:, 0] + av)
    
        strain_ele = element_strain_column(strain_ele,
                                           strain_column, ele, grid)
    
        ele += 1
    for column in range(0, grid[0]-1):
        # Calculating the original co-ordinates of each spring for horizontal
        x_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
        y_co = np.zeros(shape=(num_spring, 2))  # [spring num, element num]
        spr = 0
        for sign in range(1, 3):
            for n in range(1, (num_spring/2)+1):
                bn1 = ((-1)**sign)*(dh/2 + dh*(n-1))
                bn2 = ((-1)**sign)*(dh/2 + dh*(n-1))
                x_co[spr, 0] = el_center[ele, 1] + a1/2
                x_co[spr, 1] = el_center[ele+1, 1] - a2/2
                y_co[spr, 0] = el_center[ele, 2] + bn1
                y_co[spr, 1] = el_center[ele+1, 2] + bn2
                L0_row[spr, 0] = x_co[spr, 1] - x_co[spr, 0]    # Original L N
                L0_row[spr, 1] = y_co[spr, 1] - y_co[spr, 0]    # Original L S
                spr += 1
        dof1 = ele*3
        dof2 = (ele+1)*3
        Li_row, x_co_n, y_co_n = new_coordinates(num_spring, ele, x_co, y_co, U, dof1, dof2)
        dL_row[:, 0] = Li_row[:, 0] - L0_row[:, 0]  # x direction
        dL_row[:, 1] = Li_row[:, 1] - L0_row[:, 1]  # y direction
        strain_row[0, :] = dL_row[:, 0]/(L0_row[:, 0] + ah)
        strain_row[2, :] = dL_row[:, 1]/(L0_row[:, 0] + ah)
    
        strain_ele = element_strain_row(strain_ele, strain_row, ele, grid)
        ele += 1
    
    #@ Stress calculation per element
    stress_Pstrain = np.zeros(shape=(3, grid[0]*grid[1]))
    stress_Pstress = np.zeros(shape=(3, grid[0]*grid[1]))
    C_strain = np.array([[1-nu, nu, 0],
                         [nu, 1-nu, 0],
                         [0, 0, 1-2*nu]])
    C_stress = np.array([[1, nu, 0],
                         [nu, 1, 0],
                         [0, 0, 1-nu]])
    
    stress_Pstrain = (E/((1+nu)*(1-2*nu))) * C_strain.dot(strain_ele)
    stress_Pstress = (E/(1-nu**2)) * C_stress.dot(strain_ele)
    ## End of calculating stresses and strains
    return U, F, strain_ele, stress_Pstrain, stress_Pstress
U, F, strain_ele, stress_Pstrain, stress_Pstress = AEM_basic(E, nu, G, T, grid, num_ele,
              el_center, a1, a2, b1, b2, theta1, theta2, gap,
              num_spring, Up, dof, F)
# Draw deformed shape
#print('U=')
#print(U)
#print('F=')
#print(F)
#print('strain=')
#print(strain_ele)
#print('Pstrain=')
#print(stress_Pstrain)
#print('Pstress=')
#print(stress_Pstress)
