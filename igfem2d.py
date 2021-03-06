from scipy import  sparse

import scipy.sparse.linalg.dsolve as linsolve
import math
from quad2d import quad2d
from isoparam import *
from probstatement import *
import numpy as np
from findIntersection import *
import itertools
import sympy
from quad2d import *
import scipy.io

from libFcts import *

from main import *
from scipy import interpolate
import time

import __builtin__
sum = __builtin__.sum


# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This fuction
# returns True or False.  The algorithm is called
# "Ray Casting Method".

def point_in_poly(x,y,poly):
    n = len(poly)
    inside = False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside

def point_in_on_poly(x,y,poly):
    inside = point_in_poly(x,y,poly)
    n = len(poly)
    c = Point(x,y)
    for i in range(0,n-1):
        p1x,p1y = poly[i]
        a = Point(p1x,p1y)
        for j in range(i+1,n):
            p2x,p2y = poly[j]
            b = Point(p2x,p2y)
            point_on_line = is_on(a,b,c)
            if point_on_line == True:
                return True
    return inside

def number_nodes(m,pp,tt):
  nr_of_nodes = (m+1)**2 - 1
  T = m*m

  for e in range(0,T):
    nodes = tt[e]
    coords = pp[nodes,:]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    for k in range(nr_of_nodes+1,len(pp)):
      in_poly = point_in_square(pp[k][0],pp[k][1],x0,y0,x1,y1)
      if in_poly == True:
        tt[e] =  tt[e] + [k]

  return tt

def point_in_square(x,y,x0,y0,x1,y1):

  if (x0 <= x and x <= x1) and (y0<=y and y <= y1):# and 
    # (x != x0 and y != y0) and :
    return True
  else:
    return False

def is_on(a, b, c):
    "Return true iff point c intersects the line segment from a to b."
    # (or the degenerate case that all 3 points are coincident)
    return (collinear(a, b, c)
        and (within(a.x, c.x, b.x) if a.x != b.x else 
        within(a.y, c.y, b.y)))

def collinear(a, b, c):
    "Return true iff a, b, and c all lie on the same line."
    return (b.x - a.x) * (c.y - a.y) == (c.x - a.x) * (b.y - a.y)

def within(p, q, r):
    "Return true iff q is between p and r (inclusive)."
    return p <= q <= r or r <= q <= p

# Use Simpson Rule for integration over isoparametric elements
def simpson_rule(f):

    I = 1.0/6.0 * f(1.0/2.0,1.0/2.0) + 1.0/6.0 * f(1.0/2.0,0.0) + 1.0/6.0 * f(0.0,1.0/2.0)
    return I

def gauss_quadrat(f,ui,wi):
    J = 0
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            J = J + wi[i] * wi[j] * f(ui[i],ui[j])

    return J

def my_gauss_rule(f,ui,wi):
# Gauss Rule with 3 points: exact for polynomial of degree 2

    if POL_APPROX == 0:
        return 1.0/2.0 * f(1.0/3.0, 1.0/3.0)
    if POL_APPROX == 1:
        return simpson_rule3(f)
    if POL_APPROX == 2:
        return simpson_rule3(f)

def simpson_rule3(f):
# Gauss Rule with 4 points: exact for polynomial of degree 3
    I = -27.0/96.0 * f(1.0/3.0,1.0/3.0) + 25.0/96.0 * f(0.2,0.6) + 25.0/96.0 * f(0.2, 0.2) + 25.0/96.0 * f(0.6, 0.2)
    return I

def simpson_rule4(f):
    
    I = 0
    I = I + f(0.44594849091597, 0.44594849091597) * 0.22338158967801 * 1.0/2.0
    I = I + f(0.44594849091597, 0.10810301816807) * 0.22338158967801 * 1.0/2.0
    I = I + f(0.10810301816807, 0.44594849091597) * 0.22338158967801 * 1.0/2.0
    I = I + f(0.09157621350977, 0.09157621350977) * 0.10995174365532 * 1.0/2.0
    I = I + f(0.09157621350977, 0.81684757298046) * 0.10995174365532 * 1.0/2.0
    I = I + f(0.81684757298046, 0.09157621350977) * 0.10995174365532 * 1.0/2.0

    return I

def simpson_rule5(f):
    I = 0
    I = I + f(0.33333333333333, 0.33333333333333) * 0.22500000000000 * 1.0/2.0
    I = I + f(0.47014206410511, 0.47014206410511) * 0.13239415278851 * 1.0/2.0
    I = I + f(0.47014206410511, 0.05971587178977) * 0.13239415278851 * 1.0/2.0
    I = I + f(0.05971587178977, 0.47014206410511) * 0.13239415278851 * 1.0/2.0
    I = I + f(0.10128650732346, 0.10128650732346) * 0.12593918054483 * 1.0/2.0
    I = I + f(0.10128650732346, 0.79742698535309) * 0.12593918054483 * 1.0/2.0
    I = I + f(0.79742698535309, 0.10128650732346) * 0.12593918054483 * 1.0/2.0

    return I

def simpson_rule6(f):
    
    I = 0
    
    I = I + f(0.24928674517091, 0.24928674517091) * 0.11678627572638 * 1.0/2.0
    I = I + f(0.24928674517091, 0.50142650965818) * 0.11678627572638 * 1.0/2.0
    I = I + f(0.50142650965818, 0.24928674517091) * 0.11678627572638 * 1.0/2.0
    I = I + f(0.06308901449150, 0.06308901449150) * 0.05084490637021 * 1.0/2.0
    I = I + f(0.06308901449150, 0.87382197101700) * 0.05084490637021 * 1.0/2.0
    I = I + f(0.87382197101700, 0.06308901449150) * 0.05084490637021 * 1.0/2.0
    I = I + f(0.31035245103378, 0.63650249912140) * 0.08285107561837 * 1.0/2.0
    I = I + f(0.63650249912140, 0.05314504984482) * 0.08285107561837 * 1.0/2.0
    I = I + f(0.05314504984482, 0.31035245103378) * 0.08285107561837 * 1.0/2.0
    I = I + f(0.63650249912140, 0.31035245103378) * 0.08285107561837 * 1.0/2.0
    I = I + f(0.31035245103378, 0.05314504984482) * 0.08285107561837 * 1.0/2.0
    I = I + f(0.05314504984482, 0.63650249912140) * 0.08285107561837 * 1.0/2.0
    return I

def simpson_rule7(f):
    I = 0

    I = I + f(0.33333333333333,0.33333333333333) * -0.14957004446768 * 1.0/2.0
    I = I + f(0.26034596607904,0.26034596607904) * 0.17561525743321 * 1.0/2.0
    I = I + f(0.26034596607904,0.47930806784192) * 0.17561525743321 * 1.0/2.0
    I = I + f(0.47930806784192,0.26034596607904) * 0.17561525743321 * 1.0/2.0
    I = I + f(0.06513010290222,0.06513010290222) * 0.05334723560884 * 1.0/2.0
    I = I + f(0.06513010290222,0.86973979419557) * 0.05334723560884 * 1.0/2.0
    I = I + f(0.86973979419557,0.06513010290222) * 0.05334723560884 * 1.0/2.0
    I = I + f(0.31286549600487,0.63844418856981) * 0.07711376089026 * 1.0/2.0
    I = I + f(0.63844418856981,0.04869031542532) * 0.07711376089026 * 1.0/2.0
    I = I + f(0.04869031542532,0.31286549600487) * 0.07711376089026 * 1.0/2.0
    I = I + f(0.63844418856981,0.31286549600487) * 0.07711376089026 * 1.0/2.0
    I = I + f(0.31286549600487,0.04869031542532) * 0.07711376089026 * 1.0/2.0
    I = I + f(0.04869031542532,0.63844418856981) * 0.07711376089026 * 1.0/2.0
    
    return I
def gauss_integration_HN(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem_HN,detJ_HN):
    J = 0
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            J = J + wi[i] * wi[j] * diff_fct_triangles_HN(ui[i],ui[j],x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem_HN,detJ_HN)
    return J

def gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem,detJ):
    J = 0
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            J = J + wi[i] * wi[j] * ufunction(ui[i],ui[j],x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ)
    return J


def gauss_integration(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem,detJ):

    fct = lambda e,n: ufunction(e,n,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ)
    
    if POL_APPROX == 0:
        return 1.0/2.0 * fct(1.0/3.0, 1.0/3.0)
    if POL_APPROX == 1:
        return simpson_rule3(fct)
    if POL_APPROX == 2:
        return simpson_rule3(fct)
    
def gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem,detJ):
    J = 0
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            J = J + wi[i] * wi[j] * ufunction(ui[i],ui[j],x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ)
    return J

def found_in_FEM(point_x,point_y,pConf,tConf,UConf):
    u_conf = 0

    T = len(tConf)
    for e in range(0,T):
        nodes = tConf[e]
        nodes = np.array(nodes)
        coords = pConf[nodes,:]
        triangle_FEM = [ (coords[0,0],coords[0,1]), (coords[1,0],coords[1,1]), (coords[2,0],coords[2,1]) ]

        if point_in_on_poly(point_x,point_y,triangle_FEM):
                
            Pe = np.zeros((3,3))
            Pe[:,0] = np.ones((3,1)).transpose()
            Pe[:,1] = pConf[nodes[0:3],0]
            Pe[:,2] = pConf[nodes[0:3],1]
            C = np.linalg.inv(Pe)
            Nbasis = tribasisFct(C)
            u_conf = (UConf[tConf[e][0],0] * Nbasis[0](point_x, point_y) +
                                UConf[tConf[e][1],0] * Nbasis[1](point_x, point_y) +
                                UConf[tConf[e][2],0] * Nbasis[2](point_x, point_y) )

            return u_conf

    return 0    
    
def diff_fct_triangles_HN(x,y,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem_HN,detJ_HN):
    # convert xi, yi from (-1,1) to physical coordinates (A,B)
    xi_phys = x_trans_fct(x,y)
    yi_phys = y_trans_fct(x,y)

    u_conf = 0

    T = len(tConf)
    for e in range(0,T):
        nodes = tConf[e]
        nodes = np.array(nodes)
        coords = pConf[nodes,:]
        triangle_FEM = [ (coords[0,0],coords[0,1]), (coords[1,0],coords[1,1]), (coords[2,0],coords[2,1]) ]
        point_x = xi_phys
        point_y = yi_phys

        if point_in_on_poly(point_x,point_y,triangle_FEM):
                
            Pe = np.zeros((3,3))
            Pe[:,0] = np.ones((3,1)).transpose()
            Pe[:,1] = pConf[nodes[0:3],0]
            Pe[:,2] = pConf[nodes[0:3],1]
            C = np.linalg.inv(Pe)
            Nbasis = tribasisFct(C)
            u_conf = ( UConf[tConf[e][0],0] * Nbasis[0](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][1],0] * Nbasis[1](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][2],0] * Nbasis[2](x_trans_fct(x,y),y_trans_fct(x,y)) )

            difference_sum = (u_conf - uh_elem_HN(x,y))**2 * detJ_HN(x,y)
            return difference_sum

    return 0

def diff_fct_triangles(x,y,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ):
    # convert xi, yi from (-1,1) to physical coordinates (A,B)
    xi_phys = x_trans_fct(x,y)
    yi_phys = y_trans_fct(x,y)

    u_conf = 0

    T = len(tConf)
    for e in range(0,T):
        nodes = tConf[e]
        nodes = np.array(nodes)
        coords = pConf[nodes,:]
        triangle_FEM = [ (coords[0,0],coords[0,1]), (coords[1,0],coords[1,1]), (coords[2,0],coords[2,1]) ]
        point_x = xi_phys
        point_y = yi_phys

        if point_in_on_poly(point_x,point_y,triangle_FEM):
                
            Pe = np.zeros((3,3))
            Pe[:,0] = np.ones((3,1)).transpose()
            Pe[:,1] = pConf[nodes[0:3],0]
            Pe[:,2] = pConf[nodes[0:3],1]
            C = np.linalg.inv(Pe)
            Nbasis = tribasisFct(C)
            u_conf = ( UConf[tConf[e][0],0] * Nbasis[0](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][1],0] * Nbasis[1](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][2],0] * Nbasis[2](x_trans_fct(x,y),y_trans_fct(x,y)) )

            difference_sum = (u_conf - uh_elem(x_trans_fct(x,y),y_trans_fct(x,y)))**2 * detJ(x_trans_fct(x,y),y_trans_fct(x,y))
            return difference_sum

    return 0
    
def diff_fct(x,y,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ):
    # convert xi, yi from (-1,1) to physical coordinates (A,B)
    xi_phys = x_trans_fct(x,y)
    yi_phys = y_trans_fct(x,y)

    u_conf = 0

    T = len(tConf)
    for e in range(0,T):
        nodes = tConf[e]
        nodes = np.array(nodes)
        coords = pConf[nodes,:]
        x1 = coords[0,0]
        x2 = coords[1,0]
        y1 = coords[0,1]
        y2 = coords[2,1]
        if (x1 <= xi_phys and xi_phys <= x2) and (y1 <= yi_phys and yi_phys <= y2):
            
            Pe = np.zeros((4,4))
            Pe[:,0] = np.ones((4,1)).transpose()
            Pe[:,1] = pConf[nodes[0:4],0]
            Pe[:,2] = pConf[nodes[0:4],1]
            Pe[:,3] = pConf[nodes[0:4],0]*pConf[nodes[0:4],1]
            C = np.linalg.inv(Pe)
            Nbasis = basisFct(C)
            u_conf = ( UConf[tConf[e][0],0] * Nbasis[0](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][1],0] * Nbasis[1](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][2],0] * Nbasis[2](x_trans_fct(x,y),y_trans_fct(x,y)) +
                                UConf[tConf[e][3],0] * Nbasis[3](x_trans_fct(x,y),y_trans_fct(x,y)) )

            difference_sum = (u_conf - uh_elem(x_trans_fct(x,y),y_trans_fct(x,y)))**2 * detJ(x_trans_fct(x,y),y_trans_fct(x,y))
            return difference_sum

def ufunction(x,y,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ):
    return diff_fct_triangles(x,y,x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem,detJ)

def Nbasis68(x,y,p,nodes,x0,x1,y0,y1):
    
    PeB = np.zeros((4,4))
    PeB[:,0] = np.ones((4,1)).transpose()
    PeB[:,1] = p[nodes[0:4],0]
    PeB[:,2] = p[nodes[0:4],1]
    PeB[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    PeB[2,2] = (y0+y1)/2.0
    PeB[3,2] = (y0+y1)/2.0
    PeB[2,3] = (y0+y1)/2.0 * x1
    PeB[3,3] = (y0+y1)/2.0 * x0
    
    CB = np.linalg.inv(PeB)
    NbasisB = basisFct(CB)

    PeT = np.zeros((4,4))
    PeT[:,0] = np.ones((4,1)).transpose()
    PeT[:,1] = p[nodes[0:4],0]
    PeT[:,2] = p[nodes[0:4],1]
    PeT[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    PeT[0,2] = (y0+y1)/2.0
    PeT[1,2] = (y0+y1)/2.0
    PeT[0,3] = (y0+y1)/2.0 * x0
    PeT[1,3] = (y0+y1)/2.0 * x1
            
    CT = np.linalg.inv(PeT)
    NbasisT = basisFct(CT)
            
    if y <= (y0+y1)/2.0:
        Nbasis6 = NbasisB[2](x,y)
        Nbasis8 = NbasisB[3](x,y)
    else:
        Nbasis6 = NbasisT[1](x,y)
        Nbasis8 = NbasisT[0](x,y)
        
    return [Nbasis6, Nbasis8]
    
def Nbasis57(x,y,p,nodes,x0,x1,y0,y1):
    
    PeL = np.zeros((4,4))
    PeL[:,0] = np.ones((4,1)).transpose()
    PeL[:,1] = p[nodes[0:4],0]
    PeL[:,2] = p[nodes[0:4],1]
    PeL[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    PeL[1,1] = (x0+x1)/2.0
    PeL[2,1] = (x0+x1)/2.0
    PeL[1,3] = (x0+x1)/2.0 * y0
    PeL[2,3] = (x0+x1)/2.0 * y1
    
    CL = np.linalg.inv(PeL)
    NbasisL = basisFct(CL)

    PeR = np.zeros((4,4))
    PeR[:,0] = np.ones((4,1)).transpose()
    PeR[:,1] = p[nodes[0:4],0]
    PeR[:,2] = p[nodes[0:4],1]
    PeR[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    PeR[0,1] = (x0+x1)/2.0
    PeR[3,1] = (x0+x1)/2.0
    PeR[0,3] = (x0+x1)/2.0 * y0
    PeR[3,3] = (x0+x1)/2.0 * y1
            
    CR = np.linalg.inv(PeR)
    NbasisR = basisFct(CR)
            
    if x <= (x0+x1)/2.0:
        Nbasis5 = NbasisL[1](x,y)
        Nbasis7 = NbasisL[2](x,y)
    else:
        Nbasis5 = NbasisR[0](x,y)
        Nbasis7 = NbasisR[3](x,y)
        
    return [Nbasis5,Nbasis7]

def Nbases(Nbasis,x0,x1,y0,y1,p,nodes,N,S,E,W):

    if S == 1:
        N5 = lambda x,y: Nbasis57(x,y,p,nodes,x0,x1,y0,y1)[0]
    else:
        N5 = lambda x,y: 0.0
    
    if E == 1:
        N6 = lambda x,y: Nbasis68(x,y,p,nodes,x0,x1,y0,y1)[0]
    else:
        N6 = lambda x,y: 0.0
    if N == 1:
        N7 = lambda x,y: Nbasis57(x,y,p,nodes,x0,x1,y0,y1)[1]
    else:
        N7 = lambda x,y: 0.0

    if W == 1:
        N8 = lambda x,y: Nbasis68(x,y,p,nodes,x0,x1,y0,y1)[1]
    else:
        N8 = lambda x,y: 0.0


    N1 = lambda x,y: Nbasis[0](x,y) - 1.0/2.0 * (N5(x,y) + N8(x,y))
    N2 = lambda x,y: Nbasis[1](x,y) - 1.0/2.0 * (N5(x,y) + N6(x,y))
    N3 = lambda x,y: Nbasis[2](x,y) - 1.0/2.0 * (N6(x,y) + N7(x,y))
    N4 = lambda x,y: Nbasis[3](x,y) - 1.0/2.0 * (N8(x,y) + N7(x,y))
    
    return [N1,N2,N3,N4]

def myquad(k1,k2,ui,wi,p,t,masterNode,llist,image,lenGridPts):
    
    tbcs = []
    bbcs = []
    lbcs = []
    rbcs = []
    for i in range(0,lenGridPts):
         # bottom boundary
        if p[i,1] == 0.0 and (p[i,0] != 0.0 and p[i,0] != 1.0):
            bbcs = bbcs + [i]
        else:
            # top boundary            
            if p[i,1] == 1.0 and (p[i,0] != 0.0 and p[i,0] != 1.0):
                tbcs = tbcs + [i]
            else:
                # left boundary
                if p[i,0] == 0.0:
                   lbcs = lbcs + [i]
                else:
                    # right boundary
                    if p[i,0] == 1.0:
                        rbcs = rbcs + [i]

    # temperatures for Dirichlet BCs
    leftDirich = 1
    rightDirich = 1
    topDirich = 0
    bottomDirich = 0
    Temp_left = 0
    Temp_right = 0
    Temp_top = 0
    Temp_bottom = 0

    # Neumann BCs
    leftNeumann = 0
    rightNeumann = 0
    topNeumann = 1
    bottomNeumann = 1
    g1_left = 0
    g1_right = 0
    g1_top = 0
    g1_bottom = 0

    N = len(p)  #+ 5 # number of nodes
    T = len(t) #+ 3 # number of elements
    
    nb = 4 # number of basis functions per element|
    K_col = numpy.zeros((T * nb**2, 1))
    K_row = numpy.zeros((T * nb**2, 1))
    K_val = numpy.zeros((T * nb**2, 1))
    
    F_row = numpy.zeros((T * nb**2, 1))
    F_val = numpy.zeros((T * nb**2, 1))
    
    print 'number of elements:', T
    sstart = time.time()
    list_hanging_nodes = []
    
    
    for e in range(0,T): #800, 833
        
        nodes = t[e] 
        
        root = get_node_by_id(masterNode,llist[e])
                     
        p1,p2,p3,p4 = root.rect
        
        nodes6 = Coordinate(0,0)
        nodes7 = Coordinate(0,0)
        
        if len(root.enrichNodes) == 3:
            
            nodes6.x =  root.enrichNodes[2].x / 1000.0
            if nodes6.x != 0.0:
                nodes6.x += 0.001
            nodes6.y =  root.enrichNodes[2].y / 1000.0
            if nodes6.y != 0.0:
                nodes6.y += 0.001
            nodes6.y = 1 - nodes6.y
    
            old_node6 = Coordinate(nodes6.x, nodes6.y)
            [c_center, c_rad] = which_circle(nodes6)
            new_nodes6 = Coordinate(0.0,0.0)
            new_nodes6.x = c_center.x + c_rad * ( nodes6.x - c_center.x) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            new_nodes6.y = c_center.y + c_rad * ( nodes6.y - c_center.y) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            nodes6 = Coordinate(new_nodes6.x, new_nodes6.y)
            
        if len(root.enrichNodes) == 4:
            
            if root.enrichNodes[2].y == root.enrichNodes[3].y:
                if root.enrichNodes[2].x <= root.enrichNodes[3].x: 
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[2].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[2].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[3].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[3].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                    
                else:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[3].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[3].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[2].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[2].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
            else:
                if root.enrichNodes[2].y >= root.enrichNodes[3].y:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[2].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[2].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[3].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[3].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                else:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[3].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[3].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[2].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[2].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                    
            old_node6 = Coordinate(nodes6.x, nodes6.y)
            [c_center, c_rad] = which_circle(nodes6)
            new_nodes6 = Coordinate(0.0,0.0)
            new_nodes6.x = c_center.x + c_rad * ( nodes6.x - c_center.x) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            new_nodes6.y = c_center.y + c_rad * ( nodes6.y - c_center.y) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            nodes6 = Coordinate(new_nodes6.x, new_nodes6.y)
            
            old_node7 = Coordinate(nodes7.x, nodes7.y)
            [c_center, c_rad] = which_circle(nodes7)
            new_nodes7 = Coordinate(0.0,0.0)
            new_nodes7.x = c_center.x + c_rad * ( nodes7.x - c_center.x) / math.sqrt( math.pow(nodes7.x - c_center.x, 2) + math.pow( nodes7.y - c_center.y,2) )
            new_nodes7.y = c_center.y + c_rad * ( nodes7.y - c_center.y) / math.sqrt( math.pow(nodes7.x - c_center.x, 2) + math.pow( nodes7.y - c_center.y,2) )
            nodes7 = Coordinate(new_nodes7.x, new_nodes7.y)
            
                        
        pxVal1 = image.GetPixel(int(p1.x), int(p1.y))
        pxVal2 = image.GetPixel(int(p2.x), int(p2.y))
        pxVal3 = image.GetPixel(int(p3.x), int(p3.y))
        pxVal4 = image.GetPixel(int(p4.x), int(p4.y))
        
        pxVals = [pxVal1,pxVal2,pxVal3,pxVal4]
        # 2-column matrix containing on each row the coordinates of each of the nodes
        coords = p[nodes,:]    

        coords[0,:] = p[nodes[0]]
        coords[1,:] = p[nodes[1]]
        coords[2,:] = p[nodes[2]]
        coords[3,:] = p[nodes[3]]
        
        Pe = np.zeros((4,4))
        Pe[:,0] = np.ones((4,1)).transpose()
        Pe[:,1] = p[nodes[0:4],0]
        Pe[:,2] = p[nodes[0:4],1]
        Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

        C = np.linalg.inv(Pe)

        Nbasis = basisFct(C)
        Nx = derivX(C)
        Ny = derivY(C)

        x0 = coords[0,0]
        x1 = coords[1,0]
        y0 = coords[0,1]
        y1 = coords[2,1]

        thru_corner = False

        if len(nodes) == 5:
                enrich1 = np.array(p[nodes[4]])

                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )

                which_corner = -1
                if corner0 == True and corner1 == False and corner2 == False and corner3 == False:
                    thru_corner = True
                    which_corner = 1
                if corner0 == False and corner1 == True and corner2 == False and corner3 == False:
                    thru_corner = True
                    which_corner = 2
                if corner0 == False and corner1 == False and corner2 == True and corner3 == False:
                    thru_corner = True
                    which_corner = 3
                if corner0 == False and corner1 == False and corner2 == False and corner3 == True:
                    thru_corner = True
                    which_corner = 4   
                     
        if (len(nodes) == 4  or root.ishomog == 1) or thru_corner == True:        # elements need no enrichment
                      
                  
            Ke = np.zeros((4,4))
            Fe = np.zeros((4,1))
         
            p1,p2,p3,p4 = root.rect
            
            pxVal1 = image.GetPixel(int(p1.x), int(p1.y))
            pxVal2 = image.GetPixel(int(p2.x), int(p2.y))
            pxVal3 = image.GetPixel(int(p3.x), int(p3.y))
            pxVal4 = image.GetPixel(int(p4.x), int(p4.y))
            pxValMed = image.GetPixel(int((p1.x+p3.x)/2.0), int((p1.y+p3.y)/2.0))
            
            if ( is_in_same_bin(pxVal1,pxVal2) == True and is_in_same_bin(pxVal3,pxVal4)==True and 
                 is_in_same_bin(pxVal3,pxVal2)==True and  pxVal1 <= binBnd[1] ):
                K_cst = k1
            else: 
                if ( is_in_same_bin(pxVal1,pxVal2) == True and is_in_same_bin(pxVal3,pxVal4)==True and 
                 is_in_same_bin(pxVal3,pxVal2)==True and  pxVal1 > binBnd[1] ):
                    K_cst = k2
                elif root.ishomog == 1:
                    if pxValMed > binBnd[1]:
                        K_cst = k2
                    else:
                        K_cst = k1
        
            if thru_corner == True:
                if which_corner == 1: 
                    if pxVal1 <= binBnd[1]:
                        K_cst = k1
                    else:
                        K_cst = k2
                if which_corner == 2: 
                    if pxVal2 <= binBnd[1]:
                        K_cst = k1
                    else:
                        K_cst = k2
                if which_corner == 3: 
                    if pxVal3 <= binBnd[1]:
                        K_cst = k1
                    else:
                        K_cst = k2
                if which_corner == 4: 
                    if pxVal4 <= binBnd[1]:
                        K_cst = k1
                    else:
                        K_cst = k2                                                             

            #N,S,E,W
            [N_hn,S_hn,E_hn,W_hn] = root.nsew
            
            # Method 2 for computing HangingNodes
            #===================================================================
            # [x_fct,y_fct] = xy_fct_HN(coords[0:4,0],coords[0:4,1],N_hn,S_hn,E_hn,W_hn)
            # Jac = jacobian_mat_HN( coords[0:4,0], coords[0:4,1],N_hn,S_hn,E_hn,W_hn)
            # det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
            #===================================================================
            
            # Method 1 for computing Hanging Nodes: averaging
            if N_hn == 1: # North
                north_hn = root.hn[0]
                nhn = [north_hn.x, north_hn.y]
                n_ind = numpy.where(numpy.all(p==nhn,axis=1))
                north_node = n_ind[0][0]           
                list_hanging_nodes = list_hanging_nodes + [[ north_node, nodes[2], nodes[3] ]]

            if S_hn == 1:
                south_hn = root.hn[1]
                shn = [south_hn.x, south_hn.y]
                s_ind = numpy.where(numpy.all(p==shn,axis=1))
                south_node = s_ind[0][0]
                list_hanging_nodes = list_hanging_nodes + [[ south_node, nodes[0], nodes[1] ]]

                    
            if E_hn == 1:
                east_hn = root.hn[2]
                ehn = [east_hn.x, east_hn.y]
                e_ind = numpy.where(numpy.all(p==ehn,axis=1))
                east_node = e_ind[0][0]   
                list_hanging_nodes = list_hanging_nodes + [[ east_node, nodes[1], nodes[2] ]]

                                                 
            if W_hn == 1:
                west_hn = root.hn[3]
                whn = [west_hn.x, west_hn.y]
                w_ind = numpy.where(numpy.all(p==whn,axis=1))
                west_node = w_ind[0][0]            
                list_hanging_nodes = list_hanging_nodes + [[ west_node, nodes[0], nodes[3] ]]


            [x_fct,y_fct] = xy_fct(coords[0:4,0],coords[0:4,1])
            Jac = jacobian_mat( coords[0:4,0], coords[0:4,1])
            det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
                        
            # Method 1
            # construct the local matrix and local components of the load vector    
            for i in range(0,4):
                for j in range(0,4):
                    if i>= j:
                            Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                            Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)
                            Ke[j,i] = Ke[i,j]
                # construct the local load vector
                fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)
            
            # add the local stiffness matrix and local load vector to the global K and F
            for i in range(0,4):
                for j in range(0,4):
                        K_row = numpy.vstack([K_row,nodes[i]])
                        K_col = numpy.vstack([K_col,nodes[j]])
                        K_val = numpy.vstack([K_val,Ke[i,j]])
           
                F_row = numpy.vstack([F_row,nodes[i]])
                F_val = numpy.vstack([F_val,Fe[i]])
            
        else: # element has more than 4 nodes, it is an element that needs enrichment at these additional nodes
            enrich1 = np.array(p[nodes[4]])

            if len(nodes) == 5:
                
                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )

                elCorner = False
    
             # FALSE POSITIVIES: only corner is in a different material
             # HOMOGENEOUS element
                if (corner1 == True or corner3 == True):
                    elCorner = True
                    midInside = Point( (x0+x1)/2.0, (y0+y1)/2.0 )
                    
                    pxValMid = image.GetPixel(int(midInside.x), int(midInside.y))
                    if pxValMid > binBnd[1]:
                        K_cst = k2
                    else:
                        K_cst = k1
                    
                    Ke = np.zeros((4,4))
                    Fe = np.zeros((4,1))
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,4):
                        for j in range(0,4):
                            if i>=j:
                                Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                                Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)
                                Ke[j,i] = Ke[i,j]
                        # construct the local load vector
                        fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                        Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)
            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,4):
                        for j in range(0,4):
                                K_row = numpy.vstack([K_row,nodes[i]])
                                K_col = numpy.vstack([K_col,nodes[j]])
                                K_val = numpy.vstack([K_val,Ke[i,j]])
                        F_row = numpy.vstack([F_row,nodes[i]])
                        F_val = numpy.vstack([F_val,Fe[i]])            

             # FALSE POSITIVIES: only corner is in a different material
             # HOMOGENEOUS element
                if (corner0 == True or corner2 == True):
                    elCorner = True
                    midInside = Point( (x0+x1)/2.0, (y0+y1)/2.0 )
                                        
                    pxValMid = image.GetPixel(int(midInside.x), int(midInside.y))
                    if pxValMid > binBnd[1]:
                        K_cst = k2
                    else:
                        K_cst = k1
                    
                    Ke = np.zeros((4,4))
                    Fe = np.zeros((4,1))
            
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,4):
                        for j in range(0,4):
                            if i>= j:
                                Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                                Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)
                                Ke[j,i] = Ke[i,j]
                        # construct the local load vector
                        fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                        Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)
            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,4):
                        for j in range(0,4):
                                K_row = numpy.vstack([K_row,nodes[i]])
                                K_col = numpy.vstack([K_col,nodes[j]])
                                K_val = numpy.vstack([K_val,Ke[i,j]])
                        F_row = numpy.vstack([F_row,nodes[i]])
                        F_val = numpy.vstack([F_val,Fe[i]])                        


            else: # or (len(nodes) == 6)
                
                # enrichment nodes: enrichmentNode = [x y], x = enrichmentNode[0], y = enrichmentNode[1]
                enrich1 = np.array(p[nodes[4]])
                enrich2 = np.array(p[nodes[5]])

                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich2[0] - [x1]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich2[0] - [x0]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )
                
                # testing if interface is along a diagonal
                if (corner0 == True and corner2 == True) or (corner1 == True and corner3 == True):
                    Ke_trid1 = np.zeros((3,3))
                    Fe_trid1 = np.zeros((3,1))
                    Ke_trid2 = np.zeros((3,3))
                    Fe_trid2 = np.zeros((3,1))

                    # even diagonal: SW - NE
                    if(corner0 == True and corner2 == True):
                        nodes_trid1 = [nodes[0], nodes[1], nodes[2]]
                        nodes_trid2 = [nodes[0], nodes[2], nodes[3]]
            
                        if pxVal1 > binBnd[1] and pxVal3 <= binBnd[1]:
                            K_cst_trid2 = k2 
                            K_cst_trid1 = k1
                        else:
                            K_cst_trid2 = k1
                            K_cst_trid1 = k2
                        
                    # odd diagonal: SE - NW
                    else:
                        nodes_trid1 = [nodes[0], nodes[1], nodes[3]]
                        nodes_trid2 = [nodes[1], nodes[2], nodes[3]]

                        if pxVal2 > binBnd[1] and pxVal4 <= binBnd[1]:
                            K_cst_trid2 = k2 
                            K_cst_trid1 = k1
                        else:
                            K_cst_trid2 = k1
                            K_cst_trid1 = k2

                    coords_trid1 = p[nodes_trid1]
                    coords_trid2 = p[nodes_trid2]

                    Pe_trid1 = np.zeros((3,3))
                    Pe_trid2 = np.zeros((3,3))
    
                    Pe_trid1[:,0] = np.ones((3,1)).transpose()                    
                    Pe_trid1[:,1] = p[nodes_trid1[0:3],0]
                    Pe_trid1[:,2] = p[nodes_trid1[0:3],1]

                    Pe_trid2[:,0] = np.ones((3,1)).transpose()                    
                    Pe_trid2[:,1] = p[nodes_trid2[0:3],0]
                    Pe_trid2[:,2] = p[nodes_trid2[0:3],1]

                    C_trid1 = np.linalg.inv(Pe_trid1)
                    C_trid2 = np.linalg.inv(Pe_trid2)

                    Nbasis_trid1 = tribasisFct(C_trid1)
                    Nbasis_trid2 = tribasisFct(C_trid2)

                    Nx_trid1 = triderivX(C_trid1)
                    Ny_trid1 = triderivY(C_trid1)

                    Nx_trid2 = triderivX(C_trid2)
                    Ny_trid2 = triderivY(C_trid2)

                    if len(root.enrichNodes) < 3:
                    
                        [x_fct_1, y_fct_1] = tri_xy_fct( coords_trid1[:,0], coords_trid1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( coords_trid2[:,0], coords_trid2[:,1] )
        
                        J1 = tri_jacobian_mat( coords_trid1[:,0], coords_trid1[:,1] )
                        J2 = tri_jacobian_mat( coords_trid2[:,0], coords_trid2[:,1] )
                    
                    if len(root.enrichNodes) == 3:
                           
                        coord_enrich = coord_enrich_comp_quad_circle(p[nodess], nodes6)
                          
                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [1,2,0] 
                        else:
                            lOrd = [0,1,2] 
                          
                        vec1_x = [ coords_trid1[ lOrd[0],0], coords_trid1[lOrd[1],0], coords_trid1[lOrd[2],0], (coords_trid1[ lOrd[0],0] + coords_trid1[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid1[lOrd[0],0] + coords_trid1[lOrd[2],0])/2.0  ]
                        vec1_y = [ coords_trid1[ lOrd[0],1], coords_trid1[lOrd[1],1], coords_trid1[lOrd[2],1], (coords_trid1[ lOrd[0],1] + coords_trid1[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid1[lOrd[0],1] + coords_trid1[lOrd[2],1])/2.0  ]
                
                        [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                        J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [2,0,1] 
                        else:
                            lOrd = [1,2,0] 
                            
                        vec2_x = [ coords_trid2[ lOrd[0],0], coords_trid2[lOrd[1],0], coords_trid2[lOrd[2],0], (coords_trid2[ lOrd[0],0] + coords_trid2[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid2[lOrd[0],0] + coords_trid2[lOrd[2],0])/2.0  ]
                        vec2_y = [ coords_trid2[ lOrd[0],1], coords_trid2[lOrd[1],1], coords_trid2[lOrd[2],1], (coords_trid2[ lOrd[0],1] + coords_trid2[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid2[lOrd[0],1] + coords_trid2[lOrd[2],1])/2.0  ]
                
                        [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                        J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )

                    if len(root.enrichNodes) == 4:
                
                        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess], nodes6)
                        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess], nodes7)
                        
                        circumcenter_pt1 = circumcenter_tri(coords_trid1)
                        circumcenter_pt2 = circumcenter_tri(coords_trid2)

                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [1,2,0] 
                            pt6.x = coord_enrich2.x
                            pt6.y = coord_enrich2.y 
                        else:
                            lOrd = [0,1,2]   
                            pt7.x = coord_enrich1.x
                            pt7.y = coord_enrich1.y                      
                        
                        pt4 = Coordinate(0,0)
                        pt5 = Coordinate(0,0)
                        pt6 = Coordinate(0,0)
                        pt7 = Coordinate(0,0)
                        pt8 = Coordinate(0,0)
                        pt9 = Coordinate(0,0)
                         
                        pt4.x = 2.0/3.0 * coords_trid1[lOrd[0],0] + 1.0/3.0 * coords_trid1[lOrd[1],0]
                        pt4.y = 2.0/3.0 * coords_trid1[lOrd[0],1] + 1.0/3.0 * coords_trid1[lOrd[1],1]
                         
                        pt5.x = 1.0/3.0 * coords_trid1[lOrd[0],0] + 2.0/3.0 * coords_trid1[lOrd[1],0]
                        pt5.y = 1.0/3.0 * coords_trid1[lOrd[0],1] + 2.0/3.0 * coords_trid1[lOrd[1],1]
                         
                        pt8.x = 2.0/3.0 * coords_trid1[lOrd[2],0] + 1.0/3.0 * coords_trid1[lOrd[0],0]
                        pt8.y = 2.0/3.0 * coords_trid1[lOrd[2],1] + 1.0/3.0 * coords_trid1[lOrd[0],1]
                         
                        pt9.x = 1.0/3.0 * coords_trid1[lOrd[2],0] + 2.0/3.0 * coords_trid1[lOrd[0],0]
                        pt9.y = 1.0/3.0 * coords_trid1[lOrd[2],1] + 2.0/3.0 * coords_trid1[lOrd[0],1]
                 
                         
                        vec1_x = [coords_trid1[lOrd[0],0], 
                                  coords_trid1[lOrd[1],0], 
                                  coords_trid1[lOrd[2],0], 
                                  pt4.x, 
                                  pt5.x,
                                  pt6.x,
                                  pt7.x,
                                  pt8.x,
                                  pt9.x,
                                  circumcenter_pt1.x  
                                  ]
                        vec1_y = [coords_trid1[lOrd[0],1], 
                                  coords_trid1[lOrd[1],1], 
                                  coords_trid1[lOrd[2],1], 
                                  pt4.y, 
                                  pt5.y,
                                  pt6.y,
                                  pt7.y,
                                  pt8.y,
                                  pt9.y,
                                  circumcenter_pt1.y  
                                  ]
                        [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                        J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                        

                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [2,0,1] 
                            pt6.x = coord_enrich1.x
                            pt6.y = coord_enrich1.y  
                        else:
                            lOrd = [1,2,0] 
                            pt7.x = coord_enrich2.x
                            pt7.y = coord_enrich2.y
                            
                        pt4 = Coordinate(0,0)
                        pt5 = Coordinate(0,0)
                        pt6 = Coordinate(0,0)
                        pt7 = Coordinate(0,0)
                        pt8 = Coordinate(0,0)
                        pt9 = Coordinate(0,0)
                    
                        pt4.x = 2.0/3.0 * coords_trid2[lOrd[0],0] + 1.0/3.0 * coords_trid2[lOrd[1],0]
                        pt4.y = 2.0/3.0 * coords_trid2[lOrd[0],1] + 1.0/3.0 * coords_trid2[lOrd[1],1]
                        
                        pt5.x = 1.0/3.0 * coords_trid2[lOrd[0],0] + 2.0/3.0 * coords_trid2[lOrd[1],0]
                        pt5.y = 1.0/3.0 * coords_trid2[lOrd[0],1] + 2.0/3.0 * coords_trid2[lOrd[1],1]
                        
                        pt8.x = 2.0/3.0 * coords_trid2[lOrd[2],0] + 1.0/3.0 * coords_trid2[lOrd[0],0]
                        pt8.y = 2.0/3.0 * coords_trid2[lOrd[2],1] + 1.0/3.0 * coords_trid2[lOrd[0],1]
                        
                        pt9.x = 1.0/3.0 * coords_trid2[lOrd[2],0] + 2.0/3.0 * coords_trid2[lOrd[0],0]
                        pt9.y = 1.0/3.0 * coords_trid2[lOrd[2],1] + 2.0/3.0 * coords_trid2[lOrd[0],1]
                    
                        vec2_x = [coords_trid2[lOrd[0],0], 
                                  coords_trid2[lOrd[1],0],
                                  coords_trid2[lOrd[2],0],
                                  pt4.x, 
                                  pt5.x,
                                  pt6.x,
                                  pt7.x,
                                  pt8.x,
                                  pt9.x,
                                  circumcenter_pt2.x  
                                  ]
                        vec2_y = [coords_trid2[lOrd[0],1], 
                                  coords_trid2[lOrd[1],1], 
                                  coords_trid2[lOrd[2],1], 
                                  pt4.y, 
                                  pt5.y,
                                  pt6.y, 
                                  pt7.y,
                                  pt8.y,
                                  pt9.y,
                                  circumcenter_pt2.y 
                                  ]
                        
                        [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                        J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )  
                                    
                        
                    det_J1 = lambda e,n: determinant(J1)(e,n)
                    det_J2 = lambda e,n: determinant(J2)(e,n)
                        
                    ## FIRST TRIANGLE
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,3):
                        for j in range(0,3):
                                Kefunc_trid1 = lambda e,n: K_cst_trid1 * ( Nx_trid1[i](x_fct_1(e,n), y_fct_1(e,n)) * Nx_trid1[j](x_fct_1(e,n), y_fct_1(e,n)) + Ny_trid1[i](x_fct_1(e,n), y_fct_1(e,n)) * Ny_trid1[j](x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)
                                Ke_trid1[i,j] = 1.0/2.0 * quad2d(Kefunc_trid1,x0,x1,y0,y1,ui,wi)
                        # construct the local load vector
                        fv_trid1 = lambda x,y: rhs(x,y) * Nbasis_trid1[i](x,y)
                        Fe_trid1[i] = 1.0/2.0 * quad2d(fv_trid1,x0,x1,y0,y1,ui,wi)                      

                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,3):
                        for j in range(0,3):
                                K_row = numpy.vstack([K_row,nodes_trid1[i]])
                                K_col = numpy.vstack([K_col,nodes_trid1[j]])
                                K_val = numpy.vstack([K_val,Ke_trid1[i,j]])
                                
                        F_row = numpy.vstack([F_row,nodes_trid1[i]])
                        F_val = numpy.vstack([F_val,Fe_trid1[i]])
                        
                    ## SECOND TRIANGLE
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,3):
                        for j in range(0,3):
                                Kefunc_trid2 = lambda e,n: K_cst_trid2 * ( Nx_trid2[i](x_fct_2(e,n), y_fct_2(e,n)) * Nx_trid2[j](x_fct_2(e,n), y_fct_2(e,n)) + Ny_trid2[i](x_fct_2(e,n), y_fct_2(e,n)) * Ny_trid2[j](x_fct_2(e,n), y_fct_2(e,n)) )* det_J2(e,n)
                                Ke_trid2[i,j] = 1.0/2.0 * quad2d(Kefunc_trid2,x0,x1,y0,y1,ui,wi)                               
                                
                        # construct the local load vector
                        fv_trid2 = lambda x,y: rhs(x,y) * Nbasis_trid2[i](x,y)
                        Fe_trid2[i] = 1.0/2.0 * quad2d(fv_trid2,x0,x1,y0,y1,ui,wi)
            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,3):
                        for j in range(0,3):
                                K_row = numpy.vstack([K_row,nodes_trid2[i]])
                                K_col = numpy.vstack([K_col,nodes_trid2[j]])
                                K_val = numpy.vstack([K_val,Ke_trid2[i,j]])
                        F_row = numpy.vstack([F_row,nodes_trid2[i]])
                        F_val = numpy.vstack([F_val,Fe_trid2[i]])                        
                        
                else:

                    # the North-West corner is cut, 0-4-3, 2-5-3

                    if (
                        ((enrich1[0] == x0 and enrich2[1] == y1) or
                         (enrich2[0] == x0 and enrich1[1] == y1)) and
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))                         
                        ):

                        [Ke_NW,Fe_NW] = NW_corner(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)

                        
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_NW[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_NW[i]])

                    # the South-East corner is cut, 0-4-1, 1-5-2
                    if ( 
                        ((enrich1[1] == y0 and enrich2[0] == x1) or
                          (enrich2[1] == y0 and enrich1[0] == x1)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))                        
                        ):
                        
                        [Ke_SE,Fe_SE] = SE_corner(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)
    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_SE[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_SE[i]])
                                
                    # the North East corner is cut, 1-4-2, 2-5-3
                    if ( 
                        ((enrich1[0] == x1 and enrich2[1] == y1) or
                          (enrich2[0] == x1 and enrich1[1] == y1)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
                        ):
                        [Ke_NE,Fe_NE] = NE_corner(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)
    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_NE[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_NE[i]])
                            
                    # the South-West corner is cut, 0-4-1, and 0-5-3
                    if ( 
                        ((enrich1[1] == y0 and enrich2[0] == x0) or
                         (enrich2[1] == y0 and enrich1[0] == x0)) and
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) 
                        ):
                        [Ke_SW,Fe_SW] = SW_corner(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)
        
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_SW[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_SW[i]])
                                                       
                    pxVal14 = image.GetPixel(int( (p1.x+p4.x)/2.0),int( (p1.y+p4.y)/2.0) )
                    pxVal12 = image.GetPixel(int( (p1.x+p2.x)/2.0),int( (p1.y+p2.y)/2.0) )
                    pxVal23 = image.GetPixel(int( (p2.x+p3.x)/2.0),int( (p2.y+p3.y)/2.0) )
                    pxVal34 = image.GetPixel(int( (p3.x+p4.x)/2.0),int( (p3.y+p4.y)/2.0) ) 
                    
                    pxVals2 = [pxVal14,pxVal12,pxVal23,pxVal34]
                    # the South edge
                    if (  ((enrich1[1] == y0 and enrich2[1] == y1) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x1) ) or
                        ( (enrich2[1] == y0 and enrich1[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1 ) ) ):

                        if not(on_corners(enrich2,coords)) :
                            south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_South,Fe_South] = South_edge(p,ui,wi,k1,k2,south_nodes,root,image,nodes6,nodes7,nodes,pxVals,pxVals2)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
                                    K_row = numpy.vstack([K_row,south_nodes[i]])
                                    K_col = numpy.vstack([K_col,south_nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_South[i,j]])
                            F_row = numpy.vstack([F_row,south_nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_South[i]])
                            
                    # the North edge
                    if ( ( (enrich1[1] == y0 and enrich2[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1)) or
                        ( (enrich1[1] == y1 and enrich2[1] == y0) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x0) )  ):

                        if not(on_corners(enrich2,coords)):
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5] ]
                        else:
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] ]
                        [Ke_North,Fe_North] = North_edge(p,ui,wi,k1,k2,north_nodes,root,image,nodes6,nodes7,nodes,pxVals,pxVals2)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
                                    K_row = numpy.vstack([K_row,north_nodes[i]])
                                    K_col = numpy.vstack([K_col,north_nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_North[i,j]])
                            F_row = numpy.vstack([F_row,north_nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_North[i]])
                    
                    # the West edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1)) or
                        ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)  ) ):

                        enrich1 = np.array(p[nodes[4]])
                        enrich2 = np.array(p[nodes[5]])

                        if not(on_corners(enrich2,coords)) :
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_West,Fe_West] = West_edge(p,ui,wi,k1,k2,west_nodes,root,image,nodes6,nodes7,nodes,pxVals,pxVals2)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
                                    K_row = numpy.vstack([K_row,west_nodes[i]])
                                    K_col = numpy.vstack([K_col,west_nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_West[i,j]])
                            F_row = numpy.vstack([F_row,west_nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_West[i]])


                    # the East edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)) or
                            ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1) )  ):

                        if not(on_corners(enrich2,coords)) :

                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_East,Fe_East] = East_edge(p,ui,wi,k1,k2,east_nodes,root,image,nodes6,nodes7,nodes,pxVals,pxVals2)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
                                    K_row = numpy.vstack([K_row,east_nodes[i]])
                                    K_col = numpy.vstack([K_col,east_nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_East[i,j]])
                            F_row = numpy.vstack([F_row,east_nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_East[i]])    

                    # interface cuts the element horizontally into two quads, 0-4-3, 1-5-2 
                    if ( ((enrich1[0] == x0  and enrich2[0] == x1) or (enrich1[0] == x1 and enrich2[0] == x0)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):

                        [Ke_Horiz,Fe_Horiz] = horizontal_cut(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)
                    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_Horiz[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_Horiz[i]])            

                    # interface cuts the element vertically into two quads, 0-4-1, 3-5-2
                    if ( ((enrich1[1] == y0 and enrich2[1] == y1) or (enrich1[1] == y1 and enrich2[1] == y0 )) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):
                        [Ke_Vertical,Fe_Vertical] = vertical_cut(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals)
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
                                    K_row = numpy.vstack([K_row,nodes[i]])
                                    K_col = numpy.vstack([K_col,nodes[j]])
                                    K_val = numpy.vstack([K_val,Ke_Vertical[i,j]])
                            F_row = numpy.vstack([F_row,nodes[i]])
                            F_val = numpy.vstack([F_val,Fe_Vertical[i]])            
    # end of loop
    # BCs: a * U + b * dU/dx + c * dU/dy + d = 0
    # Dirichlet: b,c = 0, homogeneous Dirichlet: b,c = 0, d = 0
    # Neumann: a = 0, and b or c may be 0, but not both

    K = sparse.coo_matrix((numpy.array(K_val.T)[0],(numpy.array(K_row.T)[0],numpy.array(K_col.T)[0])),shape=(len(p),len(p)))
    K = K.tocsr()
    
    F_col = numpy.zeros((1,len(F_row)))[0]
    F_val = numpy.array(F_val.T)[0]
    F_row = numpy.array(F_row.T)[0]

    F = sparse.coo_matrix((F_val,(F_row,F_col)),shape=(len(p),1))
    F = F.tocsr()
    U = numpy.zeros((N,1))

    # Setting Dirichlet BCs
    # left side of the domain
    if leftDirich == 1:
        for l in lbcs:
            U[l,0] = Temp_left

    # right side of the domain
    if rightDirich == 1:
        for l in rbcs:
            U[l,0] = Temp_right

    #print "bcs", lbcs,rbcs
    #print "top/bottom",range(1,m-1),range(lbcs[-1]+1,rbcs[-1])

    # top side of the domain
    if topDirich == 1:
        for l in tbcs:
            U[l,0] = Temp_top

    # bottom side of the domain
    if bottomDirich == 1:
        for l in bbcs:
            U[l,0] = Temp_bottom

    
    F = F - K*U

    # in case the nodes have duplicate (x,y) coordinates 
    # such as in the case of false positive need for enrichment
    # those nodes need to be removed, as they are not contributing to the stiffness matrix
    # zNodes contains the vector p, and the third column is the row number
    zNodes = numpy.zeros((len(p),3))
    zNodes[:,0] = p[:,0]
    zNodes[:,1] = p[:,1]
    zNodes[:,2] = range(0,len(p))

    # starting the process of removing the row that have the first and second column duplicated
    # we do not consider the third column in the duplication criterion
    # for example: z = [ 0 0 1; 1 0 2; 2 0 3; 1 0 4; 3 0 5]
    # in this ex we remove the duplicate and obtain:
    # z = [0 0 1; 1 0 2; 2 0 3; 3 0 5]
    keyfunc = lambda kf: kf[:2]
    mypoints = []
    for kind, gind in groupby( sorted( zip(zNodes[:,0], zNodes[:,1], zNodes[:,2]), key = keyfunc), keyfunc):
        mypoints.append(list(gind)[0])

    arr_mypoints = numpy.array(mypoints)
    #unique_nodes = range(0,len(p))#arr_mypoints[:,2]
    unique_nodes = arr_mypoints[:,2]
    
    extraNodes = []
    # from all the nodes remove those corresponding to the left and right boundary
    FreeNodes = list( ((set( unique_nodes) - set(lbcs)) - set(rbcs)) - set(extraNodes))

    eend2 = time.time()
    
    # Need to reduce the Kb matrix in order to be able to use it with SpSolve
    Kbreduced = numpy.zeros((len(FreeNodes), len(FreeNodes)))
    
    eend3 = time.time()
    for i in range(0,len(FreeNodes)):
        for j in range(0,len(FreeNodes)):
            Kbreduced[i,j] = K[FreeNodes[i],FreeNodes[j]]
    Kbreduced = sparse.csr_matrix(Kbreduced)

    # solve for the numerical solution
    numericalSoln = linsolve.spsolve(Kbreduced,F[FreeNodes,0])
    
    for i in range(0, len(FreeNodes)):
        U[FreeNodes[i],0] = numericalSoln[0][i]
    
    eend = time.time()
    print 'ELAPSED TIME = ', eend-sstart

#     HANGING NODES implementation
    for i in range(0,len(list_hanging_nodes)):
        listHN = list_hanging_nodes[i]
        U[listHN[0],0] = ( U[listHN[1],0] + U[listHN[2],0] ) / 2.0

    return  sparse.csr_matrix(U)

# semi-circle def:
def f_circle_s(x,y):
  return (x-1.0)**2 + (y-0.0)**2

# circle 1 def:
def f_circle1(x,y):
  return (x-0.25)**2 + (y-0.25)**2

#circle 2 def:
def f_circle2(x,y):
  return (x-0.5)**2 + (y-0.75)**2
  
def which_circle(node):
    A = abs( f_circle_s(node.x, node.y) - Rs**2)
    B = abs( f_circle1(node.x, node.y) - R1**2)
    C = abs( f_circle2(node.x, node.y) - R2**2)
    
    my_min = min(A,B,C)
    if my_min == A:
        c_rad = Rs
        c_center = Coordinate(1.0, 0.0)
    if my_min == B:
        c_rad = R1
        c_center = Coordinate(0.25, 0.25)
    if my_min == C:
        c_rad = R2
        c_center = Coordinate(0.5, 0.75)
    return [c_center, c_rad]  
  
def computeNorm(p,t,pConf,tConf,ui,wi,k1,k2,U,UConf,masterNode,llist, p_extra, P_quad, P_cub):
    print 'compute the L-2 norm ...'

    T = len(t)
    all_elems_sum = 0

    en_arr = [pp for pp in itertools.product(ui,repeat=2)]
    en_arr = np.array(en_arr)
    e_arr = en_arr[:,0] # contains the epsilon coordinates
    n_arr = en_arr[:,1] # contains the niu coordinates

    yloc = 0.1

    Usolution = np.zeros((len(p) + len(p_extra),1))
    polygonList = []

    # COMPUTING THE L-2 NORM
    for e in range(0,T):
        
        nodes = t[e] # row of t =  node numbers of the 4 corners of element e
        nodes = np.array(nodes)
    
        root = get_node_by_id(masterNode,llist[e])

        if len(root.enrichNodes) == 3:
            nodes6 = Coordinate(0,0)
            nodes6.x =  root.enrichNodes[2].x / 1000.0
            if nodes6.x != 0.0:
                nodes6.x += 0.001
            nodes6.y =  root.enrichNodes[2].y / 1000.0
            if nodes6.y != 0.0:
                nodes6.y += 0.001
            nodes6.y = 1 - nodes6.y
    
            old_node6 = Coordinate(nodes6.x, nodes6.y)
            [c_center, c_rad] = which_circle(nodes6)
            new_nodes6 = Coordinate(0.0,0.0)
            new_nodes6.x = c_center.x + c_rad * ( nodes6.x - c_center.x) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            new_nodes6.y = c_center.y + c_rad * ( nodes6.y - c_center.y) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            nodes6 = Coordinate(new_nodes6.x, new_nodes6.y)
            
        if len(root.enrichNodes) == 4:
            
            if root.enrichNodes[2].y == root.enrichNodes[3].y:
                if root.enrichNodes[2].x <= root.enrichNodes[3].x: 
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[2].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[2].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[3].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[3].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                    
                else:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[3].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[3].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[2].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[2].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
            else:
                if root.enrichNodes[2].y >= root.enrichNodes[3].y:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[2].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[2].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[3].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[3].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                else:
                    nodes6 = Coordinate(0,0)
                    nodes6.x =  root.enrichNodes[3].x / 1000.0
                    if nodes6.x != 0.0:
                        nodes6.x += 0.001
                    nodes6.y =  root.enrichNodes[3].y / 1000.0
                    if nodes6.y != 0.0:
                        nodes6.y += 0.001
                    nodes6.y = 1 - nodes6.y
                    
                    nodes7 = Coordinate(0,0)
                    nodes7.x =  root.enrichNodes[2].x / 1000.0
                    if nodes7.x != 0.0:
                        nodes7.x += 0.001
                    nodes7.y =  root.enrichNodes[2].y / 1000.0
                    if nodes7.y != 0.0:
                        nodes7.y += 0.001
                    nodes7.y = 1 - nodes7.y
                    
            old_node6 = Coordinate(nodes6.x, nodes6.y)
            [c_center, c_rad] = which_circle(nodes6)
            new_nodes6 = Coordinate(0.0,0.0)
            new_nodes6.x = c_center.x + c_rad * ( nodes6.x - c_center.x) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            new_nodes6.y = c_center.y + c_rad * ( nodes6.y - c_center.y) / math.sqrt( math.pow(nodes6.x - c_center.x, 2) + math.pow( nodes6.y - c_center.y,2) )
            nodes6 = Coordinate(new_nodes6.x, new_nodes6.y)
            
            old_node7 = Coordinate(nodes7.x, nodes7.y)
            [c_center, c_rad] = which_circle(nodes7)
            new_nodes7 = Coordinate(0.0,0.0)
            new_nodes7.x = c_center.x + c_rad * ( nodes7.x - c_center.x) / math.sqrt( math.pow(nodes7.x - c_center.x, 2) + math.pow( nodes7.y - c_center.y,2) )
            new_nodes7.y = c_center.y + c_rad * ( nodes7.y - c_center.y) / math.sqrt( math.pow(nodes7.x - c_center.x, 2) + math.pow( nodes7.y - c_center.y,2) )
            nodes7 = Coordinate(new_nodes7.x, new_nodes7.y)
           
        # 2-column matrix containing on each row the coordinates of each of the nodes
        coords = p[nodes,:]
        
        Pe = np.zeros((4,4))
        Pe[:,0] = np.ones((4,1)).transpose()
        Pe[:,1] = p[nodes[0:4],0]
        Pe[:,2] = p[nodes[0:4],1]
        Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

        C = np.linalg.inv(Pe)

        Nbasis = basisFct(C)
        c_fun = coords[1,1]
        d_fun = coords[2,1]

        x0 = coords[0,0]
        x1 = coords[1,0]
        y0 = coords[0,1]
        y1 = coords[2,1]

        thru_corner = False

        if len(nodes) == 5:
            enrich1 = np.array(p[nodes[4]])
    
            corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
            corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
            corner2 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )
            corner3 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )
    
            which_corner = -1
            if corner0 == True and corner1 == False and corner2 == False and corner3 == False:
                thru_corner = True
                which_corner = 1
            if corner0 == False and corner1 == True and corner2 == False and corner3 == False:
                thru_corner = True
                which_corner = 2
            if corner0 == False and corner1 == False and corner2 == True and corner3 == False:
                thru_corner = True
                which_corner = 3
            if corner0 == False and corner1 == False and corner2 == False and corner3 == True:
                thru_corner = True
                which_corner = 4 

        if (len(nodes) == 4 or  root.ishomog == 1) or thru_corner == True:
            
            x_coords = coords[:,0]
            y_coords = coords[:,1]

            uh_elem = lambda x,y: ( U[t[e][0],0] * Nbasis[0](x,y) +
                    U[t[e][1],0] * Nbasis[1](x,y) +
                    U[t[e][2],0] * Nbasis[2](x,y) +
                    U[t[e][3],0] * Nbasis[3](x,y) 
                    )

            Usolution[nodes[0],0] = uh_elem(p[nodes[0],0],p[nodes[0],1])
            Usolution[nodes[1],0] = uh_elem(p[nodes[1],0],p[nodes[1],1])
            Usolution[nodes[2],0] = uh_elem(p[nodes[2],0],p[nodes[2],1])
            Usolution[nodes[3],0] = uh_elem(p[nodes[3],0],p[nodes[3],1])
            
            polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], nodes[3] ]]

                  
# canceling the norm computation
            if NORM_COMP == 1:
                 # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
                 # for transformation from the parametric element to phisycal element
                 # of the Gauss nodes ui
                [x_transform_fct,y_transform_fct] = xy_fct(x_coords,y_coords)            
     
                Jac = jacobian_mat( coords[:,0], coords[:,1] )
                detJ = lambda eps,niu: determinant(Jac)(eps,niu)
                el_sum =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct,y_transform_fct,uh_elem,detJ)
                all_elems_sum = all_elems_sum + el_sum
                            
                print 'element e',e,' norm is :', el_sum
        else: # element has more than 4 nodes, meaning it is an element that needs enrichment at these additional nodes

            enrich1 = np.array(p[nodes[4]])
            if len(nodes) == 5:

                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y1])) <= 1e-12 )

                elCorner = False

                # FALSE POSITIVIES: only corner is in a different material
                # HOMOGENEOUS element
                if (corner1 == True or corner3 == True) or (corner0 == True or corner2 == True):
                    print 'False positive: element', e
                    elCorner = True
                    uh_elem = lambda x,y: ( 
                                U[t[e][0],0] * Nbasis[0](x,y) +
                                U[t[e][1],0] * Nbasis[1](x,y) +
                                U[t[e][2],0] * Nbasis[2](x,y) +
                                U[t[e][3],0] * Nbasis[3](x,y) 
                        )
    
                    x_coords = coords[:,0]
                    y_coords = coords[:,1]
                    
                    Usolution[nodes[0],0] = uh_elem(p[nodes[0],0],p[nodes[0],1])
                    Usolution[nodes[1],0] = uh_elem(p[nodes[1],0],p[nodes[1],1])
                    Usolution[nodes[2],0] = uh_elem(p[nodes[2],0],p[nodes[2],1])
                    Usolution[nodes[3],0] = uh_elem(p[nodes[3],0],p[nodes[3],1])

                    polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], nodes[3] ]]

# canceling the norm computation
                    if NORM_COMP == 1:
                         # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
                         # for transformation from the parametric element to phisycal element
                         # of the Gauss nodes ui
                        [x_transform_fct,y_transform_fct] = xy_fct(x_coords,y_coords)            
     
                        Jac = jacobian_mat( coords[:,0], coords[:,1] )
                        detJ = lambda eps,niu: determinant(Jac)(eps,niu)
     
                        el_sum =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct,y_transform_fct,uh_elem,detJ)
     
                        all_elems_sum = all_elems_sum + el_sum;

            else:
                enrich1 = np.array(p[nodes[4]])
                enrich2 = np.array(p[nodes[5]])

                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich2[0] - [x1]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich2[0] - [x0]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )

                if (corner0 == True and corner2 == True) or (corner1 == True and corner3 == True):
                    if(corner0 == True and corner2 == True):
                        nodes_trid1 = [nodes[0], nodes[1], nodes[2]]
                        nodes_trid2 = [nodes[0], nodes[2], nodes[3]]
                        tc1 = [0,1,2]
                        tc2 = [0,2,3]
                        even_diag = True
                    else:
                        nodes_trid1 = [nodes[0], nodes[1], nodes[3]]
                        nodes_trid2 = [nodes[1], nodes[2], nodes[3]]
                        tc1 = [0,1,3]
                        tc2 = [1,2,3]
                        even_diag = False

                    coords_trid1 = p[nodes_trid1]
                    coords_trid2 = p[nodes_trid2]
    
                    Pe_trid1 = np.zeros((3,3))
                    Pe_trid2 = np.zeros((3,3))
                        
                    Pe_trid1[:,0] = np.ones((3,1)).transpose()
                    Pe_trid1[:,1] = p[nodes_trid1[0:3],0]
                    Pe_trid1[:,2] = p[nodes_trid1[0:3],1]
                        
                    Pe_trid2[:,0] = np.ones((3,1)).transpose()
                    Pe_trid2[:,1] = p[nodes_trid2[0:3],0]
                    Pe_trid2[:,2] = p[nodes_trid2[0:3],1]
                        
                    C_trid1 = np.linalg.inv(Pe_trid1)
                    C_trid2 = np.linalg.inv(Pe_trid2)
                        
                    Nbasis_trid1 = tribasisFct(C_trid1)
                    Nbasis_trid2 = tribasisFct(C_trid2)
                        
                    Nx_trid1 = triderivX(C_trid1)
                    Ny_trid1 = triderivY(C_trid1)
                        
                    Nx_trid2 = triderivX(C_trid2)
                    Ny_trid2 = triderivY(C_trid2)
    
                    uh_elem_trid1 = lambda x,y: ( 
                                U[t[e][tc1[0]],0] * Nbasis_trid1[0](x,y) +
                                U[t[e][tc1[1]],0] * Nbasis_trid1[1](x,y) +
                                U[t[e][tc1[2]],0] * Nbasis_trid1[2](x,y)
                    )
    
                    uh_elem_trid2 = lambda x,y: ( 
                                U[t[e][tc2[0]],0] * Nbasis_trid2[0](x,y) +
                                U[t[e][tc2[1]],0] * Nbasis_trid2[1](x,y) +
                                U[t[e][tc2[2]],0] * Nbasis_trid2[2](x,y)
                    )
            
                    x_coords_trid1 = coords_trid1[:,0]
                    y_coords_trid1 = coords_trid1[:,1]
    
                    x_coords_trid2 = coords_trid2[:,0]
                    y_coords_trid2 = coords_trid2[:,1]
    
                    Usolution[nodes_trid1[0],0] = uh_elem_trid1(p[nodes_trid1[0],0],p[nodes_trid1[0],1])
                    Usolution[nodes_trid1[1],0] = uh_elem_trid1(p[nodes_trid1[1],0],p[nodes_trid1[1],1])
                    Usolution[nodes_trid1[2],0] = uh_elem_trid1(p[nodes_trid1[2],0],p[nodes_trid1[2],1])
                    
    
                    Usolution[nodes_trid2[0],0] = uh_elem_trid2(p[nodes_trid2[0],0],p[nodes_trid2[0],1])
                    Usolution[nodes_trid2[1],0] = uh_elem_trid2(p[nodes_trid2[1],0],p[nodes_trid2[1],1])
                    Usolution[nodes_trid2[2],0] = uh_elem_trid2(p[nodes_trid2[2],0],p[nodes_trid2[2],1])
    
                                                
                    if len(root.enrichNodes) < 3:
                        polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], nodes_trid1[2] ]]
                        polygonList = polygonList + [[nodes_trid2[0], nodes_trid2[1], nodes_trid2[2] ]]
                        
                    if len(root.enrichNodes) == 3:
                        index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                        pt6 = index6[0][0] + len(p)
                        Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                        if even_diag == True:
                            polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], pt6 ]]
                            polygonList = polygonList + [[pt6, nodes_trid1[1], nodes_trid1[2] ]]
                                
                            polygonList = polygonList + [[nodes_trid2[0], pt6, nodes_trid2[2] ]]
                            polygonList = polygonList + [[pt6, nodes_trid2[1], nodes_trid2[2] ]]
                        else:
                            
                            polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], pt6 ]]
                            polygonList = polygonList + [[pt6, nodes_trid1[0], nodes_trid1[2] ]]
                                
                            polygonList = polygonList + [[nodes_trid2[0], pt6, nodes_trid2[1] ]]
                            polygonList = polygonList + [[pt6, nodes_trid2[1], nodes_trid2[2] ]]
                           
                    if len(root.enrichNodes) == 4:
                        index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                        pt6 = index6[0][0] + len(p)
                        Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                        
                        index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                        pt7 = index7[0][0] + len(p)
                        Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                            
                        if even_diag == True:
                            polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], pt6 ]]
                            polygonList = polygonList + [[pt6, nodes_trid1[1], pt7 ]]
                            polygonList = polygonList + [[pt7, nodes_trid1[1], nodes_trid1[2] ]]
                                
                            polygonList = polygonList + [[nodes_trid2[0], pt6, nodes_trid2[2] ]]
                            polygonList = polygonList + [[pt6, pt7, nodes_trid2[2] ]]
                            polygonList = polygonList + [[pt7, nodes_trid2[1], nodes_trid2[2] ]]
                        else:
                            
                            polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], pt6 ]]
                            polygonList = polygonList + [[pt6, nodes_trid1[0], pt7 ]]
                            polygonList = polygonList + [[pt7, nodes_trid1[0], nodes_trid1[2] ]]
                                
                            polygonList = polygonList + [[nodes_trid2[0], pt6, nodes_trid2[1] ]]
                            polygonList = polygonList + [[pt6, nodes_trid2[1], pt7 ]]
                            polygonList = polygonList + [[pt7, nodes_trid2[1], nodes_trid2[2] ]]
                              
                                 
# canceling the norm computation
                    if NORM_COMP == 1:   
                        if len(root.enrichNodes) < 3:
                            [x_transform_fct_trid1,y_transform_fct_trid1] = tri_xy_fct(x_coords_trid1,y_coords_trid1)            
                            Jac_trid1 = tri_jacobian_mat( coords_trid1[:,0], coords_trid1[:,1] )
                            [x_transform_fct_trid2,y_transform_fct_trid2] = tri_xy_fct(x_coords_trid2,y_coords_trid2)            
                            Jac_trid2 = tri_jacobian_mat( coords_trid2[:,0], coords_trid2[:,1] )
                        
                        if len(root.enrichNodes) == 3:
                            coord_enrich = coord_enrich_comp_quad_circle(p[nodess], nodes6)
                              
                            if(corner0 == True and corner2 == True): 
                            # even diagonal: SW - NE
                                lOrd = [1,2,0] 
                            else:
                                lOrd = [0,1,2] 
                              
                            vec1_x = [ coords_trid1[ lOrd[0],0], coords_trid1[lOrd[1],0], coords_trid1[lOrd[2],0], (coords_trid1[ lOrd[0],0] + coords_trid1[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid1[lOrd[0],0] + coords_trid1[lOrd[2],0])/2.0  ]
                            vec1_y = [ coords_trid1[ lOrd[0],1], coords_trid1[lOrd[1],1], coords_trid1[lOrd[2],1], (coords_trid1[ lOrd[0],1] + coords_trid1[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid1[lOrd[0],1] + coords_trid1[lOrd[2],1])/2.0  ]
                    
                            [x_transform_fct_trid1,y_transform_fct_trid1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                            Jac_trid1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                            
                            
                            if(corner0 == True and corner2 == True): 
                            # even diagonal: SW - NE
                                lOrd = [2,0,1] 
                            else:
                                lOrd = [1,2,0] 
                                
                            vec2_x = [ coords_trid2[ lOrd[0],0], coords_trid2[lOrd[1],0], coords_trid2[lOrd[2],0], (coords_trid2[ lOrd[0],0] + coords_trid2[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid2[lOrd[0],0] + coords_trid2[lOrd[2],0])/2.0  ]
                            vec2_y = [ coords_trid2[ lOrd[0],1], coords_trid2[lOrd[1],1], coords_trid2[lOrd[2],1], (coords_trid2[ lOrd[0],1] + coords_trid2[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid2[lOrd[0],1] + coords_trid2[lOrd[2],1])/2.0  ]
                    
                            [x_transform_fct_trid2,y_transform_fct_trid2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                            Jac_trid2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                            
                        if len(root.enrichNodes) == 4:
                            
                            coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess], nodes6)
                            coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess], nodes7)
                            
                            circumcenter_pt1 = circumcenter_tri(coords_trid1)
                            circumcenter_pt2 = circumcenter_tri(coords_trid2)
    
                            pt4 = Coordinate(0,0)
                            pt5 = Coordinate(0,0)
                            pt6 = Coordinate(0,0)
                            pt7 = Coordinate(0,0)
                            pt8 = Coordinate(0,0)
                            pt9 = Coordinate(0,0)
                            
                            if(corner0 == True and corner2 == True): 
                            # even diagonal: SW - NE
                                lOrd = [1,2,0] 
                                pt6.x = coord_enrich2.x
                                pt6.y = coord_enrich2.y 
                            else:
                                lOrd = [0,1,2]   
                                pt7.x = coord_enrich1.x
                                pt7.y = coord_enrich1.y                      
                                                       
                             
                            pt4.x = 2.0/3.0 * coords_trid1[lOrd[0],0] + 1.0/3.0 * coords_trid1[lOrd[1],0]
                            pt4.y = 2.0/3.0 * coords_trid1[lOrd[0],1] + 1.0/3.0 * coords_trid1[lOrd[1],1]
                             
                            pt5.x = 1.0/3.0 * coords_trid1[lOrd[0],0] + 2.0/3.0 * coords_trid1[lOrd[1],0]
                            pt5.y = 1.0/3.0 * coords_trid1[lOrd[0],1] + 2.0/3.0 * coords_trid1[lOrd[1],1]
                             
                            pt8.x = 2.0/3.0 * coords_trid1[lOrd[2],0] + 1.0/3.0 * coords_trid1[lOrd[0],0]
                            pt8.y = 2.0/3.0 * coords_trid1[lOrd[2],1] + 1.0/3.0 * coords_trid1[lOrd[0],1]
                             
                            pt9.x = 1.0/3.0 * coords_trid1[lOrd[2],0] + 2.0/3.0 * coords_trid1[lOrd[0],0]
                            pt9.y = 1.0/3.0 * coords_trid1[lOrd[2],1] + 2.0/3.0 * coords_trid1[lOrd[0],1]
                     
                             
                            vec1_x = [coords_trid1[lOrd[0],0], 
                                      coords_trid1[lOrd[1],0], 
                                      coords_trid1[lOrd[2],0], 
                                      pt4.x, 
                                      pt5.x,
                                      pt6.x,
                                      pt7.x,
                                      pt8.x,
                                      pt9.x,
                                      circumcenter_pt1.x  
                                      ]
                            vec1_y = [coords_trid1[lOrd[0],1], 
                                      coords_trid1[lOrd[1],1], 
                                      coords_trid1[lOrd[2],1], 
                                      pt4.y, 
                                      pt5.y,
                                      pt6.y,
                                      pt7.y,
                                      pt8.y,
                                      pt9.y,
                                      circumcenter_pt1.y  
                                      ]
                            [x_transform_fct_trid1,y_transform_fct_trid1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                            Jac_trid2 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                            
    
                            pt4 = Coordinate(0,0)
                            pt5 = Coordinate(0,0)
                            pt6 = Coordinate(0,0)
                            pt7 = Coordinate(0,0)
                            pt8 = Coordinate(0,0)
                            pt9 = Coordinate(0,0)
                            
                            if(corner0 == True and corner2 == True): 
                            # even diagonal: SW - NE
                                lOrd = [2,0,1] 
                                pt6.x = coord_enrich1.x
                                pt6.y = coord_enrich1.y  
                            else:
                                lOrd = [1,2,0] 
                                pt7.x = coord_enrich2.x
                                pt7.y = coord_enrich2.y
                                
                            pt4.x = 2.0/3.0 * coords_trid2[lOrd[0],0] + 1.0/3.0 * coords_trid2[lOrd[1],0]
                            pt4.y = 2.0/3.0 * coords_trid2[lOrd[0],1] + 1.0/3.0 * coords_trid2[lOrd[1],1]
                            
                            pt5.x = 1.0/3.0 * coords_trid2[lOrd[0],0] + 2.0/3.0 * coords_trid2[lOrd[1],0]
                            pt5.y = 1.0/3.0 * coords_trid2[lOrd[0],1] + 2.0/3.0 * coords_trid2[lOrd[1],1]
                            
                            pt8.x = 2.0/3.0 * coords_trid2[lOrd[2],0] + 1.0/3.0 * coords_trid2[lOrd[0],0]
                            pt8.y = 2.0/3.0 * coords_trid2[lOrd[2],1] + 1.0/3.0 * coords_trid2[lOrd[0],1]
                            
                            pt9.x = 1.0/3.0 * coords_trid2[lOrd[2],0] + 2.0/3.0 * coords_trid2[lOrd[0],0]
                            pt9.y = 1.0/3.0 * coords_trid2[lOrd[2],1] + 2.0/3.0 * coords_trid2[lOrd[0],1]
                        
                            vec2_x = [coords_trid2[lOrd[0],0], 
                                      coords_trid2[lOrd[1],0],
                                      coords_trid2[lOrd[2],0],
                                      pt4.x, 
                                      pt5.x,
                                      pt6.x,
                                      pt7.x,
                                      pt8.x,
                                      pt9.x,
                                      circumcenter_pt2.x  
                                      ]
                            vec2_y = [coords_trid2[lOrd[0],1], 
                                      coords_trid2[lOrd[1],1], 
                                      coords_trid2[lOrd[2],1], 
                                      pt4.y, 
                                      pt5.y,
                                      pt6.y, 
                                      pt7.y,
                                      pt8.y,
                                      pt9.y,
                                      circumcenter_pt2.y 
                                      ]
                            
                            [x_transform_fct_trid2,y_transform_fct_trid2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                            Jac_trid2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                              
                        detJ_trid1 = lambda eps,niu: determinant(Jac_trid1)(eps,niu)
                        el_sum_trid1 =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_trid1,y_transform_fct_trid1,uh_elem_trid1,detJ_trid1)
                        all_elems_sum = all_elems_sum + el_sum_trid1;
    
                        detJ_trid2 = lambda eps,niu: determinant(Jac_trid2)(eps,niu)
                        el_sum_trid2 =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_trid2,y_transform_fct_trid2,uh_elem_trid2,detJ_trid2)
                        all_elems_sum = all_elems_sum + el_sum_trid2
    
    
                else:
                    if ( ( (enrich1[1] == y0 and enrich2[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1)) or
                          ( (enrich1[1] == y1 and enrich2[1] == y0) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x0) )  ):
                        print 'norm computation: North edge'
    
                        if ( ( (enrich1[1] == y0 and enrich2[1] == y1) and (enrich1[0] == x0) and (enrich2[0] != x0 and enrich2[0] != x1)) or
                          ( (enrich1[1] == y1 and enrich2[1] == y0) and (enrich2[0] == x0) and (enrich1[0] != x0 and enrich1[0] != x0) )  ):
                            triangle_1 = True
                        else:
                            triangle_1 = False
                                                        
                        if not(on_corners(enrich2,coords)):
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5] ]
                        else:
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] ]
    
                        tri_nodes1 = [north_nodes[0],north_nodes[4],north_nodes[3]]
                        tri_nodes2 = [north_nodes[0],north_nodes[1],north_nodes[4]]
                        tri_nodes3 = [north_nodes[1],north_nodes[2],north_nodes[4]]
                        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
    
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
        
    
                        node_enr_4 = [north_nodes[4]]
                        coords_enr_4 = p[node_enr_4]
    
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
    
                        n1 = x4 - x0
                        n2 = x1 - x4
                        factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
                        if factor_N > EPS_FACTOR:
                            factor_N = 1
    
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_N )
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
                                                
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_N)
# canceling the norm computation  
                        if NORM_COMP == 1:
                            
                            if len(root.enrichNodes) == 2:
                                [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                
                                J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
        
                            if len(root.enrichNodes) == 3:
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                    
                                    lOrd = [1,2,0]
                                    vec3_x = [ tri_coords3[ lOrd[0],0], tri_coords3[lOrd[1],0], tri_coords3[lOrd[2],0], (tri_coords3[ lOrd[0],0] + tri_coords3[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords3[lOrd[0],0] + tri_coords3[lOrd[2],0])/2.0  ]
                                    vec3_y = [ tri_coords3[ lOrd[0],1], tri_coords3[lOrd[1],1], tri_coords3[lOrd[2],1], (tri_coords3[ lOrd[0],1] + tri_coords3[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords3[lOrd[0],1] + tri_coords3[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is curved
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [2,0,1]
                                    vec1_x = [ tri_coords1[ lOrd[0],0], tri_coords1[lOrd[1],0], tri_coords1[lOrd[2],0], (tri_coords1[ lOrd[0],0] + tri_coords1[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords1[lOrd[0],0] + tri_coords1[lOrd[2],0])/2.0  ]
                                    vec1_y = [ tri_coords1[ lOrd[0],1], tri_coords1[lOrd[1],1], tri_coords1[lOrd[2],1], (tri_coords1[ lOrd[0],1] + tri_coords1[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords1[lOrd[0],1] + tri_coords1[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                                                
                                    lOrd = [1,2,0] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )                                
                                
                            if len(root.enrichNodes) == 4:
                        
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes], nodes7)
                                
                                circumcenter_pt1 = circumcenter_tri(tri_coords1)
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                circumcenter_pt3 = circumcenter_tri(tri_coords3)
                        
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2] # local order    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                    
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                            
                                    lOrd = [1,2,0]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords3[lOrd[0],0] + 1.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords3[lOrd[0],1] + 1.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords3[lOrd[0],0] + 2.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords3[lOrd[0],1] + 2.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords3[lOrd[2],0] + 1.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords3[lOrd[2],1] + 1.0/3.0 * tri_coords3[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords3[lOrd[2],0] + 2.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords3[lOrd[2],1] + 2.0/3.0 * tri_coords3[lOrd[0],1]
                             
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                     
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                    
                                    vec3_x = [tri_coords3[ lOrd[0],0], 
                                              tri_coords3[lOrd[1],0],
                                              tri_coords3[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt3.x  
                                              ]
                                    vec3_y = [tri_coords3[ lOrd[0],1], 
                                              tri_coords3[lOrd[1],1], 
                                              tri_coords3[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt3.y  
                                              ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [2,0,1]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords1[lOrd[0],0] + 1.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords1[lOrd[0],1] + 1.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords1[lOrd[0],0] + 2.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords1[lOrd[0],1] + 2.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords1[lOrd[2],0] + 1.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords1[lOrd[2],1] + 1.0/3.0 * tri_coords1[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords1[lOrd[2],0] + 2.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords1[lOrd[2],1] + 2.0/3.0 * tri_coords1[lOrd[0],1]
                             
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                     
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                     
                                    vec1_x = [tri_coords1[lOrd[0],0], 
                                              tri_coords1[lOrd[1],0], 
                                              tri_coords1[lOrd[2],0], 
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt1.x  
                                              ]
                                    vec1_y = [tri_coords1[lOrd[0],1], 
                                              tri_coords1[lOrd[1],1], 
                                              tri_coords1[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt1.y  
                                              ]
                                    [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                                    
                                                
                                    lOrd = [1,2,0] # local order    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                    
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )  
                                                                   
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                            detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
                        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)      
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
            
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_2( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
    
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
                        
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], pt6, tri_nodes1[2]]]
                                polygonList = polygonList + [[pt6, tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], tri_nodes3[2]]]
                                
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], pt6, tri_nodes1[2]]]
                                polygonList = polygonList + [[pt6, pt7, tri_nodes1[2]]]
                                polygonList = polygonList + [[pt7, tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], pt7 ]]
                                polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[tri_nodes2[0], pt6, pt7 ]]
                                polygonList = polygonList + [[tri_nodes2[0], pt7, tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], pt7 ]]
                                polygonList = polygonList + [[pt7, tri_nodes3[1], tri_nodes3[2]]]
                            
        
                    # the South edge
                    if (  ((enrich1[1] == y0 and enrich2[1] == y1) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x1) ) or
                          ( (enrich2[1] == y0 and enrich1[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1 ) ) ):
                        print 'norm computation: South edge'
                        
                        if (  ((enrich1[1] == y0 and enrich2[1] == y1) and (enrich2[0] == x0) and (enrich1[0] != x0 and enrich1[0] != x1) ) or
                          ( (enrich2[1] == y0 and enrich1[1] == y1) and (enrich1[0] == x0) and (enrich2[0] != x0 and enrich2[0] != x1 ) ) ):
                            triangle_1 = True
                        else:
                            triangle_1 = False

                        if not(on_corners(enrich2,coords)) :
                          south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                          south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
                        tri_nodes1 = [south_nodes[0],south_nodes[4],south_nodes[3]]
                        tri_nodes2 = [south_nodes[4],south_nodes[2],south_nodes[3]]
                        tri_nodes3 = [south_nodes[4],south_nodes[1],south_nodes[2]] # the one triangle in a diff material
                        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
        
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
        
        
                        node_enr_4 = [south_nodes[4]]
                        coords_enr_4 = p[node_enr_4]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        s1 = x4 - x0
                        s2 = x1 - x4
                        factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
                        if factor_S > EPS_FACTOR:
                            factor_S = 1
    
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_S)
            
        
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_S )
        
        
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_S )
# canceling the norm computation   
                        if NORM_COMP == 1: 
                            if len(root.enrichNodes) < 3:
                                [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                
                
                                J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                
                            if len(root.enrichNodes) == 3:
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
        
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [2,0,1] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                    
                                    lOrd = [1,2,0]
                                    vec3_x = [ tri_coords3[ lOrd[0],0], tri_coords3[lOrd[1],0], tri_coords3[lOrd[2],0], (tri_coords3[ lOrd[0],0] + tri_coords3[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords3[lOrd[0],0] + tri_coords3[lOrd[2],0])/2.0  ]
                                    vec3_y = [ tri_coords3[ lOrd[0],1], tri_coords3[lOrd[1],1], tri_coords3[lOrd[2],1], (tri_coords3[ lOrd[0],1] + tri_coords3[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords3[lOrd[0],1] + tri_coords3[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [0,1,2]
                                    vec1_x = [ tri_coords1[ lOrd[0],0], tri_coords1[lOrd[1],0], tri_coords1[lOrd[2],0], (tri_coords1[ lOrd[0],0] + tri_coords1[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords1[lOrd[0],0] + tri_coords1[lOrd[2],0])/2.0  ]
                                    vec1_y = [ tri_coords1[ lOrd[0],1], tri_coords1[lOrd[1],1], tri_coords1[lOrd[2],1], (tri_coords1[ lOrd[0],1] + tri_coords1[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords1[lOrd[0],1] + tri_coords1[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                                                
                                    lOrd = [1,2,0] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                        
                            if len(root.enrichNodes) == 4:
                        
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes], nodes7)
                                
                                circumcenter_pt1 = circumcenter_tri(tri_coords1)
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                circumcenter_pt3 = circumcenter_tri(tri_coords3)
                                
                                            
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [2,0,1] # local order    
                                    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                    
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                    
                                            
                                    lOrd = [1,2,0]
                                    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords3[lOrd[0],0] + 1.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords3[lOrd[0],1] + 1.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords3[lOrd[0],0] + 2.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords3[lOrd[0],1] + 2.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords3[lOrd[2],0] + 1.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords3[lOrd[2],1] + 1.0/3.0 * tri_coords3[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords3[lOrd[2],0] + 2.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords3[lOrd[2],1] + 2.0/3.0 * tri_coords3[lOrd[0],1]
                             
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                     
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                    
                                    vec3_x = [tri_coords3[ lOrd[0],0], 
                                              tri_coords3[lOrd[1],0],
                                              tri_coords3[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt3.x  
                                              ]
                                    vec3_y = [tri_coords3[ lOrd[0],1], 
                                              tri_coords3[lOrd[1],1], 
                                              tri_coords3[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt3.y  
                                              ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [0,1,2]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords1[lOrd[0],0] + 1.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords1[lOrd[0],1] + 1.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords1[lOrd[0],0] + 2.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords1[lOrd[0],1] + 2.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords1[lOrd[2],0] + 1.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords1[lOrd[2],1] + 1.0/3.0 * tri_coords1[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords1[lOrd[2],0] + 2.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords1[lOrd[2],1] + 2.0/3.0 * tri_coords1[lOrd[0],1]
                             
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                     
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                     
                                    vec1_x = [tri_coords1[lOrd[0],0], 
                                              tri_coords1[lOrd[1],0], 
                                              tri_coords1[lOrd[2],0], 
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt1.x  
                                              ]
                                    vec1_y = [tri_coords1[lOrd[0],1], 
                                              tri_coords1[lOrd[1],1], 
                                              tri_coords1[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt1.y  
                                              ]
                                    [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                                    
                                                
                                    lOrd = [1,2,0] # local order    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                    
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )                                
                                
                                   
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
       
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
    
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_3( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
        
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
        
        
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[tri_nodes1[0], pt6, tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], tri_nodes3[2]]]


                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                            
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[tri_nodes1[0], pt6, pt7]]
                                polygonList = polygonList + [[tri_nodes1[0], pt7, tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], pt7 ]]
                                polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt6, pt7, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], pt7]]
                                polygonList = polygonList + [[pt7, tri_nodes3[1], tri_nodes3[2]]]

    
                    # the West edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1)) or
                          ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)  ) ):
                        print 'norm computation: West edge'
                        
                        
                        if ( 
                            ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1))
                             or
                            ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)  )
                             ) :
                            triangle_1 = True
                        else:
                            triangle_1 = False  
                        
                        if not(on_corners(enrich2,coords)) :
                          west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
                        tri_nodes1 = [west_nodes[0],west_nodes[1],west_nodes[4]]
                        tri_nodes2 = [west_nodes[1],west_nodes[2],west_nodes[4]]
                        tri_nodes3 = [west_nodes[4],west_nodes[2],west_nodes[3]] 
        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
    
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
        
    
                        # scaling factor
                        node_enr_4 = [west_nodes[4]]
                        coords_enr_4 = p[node_enr_4]
    
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
    
                        w1 = y1 - y4
                        w2 = y4 - y0
                        factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
                        if factor_W > EPS_FACTOR:
                            factor_W = 1
    
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[2],0] * Nbasis_tri1[2](x,y) * factor_W )
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_W )
                                                
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_W)
# canceling the norm computation 
                        if NORM_COMP == 1:
                            if len(root.enrichNodes) < 3:
                                
                                [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                
                                J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                
                            if len(root.enrichNodes) == 3:
                                
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                    
                                    lOrd = [2,0,1]
                                    vec3_x = [ tri_coords3[ lOrd[0],0], tri_coords3[lOrd[1],0], tri_coords3[lOrd[2],0], (tri_coords3[ lOrd[0],0] + tri_coords3[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords3[lOrd[0],0] + tri_coords3[lOrd[2],0])/2.0  ]
                                    vec3_y = [ tri_coords3[ lOrd[0],1], tri_coords3[lOrd[1],1], tri_coords3[lOrd[2],1], (tri_coords3[ lOrd[0],1] + tri_coords3[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords3[lOrd[0],1] + tri_coords3[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [0,1,2]
                                    vec1_x = [ tri_coords1[ lOrd[0],0], tri_coords1[lOrd[1],0], tri_coords1[lOrd[2],0], (tri_coords1[ lOrd[0],0] + tri_coords1[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords1[lOrd[0],0] + tri_coords1[lOrd[2],0])/2.0  ]
                                    vec1_y = [ tri_coords1[ lOrd[0],1], tri_coords1[lOrd[1],1], tri_coords1[lOrd[2],1], (tri_coords1[ lOrd[0],1] + tri_coords1[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords1[lOrd[0],1] + tri_coords1[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                                                
                                    lOrd = [1,2,0] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                           
                            if len(root.enrichNodes) == 4:
                        
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes], nodes7)
                                
                                circumcenter_pt1 = circumcenter_tri(tri_coords1)
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                circumcenter_pt3 = circumcenter_tri(tri_coords3)
                                
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2] # local order    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                    
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                            
                                    lOrd = [2,0,1]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords3[lOrd[0],0] + 1.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords3[lOrd[0],1] + 1.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords3[lOrd[0],0] + 2.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords3[lOrd[0],1] + 2.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords3[lOrd[2],0] + 1.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords3[lOrd[2],1] + 1.0/3.0 * tri_coords3[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords3[lOrd[2],0] + 2.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords3[lOrd[2],1] + 2.0/3.0 * tri_coords3[lOrd[0],1]
                             
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                     
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                    
                                    vec3_x = [tri_coords3[ lOrd[0],0], 
                                              tri_coords3[lOrd[1],0],
                                              tri_coords3[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt3.x  
                                              ]
                                    vec3_y = [tri_coords3[ lOrd[0],1], 
                                              tri_coords3[lOrd[1],1], 
                                              tri_coords3[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt3.y  
                                              ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    lOrd = [0,1,2]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords1[lOrd[0],0] + 1.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords1[lOrd[0],1] + 1.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords1[lOrd[0],0] + 2.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords1[lOrd[0],1] + 2.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords1[lOrd[2],0] + 1.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords1[lOrd[2],1] + 1.0/3.0 * tri_coords1[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords1[lOrd[2],0] + 2.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords1[lOrd[2],1] + 2.0/3.0 * tri_coords1[lOrd[0],1]
                             
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                     
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                     
                                    vec1_x = [tri_coords1[lOrd[0],0], 
                                              tri_coords1[lOrd[1],0], 
                                              tri_coords1[lOrd[2],0], 
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt1.x  
                                              ]
                                    vec1_y = [tri_coords1[lOrd[0],1], 
                                              tri_coords1[lOrd[1],1], 
                                              tri_coords1[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt1.y  
                                              ]
                                    [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                                                
                                    lOrd = [1,2,0] # local order    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                    
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )                                   
                                
                                
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)       
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_3( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
                        
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
        
        
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[tri_nodes1[0], pt6, tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[2]]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], tri_nodes3[2]]]
                                
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[tri_nodes1[0], pt6, pt7 ]]
                                polygonList = polygonList + [[tri_nodes1[0], pt7, tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], pt7]]
                                polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt7]]
                                polygonList = polygonList + [[pt6, pt7, tri_nodes2[0]]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[2]]]
                                polygonList = polygonList + [[pt6, pt7, tri_nodes3[2]]]
                                polygonList = polygonList + [[pt7, tri_nodes3[1], tri_nodes3[2]]]
                          
    
                    # the East edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)) or
                          ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1) )  ):
                        print 'norm computation: East edge'
    
                        if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)) or
                          ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1) )  ):
                            triangle_1 = True
                        else:
                            triangle_1 = False

                        if not(on_corners(enrich2,coords)) :
                          east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
    
                        tri_nodes1 = [east_nodes[0],east_nodes[1],east_nodes[4]]
                        tri_nodes2 = [east_nodes[0],east_nodes[4],east_nodes[3]]
                        tri_nodes3 = [east_nodes[4],east_nodes[2],east_nodes[3]] # the one triangle in a diff material
        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
        
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
        
        
                        node_enr_4 = [east_nodes[4]]
                        coords_enr_4 = p[node_enr_4]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
    
                        e1 = y1 - y4
                        e2 = y4 - y0
                        factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
                        if factor_E > EPS_FACTOR:
                            factor_E = 1
        
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[2],0] * Nbasis_tri1[2](x,y) * factor_E)
            
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E )
        
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_E )
# canceling the norm computation  
                        if NORM_COMP == 1:  
                            if len(root.enrichNodes) < 3:
                                [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                
                                J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            if len(root.enrichNodes) == 3:
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                    
                                    lOrd = [1,2,0]
                                    vec3_x = [ tri_coords3[ lOrd[0],0], tri_coords3[lOrd[1],0], tri_coords3[lOrd[2],0], (tri_coords3[ lOrd[0],0] + tri_coords3[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords3[lOrd[0],0] + tri_coords3[lOrd[2],0])/2.0  ]
                                    vec3_y = [ tri_coords3[ lOrd[0],1], tri_coords3[lOrd[1],1], tri_coords3[lOrd[2],1], (tri_coords3[ lOrd[0],1] + tri_coords3[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords3[lOrd[0],1] + tri_coords3[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    
                                    lOrd = [1,2,0]
                                    vec1_x = [ tri_coords1[ lOrd[0],0], tri_coords1[lOrd[1],0], tri_coords1[lOrd[2],0], (tri_coords1[ lOrd[0],0] + tri_coords1[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords1[lOrd[0],0] + tri_coords1[lOrd[2],0])/2.0  ]
                                    vec1_y = [ tri_coords1[ lOrd[0],1], tri_coords1[lOrd[1],1], tri_coords1[lOrd[2],1], (tri_coords1[ lOrd[0],1] + tri_coords1[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords1[lOrd[0],1] + tri_coords1[lOrd[2],1])/2.0  ]
                        
                                    [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                                                
                                    lOrd = [2,0,1] # local order    
                                    vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                    vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                        
                                
                            if len(root.enrichNodes) == 4:
                                               
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes], nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes], nodes7)
                                
                                circumcenter_pt1 = circumcenter_tri(tri_coords1)
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                circumcenter_pt3 = circumcenter_tri(tri_coords3)
                                
                                if triangle_1 == False:
                                    # triangle 3 is the one with curved edge
                                    [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                    J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                    
                                    lOrd = [0,1,2]
                                
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                    
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                            
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                
                                    lOrd = [1,2,0]
                                    
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords3[lOrd[0],0] + 1.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords3[lOrd[0],1] + 1.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords3[lOrd[0],0] + 2.0/3.0 * tri_coords3[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords3[lOrd[0],1] + 2.0/3.0 * tri_coords3[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords3[lOrd[2],0] + 1.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords3[lOrd[2],1] + 1.0/3.0 * tri_coords3[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords3[lOrd[2],0] + 2.0/3.0 * tri_coords3[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords3[lOrd[2],1] + 2.0/3.0 * tri_coords3[lOrd[0],1]
                             
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                     
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                    
                                    vec3_x = [tri_coords3[ lOrd[0],0], 
                                              tri_coords3[lOrd[1],0],
                                              tri_coords3[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt3.x  
                                              ]
                                    vec3_y = [tri_coords3[ lOrd[0],1], 
                                              tri_coords3[lOrd[1],1], 
                                              tri_coords3[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt3.y  
                                              ]
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
                                    J_tri3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
                                    
                                if triangle_1 == True:
                                    # triangle 1 is the one with curved edge
                        
                                    [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                                    J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                                    
                                    
                                    lOrd = [1,2,0]
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                     
                                    pt4.x = 2.0/3.0 * tri_coords1[lOrd[0],0] + 1.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords1[lOrd[0],1] + 1.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt5.x = 1.0/3.0 * tri_coords1[lOrd[0],0] + 2.0/3.0 * tri_coords1[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords1[lOrd[0],1] + 2.0/3.0 * tri_coords1[lOrd[1],1]
                                     
                                    pt8.x = 2.0/3.0 * tri_coords1[lOrd[2],0] + 1.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords1[lOrd[2],1] + 1.0/3.0 * tri_coords1[lOrd[0],1]
                                     
                                    pt9.x = 1.0/3.0 * tri_coords1[lOrd[2],0] + 2.0/3.0 * tri_coords1[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords1[lOrd[2],1] + 2.0/3.0 * tri_coords1[lOrd[0],1]
                             
                                    pt6.x = coord_enrich2.x
                                    pt6.y = coord_enrich2.y     
                                     
                                    pt7.x = coord_enrich1.x
                                    pt7.y = coord_enrich1.y
                                     
                                    vec1_x = [tri_coords1[lOrd[0],0], 
                                              tri_coords1[lOrd[1],0], 
                                              tri_coords1[lOrd[2],0], 
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt1.x  
                                              ]
                                    vec1_y = [tri_coords1[lOrd[0],1], 
                                              tri_coords1[lOrd[1],1], 
                                              tri_coords1[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y,
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt1.y  
                                              ]
                                    [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                                    J_tri1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                                                
                                    lOrd = [2,0,1] 
                                      
                                    pt4 = Coordinate(0,0)
                                    pt5 = Coordinate(0,0)
                                    pt6 = Coordinate(0,0)
                                    pt7 = Coordinate(0,0)
                                    pt8 = Coordinate(0,0)
                                    pt9 = Coordinate(0,0)
                                
                                    pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                    pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                    
                                    pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                    
                                    pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                    pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                            
                                    pt6.x = coord_enrich1.x
                                    pt6.y = coord_enrich1.y     
                                    
                                    pt7.x = coord_enrich2.x
                                    pt7.y = coord_enrich2.y
                                
                                    vec2_x = [tri_coords2[lOrd[0],0], 
                                              tri_coords2[lOrd[1],0],
                                              tri_coords2[lOrd[2],0],
                                              pt4.x, 
                                              pt5.x,
                                              pt6.x,
                                              pt7.x,
                                              pt8.x,
                                              pt9.x,
                                              circumcenter_pt2.x  
                                              ]
                                    vec2_y = [tri_coords2[lOrd[0],1], 
                                              tri_coords2[lOrd[1],1], 
                                              tri_coords2[lOrd[2],1], 
                                              pt4.y, 
                                              pt5.y,
                                              pt6.y, 
                                              pt7.y,
                                              pt8.y,
                                              pt9.y,
                                              circumcenter_pt2.y 
                                              ]
                                    
                                    [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                    J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                                                    
                                
                                
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                                                    
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
         
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_3( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
                        
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], tri_nodes3[2]]]
                                
                        
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                            
                            
                            if triangle_1 == True:
                                # triangle 1 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                                polygonList = polygonList + [[pt6, tri_nodes1[1], pt7]]
                                polygonList = polygonList + [[pt7, tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt6, pt7, tri_nodes2[2]]]
                                polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            else:
                                # triangle 3 is curved
                                polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                                
                                polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                                polygonList = polygonList + [[tri_nodes2[0], pt6, pt7]]
                                polygonList = polygonList + [[tri_nodes2[0], pt7, tri_nodes2[2]]]
                                
                                polygonList = polygonList + [[tri_nodes3[0], pt6, tri_nodes3[1]]]
                                polygonList = polygonList + [[pt6, tri_nodes3[1], pt7]]
                                polygonList = polygonList + [[pt7 , tri_nodes3[1], tri_nodes3[2]]]
                                
                    
                    # the North-West corner is cut, 0-4-3, 2-5-3
                    if ( 
                        ((enrich1[0] == x0 and enrich2[1] == y1) or
                        (enrich2[0] == x0 and enrich1[1] == y1)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
                        ):
                        
                        print "norm computation: NW corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (lower right corner cut)
                        tri_nodes1 = [nodes[0],nodes[1],nodes[4]]
                        tri_nodes2 = [nodes[4],nodes[1],nodes[5]]
                        tri_nodes3 = [nodes[1],nodes[2],nodes[5]]
                        tri_nodes4 = [nodes[4],nodes[5],nodes[3]] 
    
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
    
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
    
                        Pe4 = np.zeros((3,3))
                        Pe4[:,0] = np.ones((3,1)).transpose()
                        Pe4[:,1] = p[tri_nodes4[0:3],0]
                        Pe4[:,2] = p[tri_nodes4[0:3],1]
                        C4 = np.linalg.inv(Pe4)
                        Nbasis_tri4 = tribasisFct(C4)
                        Nx_tri4 = triderivX(C4)
                        Ny_tri4 = triderivY(C4)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
                        tri_coords4 = p[tri_nodes4]
            
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
    
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
    
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
    
                        w1 = y1 - y4
                        w2 = y4 - y0
                        factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
                        if factor_W > EPS_FACTOR:
                            factor_W = 1
    
                        n1 = x5 - x0
                        n2 = x1 - x5
                        factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
                        if factor_N > EPS_FACTOR:
                            factor_N = 1
    
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[2],0] * Nbasis_tri1[2](x,y) * factor_W )
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_W +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_N)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_W +
                                                U[tri_nodes4[1],0] * Nbasis_tri4[1](x,y) * factor_N )
# canceling the norm computation        
                                   
                        if NORM_COMP == 1:
                            [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                            [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                            J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            if len(root.enrichNodes) < 3:
                                
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
                
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
                            
                            if len(root.enrichNodes) == 3:
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                
                                lOrd = [1,2,0] # local order 
                                vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                        
                                [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                
                                lOrd = [2,0,1]
                                vec4_x = [ tri_coords4[ lOrd[0],0], tri_coords4[lOrd[1],0], tri_coords4[lOrd[2],0], (tri_coords4[ lOrd[0],0] + tri_coords4[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords4[lOrd[0],0] + tri_coords4[lOrd[2],0])/2.0  ]
                                vec4_y = [ tri_coords4[ lOrd[0],1], tri_coords4[lOrd[1],1], tri_coords4[lOrd[2],1], (tri_coords4[ lOrd[0],1] + tri_coords4[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords4[lOrd[0],1] + tri_coords4[lOrd[2],1])/2.0  ]
                        
                                [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
        
                            if len(root.enrichNodes) == 4:
                        
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                                
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                
                                lOrd = [1,2,0]
                                
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                        
                                pt6.x = coord_enrich2.x
                                pt6.y = coord_enrich2.y     
                                
                                pt7.x = coord_enrich1.x
                                pt7.y = coord_enrich1.y
                                
                                vec2_x = [ tri_coords2[ lOrd[0],0], 
                                          tri_coords2[lOrd[1],0],
                                          tri_coords2[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt2.x  
                                          ]
                                vec2_y = [ tri_coords2[ lOrd[0],1], 
                                          tri_coords2[lOrd[1],1], 
                                          tri_coords2[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y, 
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt2.y 
                                          ]
                        
                                [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                
                                circumcenter_pt4 = circumcenter_tri(tri_coords4)
                                lOrd = [2,0,1]
                                
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords4[lOrd[0],0] + 1.0/3.0 * tri_coords4[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords4[lOrd[0],1] + 1.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords4[lOrd[0],0] + 2.0/3.0 * tri_coords4[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords4[lOrd[0],1] + 2.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords4[lOrd[2],0] + 1.0/3.0 * tri_coords4[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords4[lOrd[2],1] + 1.0/3.0 * tri_coords4[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords4[lOrd[2],0] + 2.0/3.0 * tri_coords4[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords4[lOrd[2],1] + 2.0/3.0 * tri_coords4[lOrd[0],1]
                        
                                pt6.x = coord_enrich1.x
                                pt6.y = coord_enrich1.y     
                                
                                pt7.x = coord_enrich2.x
                                pt7.y = coord_enrich2.y
                                
                                vec4_x = [tri_coords4[lOrd[0],0], 
                                          tri_coords4[lOrd[1],0],
                                          tri_coords4[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt4.x  
                                          ]
                                vec4_y = [tri_coords4[lOrd[0],1], 
                                          tri_coords4[lOrd[1],1],
                                          tri_coords4[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y, 
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt4.y  
                                          ]
                        
                                [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
        
                                        
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                            detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
                        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)   
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
                            el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
                            
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
                            
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_4( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
                
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_4( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                            polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], pt6, tri_nodes4[2]]]
                            polygonList = polygonList + [[pt6, tri_nodes4[1], tri_nodes4[2]]]
                            
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_4( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_4( nodes7.x, nodes7.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                            polygonList = polygonList + [[pt6, tri_nodes2[1], pt7]]
                            polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], pt6, tri_nodes4[2]]]
                            polygonList = polygonList + [[pt7, pt6, tri_nodes4[2]]]
                            polygonList = polygonList + [[pt7, tri_nodes4[1], tri_nodes4[2]]]
                           
        
                    # the South-East corner is cut, 0-4-1, 1-5-2
                    if (
                        ((enrich1[1] == y0 and enrich2[0] == x1) or 
                         (enrich2[1] == y0 and enrich1[0] == x1)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
                         ):
                        
                        print 'norm computation: SE corner'
                            
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (lower right corner cut)
                        tri_nodes1 = [nodes[0],nodes[4],nodes[3]]
                        tri_nodes2 = [nodes[4],nodes[5],nodes[3]]
                        tri_nodes3 = [nodes[5],nodes[2],nodes[3]]
                        tri_nodes4 = [nodes[4],nodes[1],nodes[5]] # the one triangle in a diff material
        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
        
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
        
                        Pe4 = np.zeros((3,3))
                        Pe4[:,0] = np.ones((3,1)).transpose()
                        Pe4[:,1] = p[tri_nodes4[0:3],0]
                        Pe4[:,2] = p[tri_nodes4[0:3],1]
                        C4 = np.linalg.inv(Pe4)
                        Nbasis_tri4 = tribasisFct(C4)
                        Nx_tri4 = triderivX(C4)
                        Ny_tri4 = triderivY(C4)
        
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
                        tri_coords4 = p[tri_nodes4]
        
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
        
                        e1 = y1 - y5
                        e2 = y5 - y0
                        factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
                        if factor_E > EPS_FACTOR:
                            factor_E = 1
        
                        s1 = x4 - x0
                        s2 = x1 - x4
                        factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
                        if factor_S > EPS_FACTOR:
                            factor_S = 1
        
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_S)
            
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_S +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E ) 
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_E )
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_S +
                                                U[tri_nodes4[2],0] * Nbasis_tri4[2](x,y) * factor_E )
# canceling the norm computation
                        if NORM_COMP == 1:
                            [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                            [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                            J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            if len(root.enrichNodes) < 3:
                                
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
                            
                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )

                            if len(root.enrichNodes) == 3:

                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                
                                lOrd = [2,0,1] # local order    
                                vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                                                        
                                [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                
                                lOrd = [1,2,0]
                                vec4_x = [ tri_coords4[ lOrd[0],0], tri_coords4[lOrd[1],0], tri_coords4[lOrd[2],0], (tri_coords4[ lOrd[0],0] + tri_coords4[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords4[lOrd[0],0] + tri_coords4[lOrd[2],0])/2.0  ]
                                vec4_y = [ tri_coords4[ lOrd[0],1], tri_coords4[lOrd[1],1], tri_coords4[lOrd[2],1], (tri_coords4[ lOrd[0],1] + tri_coords4[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords4[lOrd[0],1] + tri_coords4[lOrd[2],1])/2.0  ]
                        
                                [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
                                
                            if len(root.enrichNodes) == 4:
                                
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                                                        
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                lOrd = [2,0,1] 
                                
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                        
                                pt6.x = coord_enrich1.x
                                pt6.y = coord_enrich1.y     
                                
                                pt7.x = coord_enrich2.x
                                pt7.y = coord_enrich2.y
                                
                                vec2_x = [tri_coords2[lOrd[0],0], 
                                          tri_coords2[lOrd[1],0],
                                          tri_coords2[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt2.x  
                                          ]
                                vec2_y = [tri_coords2[lOrd[0],1], 
                                          tri_coords2[lOrd[1],1], 
                                          tri_coords2[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y, 
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt2.y 
                                          ]
                        
                                [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )

                                
                                circumcenter_pt4 = circumcenter_tri(tri_coords4)
                                
                                lOrd = [1,2,0]
                                        
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords4[lOrd[0],0] + 1.0/3.0 * tri_coords4[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords4[lOrd[0],1] + 1.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords4[lOrd[0],0] + 2.0/3.0 * tri_coords4[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords4[lOrd[0],1] + 2.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords4[lOrd[2],0] + 1.0/3.0 * tri_coords4[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords4[lOrd[2],1] + 1.0/3.0 * tri_coords4[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords4[lOrd[2],0] + 2.0/3.0 * tri_coords4[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords4[lOrd[2],1] + 2.0/3.0 * tri_coords4[lOrd[0],1]
                        
                                pt6.x = coord_enrich2.x
                                pt6.y = coord_enrich2.y     
                                
                                pt7.x = coord_enrich1.x
                                pt7.y = coord_enrich1.y
                                 
                                vec4_x = [tri_coords4[lOrd[0],0], 
                                          tri_coords4[lOrd[1],0],
                                          tri_coords4[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt4.x  
                                          ]
                                vec4_y = [tri_coords4[lOrd[0],1], 
                                          tri_coords4[lOrd[1],1],
                                          tri_coords4[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y, 
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt4.y  
                                          ]
                                 
                                [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
                                                                
                                
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                            detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
                            el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
             
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_4( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
    
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]

                        
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                            polygonList = polygonList + [[pt6, tri_nodes2[1], tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt6, tri_nodes4[1], tri_nodes4[2]]]

                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_2( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_2( nodes7.x, nodes7.y )
                          
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]
                            polygonList = polygonList + [[pt6, pt7, tri_nodes2[2]]]
                            polygonList = polygonList + [[pt7, tri_nodes2[1], tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes4[1], tri_nodes4[2]]]
                            
    
                    # the North East corner is cut, 1-4-2, 2-5-3
                    if (
                        ((enrich1[0] == x1 and enrich2[1] == y1) or
                         (enrich2[0] == x1 and enrich1[1] == y1)) and
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) 
                        ):
                        
                        print "norm computation: NE corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
                        tri_nodes1 = [nodes[0],nodes[5],nodes[3]]
                        tri_nodes2 = [nodes[0],nodes[4],nodes[5]]
                        tri_nodes3 = [nodes[0],nodes[1],nodes[4]]
                        tri_nodes4 = [nodes[4],nodes[2],nodes[5]] # the one triangle in a diff material
                        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
        
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
        
                        Pe4 = np.zeros((3,3))
                        Pe4[:,0] = np.ones((3,1)).transpose()
                        Pe4[:,1] = p[tri_nodes4[0:3],0]
                        Pe4[:,2] = p[tri_nodes4[0:3],1]
                        C4 = np.linalg.inv(Pe4)
                        Nbasis_tri4 = tribasisFct(C4)
                        Nx_tri4 = triderivX(C4)
                        Ny_tri4 = triderivY(C4)
    
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
                        tri_coords4 = p[tri_nodes4]
        
                        Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri1[1]]
                        Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri2[1],Nbasis_tri2[2]]
                        Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri3[2],lambda x,y: 0]
                        Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],Nbasis_tri4[2]]
        
                                
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
        
                        e1 = y1 - y4
                        e2 = y4 - y0
                        factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
                        if factor_E > EPS_FACTOR:
                            factor_E = 1
        
                        n1 = x5 - x0
                        n2 = x1 - x5
                        factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
                        if factor_N > EPS_FACTOR:
                            factor_N = 1
        
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_N ) 

                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_E)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_E +
                                                    U[tri_nodes4[2],0] * Nbasis_tri4[2](x,y) * factor_N )
# canceling the norm computation        
                        if NORM_COMP == 1:
                            [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                            [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )

                            J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                            J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                            
                            if len(root.enrichNodes) < 3:
                                
                                [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                                [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
                

                                J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                                J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
        
                            if len(root.enrichNodes) == 3:
        
#                                coord_enrich = coord_enrich_computation(root.enrichNodes[2])
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                         
                                lOrd = [0,1,2] # local order 
                                vec2_x = [ tri_coords2[ lOrd[0],0], tri_coords2[lOrd[1],0], tri_coords2[lOrd[2],0], (tri_coords2[ lOrd[0],0] + tri_coords2[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords2[lOrd[0],0] + tri_coords2[lOrd[2],0])/2.0  ]
                                vec2_y = [ tri_coords2[ lOrd[0],1], tri_coords2[lOrd[1],1], tri_coords2[lOrd[2],1], (tri_coords2[ lOrd[0],1] + tri_coords2[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords2[lOrd[0],1] + tri_coords2[lOrd[2],1])/2.0  ]
                        
                                [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                                
                                lOrd = [1,2,0]
                                vec4_x = [ tri_coords4[ lOrd[0],0], tri_coords4[lOrd[1],0], tri_coords4[lOrd[2],0], (tri_coords4[ lOrd[0],0] + tri_coords4[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords4[lOrd[0],0] + tri_coords4[lOrd[2],0])/2.0  ]
                                vec4_y = [ tri_coords4[ lOrd[0],1], tri_coords4[lOrd[1],1], tri_coords4[lOrd[2],1], (tri_coords4[ lOrd[0],1] + tri_coords4[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords4[lOrd[0],1] + tri_coords4[lOrd[2],1])/2.0  ]
                        
                                [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
        
                            if len(root.enrichNodes) == 4:
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                                
                                circumcenter_pt2 = circumcenter_tri(tri_coords2)
                                
                                lOrd = [0,1,2] 
                        
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords2[lOrd[0],0] + 1.0/3.0 * tri_coords2[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords2[lOrd[0],1] + 1.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords2[lOrd[0],0] + 2.0/3.0 * tri_coords2[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords2[lOrd[0],1] + 2.0/3.0 * tri_coords2[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords2[lOrd[2],0] + 1.0/3.0 * tri_coords2[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords2[lOrd[2],1] + 1.0/3.0 * tri_coords2[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords2[lOrd[2],0] + 2.0/3.0 * tri_coords2[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords2[lOrd[2],1] + 2.0/3.0 * tri_coords2[lOrd[0],1]
                        
                                pt6.x = coord_enrich1.x
                                pt6.y = coord_enrich1.y     
                                
                                pt7.x = coord_enrich2.x
                                pt7.y = coord_enrich2.y
                                
                                vec2_x = [tri_coords2[lOrd[0],0], 
                                          tri_coords2[lOrd[1],0],
                                          tri_coords2[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt2.x  
                                          ]
                                vec2_y = [tri_coords2[lOrd[0],1], 
                                          tri_coords2[lOrd[1],1], 
                                          tri_coords2[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y, 
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt2.y 
                                          ]
                                                                
                                [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
                                J_tri2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                
                                circumcenter_pt4 = circumcenter_tri(tri_coords4)
                                lOrd = [1,2,0]
                                
                                
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                
                                pt4.x = 2.0/3.0 * tri_coords4[lOrd[0],0] + 1.0/3.0 * tri_coords4[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords4[lOrd[0],1] + 1.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt5.x = 1.0/3.0 * tri_coords4[lOrd[0],0] + 2.0/3.0 * tri_coords4[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords4[lOrd[0],1] + 2.0/3.0 * tri_coords4[lOrd[1],1]
                                
                                pt8.x = 2.0/3.0 * tri_coords4[lOrd[2],0] + 1.0/3.0 * tri_coords4[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords4[lOrd[2],1] + 1.0/3.0 * tri_coords4[lOrd[0],1]
                                
                                pt9.x = 1.0/3.0 * tri_coords4[lOrd[2],0] + 2.0/3.0 * tri_coords4[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords4[lOrd[2],1] + 2.0/3.0 * tri_coords4[lOrd[0],1]
                        
                                pt6.x = coord_enrich2.x
                                pt6.y = coord_enrich2.y     
                                
                                pt7.x = coord_enrich1.x
                                pt7.y = coord_enrich1.y
                                
                                vec4_x = [tri_coords4[lOrd[0],0], 
                                            tri_coords4[lOrd[1],0],
                                            tri_coords4[lOrd[2],0],
                                            pt4.x, 
                                            pt5.x,
                                            pt6.x,
                                            pt7.x,
                                            pt8.x,
                                            pt9.x,
                                            circumcenter_pt4.x  
                                           ]
                                vec4_y = [tri_coords4[lOrd[0],1], 
                                           tri_coords4[lOrd[1],1],
                                           tri_coords4[lOrd[2],1], 
                                           pt4.y, 
                                           pt5.y,
                                           pt6.y, 
                                           pt7.y,
                                           pt8.y,
                                           pt9.y,
                                           circumcenter_pt4.y  
                                           ]
                                        
                                [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
                                J_tri4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
                                        
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                            detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
                            el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
             
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
                            
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_3( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_4( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
    
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_4( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                            polygonList = polygonList + [[tri_nodes2[0], pt6, tri_nodes2[2]]]

                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt6, tri_nodes4[1], tri_nodes4[2]]]
                        
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_4( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_4( nodes7.x, nodes7.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]

                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], pt6]]
                            polygonList = polygonList + [[tri_nodes2[0], pt7, pt6]]                            
                            polygonList = polygonList + [[tri_nodes2[0], pt7, tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes4[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes4[1], tri_nodes4[2]]]
                       
        
                    # the South-West corner is cut, 0-4-1, and 0-5-3
                    if (
                        ((enrich1[1] == y0 and enrich2[0] == x0) or
                         (enrich2[1] == y0 and enrich1[0] == x0)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
                        ):
                        
                        print "norm computation: SW corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
                        tri_nodes1 = [nodes[0],nodes[4],nodes[5]]
                        tri_nodes2 = [nodes[5],nodes[2],nodes[3]]
                        tri_nodes3 = [nodes[4],nodes[2],nodes[5]]
                        tri_nodes4 = [nodes[4],nodes[1],nodes[2]] # the one triangle in a diff material
        
                        Pe1 = np.zeros((3,3))
                        Pe1[:,0] = np.ones((3,1)).transpose()
                        Pe1[:,1] = p[tri_nodes1[0:3],0]
                        Pe1[:,2] = p[tri_nodes1[0:3],1]
                        C1 = np.linalg.inv(Pe1)
                        Nbasis_tri1 = tribasisFct(C1)
                        Nx_tri1 = triderivX(C1)
                        Ny_tri1 = triderivY(C1)
    
                        Pe2 = np.zeros((3,3))
                        Pe2[:,0] = np.ones((3,1)).transpose()
                        Pe2[:,1] = p[tri_nodes2[0:3],0]
                        Pe2[:,2] = p[tri_nodes2[0:3],1]
                        C2 = np.linalg.inv(Pe2)
                        Nbasis_tri2 = tribasisFct(C2)
                        Nx_tri2 = triderivX(C2)
                        Ny_tri2 = triderivY(C2)
        
                        Pe3 = np.zeros((3,3))
                        Pe3[:,0] = np.ones((3,1)).transpose()
                        Pe3[:,1] = p[tri_nodes3[0:3],0]
                        Pe3[:,2] = p[tri_nodes3[0:3],1]
                        C3 = np.linalg.inv(Pe3)
                        Nbasis_tri3 = tribasisFct(C3)
                        Nx_tri3 = triderivX(C3)
                        Ny_tri3 = triderivY(C3)
        
                        Pe4 = np.zeros((3,3))
                        Pe4[:,0] = np.ones((3,1)).transpose()
                        Pe4[:,1] = p[tri_nodes4[0:3],0]
                        Pe4[:,2] = p[tri_nodes4[0:3],1]
                        C4 = np.linalg.inv(Pe4)
                        Nbasis_tri4 = tribasisFct(C4)
                        Nx_tri4 = triderivX(C4)
                        Ny_tri4 = triderivY(C4)
        
                        tri_coords1 = p[tri_nodes1]
                        tri_coords2 = p[tri_nodes2]
                        tri_coords3 = p[tri_nodes3]
                        tri_coords4 = p[tri_nodes4]
        
                        Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri1[1],Nbasis_tri1[2]]
                        Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri2[0]]
                        Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri3[0],Nbasis_tri3[2]]
                        Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],lambda x,y: 0]
        
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
        
                        w1 = y1 - y5
                        w2 = y5 - y0
                        factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
                        if factor_W > EPS_FACTOR:
                            factor_W = 1
        
                        s1 = x4 - x0
                        s2 = x1 - x4
                        factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
                        if factor_S > EPS_FACTOR:
                            factor_S = 1
        
                        uh_elem_1 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_S+
                                                U[tri_nodes1[2],0] * Nbasis_tri1[2](x,y) * factor_W )
        
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_W )

                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_S +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_W )

                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_S )
# canceling the norm computation        
                        if NORM_COMP == 1:
                            
                            [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                            [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
                            
                            J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                            J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
                            
                            if len(root.enrichNodes) < 3:
                                [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                                [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                            
                                J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                                J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                    
                            if len(root.enrichNodes) == 3:

                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                   
                                lOrd = [0,1,2] # local order 
                                vec1_x = [ tri_coords1[ lOrd[0],0], tri_coords1[lOrd[1],0], tri_coords1[lOrd[2],0], (tri_coords1[ lOrd[0],0] + tri_coords1[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords1[lOrd[0],0] + tri_coords1[lOrd[2],0])/2.0  ]
                                vec1_y = [ tri_coords1[ lOrd[0],1], tri_coords1[lOrd[1],1], tri_coords1[lOrd[2],1], (tri_coords1[ lOrd[0],1] + tri_coords1[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords1[lOrd[0],1] + tri_coords1[lOrd[2],1])/2.0  ]
                        
                                [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                                J_tri1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                                
                                lOrd = [1,2,0]
                                vec3_x = [ tri_coords3[ lOrd[0],0], tri_coords3[lOrd[1],0], tri_coords3[lOrd[2],0], (tri_coords3[ lOrd[0],0] + tri_coords3[lOrd[1],0])/2.0, coord_enrich.x, (tri_coords3[lOrd[0],0] + tri_coords3[lOrd[2],0])/2.0  ]
                                vec3_y = [ tri_coords3[ lOrd[0],1], tri_coords3[lOrd[1],1], tri_coords3[lOrd[2],1], (tri_coords3[ lOrd[0],1] + tri_coords3[lOrd[1],1])/2.0, coord_enrich.y, (tri_coords3[lOrd[0],1] + tri_coords3[lOrd[2],1])/2.0  ]
                        
                                [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
                                J_tri3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
                            
                            if len(root.enrichNodes) == 4:
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                                                        
                                circumcenter_pt1 = circumcenter_tri(tri_coords1)
                                lOrd = [0,1,2] # local order 
                                
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                 
                                pt4.x = 2.0/3.0 * tri_coords1[lOrd[0],0] + 1.0/3.0 * tri_coords1[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords1[lOrd[0],1] + 1.0/3.0 * tri_coords1[lOrd[1],1]
                                 
                                pt5.x = 1.0/3.0 * tri_coords1[lOrd[0],0] + 2.0/3.0 * tri_coords1[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords1[lOrd[0],1] + 2.0/3.0 * tri_coords1[lOrd[1],1]
                                 
                                pt8.x = 2.0/3.0 * tri_coords1[lOrd[2],0] + 1.0/3.0 * tri_coords1[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords1[lOrd[2],1] + 1.0/3.0 * tri_coords1[lOrd[0],1]
                                 
                                pt9.x = 1.0/3.0 * tri_coords1[lOrd[2],0] + 2.0/3.0 * tri_coords1[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords1[lOrd[2],1] + 2.0/3.0 * tri_coords1[lOrd[0],1]
                         
                                pt6.x = coord_enrich1.x
                                pt6.y = coord_enrich1.y     
                                 
                                pt7.x = coord_enrich2.x
                                pt7.y = coord_enrich2.y
                                 
                                vec1_x = [tri_coords1[lOrd[0],0], 
                                          tri_coords1[lOrd[1],0], 
                                          tri_coords1[lOrd[2],0], 
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt1.x  
                                          ]
                                vec1_y = [tri_coords1[lOrd[0],1], 
                                          tri_coords1[lOrd[1],1], 
                                          tri_coords1[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y,
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt1.y  
                                          ]
                        
                        
                                [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
                                J_tri1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                                
                                circumcenter_pt3 = circumcenter_tri(tri_coords3)
                                lOrd = [1,2,0]
                        
                                pt4 = Coordinate(0,0)
                                pt5 = Coordinate(0,0)
                                pt6 = Coordinate(0,0)
                                pt7 = Coordinate(0,0)
                                pt8 = Coordinate(0,0)
                                pt9 = Coordinate(0,0)
                                 
                                pt4.x = 2.0/3.0 * tri_coords3[lOrd[0],0] + 1.0/3.0 * tri_coords3[lOrd[1],0]
                                pt4.y = 2.0/3.0 * tri_coords3[lOrd[0],1] + 1.0/3.0 * tri_coords3[lOrd[1],1]
                                 
                                pt5.x = 1.0/3.0 * tri_coords3[lOrd[0],0] + 2.0/3.0 * tri_coords3[lOrd[1],0]
                                pt5.y = 1.0/3.0 * tri_coords3[lOrd[0],1] + 2.0/3.0 * tri_coords3[lOrd[1],1]
                                 
                                pt8.x = 2.0/3.0 * tri_coords3[lOrd[2],0] + 1.0/3.0 * tri_coords3[lOrd[0],0]
                                pt8.y = 2.0/3.0 * tri_coords3[lOrd[2],1] + 1.0/3.0 * tri_coords3[lOrd[0],1]
                                 
                                pt9.x = 1.0/3.0 * tri_coords3[lOrd[2],0] + 2.0/3.0 * tri_coords3[lOrd[0],0]
                                pt9.y = 1.0/3.0 * tri_coords3[lOrd[2],1] + 2.0/3.0 * tri_coords3[lOrd[0],1]
                         
                                pt6.x = coord_enrich2.x
                                pt6.y = coord_enrich2.y     
                                 
                                pt7.x = coord_enrich1.x
                                pt7.y = coord_enrich1.y
                                
                                vec3_x = [tri_coords3[lOrd[0],0], 
                                          tri_coords3[lOrd[1],0],
                                          tri_coords3[lOrd[2],0],
                                          pt4.x, 
                                          pt5.x,
                                          pt6.x,
                                          pt7.x,
                                          pt8.x,
                                          pt9.x,
                                          circumcenter_pt3.x  
                                          ]
                                vec3_y = [tri_coords3[lOrd[0],1], 
                                          tri_coords3[lOrd[1],1], 
                                          tri_coords3[lOrd[2],1], 
                                          pt4.y, 
                                          pt5.y,
                                          pt6.y,
                                          pt7.y,
                                          pt8.y,
                                          pt9.y,
                                          circumcenter_pt3.y  
                                          ]
                                
                                [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
                                J_tri3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
                                
                                
                            detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                            detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                            detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                            detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
                            
                        
                            el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                            el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                            el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)                            
                            el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
         
                            all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
                        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_4( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_4( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_2( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_1( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_1( p[nodes[5],0], p[nodes[5],1]  )
                        
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_1( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                            polygonList = polygonList + [[tri_nodes1[0], pt6, tri_nodes1[2]]]
                            
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]


                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], pt6]]
                            polygonList = polygonList + [[pt6, tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
                         
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_1( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_1( nodes7.x, nodes7.y )
                            
                            polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], pt6]]
                            polygonList = polygonList + [[tri_nodes1[0], pt6, pt7]]
                            polygonList = polygonList + [[tri_nodes1[0], pt7, tri_nodes1[2]]]
                            
                            polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]

                            polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes3[1], pt6]]
                            polygonList = polygonList + [[pt7, tri_nodes3[1], tri_nodes3[2]]]
                            
                            polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
                           
        
                    # interface cuts the element horizontally into two quads, 0-4-3, 1-5-2 
                    if ( ((enrich1[0] == x0  and enrich2[0] == x1) or 
                        (enrich1[0] == x1  and enrich2[0] == x0)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):
        
                        # nodes on the top and bottom side of the interface
                        top_nodes = [nodes[4], nodes[5], nodes[2],nodes[3]]
                        bottom_nodes = [nodes[0],nodes[1],nodes[5],nodes[4]]
    
                        if (enrich1[0] == x1  and enrich2[0] == x0):
                            top_nodes = [nodes[5], nodes[4], nodes[2],nodes[3]]
                            bottom_nodes = [nodes[0],nodes[1],nodes[4],nodes[5]]
    
        
                        top_coords = p[top_nodes,:]
                        bottom_coords = p[bottom_nodes,:]
    
                        # build the shape functions at the enrichment nodes
                        Pe_enr_top = np.zeros((4,4))
                        Pe_enr_top[:,0] = np.ones((4,1)).transpose()
                        Pe_enr_top[:,1] = p[top_nodes[0:4],0]
                        Pe_enr_top[:,2] = p[top_nodes[0:4],1]
                        Pe_enr_top[:,3] = p[top_nodes[0:4],0]*p[top_nodes[0:4],1]
                        C_enr_top = np.linalg.inv(Pe_enr_top)
        
                        # left enrichment shape function and its derivative wrt x & y
                        Nbasis_enr_top = basisFct(C_enr_top)
        
                        Pe_enr_bottom = np.zeros((4,4))
                        Pe_enr_bottom[:,0] = np.ones((4,1)).transpose()
                        Pe_enr_bottom[:,1] = p[bottom_nodes[0:4],0]
                        Pe_enr_bottom[:,2] = p[bottom_nodes[0:4],1]
                        Pe_enr_bottom[:,3] = p[bottom_nodes[0:4],0]*p[bottom_nodes[0:4],1]
                        C_enr_bottom = np.linalg.inv(Pe_enr_bottom)
        
                        # bottom enrichment shape function and its derivatives wrt x & y
                        Nbasis_enr_bottom = basisFct(C_enr_bottom)
        
                        # shape functions at enrichment nodes
                        psi_left_B = lambda x,y: Nbasis_enr_bottom[3](x,y)
                        psi_left_T = lambda x,y: Nbasis_enr_top[0](x,y)
                        psi_right_B = lambda x,y: Nbasis_enr_bottom[2](x,y)
                        psi_right_T = lambda x,y: Nbasis_enr_top[1](x,y)
        
                        [x_fct_T, y_fct_T] = xy_fct( top_coords[:,0], top_coords[:,1] )
                        [x_fct_B, y_fct_B] = xy_fct( bottom_coords[:,0], bottom_coords[:,1] )
        
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
        
                        e1 = y1 - y5
                        e2 = y5 - y0
                        factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
                        if factor_E > EPS_FACTOR:
                            factor_E = 1
        
                        w1 = y1 - y4
                        w2 = y4 - y0
                        factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
                        if factor_W > EPS_FACTOR:
                                factor_W = 1
        
                        # on the top side of the interface
                        uh_elem_T = lambda x,y: (
                                                U[top_nodes[0],0] * psi_left_T(x,y) * factor_W +
                                                U[top_nodes[1],0] * psi_right_T(x,y) * factor_E +
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) )
         
                        top_coords = p[top_nodes,:]
                        x_coords_T = top_coords[:,0]
                        y_coords_T = top_coords[:,1]

                        # on the bottom of the interface
                        uh_elem_B = lambda x,y: (
                                                U[bottom_nodes[3],0] * psi_left_B(x,y) * factor_W +
                                                U[bottom_nodes[2],0] * psi_right_B(x,y) * factor_E +
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) )
        
        
                        bottom_coords = p[bottom_nodes,:]
                        x_coords_B = bottom_coords[:,0]
                        y_coords_B = bottom_coords[:,1]
                        
                        
# canceling the norm computation            
                        if NORM_COMP == 1:

                            if len(root.enrichNodes) < 3:
                                # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
                                # for transformation from the parametric element to phisycal element
                                # of the Gauss nodes ui
                                [x_transform_fct_T,y_transform_fct_T] = xy_fct(x_coords_T,y_coords_T)            
                                [x_transform_fct_B,y_transform_fct_B] = xy_fct(x_coords_B,y_coords_B)            
                                # computing the Jacobian and the determinant of the left and right children of the parent element
                                J_top = jacobian_mat( top_coords[:,0], top_coords[:,1] )
                                J_bottom = jacobian_mat( bottom_coords[:,0], bottom_coords[:,1] )
                              
                            if len(root.enrichNodes) == 3:
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                
                                lOrd = [0,1,2,3] # local order    
                                vecB_x = [ bottom_coords[ lOrd[0],0], 
                                          (bottom_coords[ lOrd[0],0] + bottom_coords[lOrd[1],0])/2.0,
                                          bottom_coords[lOrd[1],0],
                                          (bottom_coords[lOrd[1],0] + bottom_coords[lOrd[2],0])/2.0,
                                          bottom_coords[lOrd[2],0], 
                                          coord_enrich.x, 
                                          bottom_coords[lOrd[3],0],
                                          (bottom_coords[lOrd[3],0] + bottom_coords[lOrd[0],0])/2.0
                                         ]
                                vecB_y = [ bottom_coords[ lOrd[0],1], 
                                          (bottom_coords[ lOrd[0],1] + bottom_coords[lOrd[1],1])/2.0,
                                          bottom_coords[lOrd[1],1],
                                          (bottom_coords[lOrd[1],1] + bottom_coords[lOrd[2],1])/2.0,
                                          bottom_coords[lOrd[2],1], 
                                          coord_enrich.y, 
                                          bottom_coords[lOrd[3],1],
                                          (bottom_coords[lOrd[3],1] + bottom_coords[lOrd[0],1])/2.0
                                         ]
                        
                                [x_transform_fct_B, y_transform_fct_B] = quad_xy_fct_bi_quadratic( vecB_x, vecB_y )
                                J_bottom = quad_jacobian_mat_bi_quadratic( vecB_x, vecB_y )
                                
                                lOrd = [2,3,0,1]
                                vecT_x = [ top_coords[ lOrd[0],0], 
                                          (top_coords[ lOrd[0],0] + top_coords[lOrd[1],0])/2.0,
                                          top_coords[lOrd[1],0],
                                          (top_coords[lOrd[1],0] + top_coords[lOrd[2],0])/2.0,
                                          top_coords[lOrd[2],0], 
                                          coord_enrich.x, 
                                          top_coords[lOrd[3],0],
                                          (top_coords[lOrd[3],0] + top_coords[lOrd[0],0])/2.0
                                         ]
                                vecT_y = [ top_coords[ lOrd[0],1], 
                                          (top_coords[ lOrd[0],1] + top_coords[lOrd[1],1])/2.0,
                                          top_coords[lOrd[1],1],
                                          (top_coords[lOrd[1],1] + top_coords[lOrd[2],1])/2.0,
                                          top_coords[lOrd[2],1], 
                                          coord_enrich.y, 
                                          top_coords[lOrd[3],1],
                                          (top_coords[lOrd[3],1] + top_coords[lOrd[0],1])/2.0
                                         ]
                        
                                [x_transform_fct_T, y_transform_fct_T] = quad_xy_fct_bi_quadratic( vecT_x, vecT_y )
                                J_top = quad_jacobian_mat_bi_quadratic( vecT_x, vecT_y )
                            
                            if len(root.enrichNodes) == 4:
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                        
                                lOrd = [0,1,2,3] # local order    
                                
                                vecB_x = [ bottom_coords[ lOrd[0],0], 
                                          2.0/3.0 * bottom_coords[lOrd[0],0] + 1.0/3.0 * bottom_coords[lOrd[1],0],
                                          1.0/3.0 * bottom_coords[lOrd[0],0] + 2.0/3.0 * bottom_coords[lOrd[1],0],
                                          bottom_coords[lOrd[1],0],
                                          2.0/3.0 * bottom_coords[lOrd[1],0] + 1.0/3.0 * bottom_coords[lOrd[2],0],
                                          1.0/3.0 * bottom_coords[lOrd[1],0] + 2.0/3.0 * bottom_coords[lOrd[2],0],
                                          bottom_coords[lOrd[2],0], 
                                          coord_enrich2.x,
                                          coord_enrich1.x,
                                          bottom_coords[lOrd[3],0],
                                          2.0/3.0 * bottom_coords[lOrd[3],0] + 1.0/3.0 * bottom_coords[lOrd[0],0],
                                          1.0/3.0 * bottom_coords[lOrd[3],0] + 2.0/3.0 * bottom_coords[lOrd[0],0]
                                         ]
                                vecB_y = [ bottom_coords[ lOrd[0],1], 
                                          2.0/3.0 * bottom_coords[lOrd[0],1] + 1.0/3.0 * bottom_coords[lOrd[1],1],
                                          1.0/3.0 * bottom_coords[lOrd[0],1] + 2.0/3.0 * bottom_coords[lOrd[1],1],
                                          bottom_coords[lOrd[1],1],
                                          2.0/3.0 * bottom_coords[lOrd[1],1] + 1.0/3.0 * bottom_coords[lOrd[2],1],
                                          1.0/3.0 * bottom_coords[lOrd[1],1] + 2.0/3.0 * bottom_coords[lOrd[2],1],
                                          bottom_coords[lOrd[2],1], 
                                          coord_enrich2.y,
                                          coord_enrich1.y,
                                          bottom_coords[lOrd[3],1],
                                          2.0/3.0 * bottom_coords[lOrd[3],1] + 1.0/3.0 * bottom_coords[lOrd[0],1],
                                          1.0/3.0 * bottom_coords[lOrd[3],1] + 2.0/3.0 * bottom_coords[lOrd[0],1]
                                         ]
                        
                                [x_transform_fct_B, y_transform_fct_B] = quad_xy_fct_bi_cubic( vecB_x, vecB_y )
                                J_bottom = quad_jacobian_mat_bi_cubic( vecB_x, vecB_y )
                                
                                lOrd = [2,3,0,1]
                                vecT_x = [ top_coords[ lOrd[0],0], 
                                          2.0/3.0 * top_coords[lOrd[0],0] + 1.0/3.0 * top_coords[lOrd[1],0],
                                          1.0/3.0 * top_coords[lOrd[0],0] + 2.0/3.0 * top_coords[lOrd[1],0],
                                          top_coords[lOrd[1],0],
                                          2.0/3.0 * top_coords[lOrd[1],0] + 1.0/3.0 * top_coords[lOrd[2],0],
                                          1.0/3.0 * top_coords[lOrd[1],0] + 2.0/3.0 * top_coords[lOrd[2],0],
                                          top_coords[lOrd[2],0], 
                                          coord_enrich1.x,
                                          coord_enrich2.x,
                                          top_coords[lOrd[3],0],
                                          2.0/3.0 * top_coords[lOrd[3],0] + 1.0/3.0 * top_coords[lOrd[0],0],
                                          1.0/3.0 * top_coords[lOrd[3],0] + 2.0/3.0 * top_coords[lOrd[0],0]
                                         ]
                                vecT_y = [ top_coords[ lOrd[0],1], 
                                          2.0/3.0 * top_coords[lOrd[0],1] + 1.0/3.0 * top_coords[lOrd[1],1],
                                          1.0/3.0 * top_coords[lOrd[0],1] + 2.0/3.0 * top_coords[lOrd[1],1],
                                          top_coords[lOrd[1],1],
                                          2.0/3.0 * top_coords[lOrd[1],1] + 1.0/3.0 * top_coords[lOrd[2],1],
                                          1.0/3.0 * top_coords[lOrd[1],1] + 2.0/3.0 * top_coords[lOrd[2],1],
                                          top_coords[lOrd[2],1], 
                                          coord_enrich1.y,
                                          coord_enrich2.y,
                                          top_coords[lOrd[3],1],
                                          2.0/3.0 * top_coords[lOrd[3],1] + 1.0/3.0 * top_coords[lOrd[0],1],
                                          1.0/3.0 * top_coords[lOrd[3],1] + 2.0/3.0 * top_coords[lOrd[0],1]
                                        ]
                        
                                [x_transform_fct_T, y_transform_fct_T] = quad_xy_fct_bi_quadratic( vecT_x, vecT_y )
                                J_top = quad_jacobian_mat_bi_quadratic( vecT_x, vecT_y )        
                                
                                
                            
                            detJ_top = lambda e,n: determinant(J_top)(e,n)
                            detJ_bottom = lambda e,n: determinant(J_bottom)(e,n)   
                            
                            el_sum_T =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct_T,y_transform_fct_T,uh_elem_T,detJ_top)
                            el_sum_B =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct_B,y_transform_fct_B,uh_elem_B,detJ_bottom)
             
                            all_elems_sum = all_elems_sum + el_sum_T + el_sum_B
                            
                        Usolution[nodes[0],0] = uh_elem_B( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_B( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_T( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_T( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_B( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_B( p[nodes[5],0], p[nodes[5],1]  )
    
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[top_nodes[0],top_nodes[1],top_nodes[2],top_nodes[3]]]
                            polygonList = polygonList + [[bottom_nodes[0],bottom_nodes[1],bottom_nodes[2],bottom_nodes[3]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_B( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[top_nodes[0],pt6,top_nodes[3]]]
                            polygonList = polygonList + [[pt6,top_nodes[2],top_nodes[3]]]
                            polygonList = polygonList + [[pt6,top_nodes[1],top_nodes[2]]]
                            
                            polygonList = polygonList + [[bottom_nodes[0],pt6,bottom_nodes[3]]]
                            polygonList = polygonList + [[bottom_nodes[0],bottom_nodes[1],pt6]]
                            polygonList = polygonList + [[pt6,bottom_nodes[1],bottom_nodes[2]]]
    
    
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_B( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_B( nodes7.x, nodes7.y )
    
                            polygonList = polygonList + [[top_nodes[0],pt6,top_nodes[3]]]
                            polygonList = polygonList + [[pt6, pt7, top_nodes[2], top_nodes[3]]]
                            polygonList = polygonList + [[pt7,top_nodes[1],top_nodes[2]]]
                            
                            polygonList = polygonList + [[bottom_nodes[0],pt6,bottom_nodes[3]]]
                            polygonList = polygonList + [[bottom_nodes[0], bottom_nodes[1], pt7, pt6]]
                            polygonList = polygonList + [[pt7,bottom_nodes[1],bottom_nodes[2]]]
                                    
                    # interface cuts the element vertically into two quads, 0-4-1, 3-5-2
                    if ( ((enrich1[1] == y0 and enrich2[1] == y1) or 
                        (enrich1[1] == y1 and enrich2[1] == y0)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):
                        print "norm computation: vertical slide: quad-quad"
                        # nodes on the left side of the interface
                        left_nodes = [nodes[0],nodes[4],nodes[5],nodes[3]]
                        # nodes on the right side of the interface
                        right_nodes = [nodes[4],nodes[1],nodes[2],nodes[5]]
        
                        if (enrich1[1] == y1 and enrich2[1] == y0):
                            left_nodes = [nodes[0],nodes[5],nodes[4],nodes[3]]
                            right_nodes = [nodes[5],nodes[1],nodes[2],nodes[4]]
        
                        # coordinates of the left sub-element or sub-element number 1
                        left_coords = p[left_nodes,:]
                        # coordinates of the right sub-element or the sub-element number 2
                        right_coords = p[right_nodes,:]
        
                        # build the shape functions at the enrichment nodes
                        Pe_enr_left = np.zeros((4,4))
                        Pe_enr_left[:,0] = np.ones((4,1)).transpose()
                        Pe_enr_left[:,1] = p[left_nodes[0:4],0]
                        Pe_enr_left[:,2] = p[left_nodes[0:4],1]
                        Pe_enr_left[:,3] = p[left_nodes[0:4],0]*p[left_nodes[0:4],1]
        
                        C_enr_left = np.linalg.inv(Pe_enr_left)
        
                        # left enrichment shape function
                        Nbasis_enr_left = basisFct(C_enr_left)
                    
                        Pe_enr_right = np.zeros((4,4))
                        Pe_enr_right[:,0] = np.ones((4,1)).transpose()
                        Pe_enr_right[:,1] = p[right_nodes[0:4],0]
                        Pe_enr_right[:,2] = p[right_nodes[0:4],1]
                        Pe_enr_right[:,3] = p[right_nodes[0:4],0]*p[right_nodes[0:4],1]
                 
                        C_enr_right = np.linalg.inv(Pe_enr_right)
                      
                        # right enrichment shape function
                        Nbasis_enr_right = basisFct(C_enr_right)
             
                        # shape functions at enrichment nodes
                        psi_btm_L = lambda x,y: Nbasis_enr_left[1](x,y)
                        psi_upr_L = lambda x,y: Nbasis_enr_left[2](x,y)
                        psi_btm_R = lambda x,y: Nbasis_enr_right[0](x,y)
                        psi_upr_R = lambda x,y: Nbasis_enr_right[3](x,y)
                      
                        # scaling factor
                        s_factor = 1
                    
                        node_enr_4 = [nodes[4]]
                        node_enr_5 = [nodes[5]]
                        coords_enr_4 = p[node_enr_4]
                        coords_enr_5 = p[node_enr_5]
        
                        x4 = coords_enr_4[0,0]
                        y4 = coords_enr_4[0,1]
        
                        x5 = coords_enr_5[0,0]
                        y5 = coords_enr_5[0,1]
        
                        s1 = x4 - x0
                        s2 = x1 - x4
                        factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
                        if factor_S > EPS_FACTOR:
                            factor_S = 1
        
                        n1 = x5 - x0
                        n2 = x1 - x5
                        factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
                        if factor_N > EPS_FACTOR:
                            factor_N = 1
        
                        # on the left side of the interface
                        uh_elem_L = lambda x,y: (
                                                U[left_nodes[1],0] * psi_btm_L(x,y) * factor_S +
                                                U[left_nodes[2],0] * psi_upr_L(x,y) * factor_N +
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) )    
         
                        left_coords = p[left_nodes,:]
                        x_coords_L = left_coords[:,0]
                        y_coords_L = left_coords[:,1]
                                        
                        # on the right side of the interface
                        uh_elem_R = lambda x,y: (
                                                U[right_nodes[0],0] * psi_btm_R(x,y) * factor_S +
                                                U[right_nodes[3],0] * psi_upr_R(x,y) * factor_N +
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) )
        
                        right_coords = p[right_nodes,:]
                        x_coords_R = right_coords[:,0]
                        y_coords_R = right_coords[:,1]

# canceling the norm computation            
                        if NORM_COMP == 1:
                            if len(root.enrichNodes) < 3:
                                # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
                                # for transformation from the parametric element to phisycal element
                                # of the Gauss nodes ui
                                [x_transform_fct_L,y_transform_fct_L] = xy_fct(left_coords[:,0], left_coords[:,1])            
                                [x_transform_fct_R,y_transform_fct_R] = xy_fct(right_coords[:,0], right_coords[:,1])            
         
                                # computing the Jacobian and the determinant of the left and right children of the parent element
                                J_left = jacobian_mat( left_coords[:,0], left_coords[:,1] )
                                J_right = jacobian_mat( right_coords[:,0], right_coords[:,1] )
                                                    
                            if len(root.enrichNodes) == 3:
                                
                                coord_enrich = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                
                                lOrd = [3,0,1,2] # local order    
                                vecL_x = [ left_coords[ lOrd[0],0], 
                                          (left_coords[ lOrd[0],0] + left_coords[lOrd[1],0])/2.0,
                                          left_coords[lOrd[1],0],
                                          (left_coords[lOrd[1],0] + left_coords[lOrd[2],0])/2.0,
                                          left_coords[lOrd[2],0], 
                                          coord_enrich.x, 
                                          left_coords[lOrd[3],0],
                                          (left_coords[lOrd[3],0] + left_coords[lOrd[0],0])/2.0
                                         ]
                                vecL_y = [ left_coords[ lOrd[0],1], 
                                          (left_coords[ lOrd[0],1] + left_coords[lOrd[1],1])/2.0,
                                          left_coords[lOrd[1],1],
                                          (left_coords[lOrd[1],1] + left_coords[lOrd[2],1])/2.0,
                                          left_coords[lOrd[2],1], 
                                          coord_enrich.y, 
                                          left_coords[lOrd[3],1],
                                          (left_coords[lOrd[3],1] + left_coords[lOrd[0],1])/2.0
                                         ]
                        
                                
                                [x_transform_fct_L, y_transform_fct_L] = quad_xy_fct_bi_quadratic( vecL_x, vecL_y )
                                J_left = quad_jacobian_mat_bi_quadratic( vecL_x, vecL_y )
                                
                                lOrd = [1,2,3,0]
                                vecR_x = [ right_coords[ lOrd[0],0], 
                                          (right_coords[ lOrd[0],0] + right_coords[lOrd[1],0])/2.0,
                                          right_coords[lOrd[1],0],
                                          (right_coords[lOrd[1],0] + right_coords[lOrd[2],0])/2.0,
                                          right_coords[lOrd[2],0], 
                                          coord_enrich.x, 
                                          right_coords[lOrd[3],0],
                                          (right_coords[lOrd[3],0] + right_coords[lOrd[0],0])/2.0
                                         ]
                                vecR_y = [ right_coords[ lOrd[0],1], 
                                          (right_coords[ lOrd[0],1] + right_coords[lOrd[1],1])/2.0,
                                          right_coords[lOrd[1],1],
                                          (right_coords[lOrd[1],1] + right_coords[lOrd[2],1])/2.0,
                                          right_coords[lOrd[2],1], 
                                          coord_enrich.y, 
                                          right_coords[lOrd[3],1],
                                          (right_coords[lOrd[3],1] + right_coords[lOrd[0],1])/2.0
                                         ]
                                
                                [x_transform_fct_R, y_transform_fct_R] = quad_xy_fct_bi_quadratic( vecR_x, vecR_y )
                                J_right = quad_jacobian_mat_bi_quadratic( vecR_x, vecR_y )

                            if len(root.enrichNodes) == 4:
                                
                                coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
                                coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
                                
                                
                                lOrd = [3,0,1,2] # local order    
                                vecL_x = [ left_coords[lOrd[0],0], 
                                          2.0/3.0 * left_coords[lOrd[0],0] + 1.0/3.0 * left_coords[lOrd[1],0],
                                          1.0/3.0 * left_coords[lOrd[0],0] + 2.0/3.0 * left_coords[lOrd[1],0],
                                          left_coords[lOrd[1],0],
                                          2.0/3.0 * left_coords[lOrd[1],0] + 1.0/3.0 * left_coords[lOrd[2],0],
                                          1.0/3.0 * left_coords[lOrd[1],0] + 2.0/3.0 * left_coords[lOrd[2],0],
                                          left_coords[lOrd[2],0], 
                                          coord_enrich1.x,
                                          coord_enrich2.x,
                                          left_coords[lOrd[3],0],
                                          2.0/3.0 * left_coords[lOrd[3],0] + 1.0/3.0 * left_coords[lOrd[0],0],
                                          1.0/3.0 * left_coords[lOrd[3],0] + 2.0/3.0 * left_coords[lOrd[0],0]
                                         ]
                                vecL_y = [ left_coords[lOrd[0],1], 
                                          2.0/3.0 * left_coords[lOrd[0],1] + 1.0/3.0 * left_coords[lOrd[1],1],
                                          1.0/3.0 * left_coords[lOrd[0],1] + 2.0/3.0 * left_coords[lOrd[1],1],
                                          left_coords[lOrd[1],1],
                                          2.0/3.0 * left_coords[lOrd[1],1] + 1.0/3.0 * left_coords[lOrd[2],1],
                                          1.0/3.0 * left_coords[lOrd[1],1] + 2.0/3.0 * left_coords[lOrd[2],1],
                                          left_coords[lOrd[2],1], 
                                          coord_enrich1.y,
                                          coord_enrich2.y,
                                          left_coords[lOrd[3],1],
                                          2.0/3.0 * left_coords[lOrd[3],1] + 1.0/3.0 * left_coords[lOrd[0],1],
                                          1.0/3.0 * left_coords[lOrd[3],1] + 2.0/3.0 * left_coords[lOrd[0],1]
                                        ]
                        
                                [x_transform_fct_L, y_transform_fct_L] = quad_xy_fct_bi_cubic( vecL_x, vecL_y )
                                J_left = quad_jacobian_mat_bi_cubic( vecL_x, vecL_y )
                                
                                lOrd = [1,2,3,0]
                        
                                vecR_x = [ right_coords[ lOrd[0],0], 
                                          2.0/3.0 * right_coords[lOrd[0],0] + 1.0/3.0 * right_coords[lOrd[1],0],
                                          1.0/3.0 * right_coords[lOrd[0],0] + 2.0/3.0 * right_coords[lOrd[1],0],
                                          right_coords[lOrd[1],0],
                                          2.0/3.0 * right_coords[lOrd[1],0] + 1.0/3.0 * right_coords[lOrd[2],0],
                                          1.0/3.0 * right_coords[lOrd[1],0] + 2.0/3.0 * right_coords[lOrd[2],0],
                                          right_coords[lOrd[2],0], 
                                          coord_enrich2.x,
                                          coord_enrich1.x,
                                          right_coords[lOrd[3],0],
                                          2.0/3.0 * right_coords[lOrd[3],0] + 1.0/3.0 * right_coords[lOrd[0],0],
                                          1.0/3.0 * right_coords[lOrd[3],0] + 2.0/3.0 * right_coords[lOrd[0],0]
                                        ]
                                vecR_y = [ right_coords[ lOrd[0],1], 
                                          2.0/3.0 * right_coords[lOrd[0],1] + 1.0/3.0 * right_coords[lOrd[1],1],
                                          1.0/3.0 * right_coords[lOrd[0],1] + 2.0/3.0 * right_coords[lOrd[1],1],
                                          right_coords[lOrd[1],1],
                                          2.0/3.0 * right_coords[lOrd[1],1] + 1.0/3.0 * right_coords[lOrd[2],1],
                                          1.0/3.0 * right_coords[lOrd[1],1] + 2.0/3.0 * right_coords[lOrd[2],1],
                                          right_coords[lOrd[2],1], 
                                          coord_enrich2.y,
                                          coord_enrich1.y,
                                          right_coords[lOrd[3],1],
                                          2.0/3.0 * right_coords[lOrd[3],1] + 1.0/3.0 * right_coords[lOrd[0],1],
                                          1.0/3.0 * right_coords[lOrd[3],1] + 2.0/3.0 * right_coords[lOrd[0],1]
                                          ]
                                
                                [x_transform_fct_R, y_transform_fct_R] = quad_xy_fct_bi_cubic( vecR_x, vecR_y )
                                J_right = quad_jacobian_mat_bi_cubic( vecR_x, vecR_y )

                                
                                
                            detJ_left = lambda eps,niu: determinant(J_left)(eps,niu)
                            detJ_right = lambda eps,niu: determinant(J_right)(eps,niu)
                        
                            el_sum_L =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct_L,y_transform_fct_L,uh_elem_L,detJ_left)
                            el_sum_R =  gauss_integration_quad(ui,wi,UConf,pConf,tConf,x_transform_fct_R,y_transform_fct_R,uh_elem_R,detJ_right)
    
                            all_elems_sum = all_elems_sum + el_sum_R + el_sum_L;
                            
                        Usolution[nodes[0],0] = uh_elem_L( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_R( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_R( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_L( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_L( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_L( p[nodes[5],0], p[nodes[5],1]  )
    
                        if len(root.enrichNodes) < 3:
                            polygonList = polygonList + [[left_nodes[0],left_nodes[1],left_nodes[2],left_nodes[3]]]
                            polygonList = polygonList + [[right_nodes[0],right_nodes[1],right_nodes[2],right_nodes[3]]]
    
    
                        if len(root.enrichNodes) == 3:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_L( nodes6.x, nodes6.y )
                            
                            polygonList = polygonList + [[left_nodes[0],left_nodes[1],pt6]]
                            polygonList = polygonList + [[left_nodes[0],pt6,left_nodes[3]]]
                            polygonList = polygonList + [[pt6,left_nodes[2],left_nodes[3]]]
                                                                            
                            polygonList = polygonList + [[right_nodes[0],right_nodes[1],pt6]]
                            polygonList = polygonList + [[right_nodes[1],pt6,right_nodes[2]]]
                            polygonList = polygonList + [[pt6,right_nodes[2],right_nodes[3]]]
    
                        if len(root.enrichNodes) == 4:
                            index6 = numpy.where(numpy.all(p_extra==[old_node6.x, old_node6.y],axis=1))
                            pt6 = index6[0][0] + len(p)
                            Usolution[pt6] = uh_elem_L( nodes6.x, nodes6.y )
                            
                            index7 = numpy.where(numpy.all(p_extra==[old_node7.x, old_node7.y],axis=1))
                            pt7 = index7[0][0] + len(p)
                            Usolution[pt7] = uh_elem_L( nodes7.x, nodes7.y )
                            
                            polygonList = polygonList + [[left_nodes[0],left_nodes[1],pt6]]
                            polygonList = polygonList + [[left_nodes[0],pt6,pt7,left_nodes[3]]]
                            polygonList = polygonList + [[pt7,left_nodes[2],left_nodes[3]]]
                                                                            
                            polygonList = polygonList + [[right_nodes[0],right_nodes[1],pt6]]
                            polygonList = polygonList + [[pt7,right_nodes[2],right_nodes[3]]]
                            polygonList = polygonList + [[pt6,right_nodes[1],right_nodes[2],pt7]]
    
    
    print all_elems_sum
    print 'norm is: ', math.sqrt(all_elems_sum)

    print 'Writing VTK file...' 
    print_vtk_file(p,Usolution,polygonList,p_extra)
    print ' Done.'

    return  math.sqrt(all_elems_sum)

def print_vtk_file(p,Usolution,plist,p_extra):
    
    P = len(p)
    filename = 'dataset' + str(P+len(p_extra)) + 'points.vtk'
    target = open(filename,'w')
    target.write('# vtk DataFile Version 3.1 \n')
    target.write('Circle example \n')
    target.write('ASCII \n')
    target.write('DATASET POLYDATA \n')
    str1 = 'POINTS ' +  str(P+len(p_extra)) + ' FLOAT \n'
    target.write(str1)

    for i in range(0,P):
        stri = str(p[i,0]) + '  ' + str(p[i,1]) + '  ' + str(Usolution[i,0]) + ' \n'
        target.write(stri)
    
    for i in range(0,len(p_extra)):
        stri = str(p_extra[i,0]) + '  ' + str(p_extra[i,1]) + '  ' + str(Usolution[i+P,0]) + ' \n'
        target.write(stri)
    
    NPolyg = len( sum (plist, [] ) ) + len(plist)
    str2 = '\nPOLYGONS  ' + str(len(plist)) + '   ' + str(NPolyg) + ' \n'
    target.write(str2)
    
    strk = ''
    for j in range( 0, len(plist) ):
        for k in range( 0, len(plist[j]) ):
            strk = strk + str(plist[j][k]) + '  ' 
        target.write( str(len(plist[j])) + '   ' + strk + ' \n')
        strk = ''

    str3 = '\nPOINT_DATA ' + str(P+len(p_extra)) + ' \n'
    target.write(str3)
    target.write('SCALARS Temperature FLOAT \n')
    target.write('LOOKUP_TABLE default \n')
    for z in range(0,len(Usolution)):
        strz = str(Usolution[z,0])
        target.write(strz + ' \n')

    target.close()

def coord_enrich_comp_quad_circle(p,nodes6):
    
    node4 = Coordinate(p[4,0], p[4,1])
    node5 = Coordinate(p[5,0], p[5,1])
    midPoint = Coordinate((node4.x + node5.x)/2.0, (node4.y + node5.y)/2.0)
    
    if node4.x != node5.x:
         xx = [node4.x, node5.x, nodes6.x]
         x_n = [node4.x, node5.x, nodes6.x]
         yy = [0,0,0]
         y_n = [node4.y, node5.y, nodes6.y]
         xx.sort()
         for i in range(0,3):
             indx = x_n.index(xx[i])
             yy[i] = y_n[indx]  
         
         # quadratic interpolation
         f = interpolate.interp1d(np.array(xx), np.array(yy))
         coord_enrich = Coordinate(0,0)
         coord_enrich.x = midPoint.x
         coord_enrich.y = f(midPoint.x)
 
    else:
         xx = [0,0,0]
         x_n = [node4.x, node5.x, nodes6.x]
         yy = [node4.y, node5.y, nodes6.y]
         y_n = [node4.y, node5.y, nodes6.y]
         yy.sort()
         for i in range(0,3):
             indx = y_n.index(yy[i])
             xx[i] = x_n[indx]  
                      
         f = interpolate.interp1d(np.array(yy), np.array(xx))
         coord_enrich = Coordinate(0,0)
         coord_enrich.x = f(midPoint.y)
         coord_enrich.y = midPoint.y

    new_coord_enrich1 = Coordinate(0, 0)
    [c_center, c_rad] = which_circle(coord_enrich)
    new_coord_enrich1.x = c_center.x + c_rad * ( coord_enrich.x - c_center.x) / math.sqrt( math.pow(coord_enrich.x - c_center.x, 2) + math.pow(coord_enrich.y - c_center.y,2) )
    new_coord_enrich1.y = c_center.y + c_rad * ( coord_enrich.y - c_center.y) / math.sqrt( math.pow(coord_enrich.x - c_center.x, 2) + math.pow( coord_enrich.y - c_center.y,2) )
           
    return new_coord_enrich1

    return coord_enrich

def coord_enrich_comp_quad(root):
    
    midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
    
    if root.enrichNodes[0].x != root.enrichNodes[1].x:
         xx = [root.enrichNodes[0].x, root.enrichNodes[1].x, root.enrichNodes[2].x]
         x_n = [root.enrichNodes[0].x, root.enrichNodes[1].x, root.enrichNodes[2].x]
         yy = [0,0,0]
         y_n = [root.enrichNodes[0].y, root.enrichNodes[1].y, root.enrichNodes[2].y]
         xx.sort()
         for i in range(0,3):
             indx = x_n.index(xx[i])
             yy[i] = y_n[indx]  
         
         # quadratic interpolation
         f = interpolate.interp1d(np.array(xx), np.array(yy))
         coord_enrich = Coordinate(0,0)
         coord_enrich.x = midPoint.x
         coord_enrich.y = f(midPoint.x)
 
    else:
         xx = [0,0,0]
         x_n = [root.enrichNodes[0].x, root.enrichNodes[1].x, root.enrichNodes[2].x]
         yy = [root.enrichNodes[0].y, root.enrichNodes[1].y, root.enrichNodes[2].y]
         y_n = [root.enrichNodes[0].y, root.enrichNodes[1].y, root.enrichNodes[2].y]
         yy.sort()
         for i in range(0,3):
             indx = y_n.index(yy[i])
             xx[i] = x_n[indx]  
             
         
         f = interpolate.interp1d(np.array(yy), np.array(xx))
         coord_enrich = Coordinate(0,0)
         coord_enrich.x = f(midPoint.y)
         coord_enrich.y = midPoint.y

    coord_enrich.x =  coord_enrich.x / 1000.0
    if coord_enrich.x != 0.0:
        coord_enrich.x += 0.001
    coord_enrich.y =  coord_enrich.y / 1000.0
    if coord_enrich.y != 0.0:
        coord_enrich.y += 0.001
    coord_enrich.y = 1 - coord_enrich.y
        

    new_coord_enrich1 = Coordinate(0, 0)
    [c_center, c_rad] = which_circle(coord_enrich)
    new_coord_enrich1.x = c_center.x + c_rad * ( coord_enrich.x - c_center.x) / math.sqrt( math.pow(coord_enrich.x - c_center.x, 2) + math.pow(coord_enrich.y - c_center.y,2) )
    new_coord_enrich1.y = c_center.y + c_rad * ( coord_enrich.y - c_center.y) / math.sqrt( math.pow(coord_enrich.x - c_center.x, 2) + math.pow( coord_enrich.y - c_center.y,2) )
           
    return new_coord_enrich1

    return coord_enrich

def coord_enrich_computation(E):        
        
    coord_enrich1 = Coordinate(E.x, E.y)
    coord_enrich1.x =  coord_enrich1.x / 1000.0
    if coord_enrich1.x != 0.0:
        coord_enrich1.x += 0.001
    coord_enrich1.y =  coord_enrich1.y / 1000.0
    if coord_enrich1.y != 0.0:
        coord_enrich1.y += 0.001
    coord_enrich1.y = 1 - coord_enrich1.y
    
    new_coord_enrich1 = Coordinate(0, 0)
    [c_center, c_rad] = which_circle(coord_enrich1)
    new_coord_enrich1.x = c_center.x + c_rad * ( coord_enrich1.x - c_center.x) / math.sqrt( math.pow(coord_enrich1.x - c_center.x, 2) + math.pow(coord_enrich1.y - c_center.y,2) )
    new_coord_enrich1.y = c_center.y + c_rad * ( coord_enrich1.y - c_center.y) / math.sqrt( math.pow(coord_enrich1.x - c_center.x, 2) + math.pow( coord_enrich1.y - c_center.y,2) )
           
    return new_coord_enrich1

def circumcenter_tri(coords):
# computing the circumcenter of the triangle defined by coordinates a, b, and c

    a = Coordinate(coords[0,0], coords[0,1])
    b = Coordinate(coords[1,0], coords[1,1])
    c = Coordinate(coords[2,0], coords[2,1])
        
    U =  Coordinate(0,0)
    
    d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y))

    U.x = ( ( (a.x * a.x + a.y * a.y) * (b.y - c.y) + 
              (b.x * b.x + b.y * b.y) * (c.y - a.y) + 
              (c.x * c.x + c.y * c.y) * (a.y - b.y)) 
           / d  )
    
    U.y = ( ( (a.x * a.x + a.y * a.y) * (c.x - b.x) + 
              (b.x * b.x  + b.y * b.y) * (a.x - c.x) + 
              (c.x * c.x + c.y * c.y) * (b.x - a.x) )
           /d )
    
    return U
    
def NW_corner(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7,pxVals):
    K = numpy.zeros((6,6))
    Fe = np.zeros((6,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper left corner cut)
    nodes1 = [nodess[0],nodess[1],nodess[4]]
    nodes2 = [nodess[4],nodess[1],nodess[5]]
    nodes3 = [nodess[1],nodess[2],nodess[5]]
    nodes4 = [nodess[4],nodess[5],nodess[3]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]
    coords4 = p[nodes4]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]
    
    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    #if NW
    if ( is_in_same_bin(pxVal1,pxVal4) == False and pxVal1 > binBnd[1] and
        ( is_in_same_bin(pxVal4,pxVal2)==True and is_in_same_bin(pxVal2,pxVal3)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]
        
    #ISOPARAMETRIC linear, quadratic, cubic
    [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
    J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )

    [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
    if len(root.enrichNodes) == 2:
        
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        
        [x_fct_4, y_fct_4] = tri_xy_fct( coords4[:,0], coords4[:,1] )
        J4 = tri_jacobian_mat( coords4[:,0], coords4[:,1] )
           
    if len(root.enrichNodes) == 3:

        coord_enrich = coord_enrich_comp_quad_circle(p[nodess],nodes6)
        
        lOrd = [1,2,0] # local order 
        
        vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
        vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]

        [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
        
        lOrd = [2,0,1]
        vec4_x = [ coords4[ lOrd[0],0], coords4[lOrd[1],0], coords4[lOrd[2],0], (coords4[ lOrd[0],0] + coords4[lOrd[1],0])/2.0, coord_enrich.x, (coords4[lOrd[0],0] + coords4[lOrd[2],0])/2.0  ]
        vec4_y = [ coords4[ lOrd[0],1], coords4[lOrd[1],1], coords4[lOrd[2],1], (coords4[ lOrd[0],1] + coords4[lOrd[1],1])/2.0, coord_enrich.y, (coords4[lOrd[0],1] + coords4[lOrd[2],1])/2.0  ]

        [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
        

    if len(root.enrichNodes) == 4:
     
        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess],nodes7)
        
        circumcenter_pt2 = circumcenter_tri(coords2)
        
        lOrd = [1,2,0]
        
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]

        pt6.x = coord_enrich2.x
        pt6.y = coord_enrich2.y     
        
        pt7.x = coord_enrich1.x
        pt7.y = coord_enrich1.y
        
        vec2_x = [ coords2[ lOrd[0],0], 
                  coords2[lOrd[1],0],
                  coords2[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt2.x  
                  ]
        vec2_y = [ coords2[ lOrd[0],1], 
                  coords2[lOrd[1],1], 
                  coords2[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y, 
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt2.y 
                  ]

        
        [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
        
        circumcenter_pt4 = circumcenter_tri(coords4)
        lOrd = [2,0,1]
        
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords4[lOrd[0],0] + 1.0/3.0 * coords4[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords4[lOrd[0],1] + 1.0/3.0 * coords4[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords4[lOrd[0],0] + 2.0/3.0 * coords4[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords4[lOrd[0],1] + 2.0/3.0 * coords4[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords4[lOrd[2],0] + 1.0/3.0 * coords4[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords4[lOrd[2],1] + 1.0/3.0 * coords4[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords4[lOrd[2],0] + 2.0/3.0 * coords4[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords4[lOrd[2],1] + 2.0/3.0 * coords4[lOrd[0],1]

        pt6.x = coord_enrich1.x
        pt6.y = coord_enrich1.y     
        
        pt7.x = coord_enrich2.x
        pt7.y = coord_enrich2.y
        
        vec4_x = [coords4[lOrd[0],0], 
                  coords4[lOrd[1],0],
                  coords4[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt4.x  
                  ]
        vec4_y = [coords4[lOrd[0],1], 
                  coords4[lOrd[1],1],
                  coords4[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y, 
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt4.y  
                  ]
        
        [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
        
        
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)
    det_J4 = lambda e,n: determinant(J4)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    # triangle IV
    Pe4 = numpy.zeros((3,3))
    Pe4[:,0] = numpy.ones((3,1)).transpose()
    Pe4[:,1] = p[nodes4[0:3],0]
    Pe4[:,2] = p[nodes4[0:3],1]
    C4 = numpy.linalg.inv(Pe4)
    Nbasis4 = tribasisFct(C4)
    Nx4 = triderivX(C4)
    Ny4 = triderivY(C4)
        
    node_enr_4 = [nodess[4]]
    node_enr_5 = [nodess[5]]
    coords_enr_4 = p[node_enr_4]
    coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    w1 = y1 - y4
    w2 = y4 - y0
    factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
    if factor_W > EPS_FACTOR:
        factor_W = 1

    n1 = x5 - x0
    n2 = x1 - x5
    factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
    if factor_N > EPS_FACTOR:
        factor_N = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis1[2](x,y), lambda x,y: 0]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis2[0](x,y), lambda x,y: factor_N * Nbasis2[2](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: 0, lambda x,y: factor_N * Nbasis3[2](x,y)]
    Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis4[0](x,y), lambda x,y: factor_N * Nbasis4[1](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx1[2](x,y), lambda x,y: 0]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx2[0](x,y), lambda x,y: factor_N * Nx2[2](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: 0, lambda x,y: factor_N * Nx3[2](x,y)]
    Nx_4 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx4[0](x,y), lambda x,y: factor_N * Nx4[1](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny1[2](x,y),lambda x,y: 0]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny2[0](x,y), lambda x,y: factor_N * Ny2[2](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: 0, lambda x,y: factor_N * Ny3[2](x,y)]
    Ny_4 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny4[0](x,y), lambda x,y: factor_N * Ny4[1](x,y)]

    for i in range(0,6):
        for j in range(0,6):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            Kefunc4 = lambda e,n: K_cst[3] * (
                                Nx_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Nx_4[j]( x_fct_4(e,n), y_fct_4(e,n)) +
                                Ny_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Ny_4[j]( x_fct_4(e,n), y_fct_4(e,n)) ) * det_J4(e,n)
            

            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)
            integral4 = my_gauss_rule(Kefunc4,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
            fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

            Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) + my_gauss_rule(fv4,ui,wi)

    return [K,Fe]


def SW_corner(p,ui,wi,k1,k2,nodess, root, image,nodes6,nodes7,pxVals):
    K = numpy.zeros((6,6))
    Fe = np.zeros((6,1))


    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (lower left corner cut)
    nodes1 = [nodess[0],nodess[4],nodess[5]]
    nodes2 = [nodess[5],nodess[2],nodess[3]]
    nodes3 = [nodess[4],nodess[2],nodess[5]]
    nodes4 = [nodess[4],nodess[1],nodess[2]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]
    coords4 = p[nodes4]
    

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    #if SW
    if ( is_in_same_bin(pxVal3,pxVal4) == False and pxVal4 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True and is_in_same_bin(pxVal2,pxVal3)) ):
        K_cst = [k2,k1,k1,k1]
    else:
        K_cst = [k1,k2,k2,k2]
        
             
    #ISOPARAMETRIC linear, quadratic, cubic triangles         
    [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
    J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )

    [x_fct_4, y_fct_4] = tri_xy_fct( coords4[:,0], coords4[:,1] )
    J4 = tri_jacobian_mat( coords4[:,0], coords4[:,1] )
            
    if len(root.enrichNodes) == 2:
        
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
           
    if len(root.enrichNodes) == 3:
        
        coord_enrich = coord_enrich_comp_quad_circle(p[nodess], nodes6)
              
        lOrd = [0,1,2] # local order 
          
        vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
        vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

        [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
        J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
        
        lOrd = [1,2,0]
        vec3_x = [ coords3[ lOrd[0],0], coords3[lOrd[1],0], coords3[lOrd[2],0], (coords3[ lOrd[0],0] + coords3[lOrd[1],0])/2.0, coord_enrich.x, (coords3[lOrd[0],0] + coords3[lOrd[2],0])/2.0  ]
        vec3_y = [ coords3[ lOrd[0],1], coords3[lOrd[1],1], coords3[lOrd[2],1], (coords3[ lOrd[0],1] + coords3[lOrd[1],1])/2.0, coord_enrich.y, (coords3[lOrd[0],1] + coords3[lOrd[2],1])/2.0  ]

        [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
        J3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
        
        
    if len(root.enrichNodes) == 4:
        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess],nodes7)
        
        circumcenter_pt1 = circumcenter_tri(coords1)
        lOrd = [0,1,2] # local order 
        
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
         
        pt4.x = 2.0/3.0 * coords1[lOrd[0],0] + 1.0/3.0 * coords1[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords1[lOrd[0],1] + 1.0/3.0 * coords1[lOrd[1],1]
         
        pt5.x = 1.0/3.0 * coords1[lOrd[0],0] + 2.0/3.0 * coords1[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords1[lOrd[0],1] + 2.0/3.0 * coords1[lOrd[1],1]
         
        pt8.x = 2.0/3.0 * coords1[lOrd[2],0] + 1.0/3.0 * coords1[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords1[lOrd[2],1] + 1.0/3.0 * coords1[lOrd[0],1]
         
        pt9.x = 1.0/3.0 * coords1[lOrd[2],0] + 2.0/3.0 * coords1[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords1[lOrd[2],1] + 2.0/3.0 * coords1[lOrd[0],1]
 
        pt6.x = coord_enrich1.x
        pt6.y = coord_enrich1.y     
         
        pt7.x = coord_enrich2.x
        pt7.y = coord_enrich2.y
         
        vec1_x = [coords1[lOrd[0],0], 
                  coords1[lOrd[1],0], 
                  coords1[lOrd[2],0], 
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt1.x  
                  ]
        vec1_y = [coords1[lOrd[0],1], 
                  coords1[lOrd[1],1], 
                  coords1[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y,
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt1.y  
                  ]
 

        [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
        J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
        
        circumcenter_pt3 = circumcenter_tri(coords3)
        lOrd = [1,2,0]

        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
         
        pt4.x = 2.0/3.0 * coords3[lOrd[0],0] + 1.0/3.0 * coords3[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords3[lOrd[0],1] + 1.0/3.0 * coords3[lOrd[1],1]
         
        pt5.x = 1.0/3.0 * coords3[lOrd[0],0] + 2.0/3.0 * coords3[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords3[lOrd[0],1] + 2.0/3.0 * coords3[lOrd[1],1]
         
        pt8.x = 2.0/3.0 * coords3[lOrd[2],0] + 1.0/3.0 * coords3[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords3[lOrd[2],1] + 1.0/3.0 * coords3[lOrd[0],1]
         
        pt9.x = 1.0/3.0 * coords3[lOrd[2],0] + 2.0/3.0 * coords3[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords3[lOrd[2],1] + 2.0/3.0 * coords3[lOrd[0],1]
 
        pt6.x = coord_enrich2.x
        pt6.y = coord_enrich2.y     
         
        pt7.x = coord_enrich1.x
        pt7.y = coord_enrich1.y
        
        vec3_x = [coords3[ lOrd[0],0], 
                  coords3[lOrd[1],0],
                  coords3[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt3.x  
                  ]
        vec3_y = [coords3[ lOrd[0],1], 
                  coords3[lOrd[1],1], 
                  coords3[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y,
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt3.y  
                  ]
        
        [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
        J3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
        
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)
    det_J4 = lambda e,n: determinant(J4)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    # triangle IV
    Pe4 = numpy.zeros((3,3))
    Pe4[:,0] = numpy.ones((3,1)).transpose()
    Pe4[:,1] = p[nodes4[0:3],0]
    Pe4[:,2] = p[nodes4[0:3],1]
    C4 = numpy.linalg.inv(Pe4)
    Nbasis4 = tribasisFct(C4)
    Nx4 = triderivX(C4)
    Ny4 = triderivY(C4)

    node_enr_4 = [nodess[4]]
    node_enr_5 = [nodess[5]]
    coords_enr_4 = p[node_enr_4]
    coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    w1 = y1 - y5
    w2 = y5 - y0
    factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
    if factor_W > EPS_FACTOR:
        factor_W = 1

    s1 = x4 - x0
    s2 = x1 - x4
    factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
    if factor_S > EPS_FACTOR:
        factor_S = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis1[1](x,y), lambda x,y: factor_W * Nbasis1[2](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: 0, lambda x,y: factor_W * Nbasis2[0](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis3[0](x,y), lambda x,y: factor_W * Nbasis3[2](x,y)]
    Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis4[0](x,y), lambda x,y: 0]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx1[1](x,y), lambda x,y: factor_W * Nx1[2](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: 0, lambda x,y: factor_W * Nx2[0](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx3[0](x,y), lambda x,y: factor_W * Nx3[2](x,y)]
    Nx_4 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx4[0](x,y), lambda x,y: 0]

    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny1[1](x,y), lambda x,y: factor_W * Ny1[2](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: 0, lambda x,y: factor_W * Ny2[0](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny3[0](x,y), lambda x,y: factor_W * Ny3[2](x,y)]
    Ny_4 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny4[0](x,y), lambda x,y: 0]

    for i in range(0,6):
        for j in range(0,6):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            Kefunc4 = lambda e,n: K_cst[3] * (
                                Nx_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Nx_4[j]( x_fct_4(e,n), y_fct_4(e,n)) +
                                Ny_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Ny_4[j]( x_fct_4(e,n), y_fct_4(e,n)) ) * det_J4(e,n)
            
            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)
            integral4 = my_gauss_rule(Kefunc4,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
        # construct the local load vector
        fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
        fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

        Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) + my_gauss_rule(fv4,ui,wi)

    #NEUMANN BCS are zero - code not inserted here
    return [K,Fe]

def NE_corner(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7,pxVals):
    K = numpy.zeros((6,6))
    Fe = np.zeros((6,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
    nodes1 = [nodess[0],nodess[5],nodess[3]]
    nodes2 = [nodess[0],nodess[4],nodess[5]]
    nodes3 = [nodess[0],nodess[1],nodess[4]]
    nodes4 = [nodess[4],nodess[2],nodess[5]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]
    coords4 = p[nodes4]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    #if NE
    if ( is_in_same_bin(pxVal1,pxVal2) == False and pxVal2 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal3)==True and is_in_same_bin(pxVal4,pxVal3)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]      
       
    #ISOPARAMETRIC linear, quadratic, cubic triangles   
    [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
    J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )

    [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
    if len(root.enrichNodes) == 2:
        
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        
        [x_fct_4, y_fct_4] = tri_xy_fct( coords4[:,0], coords4[:,1] )
        J4 = tri_jacobian_mat( coords4[:,0], coords4[:,1] )
           
    if len(root.enrichNodes) == 3:
        coord_enrich = coord_enrich_comp_quad_circle(p[nodess], nodes6)
        
        lOrd = [0,1,2] # local order 
        vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
        vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]

        [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
        
        lOrd = [1,2,0]
        vec4_x = [ coords4[ lOrd[0],0], coords4[lOrd[1],0], coords4[lOrd[2],0], (coords4[ lOrd[0],0] + coords4[lOrd[1],0])/2.0, coord_enrich.x, (coords4[lOrd[0],0] + coords4[lOrd[2],0])/2.0  ]
        vec4_y = [ coords4[ lOrd[0],1], coords4[lOrd[1],1], coords4[lOrd[2],1], (coords4[ lOrd[0],1] + coords4[lOrd[1],1])/2.0, coord_enrich.y, (coords4[lOrd[0],1] + coords4[lOrd[2],1])/2.0  ]

        [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
        
        
    if len(root.enrichNodes) == 4:
        
        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess],nodes7)
        
        circumcenter_pt2 = circumcenter_tri(coords2)
        
        lOrd = [0,1,2] 

        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]

        pt6.x = coord_enrich1.x
        pt6.y = coord_enrich1.y     
        
        pt7.x = coord_enrich2.x
        pt7.y = coord_enrich2.y
        
        vec2_x = [coords2[lOrd[0],0], 
                  coords2[lOrd[1],0],
                  coords2[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt2.x  
                  ]
        vec2_y = [coords2[lOrd[0],1], 
                  coords2[lOrd[1],1], 
                  coords2[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y, 
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt2.y 
                  ]
        
        
        [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
        
        circumcenter_pt4 = circumcenter_tri(coords4)
        lOrd = [1,2,0]
        
        
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords4[lOrd[0],0] + 1.0/3.0 * coords4[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords4[lOrd[0],1] + 1.0/3.0 * coords4[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords4[lOrd[0],0] + 2.0/3.0 * coords4[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords4[lOrd[0],1] + 2.0/3.0 * coords4[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords4[lOrd[2],0] + 1.0/3.0 * coords4[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords4[lOrd[2],1] + 1.0/3.0 * coords4[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords4[lOrd[2],0] + 2.0/3.0 * coords4[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords4[lOrd[2],1] + 2.0/3.0 * coords4[lOrd[0],1]

        pt6.x = coord_enrich2.x
        pt6.y = coord_enrich2.y     
        
        pt7.x = coord_enrich1.x
        pt7.y = coord_enrich1.y
        
        vec4_x = [coords4[lOrd[0],0], 
                    coords4[lOrd[1],0],
                    coords4[lOrd[2],0],
                    pt4.x, 
                    pt5.x,
                    pt6.x,
                    pt7.x,
                    pt8.x,
                    pt9.x,
                    circumcenter_pt4.x  
                   ]
        vec4_y = [coords4[lOrd[0],1], 
                   coords4[lOrd[1],1],
                   coords4[lOrd[2],1], 
                   pt4.y, 
                   pt5.y,
                   pt6.y, 
                   pt7.y,
                   pt8.y,
                   pt9.y,
                   circumcenter_pt4.y  
                   ]
                
        [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
        
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)
    det_J4 = lambda e,n: determinant(J4)(e,n)

    
    
    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    # triangle IV
    Pe4 = numpy.zeros((3,3))
    Pe4[:,0] = numpy.ones((3,1)).transpose()
    Pe4[:,1] = p[nodes4[0:3],0]
    Pe4[:,2] = p[nodes4[0:3],1]
    C4 = numpy.linalg.inv(Pe4)
    Nbasis4 = tribasisFct(C4)
    Nx4 = triderivX(C4)
    Ny4 = triderivY(C4)
        
    node_enr_4 = [nodess[4]]
    node_enr_5 = [nodess[5]]
    coords_enr_4 = p[node_enr_4]
    coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    e1 = y1 - y4
    e2 = y4 - y0
    factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
    if factor_E > EPS_FACTOR:
        factor_E = 1

    n1 = x5 - x0
    n2 = x1 - x5
    factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
    if factor_N > EPS_FACTOR:
        factor_N = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: 0, lambda x,y: factor_N * Nbasis1[1](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis2[1](x,y), lambda x,y: factor_N * Nbasis2[2](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis3[2](x,y), lambda x,y: 0]
    Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis4[0](x,y), lambda x,y: factor_N * Nbasis4[2](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: 0, lambda x,y: factor_N * Nx1[1](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx2[1](x,y), lambda x,y: factor_N * Nx2[2](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx3[2](x,y), lambda x,y: 0]
    Nx_4 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx4[0](x,y), lambda x,y: factor_N * Nx4[2](x,y)]

    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: 0, lambda x,y: factor_N * Ny1[1](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny2[1](x,y), lambda x,y: factor_N * Ny2[2](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny3[2](x,y), lambda x,y: 0]
    Ny_4 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny4[0](x,y), lambda x,y: factor_N * Ny4[2](x,y)]

    for i in range(0,6):
        for j in range(0,6):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(x_fct_2(e,n), y_fct_2(e,n))

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            Kefunc4 = lambda e,n: K_cst[3] * (
                                Nx_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Nx_4[j]( x_fct_4(e,n), y_fct_4(e,n)) +
                                Ny_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Ny_4[j]( x_fct_4(e,n), y_fct_4(e,n)) ) * det_J4(e,n)
            
            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)
            integral4 = my_gauss_rule(Kefunc4,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
            fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

            Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) + my_gauss_rule(fv4,ui,wi)

    #NEUMANN BCS are zero - code not inserted here
    return [K,Fe]

def SE_corner(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7,pxVals):
#def SE_corner(p,ui,wi,k1,k2,nodess,UConf,pConf,tConf):
    K = numpy.zeros((6,6))
    Fe = np.zeros((6,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (lower right corner cut)
    nodes1 = [nodess[0],nodess[4],nodess[3]]
    nodes2 = [nodess[4],nodess[5],nodess[3]]
    nodes3 = [nodess[5],nodess[2],nodess[3]]
    nodes4 = [nodess[4],nodess[1],nodess[5]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]
    coords4 = p[nodes4]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    #if SE
    if ( is_in_same_bin(pxVal3,pxVal4) == False and pxVal3 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True and is_in_same_bin(pxVal2,pxVal4)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]

    #ISOPARAMETRIC linear, quadratic, cubic triangles    
    [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
    J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )

    [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
    if len(root.enrichNodes) == 2:
        
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
           
        [x_fct_4, y_fct_4] = tri_xy_fct( coords4[:,0], coords4[:,1] )
        J4 = tri_jacobian_mat( coords4[:,0], coords4[:,1] )
        
    if len(root.enrichNodes) == 3:
        coord_enrich = coord_enrich_comp_quad_circle(p[nodess], nodes6)
        
        lOrd = [2,0,1] # local order    
        vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
        vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]

        [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
        
        lOrd = [1,2,0]
        vec4_x = [ coords4[ lOrd[0],0], coords4[lOrd[1],0], coords4[lOrd[2],0], (coords4[ lOrd[0],0] + coords4[lOrd[1],0])/2.0, coord_enrich.x, (coords4[lOrd[0],0] + coords4[lOrd[2],0])/2.0  ]
        vec4_y = [ coords4[ lOrd[0],1], coords4[lOrd[1],1], coords4[lOrd[2],1], (coords4[ lOrd[0],1] + coords4[lOrd[1],1])/2.0, coord_enrich.y, (coords4[lOrd[0],1] + coords4[lOrd[2],1])/2.0  ]

        [x_fct_4, y_fct_4] = tri_xy_fct_quadratic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_quadratic( vec4_x, vec4_y )
        
    if len(root.enrichNodes) == 4:

        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodess],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodess],nodes7)
        
        circumcenter_pt2 = circumcenter_tri(coords2)
        lOrd = [2,0,1] 
        
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]

        pt6.x = coord_enrich1.x
        pt6.y = coord_enrich1.y     
        
        pt7.x = coord_enrich2.x
        pt7.y = coord_enrich2.y
        
        vec2_x = [coords2[lOrd[0],0], 
                  coords2[lOrd[1],0],
                  coords2[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt2.x  
                  ]
        vec2_y = [coords2[lOrd[0],1], 
                  coords2[lOrd[1],1], 
                  coords2[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y, 
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt2.y 
                  ]

        [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
        J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                   
        circumcenter_pt4 = circumcenter_tri(coords4)
        
        lOrd = [1,2,0]
                
        pt4 = Coordinate(0,0)
        pt5 = Coordinate(0,0)
        pt6 = Coordinate(0,0)
        pt7 = Coordinate(0,0)
        pt8 = Coordinate(0,0)
        pt9 = Coordinate(0,0)
        
        pt4.x = 2.0/3.0 * coords4[lOrd[0],0] + 1.0/3.0 * coords4[lOrd[1],0]
        pt4.y = 2.0/3.0 * coords4[lOrd[0],1] + 1.0/3.0 * coords4[lOrd[1],1]
        
        pt5.x = 1.0/3.0 * coords4[lOrd[0],0] + 2.0/3.0 * coords4[lOrd[1],0]
        pt5.y = 1.0/3.0 * coords4[lOrd[0],1] + 2.0/3.0 * coords4[lOrd[1],1]
        
        pt8.x = 2.0/3.0 * coords4[lOrd[2],0] + 1.0/3.0 * coords4[lOrd[0],0]
        pt8.y = 2.0/3.0 * coords4[lOrd[2],1] + 1.0/3.0 * coords4[lOrd[0],1]
        
        pt9.x = 1.0/3.0 * coords4[lOrd[2],0] + 2.0/3.0 * coords4[lOrd[0],0]
        pt9.y = 1.0/3.0 * coords4[lOrd[2],1] + 2.0/3.0 * coords4[lOrd[0],1]

        pt6.x = coord_enrich2.x
        pt6.y = coord_enrich2.y     
        
        pt7.x = coord_enrich1.x
        pt7.y = coord_enrich1.y
         
        vec4_x = [coords4[lOrd[0],0], 
                  coords4[lOrd[1],0],
                  coords4[lOrd[2],0],
                  pt4.x, 
                  pt5.x,
                  pt6.x,
                  pt7.x,
                  pt8.x,
                  pt9.x,
                  circumcenter_pt4.x  
                  ]
        vec4_y = [coords4[lOrd[0],1], 
                  coords4[lOrd[1],1],
                  coords4[lOrd[2],1], 
                  pt4.y, 
                  pt5.y,
                  pt6.y, 
                  pt7.y,
                  pt8.y,
                  pt9.y,
                  circumcenter_pt4.y  
                  ]
                 
        [x_fct_4, y_fct_4] = tri_xy_fct_cubic( vec4_x, vec4_y )
        J4 = tri_jacobian_mat_cubic( vec4_x, vec4_y )
        
        
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)
    det_J4 = lambda e,n: determinant(J4)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    # triangle IV
    Pe4 = numpy.zeros((3,3))
    Pe4[:,0] = numpy.ones((3,1)).transpose()
    Pe4[:,1] = p[nodes4[0:3],0]
    Pe4[:,2] = p[nodes4[0:3],1]
    C4 = numpy.linalg.inv(Pe4)
    Nbasis4 = tribasisFct(C4)
    Nx4 = triderivX(C4)
    Ny4 = triderivY(C4)
        
    node_enr_4 = [nodess[4]]
    coords_enr_4 = p[node_enr_4]
    node_enr_5 = [nodess[5]]
    coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    e1 = y1 - y5
    e2 = y5 - y0
    factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
    if factor_E > EPS_FACTOR:
        factor_E = 1

    s1 = x4 - x0
    s2 = x1 - x4
    factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
    if factor_S > EPS_FACTOR:
        factor_S = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis1[1](x,y), lambda x,y: 0]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis2[0](x,y), lambda x,y: factor_E * Nbasis2[1](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: 0, lambda x,y: factor_E * Nbasis3[0](x,y)]
    Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis4[0](x,y), lambda x,y: factor_E * Nbasis4[2](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx1[1](x,y), lambda x,y: 0]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx2[0](x,y), lambda x,y: factor_E * Nx2[1](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: 0, lambda x,y: factor_E * Nx3[0](x,y)]
    Nx_4 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx4[0](x,y), lambda x,y: factor_E * Nx4[2](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny1[1](x,y), lambda x,y: 0]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny2[0](x,y), lambda x,y: factor_E * Ny2[1](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: 0, lambda x,y: factor_E * Ny3[0](x,y)]
    Ny_4 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny4[0](x,y), lambda x,y: factor_E * Ny4[2](x,y)]

    for i in range(0,6):
        for j in range(0,6):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            Kefunc4 = lambda e,n: K_cst[3] * (
                                Nx_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Nx_4[j]( x_fct_4(e,n), y_fct_4(e,n)) +
                                Ny_4[i]( x_fct_4(e,n), y_fct_4(e,n)) *
                                Ny_4[j]( x_fct_4(e,n), y_fct_4(e,n)) ) * det_J4(e,n)
            
            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)
            integral4 = my_gauss_rule(Kefunc4,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 + integral4

    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
        fv4 = lambda e,n:  rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

        Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) + my_gauss_rule(fv4,ui,wi)

    return [K,Fe]


def East_edge(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7,full_nodes,pxVals,pxVals2):
    K = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    nodes1 = [nodess[0],nodess[1],nodess[4]]
    nodes2 = [nodess[0],nodess[4],nodess[3]]
    nodes3 = [nodess[4],nodess[2],nodess[3]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    pxVal14 = pxVals2[0]
    pxVal12 = pxVals2[1]
    pxVal23 = pxVals2[2]
    pxVal34 = pxVals2[3]
    
    #if SE - East
    if ( is_in_same_bin(pxVal3,pxVal2) == False and pxVal3 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True ) and
        is_in_same_bin(pxVal12,pxVal14) == True):
        K_cst = [k2,k1,k1]
    else: 
        if (  is_in_same_bin(pxVal3,pxVal2) == False and pxVal3 <= binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True ) and
        is_in_same_bin(pxVal12,pxVal14) == True):
            K_cst = [k1,k2,k2]

    # if NE - East
    if ( is_in_same_bin(pxVal3,pxVal2) == False and pxVal2 > binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal4)==True ) and
        is_in_same_bin(pxVal14,pxVal34) == True ):
        K_cst = [k1,k1,k2]
    else:
        if( is_in_same_bin(pxVal3,pxVal2) == False and pxVal2 <= binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal4)==True ) and
        is_in_same_bin(pxVal14,pxVal34) == True ):
             K_cst = [k2,k2,k1]

    if len(root.enrichNodes) == 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        coord_enrich = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
            
            lOrd = [1,2,0]
            vec3_x = [ coords3[ lOrd[0],0], coords3[lOrd[1],0], coords3[lOrd[2],0], (coords3[ lOrd[0],0] + coords3[lOrd[1],0])/2.0, coord_enrich.x, (coords3[lOrd[0],0] + coords3[lOrd[2],0])/2.0  ]
            vec3_y = [ coords3[ lOrd[0],1], coords3[lOrd[1],1], coords3[lOrd[2],1], (coords3[ lOrd[0],1] + coords3[lOrd[1],1])/2.0, coord_enrich.y, (coords3[lOrd[0],1] + coords3[lOrd[2],1])/2.0  ]

            [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            
            lOrd = [1,2,0]
            vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
            vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

            [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
            lOrd = [2,0,1] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )

        
    if len(root.enrichNodes) == 4:

        coord_enrich1 = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[full_nodes], nodes7)
                                
        circumcenter_pt1 = circumcenter_tri(coords1)
        circumcenter_pt2 = circumcenter_tri(coords2)
        circumcenter_pt3 = circumcenter_tri(coords3)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2]
        
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
            
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
        
            lOrd = [1,2,0]
            
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords3[lOrd[0],0] + 1.0/3.0 * coords3[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords3[lOrd[0],1] + 1.0/3.0 * coords3[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords3[lOrd[0],0] + 2.0/3.0 * coords3[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords3[lOrd[0],1] + 2.0/3.0 * coords3[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords3[lOrd[2],0] + 1.0/3.0 * coords3[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords3[lOrd[2],1] + 1.0/3.0 * coords3[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords3[lOrd[2],0] + 2.0/3.0 * coords3[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords3[lOrd[2],1] + 2.0/3.0 * coords3[lOrd[0],1]
     
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
             
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
            
            vec3_x = [coords3[ lOrd[0],0], 
                      coords3[lOrd[1],0],
                      coords3[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt3.x  
                      ]
            vec3_y = [coords3[ lOrd[0],1], 
                      coords3[lOrd[1],1], 
                      coords3[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt3.y  
                      ]

            [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            
            lOrd = [1,2,0]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords1[lOrd[0],0] + 1.0/3.0 * coords1[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords1[lOrd[0],1] + 1.0/3.0 * coords1[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords1[lOrd[0],0] + 2.0/3.0 * coords1[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords1[lOrd[0],1] + 2.0/3.0 * coords1[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords1[lOrd[2],0] + 1.0/3.0 * coords1[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords1[lOrd[2],1] + 1.0/3.0 * coords1[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords1[lOrd[2],0] + 2.0/3.0 * coords1[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords1[lOrd[2],1] + 2.0/3.0 * coords1[lOrd[0],1]
     
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
             
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
             
            vec1_x = [coords1[lOrd[0],0], 
                      coords1[lOrd[1],0], 
                      coords1[lOrd[2],0], 
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt1.x  
                      ]
            vec1_y = [coords1[lOrd[0],1], 
                      coords1[lOrd[1],1], 
                      coords1[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt1.y  
                      ]
            [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                        
            lOrd = [2,0,1] 
              
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
            
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
            
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)

    node_enr_4 = [nodess[4]]
    coords_enr_4 = p[node_enr_4]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    e1 = y1 - y4
    e2 = y4 - y0
    factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
    if factor_E > EPS_FACTOR:
        factor_E = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis1[2](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis2[1](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_E * Nbasis3[0](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx1[1](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx2[0](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_E * Nx3[0](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny1[2](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny2[1](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_E * Ny3[0](x,y)]

    for i in range(0,5):
        for j in range(0,5):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 

    # construct the local matrix and local components of the load vector
    for i in range(0,5):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

        Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) 

        #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]
    

def South_edge(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7, full_nodes, pxVals, pxVals2):
    
    Ke = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    nodes1 = [nodess[0],nodess[4],nodess[3]]
    nodes2 = [nodess[4],nodess[2],nodess[3]]
    nodes3 = [nodess[4],nodess[1],nodess[2]]
#     nodes3 = [nodess[1],nodess[4],nodess[2]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    pxVal14 = pxVals2[0]
    pxVal12 = pxVals2[1]
    pxVal23 = pxVals2[2]
    pxVal34 = pxVals2[3]
    
    #if SE-South
    if ( is_in_same_bin(pxVal3,pxVal4) == False and pxVal3 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal4) == True ) and
        is_in_same_bin(pxVal14,pxVal12) == True):
        K_cst = [k1,k1,k2]
    else: 
        if (  is_in_same_bin(pxVal3,pxVal4) == False and pxVal3 <= binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal4) == True ) and
        is_in_same_bin(pxVal14,pxVal12) == True):
            K_cst = [k2,k2,k1]

    # if SW-South
    if ( is_in_same_bin(pxVal3,pxVal4) == False and pxVal4 > binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal2) == True ) and
        is_in_same_bin(pxVal12,pxVal23) == True):
        K_cst = [k2,k1,k1]
    else:
        if( is_in_same_bin(pxVal3,pxVal4) == False and pxVal4 <= binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal2)==True ) and
        is_in_same_bin(pxVal12,pxVal23) == True):
             K_cst = [k1,k2,k2]
    
    if len(root.enrichNodes) >= 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        coord_enrich = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [2,0,1] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
            
            lOrd = [1,2,0]
            vec3_x = [ coords3[ lOrd[0],0], coords3[lOrd[1],0], coords3[lOrd[2],0], (coords3[ lOrd[0],0] + coords3[lOrd[1],0])/2.0, coord_enrich.x, (coords3[lOrd[0],0] + coords3[lOrd[2],0])/2.0  ]
            vec3_y = [ coords3[ lOrd[0],1], coords3[lOrd[1],1], coords3[lOrd[2],1], (coords3[ lOrd[0],1] + coords3[lOrd[1],1])/2.0, coord_enrich.y, (coords3[lOrd[0],1] + coords3[lOrd[2],1])/2.0  ]

            [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [0,1,2]
            vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
            vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

            [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
            lOrd = [1,2,0] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )

    if len(root.enrichNodes) == 4:

        coord_enrich1 = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[full_nodes], nodes7)
        
        circumcenter_pt1 = circumcenter_tri(coords1)
        circumcenter_pt2 = circumcenter_tri(coords2)
        circumcenter_pt3 = circumcenter_tri(coords3)
        
                    
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [2,0,1] # local order    
            
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
            
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                                
            lOrd = [1,2,0]
            
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords3[lOrd[0],0] + 1.0/3.0 * coords3[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords3[lOrd[0],1] + 1.0/3.0 * coords3[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords3[lOrd[0],0] + 2.0/3.0 * coords3[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords3[lOrd[0],1] + 2.0/3.0 * coords3[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords3[lOrd[2],0] + 1.0/3.0 * coords3[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords3[lOrd[2],1] + 1.0/3.0 * coords3[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords3[lOrd[2],0] + 2.0/3.0 * coords3[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords3[lOrd[2],1] + 2.0/3.0 * coords3[lOrd[0],1]
     
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
             
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
            
            vec3_x = [coords3[ lOrd[0],0], 
                      coords3[lOrd[1],0],
                      coords3[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt3.x  
                      ]
            vec3_y = [coords3[ lOrd[0],1], 
                      coords3[lOrd[1],1], 
                      coords3[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt3.y  
                      ]

            [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [0,1,2]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords1[lOrd[0],0] + 1.0/3.0 * coords1[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords1[lOrd[0],1] + 1.0/3.0 * coords1[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords1[lOrd[0],0] + 2.0/3.0 * coords1[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords1[lOrd[0],1] + 2.0/3.0 * coords1[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords1[lOrd[2],0] + 1.0/3.0 * coords1[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords1[lOrd[2],1] + 1.0/3.0 * coords1[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords1[lOrd[2],0] + 2.0/3.0 * coords1[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords1[lOrd[2],1] + 2.0/3.0 * coords1[lOrd[0],1]
     
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
             
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
             
            vec1_x = [coords1[lOrd[0],0], 
                      coords1[lOrd[1],0], 
                      coords1[lOrd[2],0], 
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt1.x  
                      ]
            vec1_y = [coords1[lOrd[0],1], 
                      coords1[lOrd[1],1], 
                      coords1[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt1.y  
                      ]
            [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )           
                        
            lOrd = [1,2,0] # local order    
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
            
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
            
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)
    
    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)

    node_enr_4 = [nodess[4]]
    coords_enr_4 = p[node_enr_4]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    s1 = x4 - x0
    s2 = x1 - x4
    factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
    if factor_S > EPS_FACTOR:
        factor_S = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis1[1](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis2[0](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis3[0](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx1[1](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx2[0](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx3[0](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny1[1](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny2[0](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny3[0](x,y)]
   
    
    for i in range(0,5):
        for j in range(0,5):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)

            Ke[i,j] = integral1 + integral2 + integral3 

    # construct the local matrix and local components of the load vector
    for i in range(0,5):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

        Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) 

    return [Ke,Fe]

def North_edge(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7, full_nodes,pxVals,pxVals2):
    K = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

        
    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    nodes1 = [nodess[0],nodess[4],nodess[3]]
    nodes2 = [nodess[0],nodess[1],nodess[4]]
    nodes3 = [nodess[1],nodess[2],nodess[4]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]
    
    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    pxVal14 = pxVals2[0]
    pxVal12 = pxVals2[1]
    pxVal23 = pxVals2[2]
    pxVal34 = pxVals2[3]
    
    #if NW - North
    if ( is_in_same_bin(pxVal1,pxVal2) == False and pxVal1 > binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal2) == True ) and
        is_in_same_bin(pxVal34, pxVal23) == True):
        K_cst = [k2,k1,k1]
    else: 
        if (  is_in_same_bin(pxVal1,pxVal2) == False and pxVal1 <= binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal2) == True ) and
        is_in_same_bin(pxVal34, pxVal23) == True):
            K_cst = [k1,k2,k2]

    # if NE - North 
    if ( is_in_same_bin(pxVal1,pxVal2) == False and pxVal2 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal4) == True ) and 
        is_in_same_bin(pxVal14,pxVal34) == True):
        K_cst = [k1,k1,k2]
    else:
        if( is_in_same_bin(pxVal1,pxVal2) == False and pxVal2 <= binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal4) == True ) and 
        is_in_same_bin(pxVal14,pxVal34) == True ):
             K_cst = [k2,k2,k1]
    
            
    if len(root.enrichNodes) >= 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:
        coord_enrich = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
            
            lOrd = [1,2,0]
            vec3_x = [ coords3[ lOrd[0],0], coords3[lOrd[1],0], coords3[lOrd[2],0], (coords3[ lOrd[0],0] + coords3[lOrd[1],0])/2.0, coord_enrich.x, (coords3[lOrd[0],0] + coords3[lOrd[2],0])/2.0  ]
            vec3_y = [ coords3[ lOrd[0],1], coords3[lOrd[1],1], coords3[lOrd[2],1], (coords3[ lOrd[0],1] + coords3[lOrd[1],1])/2.0, coord_enrich.y, (coords3[lOrd[0],1] + coords3[lOrd[2],1])/2.0  ]

            [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [2,0,1]
            vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
            vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

            [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
            lOrd = [1,2,0] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
            
    if len(root.enrichNodes) == 4:

        coord_enrich1 = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[full_nodes], nodes7)
        
        circumcenter_pt1 = circumcenter_tri(coords1)
        circumcenter_pt2 = circumcenter_tri(coords2)
        circumcenter_pt3 = circumcenter_tri(coords3)

        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2] # local order    
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
            
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                    
            lOrd = [1,2,0]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords3[lOrd[0],0] + 1.0/3.0 * coords3[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords3[lOrd[0],1] + 1.0/3.0 * coords3[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords3[lOrd[0],0] + 2.0/3.0 * coords3[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords3[lOrd[0],1] + 2.0/3.0 * coords3[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords3[lOrd[2],0] + 1.0/3.0 * coords3[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords3[lOrd[2],1] + 1.0/3.0 * coords3[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords3[lOrd[2],0] + 2.0/3.0 * coords3[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords3[lOrd[2],1] + 2.0/3.0 * coords3[lOrd[0],1]
     
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
             
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
            
            vec3_x = [coords3[ lOrd[0],0], 
                      coords3[lOrd[1],0],
                      coords3[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt3.x  
                      ]
            vec3_y = [coords3[ lOrd[0],1], 
                      coords3[lOrd[1],1], 
                      coords3[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt3.y  
                      ]

            [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [2,0,1]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords1[lOrd[0],0] + 1.0/3.0 * coords1[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords1[lOrd[0],1] + 1.0/3.0 * coords1[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords1[lOrd[0],0] + 2.0/3.0 * coords1[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords1[lOrd[0],1] + 2.0/3.0 * coords1[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords1[lOrd[2],0] + 1.0/3.0 * coords1[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords1[lOrd[2],1] + 1.0/3.0 * coords1[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords1[lOrd[2],0] + 2.0/3.0 * coords1[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords1[lOrd[2],1] + 2.0/3.0 * coords1[lOrd[0],1]
     
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
             
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
             
            vec1_x = [coords1[lOrd[0],0], 
                      coords1[lOrd[1],0], 
                      coords1[lOrd[2],0], 
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt1.x  
                      ]
            vec1_y = [coords1[lOrd[0],1], 
                      coords1[lOrd[1],1], 
                      coords1[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt1.y  
                      ]
            [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
            
                        
            lOrd = [1,2,0] # local order    
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
            
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )  
            
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    node_enr_4 = [nodess[4]]
    coords_enr_4 = p[node_enr_4]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    n1 = x4 - x0
    n2 = x1 - x4
    factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
    if factor_N > EPS_FACTOR:
        factor_N = 1

    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_N * Nbasis1[1](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_N * Nbasis2[2](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_N * Nbasis3[2](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_N * Nx1[1](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_N * Nx2[2](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_N * Nx3[2](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_N * Ny1[1](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_N * Ny2[2](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_N * Ny3[2](x,y)]

    for i in range(0,5):
        for j in range(0,5):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 

        # construct the local matrix and local components of the load vector
    for i in range(0,5):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

            Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) 

            #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]



def West_edge(p,ui,wi,k1,k2,nodess,root,image,nodes6,nodes7,full_nodes, pxVals, pxVals2):
    K = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper left corner cut)
    nodes1 = [nodess[0],nodess[1],nodess[4]]
    nodes2 = [nodess[1],nodess[2],nodess[4]]
    nodes3 = [nodess[4],nodess[2],nodess[3]]

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    pxVal14 = pxVals2[0]
    pxVal12 = pxVals2[1]
    pxVal23 = pxVals2[2]
    pxVal34 = pxVals2[3]
    
    #if NW - West
    if ( is_in_same_bin(pxVal1,pxVal4) == False and pxVal1 > binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal4)==True ) and
        is_in_same_bin(pxVal34,pxVal23) == True):
        K_cst = [k1,k1,k2]
    else: 
        if (  is_in_same_bin(pxVal1,pxVal4) == False and pxVal1 <= binBnd[1] and
        ( is_in_same_bin(pxVal3,pxVal4)==True ) and
        is_in_same_bin(pxVal34,pxVal23) == True):
            K_cst = [k2,k2,k1]

    # if SW - West
    if ( is_in_same_bin(pxVal1,pxVal4) == False and pxVal4 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True ) and
        is_in_same_bin(pxVal12, pxVal23) == True ):
        K_cst = [k2,k1,k1]
    else:
        if( is_in_same_bin(pxVal1,pxVal4) == False and pxVal4 <= binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2) == True ) and
        is_in_same_bin(pxVal12, pxVal23) == True):
             K_cst = [k1,k2,k2] 
        
    if len(root.enrichNodes) >= 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        coord_enrich = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
            
            lOrd = [2,0,1]
            vec3_x = [ coords3[ lOrd[0],0], coords3[lOrd[1],0], coords3[lOrd[2],0], (coords3[ lOrd[0],0] + coords3[lOrd[1],0])/2.0, coord_enrich.x, (coords3[lOrd[0],0] + coords3[lOrd[2],0])/2.0  ]
            vec3_y = [ coords3[ lOrd[0],1], coords3[lOrd[1],1], coords3[lOrd[2],1], (coords3[ lOrd[0],1] + coords3[lOrd[1],1])/2.0, coord_enrich.y, (coords3[lOrd[0],1] + coords3[lOrd[2],1])/2.0  ]

            [x_fct_3, y_fct_3] = tri_xy_fct_quadratic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_quadratic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [0,1,2]
            vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
            vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

            [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
            lOrd = [1,2,0] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
   
    if len(root.enrichNodes) == 4:

        coord_enrich1 = coord_enrich_comp_quad_circle(p[full_nodes], nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[full_nodes], nodes7)
        
        circumcenter_pt1 = circumcenter_tri(coords1)
        circumcenter_pt2 = circumcenter_tri(coords2)
        circumcenter_pt3 = circumcenter_tri(coords3)
        
        if K_cst[0] == K_cst[1]:
            # triangle 3 is the one with curved edge
            [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
            J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
            
            lOrd = [0,1,2] # local order    
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
            
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )
                    
            lOrd = [2,0,1]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords3[lOrd[0],0] + 1.0/3.0 * coords3[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords3[lOrd[0],1] + 1.0/3.0 * coords3[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords3[lOrd[0],0] + 2.0/3.0 * coords3[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords3[lOrd[0],1] + 2.0/3.0 * coords3[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords3[lOrd[2],0] + 1.0/3.0 * coords3[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords3[lOrd[2],1] + 1.0/3.0 * coords3[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords3[lOrd[2],0] + 2.0/3.0 * coords3[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords3[lOrd[2],1] + 2.0/3.0 * coords3[lOrd[0],1]
     
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
             
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
            
            vec3_x = [coords3[ lOrd[0],0], 
                      coords3[lOrd[1],0],
                      coords3[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt3.x  
                      ]
            vec3_y = [coords3[ lOrd[0],1], 
                      coords3[lOrd[1],1], 
                      coords3[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt3.y  
                      ]

            [x_fct_3, y_fct_3] = tri_xy_fct_cubic( vec3_x, vec3_y )
            J3 = tri_jacobian_mat_cubic( vec3_x, vec3_y )
            
        if K_cst[1] == K_cst[2]:
            # triangle 1 is the one with curved edge

            [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
            J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )
            
            lOrd = [0,1,2]
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
             
            pt4.x = 2.0/3.0 * coords1[lOrd[0],0] + 1.0/3.0 * coords1[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords1[lOrd[0],1] + 1.0/3.0 * coords1[lOrd[1],1]
             
            pt5.x = 1.0/3.0 * coords1[lOrd[0],0] + 2.0/3.0 * coords1[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords1[lOrd[0],1] + 2.0/3.0 * coords1[lOrd[1],1]
             
            pt8.x = 2.0/3.0 * coords1[lOrd[2],0] + 1.0/3.0 * coords1[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords1[lOrd[2],1] + 1.0/3.0 * coords1[lOrd[0],1]
             
            pt9.x = 1.0/3.0 * coords1[lOrd[2],0] + 2.0/3.0 * coords1[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords1[lOrd[2],1] + 2.0/3.0 * coords1[lOrd[0],1]
     
            pt6.x = coord_enrich1.x
            pt6.y = coord_enrich1.y     
             
            pt7.x = coord_enrich2.x
            pt7.y = coord_enrich2.y
             
            vec1_x = [coords1[lOrd[0],0], 
                      coords1[lOrd[1],0], 
                      coords1[lOrd[2],0], 
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt1.x  
                      ]
            vec1_y = [coords1[lOrd[0],1], 
                      coords1[lOrd[1],1], 
                      coords1[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y,
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt1.y  
                      ]
            [x_fct_1, y_fct_1] = tri_xy_fct_cubic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_cubic( vec1_x, vec1_y )
                        
            lOrd = [1,2,0] # local order    
            pt4 = Coordinate(0,0)
            pt5 = Coordinate(0,0)
            pt6 = Coordinate(0,0)
            pt7 = Coordinate(0,0)
            pt8 = Coordinate(0,0)
            pt9 = Coordinate(0,0)
        
            pt4.x = 2.0/3.0 * coords2[lOrd[0],0] + 1.0/3.0 * coords2[lOrd[1],0]
            pt4.y = 2.0/3.0 * coords2[lOrd[0],1] + 1.0/3.0 * coords2[lOrd[1],1]
            
            pt5.x = 1.0/3.0 * coords2[lOrd[0],0] + 2.0/3.0 * coords2[lOrd[1],0]
            pt5.y = 1.0/3.0 * coords2[lOrd[0],1] + 2.0/3.0 * coords2[lOrd[1],1]
            
            pt8.x = 2.0/3.0 * coords2[lOrd[2],0] + 1.0/3.0 * coords2[lOrd[0],0]
            pt8.y = 2.0/3.0 * coords2[lOrd[2],1] + 1.0/3.0 * coords2[lOrd[0],1]
            
            pt9.x = 1.0/3.0 * coords2[lOrd[2],0] + 2.0/3.0 * coords2[lOrd[0],0]
            pt9.y = 1.0/3.0 * coords2[lOrd[2],1] + 2.0/3.0 * coords2[lOrd[0],1]
    
            pt6.x = coord_enrich2.x
            pt6.y = coord_enrich2.y     
            
            pt7.x = coord_enrich1.x
            pt7.y = coord_enrich1.y
        
            vec2_x = [coords2[lOrd[0],0], 
                      coords2[lOrd[1],0],
                      coords2[lOrd[2],0],
                      pt4.x, 
                      pt5.x,
                      pt6.x,
                      pt7.x,
                      pt8.x,
                      pt9.x,
                      circumcenter_pt2.x  
                      ]
            vec2_y = [coords2[lOrd[0],1], 
                      coords2[lOrd[1],1], 
                      coords2[lOrd[2],1], 
                      pt4.y, 
                      pt5.y,
                      pt6.y, 
                      pt7.y,
                      pt8.y,
                      pt9.y,
                      circumcenter_pt2.y 
                      ]
            
            [x_fct_2, y_fct_2] = tri_xy_fct_cubic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_cubic( vec2_x, vec2_y )   
                    
    det_J1 = lambda e,n: determinant(J1)(e,n)
    det_J2 = lambda e,n: determinant(J2)(e,n)
    det_J3 = lambda e,n: determinant(J3)(e,n)

    # parent element
    Pe = numpy.zeros((4,4))
    Pe[:,0] = numpy.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]
    C = numpy.linalg.inv(Pe)
    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    # triangle I
    Pe1 = numpy.zeros((3,3))
    Pe1[:,0] = numpy.ones((3,1)).transpose()
    Pe1[:,1] = p[nodes1,0]
    Pe1[:,2] = p[nodes1,1]
    C1 = numpy.linalg.inv(Pe1)
    Nbasis1 = tribasisFct(C1)
    Nx1 = triderivX(C1)
    Ny1 = triderivY(C1)

    # triangle II
    Pe2 = numpy.zeros((3,3))
    Pe2[:,0] = numpy.ones((3,1)).transpose()
    Pe2[:,1] = p[nodes2[0:3],0]
    Pe2[:,2] = p[nodes2[0:3],1]
    C2 = numpy.linalg.inv(Pe2)
    Nbasis2 = tribasisFct(C2)
    Nx2 = triderivX(C2)
    Ny2 = triderivY(C2)

    # triangle III
    Pe3 = numpy.zeros((3,3))
    Pe3[:,0] = numpy.ones((3,1)).transpose()
    Pe3[:,1] = p[nodes3[0:3],0]
    Pe3[:,2] = p[nodes3[0:3],1]
    C3 = numpy.linalg.inv(Pe3)
    Nbasis3 = tribasisFct(C3)
    Nx3 = triderivX(C3)
    Ny3 = triderivY(C3)
    
    node_enr_4 = [nodess[4]]
    coords_enr_4 = p[node_enr_4]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    w1 = y1 - y4
    w2 = y4 - y0
    factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
    if factor_W > EPS_FACTOR:
        factor_W = 1


    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis1[2](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis2[2](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * Nbasis3[0](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx1[2](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx2[2](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * Nx3[0](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny1[2](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny2[2](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * Ny3[0](x,y)]

    for i in range(0,5):
        for j in range(0,5):
            Kefunc1 = lambda e,n: K_cst[0] * (
                                Nx_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Nx_1[j]( x_fct_1(e,n), y_fct_1(e,n)) + 
                                Ny_1[i]( x_fct_1(e,n), y_fct_1(e,n)) *
                                Ny_1[j]( x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)

            Kefunc2 = lambda e,n: K_cst[1] * (
                                Nx_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Nx_2[j]( x_fct_2(e,n), y_fct_2(e,n)) +
                                Ny_2[i]( x_fct_2(e,n), y_fct_2(e,n)) *
                                Ny_2[j]( x_fct_2(e,n), y_fct_2(e,n)) ) * det_J2(e,n)

            Kefunc3 = lambda e,n: K_cst[2] * (
                                Nx_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Nx_3[j]( x_fct_3(e,n), y_fct_3(e,n)) +
                                Ny_3[i]( x_fct_3(e,n), y_fct_3(e,n)) *
                                Ny_3[j]( x_fct_3(e,n), y_fct_3(e,n)) ) * det_J3(e,n)

            integral1 = my_gauss_rule(Kefunc1,ui,wi)
            integral2 = my_gauss_rule(Kefunc2,ui,wi)
            integral3 = my_gauss_rule(Kefunc3,ui,wi)

            K[i,j] = integral1 + integral2 + integral3 

        # construct the local matrix and local components of the load vector
    for i in range(0,5):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

            Fe[i] = my_gauss_rule(fv1,ui,wi) + my_gauss_rule(fv2,ui,wi) + my_gauss_rule(fv3,ui,wi) 

            #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]

def horizontal_cut(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals):

    enrich1 = np.array(p[nodes[4]])
    enrich2 = np.array(p[nodes[5]])

    coords = p[nodes]    
    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    # nodes on the top and bottom side of the interface
    bottom_nodes = [nodes[0],nodes[1],nodes[5],nodes[4]]
    top_nodes = [nodes[4], nodes[5], nodes[2],nodes[3]]

    if (enrich1[0] == x1  and enrich2[0] == x0):
        top_nodes = [nodes[5], nodes[4], nodes[2],nodes[3]]
        bottom_nodes = [nodes[0],nodes[1],nodes[4],nodes[5]]

    bottom_coords = p[bottom_nodes,:]
    top_coords = p[top_nodes,:]


    Pe = np.zeros((4,4))
    Pe[:,0] = np.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    C = np.linalg.inv(Pe)

    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    if ( is_in_same_bin(pxVal1,pxVal2) == True and pxVal1 < binBnd[1] and
         is_in_same_bin(pxVal3,pxVal4) == True and
        (is_in_same_bin(pxVal1,pxVal4) == False and is_in_same_bin(pxVal2,pxVal3) == False) ):
        K_cst = [k2,k1]
    else:
        K_cst = [k1,k2]

    # build the shape functions at the enrichment nodes
    Pe_enr_top = np.zeros((4,4))
    Pe_enr_top[:,0] = np.ones((4,1)).transpose()
    Pe_enr_top[:,1] = p[top_nodes[0:4],0]
    Pe_enr_top[:,2] = p[top_nodes[0:4],1]
    Pe_enr_top[:,3] = p[top_nodes[0:4],0]*p[top_nodes[0:4],1]
    
    C_enr_top = np.linalg.inv(Pe_enr_top)
    
    # top enrichment shape function and its derivative wrt x & y
    Nbasis_enr_top = basisFct(C_enr_top)
    Nx_enr_top = derivX(C_enr_top)
    Ny_enr_top = derivY(C_enr_top)
    
    Pe_enr_bottom = np.zeros((4,4))
    Pe_enr_bottom[:,0] = np.ones((4,1)).transpose()
    Pe_enr_bottom[:,1] = p[bottom_nodes[0:4],0]
    Pe_enr_bottom[:,2] = p[bottom_nodes[0:4],1]
    Pe_enr_bottom[:,3] = p[bottom_nodes[0:4],0]*p[bottom_nodes[0:4],1]
    
    C_enr_bottom = np.linalg.inv(Pe_enr_bottom)            

    # bottom enrichment shape function and its derivatives wrt x & y
    Nbasis_enr_bottom = basisFct(C_enr_bottom)
    Nx_enr_bottom = derivX(C_enr_bottom)
    Ny_enr_bottom = derivY(C_enr_bottom)
    
    # shape functions at enrichment nodes
    psi_left_B = lambda x,y: Nbasis_enr_bottom[3](x,y)  
    psi_left_T = lambda x,y: Nbasis_enr_top[0](x,y) 
    psi_right_B = lambda x,y: Nbasis_enr_bottom[2](x,y)  
    psi_right_T = lambda x,y: Nbasis_enr_top[1](x,y) 
    
    # derivaties wrt x & y of shape functions at enrichment nodes
    # wrt x
    psi_x_left_B = lambda x,y: Nx_enr_bottom[3](x,y) 
    psi_x_left_T = lambda x,y: Nx_enr_top[0](x,y) 
    psi_x_right_B = lambda x,y: Nx_enr_bottom[2](x,y) 
    psi_x_right_T = lambda x,y: Nx_enr_top[1](x,y) 
    # wrt y
    psi_y_left_B = lambda x,y: Ny_enr_bottom[3](x,y) 
    psi_y_left_T = lambda x,y: Ny_enr_top[0](x,y) 
    psi_y_right_B = lambda x,y: Ny_enr_bottom[2](x,y) 
    psi_y_right_T = lambda x,y: Ny_enr_top[1](x,y) 
    
    node_enr_4 = [nodes[4]]
    node_enr_5 = [nodes[5]]
    coords_enr_4 = p[node_enr_4]
    coords_enr_5 = p[node_enr_5]

    if (enrich1[0] == x1  and enrich2[0] == x0):
        node_enr_4 = [nodes[5]]
        node_enr_5 = [nodes[4]]
        coords_enr_4 = p[node_enr_4]
        coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    e1 = y1 - y5
    e2 = y5 - y0
    factor_E = ( 2 * min(e1,e2) / (e1 + e2)) ** 2 
    if factor_E > EPS_FACTOR:
        factor_E = 1

    w1 = y1 - y4
    w2 = y4 - y0
    factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
    if factor_W > EPS_FACTOR:
        factor_W = 1

    Nbasis_bottom = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * psi_left_B(x,y), lambda x,y: factor_E * psi_right_B(x,y)]
    Nbasis_top = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_W * psi_left_T(x,y), lambda x,y: factor_E * psi_right_T(x,y)]
    Nx_top = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * psi_x_left_T(x,y), lambda x,y: factor_E * psi_x_right_T(x,y)]
    Nx_bottom = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_W * psi_x_left_B(x,y), lambda x,y: factor_E * psi_x_right_B(x,y)]
    Ny_top = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * psi_y_left_T(x,y), lambda x,y: factor_E * psi_y_right_T(x,y)]
    Ny_bottom = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_W * psi_y_left_B(x,y), lambda x,y: factor_E * psi_y_right_B(x,y)] 
    
    
    # getting x and y coordinates transformed into the [-1,1] and [-1,1] intervals
    # of the isoparametric element
    # x,y = inner products of basis functions of isoparametric elements and actual
    # coordinates of the corners of the child element
    # x = [N1_e, N2_e, N3_e, N4_e ] * [x0, x1, x2, x3 ]'
    # y = [N1_e, N2_e, N3_e, N4_e ] * [y0, y1, y2, y3 ]'
    
    if len(root.enrichNodes) == 2:
        
        [x_fct_T, y_fct_T] = xy_fct( top_coords[:,0], top_coords[:,1] )
        [x_fct_B, y_fct_B] = xy_fct( bottom_coords[:,0], bottom_coords[:,1] )
    
        # computing the Jacobian and the determinant of the left and right children of the parent element
        J_top = jacobian_mat( top_coords[:,0], top_coords[:,1] )
        J_bottom = jacobian_mat( bottom_coords[:,0], bottom_coords[:,1] )
        
    if len(root.enrichNodes) == 3:
       
        coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
        
        lOrd = [0,1,2,3] # local order    
        vecB_x = [ bottom_coords[ lOrd[0],0], 
                  (bottom_coords[ lOrd[0],0] + bottom_coords[lOrd[1],0])/2.0,
                  bottom_coords[lOrd[1],0],
                  (bottom_coords[lOrd[1],0] + bottom_coords[lOrd[2],0])/2.0,
                  bottom_coords[lOrd[2],0], 
                  coord_enrich.x, 
                  bottom_coords[lOrd[3],0],
                  (bottom_coords[lOrd[3],0] + bottom_coords[lOrd[0],0])/2.0
                 ]
        vecB_y = [ bottom_coords[ lOrd[0],1], 
                  (bottom_coords[ lOrd[0],1] + bottom_coords[lOrd[1],1])/2.0,
                  bottom_coords[lOrd[1],1],
                  (bottom_coords[lOrd[1],1] + bottom_coords[lOrd[2],1])/2.0,
                  bottom_coords[lOrd[2],1], 
                  coord_enrich.y, 
                  bottom_coords[lOrd[3],1],
                  (bottom_coords[lOrd[3],1] + bottom_coords[lOrd[0],1])/2.0
                 ]

        [x_fct_B, y_fct_B] = quad_xy_fct_bi_quadratic( vecB_x, vecB_y )
        J_bottom = quad_jacobian_mat_bi_quadratic( vecB_x, vecB_y )
        
        lOrd = [2,3,0,1]
        vecT_x = [ top_coords[ lOrd[0],0], 
                  (top_coords[ lOrd[0],0] + top_coords[lOrd[1],0])/2.0,
                  top_coords[lOrd[1],0],
                  (top_coords[lOrd[1],0] + top_coords[lOrd[2],0])/2.0,
                  top_coords[lOrd[2],0], 
                  coord_enrich.x, 
                  top_coords[lOrd[3],0],
                  (top_coords[lOrd[3],0] + top_coords[lOrd[0],0])/2.0
                 ]
        vecT_y = [ top_coords[ lOrd[0],1], 
                  (top_coords[ lOrd[0],1] + top_coords[lOrd[1],1])/2.0,
                  top_coords[lOrd[1],1],
                  (top_coords[lOrd[1],1] + top_coords[lOrd[2],1])/2.0,
                  top_coords[lOrd[2],1], 
                  coord_enrich.y, 
                  top_coords[lOrd[3],1],
                  (top_coords[lOrd[3],1] + top_coords[lOrd[0],1])/2.0
                 ]

        [x_fct_T, y_fct_T] = quad_xy_fct_bi_quadratic( vecT_x, vecT_y )
        J_top = quad_jacobian_mat_bi_quadratic( vecT_x, vecT_y )
        
    if len(root.enrichNodes) == 4:
        
        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
        
        lOrd = [0,1,2,3] # local order    
        
        vecB_x = [ bottom_coords[ lOrd[0],0], 
                  2.0/3.0 * bottom_coords[lOrd[0],0] + 1.0/3.0 * bottom_coords[lOrd[1],0],
                  1.0/3.0 * bottom_coords[lOrd[0],0] + 2.0/3.0 * bottom_coords[lOrd[1],0],
                  bottom_coords[lOrd[1],0],
                  2.0/3.0 * bottom_coords[lOrd[1],0] + 1.0/3.0 * bottom_coords[lOrd[2],0],
                  1.0/3.0 * bottom_coords[lOrd[1],0] + 2.0/3.0 * bottom_coords[lOrd[2],0],
                  bottom_coords[lOrd[2],0], 
                  coord_enrich2.x,
                  coord_enrich1.x,
                  bottom_coords[lOrd[3],0],
                  2.0/3.0 * bottom_coords[lOrd[3],0] + 1.0/3.0 * bottom_coords[lOrd[0],0],
                  1.0/3.0 * bottom_coords[lOrd[3],0] + 2.0/3.0 * bottom_coords[lOrd[0],0]
                 ]
        vecB_y = [ bottom_coords[ lOrd[0],1], 
                  2.0/3.0 * bottom_coords[lOrd[0],1] + 1.0/3.0 * bottom_coords[lOrd[1],1],
                  1.0/3.0 * bottom_coords[lOrd[0],1] + 2.0/3.0 * bottom_coords[lOrd[1],1],
                  bottom_coords[lOrd[1],1],
                  2.0/3.0 * bottom_coords[lOrd[1],1] + 1.0/3.0 * bottom_coords[lOrd[2],1],
                  1.0/3.0 * bottom_coords[lOrd[1],1] + 2.0/3.0 * bottom_coords[lOrd[2],1],
                  bottom_coords[lOrd[2],1], 
                  coord_enrich2.y,
                  coord_enrich1.y,
                  bottom_coords[lOrd[3],1],
                  2.0/3.0 * bottom_coords[lOrd[3],1] + 1.0/3.0 * bottom_coords[lOrd[0],1],
                  1.0/3.0 * bottom_coords[lOrd[3],1] + 2.0/3.0 * bottom_coords[lOrd[0],1]
                 ]

        [x_fct_B, y_fct_B] = quad_xy_fct_bi_cubic( vecB_x, vecB_y )
        J_bottom = quad_jacobian_mat_bi_cubic( vecB_x, vecB_y )
        
        lOrd = [2,3,0,1]
        vecT_x = [ top_coords[ lOrd[0],0], 
                  2.0/3.0 * top_coords[lOrd[0],0] + 1.0/3.0 * top_coords[lOrd[1],0],
                  1.0/3.0 * top_coords[lOrd[0],0] + 2.0/3.0 * top_coords[lOrd[1],0],
                  top_coords[lOrd[1],0],
                  2.0/3.0 * top_coords[lOrd[1],0] + 1.0/3.0 * top_coords[lOrd[2],0],
                  1.0/3.0 * top_coords[lOrd[1],0] + 2.0/3.0 * top_coords[lOrd[2],0],
                  top_coords[lOrd[2],0], 
                  coord_enrich1.x,
                  coord_enrich2.x,
                  top_coords[lOrd[3],0],
                  2.0/3.0 * top_coords[lOrd[3],0] + 1.0/3.0 * top_coords[lOrd[0],0],
                  1.0/3.0 * top_coords[lOrd[3],0] + 2.0/3.0 * top_coords[lOrd[0],0]
                 ]
        vecT_y = [ top_coords[ lOrd[0],1], 
                  2.0/3.0 * top_coords[lOrd[0],1] + 1.0/3.0 * top_coords[lOrd[1],1],
                  1.0/3.0 * top_coords[lOrd[0],1] + 2.0/3.0 * top_coords[lOrd[1],1],
                  top_coords[lOrd[1],1],
                  2.0/3.0 * top_coords[lOrd[1],1] + 1.0/3.0 * top_coords[lOrd[2],1],
                  1.0/3.0 * top_coords[lOrd[1],1] + 2.0/3.0 * top_coords[lOrd[2],1],
                  top_coords[lOrd[2],1], 
                  coord_enrich1.y,
                  coord_enrich2.y,
                  top_coords[lOrd[3],1],
                  2.0/3.0 * top_coords[lOrd[3],1] + 1.0/3.0 * top_coords[lOrd[0],1],
                  1.0/3.0 * top_coords[lOrd[3],1] + 2.0/3.0 * top_coords[lOrd[0],1]
                ]

        [x_fct_T, y_fct_T] = quad_xy_fct_bi_quadratic( vecT_x, vecT_y )
        J_top = quad_jacobian_mat_bi_quadratic( vecT_x, vecT_y )        
        
        
    detJ_top = lambda e,n: determinant(J_top)(e,n)
    detJ_bottom = lambda e,n: determinant(J_bottom)(e,n)
    
    Ke_Horiz = np.zeros((6,6))
    Fe_Horiz = np.zeros((6,1))
    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        for j in range(0,6):
            if i>=j:
                Kefunc1 = lambda e,n: K_cst[0] * ( Nx_top[i](x_fct_T(e,n),y_fct_T(e,n)) * 
                                                    Nx_top[j](x_fct_T(e,n),y_fct_T(e,n)) + 
                                                    Ny_top[i](x_fct_T(e,n),y_fct_T(e,n)) * 
                                                    Ny_top[j](x_fct_T(e,n),y_fct_T(e,n)) ) * detJ_top(e,n) 
                Kefunc2 = lambda e,n: K_cst[1] * ( Nx_bottom[i](x_fct_B(e,n),y_fct_B(e,n)) * 
                                                    Nx_bottom[j](x_fct_B(e,n),y_fct_B(e,n)) + 
                                                    Ny_bottom[i](x_fct_B(e,n),y_fct_B(e,n)) * 
                                                    Ny_bottom[j](x_fct_B(e,n),y_fct_B(e,n)) ) * detJ_bottom(e,n)
                top_integral = quad2d(Kefunc1,-1,1,-1,1,ui,wi) 
                bottom_integral = quad2d(Kefunc2,-1,1,-1,1,ui,wi)

                Ke_Horiz[i,j] = top_integral + bottom_integral
                Ke_Horiz[j,i] = Ke_Horiz[i,j]
        # construct the local load vector
        fv1 = lambda e,n: rhs(e,n) * Nbasis_top[i](x_fct_T(e,n),y_fct_T(e,n)) * detJ_top(e,n)
        fv2 = lambda e,n: rhs(e,n) * Nbasis_bottom[i](x_fct_B(e,n),y_fct_B(e,n)) * detJ_bottom(e,n)

        Fe_Horiz[i] = quad2d(fv1,-1,1,-1,1,ui,wi) + quad2d(fv2,-1,1,-1,1,ui,wi)

    return [Ke_Horiz,Fe_Horiz]


def vertical_cut(p,ui,wi,k1,k2,nodes,root,image,nodes6,nodes7,pxVals):
    enrich1 = np.array(p[nodes[4]])
    enrich2 = np.array(p[nodes[5]])

    coords = p[nodes]    
    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

    # nodes on the left side of the interface
    left_nodes = [nodes[0],nodes[4],nodes[5],nodes[3]]
    # nodes on the right side of the interface
    right_nodes = [nodes[4],nodes[1],nodes[2],nodes[5]]

    if (enrich1[1] == y1 and enrich2[1] == y0):
        left_nodes = [nodes[0],nodes[5],nodes[4],nodes[3]]
        right_nodes = [nodes[5],nodes[1],nodes[2],nodes[4]]

    # coordinates of the left sub-element or sub-element number 1
    left_coords = p[left_nodes,:]                                                
    # coordinates of the right sub-element or the sub-element number 2
    right_coords = p[right_nodes,:]                                                                
    
    Pe = np.zeros((4,4))
    Pe[:,0] = np.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    C = np.linalg.inv(Pe)

    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

    pxVal1 = pxVals[0]
    pxVal2 = pxVals[1]
    pxVal3 = pxVals[2]
    pxVal4 = pxVals[3]
    
    
    if ( is_in_same_bin(pxVal1,pxVal4) == True and is_in_same_bin(pxVal2,pxVal3) == True and pxVal4 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)== False and is_in_same_bin(pxVal4,pxVal3) == False) ):
        K_cst = [k2,k1]
    else:
        K_cst = [k1,k2]
        
    # build the shape functions at the enrichment nodes
    Pe_enr_left = np.zeros((4,4))
    Pe_enr_left[:,0] = np.ones((4,1)).transpose()
    Pe_enr_left[:,1] = p[left_nodes[0:4],0]
    Pe_enr_left[:,2] = p[left_nodes[0:4],1]
    Pe_enr_left[:,3] = p[left_nodes[0:4],0]*p[left_nodes[0:4],1]
    
    C_enr_left = np.linalg.inv(Pe_enr_left)
    
    # left enrichment shape function and its derivative wrt x & y
    Nbasis_enr_left = basisFct(C_enr_left)
    Nx_enr_left = derivX(C_enr_left)
    Ny_enr_left = derivY(C_enr_left)
    
    Pe_enr_right = np.zeros((4,4))
    Pe_enr_right[:,0] = np.ones((4,1)).transpose()
    
    Pe_enr_right[:,1] = p[right_nodes[0:4],0]
    Pe_enr_right[:,2] = p[right_nodes[0:4],1]
    Pe_enr_right[:,3] = p[right_nodes[0:4],0]*p[right_nodes[0:4],1]
    
    C_enr_right = np.linalg.inv(Pe_enr_right)            
    
    # right enrichment shape function and its derivatives wrt x & y
    Nbasis_enr_right = basisFct(C_enr_right)
    Nx_enr_right = derivX(C_enr_right)
    Ny_enr_right = derivY(C_enr_right)
    
    # shape functions at enrichment nodes
    psi_btm_L = lambda x,y: Nbasis_enr_left[1](x,y)  
    psi_btm_R = lambda x,y: Nbasis_enr_right[0](x,y) 
    psi_upr_L = lambda x,y: Nbasis_enr_left[2](x,y)  
    psi_upr_R = lambda x,y: Nbasis_enr_right[3](x,y) 
    
    # derivaties wrt x & y of shape functions at enrichment nodes
    # wrt x
    psi_x_btm_L = lambda x,y: Nx_enr_left[1](x,y) 
    psi_x_btm_R = lambda x,y: Nx_enr_right[0](x,y) 
    psi_x_upr_L = lambda x,y: Nx_enr_left[2](x,y) 
    psi_x_upr_R = lambda x,y: Nx_enr_right[3](x,y) 
    # wrt y
    psi_y_btm_L = lambda x,y: Ny_enr_left[1](x,y) 
    psi_y_btm_R = lambda x,y: Ny_enr_right[0](x,y) 
    psi_y_upr_L = lambda x,y: Ny_enr_left[2](x,y) 
    psi_y_upr_R = lambda x,y: Ny_enr_right[3](x,y) 
    
    node_enr_4 = [nodes[4]]
    node_enr_5 = [nodes[5]]
    coords_enr_4 = p[node_enr_4]
    coords_enr_5 = p[node_enr_5]

    if (enrich1[1] == y1 and enrich2[1] == y0):
        node_enr_4 = [nodes[5]]
        node_enr_5 = [nodes[4]]
        coords_enr_4 = p[node_enr_4]
        coords_enr_5 = p[node_enr_5]

    x4 = coords_enr_4[0,0]
    y4 = coords_enr_4[0,1]

    x5 = coords_enr_5[0,0]
    y5 = coords_enr_5[0,1]

    s1 = x4 - x0
    s2 = x1 - x4
    factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
    if factor_S > EPS_FACTOR:
        factor_S = 1

    n1 = x5 - x0
    n2 = x1 - x5
    factor_N = ( 2 * min(n1,n2) / (n1 + n2)) ** 2 
    if factor_N > EPS_FACTOR:
        factor_N = 1

    Nbasis_left = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * psi_btm_L(x,y), lambda x,y: factor_N * psi_upr_L(x,y)]
    Nbasis_right = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: factor_S * psi_btm_R(x,y), lambda x,y: factor_N * psi_upr_R(x,y)]
    Nx_left = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * psi_x_btm_L(x,y), lambda x,y: factor_N * psi_x_upr_L(x,y)]
    Nx_right = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * psi_x_btm_R(x,y), lambda x,y: factor_N * psi_x_upr_R(x,y)]
    Ny_left = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * psi_y_btm_L(x,y), lambda x,y: factor_N * psi_y_upr_L(x,y)]
    Ny_right = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * psi_y_btm_R(x,y), lambda x,y: factor_N * psi_y_upr_R(x,y)] 
    
    
    # getting x and y coordinates transformed into the [-1,1] and [-1,1] intervals
    # of the isoparametric element
    # x,y = inner products of basis functions of isoparametric elements and actual
    # coordinates of the corners of the child element
    # x = [N1_e, N2_e, N3_e, N4_e ] * [x0, x1, x2, x3 ]'
    # y = [N1_e, N2_e, N3_e, N4_e ] * [y0, y1, y2, y3 ]'
    

    if len(root.enrichNodes) == 2:
        
        [x_fct_L, y_fct_L] = xy_fct( left_coords[:,0], left_coords[:,1] )
        [x_fct_R, y_fct_R] = xy_fct( right_coords[:,0], right_coords[:,1] )
        
        # computing the Jacobian and the determinant of the left and right children of the parent element
        J_left = jacobian_mat( left_coords[:,0], left_coords[:,1] )
        J_right = jacobian_mat( right_coords[:,0], right_coords[:,1] )
        
    if len(root.enrichNodes) == 3:
       
        coord_enrich = coord_enrich_comp_quad_circle(p[nodes], nodes6)
        
        lOrd = [3,0,1,2] # local order    
        vecL_x = [ left_coords[ lOrd[0],0], 
                  (left_coords[ lOrd[0],0] + left_coords[lOrd[1],0])/2.0,
                  left_coords[lOrd[1],0],
                  (left_coords[lOrd[1],0] + left_coords[lOrd[2],0])/2.0,
                  left_coords[lOrd[2],0], 
                  coord_enrich.x, 
                  left_coords[lOrd[3],0],
                  (left_coords[lOrd[3],0] + left_coords[lOrd[0],0])/2.0
                 ]
        vecL_y = [ left_coords[ lOrd[0],1], 
                  (left_coords[ lOrd[0],1] + left_coords[lOrd[1],1])/2.0,
                  left_coords[lOrd[1],1],
                  (left_coords[lOrd[1],1] + left_coords[lOrd[2],1])/2.0,
                  left_coords[lOrd[2],1], 
                  coord_enrich.y, 
                  left_coords[lOrd[3],1],
                  (left_coords[lOrd[3],1] + left_coords[lOrd[0],1])/2.0
                 ]

        
        [x_fct_L, y_fct_L] = quad_xy_fct_bi_quadratic( vecL_x, vecL_y )
        J_left = quad_jacobian_mat_bi_quadratic( vecL_x, vecL_y )
        
        lOrd = [1,2,3,0]
        vecR_x = [ right_coords[ lOrd[0],0], 
                  (right_coords[ lOrd[0],0] + right_coords[lOrd[1],0])/2.0,
                  right_coords[lOrd[1],0],
                  (right_coords[lOrd[1],0] + right_coords[lOrd[2],0])/2.0,
                  right_coords[lOrd[2],0], 
                  coord_enrich.x, 
                  right_coords[lOrd[3],0],
                  (right_coords[lOrd[3],0] + right_coords[lOrd[0],0])/2.0
                 ]
        vecR_y = [ right_coords[ lOrd[0],1], 
                  (right_coords[ lOrd[0],1] + right_coords[lOrd[1],1])/2.0,
                  right_coords[lOrd[1],1],
                  (right_coords[lOrd[1],1] + right_coords[lOrd[2],1])/2.0,
                  right_coords[lOrd[2],1], 
                  coord_enrich.y, 
                  right_coords[lOrd[3],1],
                  (right_coords[lOrd[3],1] + right_coords[lOrd[0],1])/2.0
                 ]
        
        [x_fct_R, y_fct_R] = quad_xy_fct_bi_quadratic( vecR_x, vecR_y )
        J_right = quad_jacobian_mat_bi_quadratic( vecR_x, vecR_y )

    if len(root.enrichNodes) == 4:
        
        coord_enrich1 = coord_enrich_comp_quad_circle(p[nodes],nodes6)
        coord_enrich2 = coord_enrich_comp_quad_circle(p[nodes],nodes7)
        
        lOrd = [3,0,1,2] # local order    
        vecL_x = [ left_coords[lOrd[0],0], 
                  2.0/3.0 * left_coords[lOrd[0],0] + 1.0/3.0 * left_coords[lOrd[1],0],
                  1.0/3.0 * left_coords[lOrd[0],0] + 2.0/3.0 * left_coords[lOrd[1],0],
                  left_coords[lOrd[1],0],
                  2.0/3.0 * left_coords[lOrd[1],0] + 1.0/3.0 * left_coords[lOrd[2],0],
                  1.0/3.0 * left_coords[lOrd[1],0] + 2.0/3.0 * left_coords[lOrd[2],0],
                  left_coords[lOrd[2],0], 
                  coord_enrich1.x,
                  coord_enrich2.x,
                  left_coords[lOrd[3],0],
                  2.0/3.0 * left_coords[lOrd[3],0] + 1.0/3.0 * left_coords[lOrd[0],0],
                  1.0/3.0 * left_coords[lOrd[3],0] + 2.0/3.0 * left_coords[lOrd[0],0]
                 ]
        vecL_y = [ left_coords[lOrd[0],1], 
                  2.0/3.0 * left_coords[lOrd[0],1] + 1.0/3.0 * left_coords[lOrd[1],1],
                  1.0/3.0 * left_coords[lOrd[0],1] + 2.0/3.0 * left_coords[lOrd[1],1],
                  left_coords[lOrd[1],1],
                  2.0/3.0 * left_coords[lOrd[1],1] + 1.0/3.0 * left_coords[lOrd[2],1],
                  1.0/3.0 * left_coords[lOrd[1],1] + 2.0/3.0 * left_coords[lOrd[2],1],
                  left_coords[lOrd[2],1], 
                  coord_enrich1.y,
                  coord_enrich2.y,
                  left_coords[lOrd[3],1],
                  2.0/3.0 * left_coords[lOrd[3],1] + 1.0/3.0 * left_coords[lOrd[0],1],
                  1.0/3.0 * left_coords[lOrd[3],1] + 2.0/3.0 * left_coords[lOrd[0],1]
                ]

        [x_fct_L, y_fct_L] = quad_xy_fct_bi_cubic( vecL_x, vecL_y )
        J_left = quad_jacobian_mat_bi_cubic( vecL_x, vecL_y )
        
        lOrd = [1,2,3,0]

        vecR_x = [ right_coords[ lOrd[0],0], 
                  2.0/3.0 * right_coords[lOrd[0],0] + 1.0/3.0 * right_coords[lOrd[1],0],
                  1.0/3.0 * right_coords[lOrd[0],0] + 2.0/3.0 * right_coords[lOrd[1],0],
                  right_coords[lOrd[1],0],
                  2.0/3.0 * right_coords[lOrd[1],0] + 1.0/3.0 * right_coords[lOrd[2],0],
                  1.0/3.0 * right_coords[lOrd[1],0] + 2.0/3.0 * right_coords[lOrd[2],0],
                  right_coords[lOrd[2],0], 
                  coord_enrich2.x,
                  coord_enrich1.x,
                  right_coords[lOrd[3],0],
                  2.0/3.0 * right_coords[lOrd[3],0] + 1.0/3.0 * right_coords[lOrd[0],0],
                  1.0/3.0 * right_coords[lOrd[3],0] + 2.0/3.0 * right_coords[lOrd[0],0]
                ]
        vecR_y = [ right_coords[ lOrd[0],1], 
                  2.0/3.0 * right_coords[lOrd[0],1] + 1.0/3.0 * right_coords[lOrd[1],1],
                  1.0/3.0 * right_coords[lOrd[0],1] + 2.0/3.0 * right_coords[lOrd[1],1],
                  right_coords[lOrd[1],1],
                  2.0/3.0 * right_coords[lOrd[1],1] + 1.0/3.0 * right_coords[lOrd[2],1],
                  1.0/3.0 * right_coords[lOrd[1],1] + 2.0/3.0 * right_coords[lOrd[2],1],
                  right_coords[lOrd[2],1], 
                  coord_enrich2.y,
                  coord_enrich1.y,
                  right_coords[lOrd[3],1],
                  2.0/3.0 * right_coords[lOrd[3],1] + 1.0/3.0 * right_coords[lOrd[0],1],
                  1.0/3.0 * right_coords[lOrd[3],1] + 2.0/3.0 * right_coords[lOrd[0],1]
                  ]
        
        [x_fct_R, y_fct_R] = quad_xy_fct_bi_cubic( vecR_x, vecR_y )
        J_right = quad_jacobian_mat_bi_cubic( vecR_x, vecR_y )

    detJ_left = lambda e,n: determinant(J_left)(e,n)
    detJ_right = lambda e,n: determinant(J_right)(e,n)
    
    
    Ke_Vertical = np.zeros((6,6))
    Fe_Vertical = np.zeros((6,1))
    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        for j in range(0,6):
            if i>=j :
                Kefunc1 = lambda e,n: K_cst[0] * ( Nx_left[i](x_fct_L(e,n),y_fct_L(e,n)) * 
                                                Nx_left[j](x_fct_L(e,n),y_fct_L(e,n)) + 
                                                Ny_left[i](x_fct_L(e,n),y_fct_L(e,n)) * 
                                                Ny_left[j](x_fct_L(e,n),y_fct_L(e,n)) ) * detJ_left(e,n) 
                Kefunc2 = lambda e,n: K_cst[1] * ( Nx_right[i](x_fct_R(e,n),y_fct_R(e,n)) * 
                                                Nx_right[j](x_fct_R(e,n),y_fct_R(e,n)) + 
                                                Ny_right[i](x_fct_R(e,n),y_fct_R(e,n)) * 
                                                Ny_right[j](x_fct_R(e,n),y_fct_R(e,n)) ) * detJ_right(e,n)
                left_integral = quad2d(Kefunc1,-1,1,-1,1,ui,wi) 
                right_integral = quad2d(Kefunc2,-1,1,-1,1,ui,wi)
    
                Ke_Vertical[i,j] = left_integral + right_integral
                Ke_Vertical[j,i] = Ke_Vertical[i,j]
                

        # construct the local load vector
        fv1 = lambda e,n: rhs(e,n) * Nbasis_left[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
        fv2 = lambda e,n: rhs(e,n) * Nbasis_right[i](x_fct_R(e,n),y_fct_R(e,n)) * detJ_right(e,n)
    
        Fe_Vertical[i] = quad2d(fv1,-1,1,-1,1,ui,wi) + quad2d(fv2,-1,1,-1,1,ui,wi)
    
    ## If Neumann BCs are not zero, uncomment this code accordingly
    #    # Neumann BCs: n dot grad(u) = g1
    #    # left side boundary
    #    if leftNeumann == 1 and  (nodes[i] in lbcs):
    #        int_neumn_left = lambda e,n: g1_left * Nbasis_left[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
    #        Fe_Vertical[i] = Fe_Vertical[i] + quad2d(int_neumn_left,-1,1,-1,1,ui,wi)
    #
    #    # right side boundary
    #    if rightNeumann == 1 and (nodes[i] in rbcs):
    #        int_neumn_right = lambda e,n: g1_right * Nbasis_right[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
    #        Fe_Vertical[i] = Fe_Vertical[i] + quad2d(int_neumn_right,-1,1,-1,1,ui,wi)
    #            
    #    # top boundary
    #    if topNeumann == 1 and (nodes[i] in tbcs):
    #        int_neumn_top = lambda e,n: g1_top * Nbasis_right[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
    #        Fe_Vertical[i] = Fe_Vertical[i] + quad2d(int_neumn_top,-1,1,-1,1,ui,wi)
    #
    #    # bottom boundary
    #    if bottomNeumann == 1 and (nodes[i] in bbcs):
    #        int_neumn_bottom = lambda e,n: g1_bottom * Nbasis_left[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
    #        Fe_Vertical[i] = Fe_Vertical[i] + quad2d(int_neumn_bottom,-1,1,-1,1,ui,wi)
                        
    
    return [Ke_Vertical,Fe_Vertical]

def on_corners(enrich,coords):
    x0 = coords[0,0]
    y0 = coords[0,1]
    x1 = coords[1,0]
    y1 = coords[1,1]
    x2 = coords[2,0]
    y2 = coords[2,1]
    x3 = coords[3,0]
    y3 = coords[3,1]
    
    if enrich[0] == x0 and enrich[1] == y0:
        return True

    if enrich[0] == x1 and enrich[1] == y1:
        return True
    
    if enrich[0] == x2 and enrich[1] == y2:
        return True

    if enrich[0] == x3 and enrich[1] == y3:
        return True

    return False

def on_corners1(enrich,x0,y0,x1,y1):
    if enrich[0] == x0 and enrich[1] == y0:
        return True

    if enrich[0] == x1 and enrich[1] == y0:
        return True
    
    if enrich[0] == x0 and enrich[1] == y1:
        return True

    if enrich[0] == x1 and enrich[1] == y1:
        return True

    return False



def local_quad(ui,wi,f):
# http://www.ece.ualberta.ca/~knight/ece632/fea/triangles/num_ex.html
    I = 0 
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            I = I + wi[i] * wi[j] * f(ui[i],ui[j])

    return I * 1.0/2.0 

def in_circle(center_x, center_y, radius, x, y):
    square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
    return square_dist <= radius ** 2

     
def in_same_domain(nd1,nd2,nd3):
#checking to see if midpoint along edges created by nd1, nd2, and nd3 are
# in the same domain

    midptH = Point(0.0,0.0)
    midptV = Point(0.0,0.0)

    if nd1.x == nd2.x:
        midptV.x = nd1.x
        midptV.y = (nd1.y + nd2.y) / 2.0
    
    if nd1.x == nd3.x:
        midptV.x = nd1.x
        midptV.y = (nd1.y + nd3.y) / 2.0

    if nd2.x == nd3.x:
        midptV.x = nd2.x
        midptV.y = (nd2.y + nd3.y) / 2.0

    if nd1.y == nd2.y:
        midptH.x = (nd1.x + nd2.x) / 2.0
        midptH.y = nd1.y

    if nd1.y == nd3.y:
        midptH.x = (nd1.x + nd3.x) / 2.0
        midptH.y = nd1.y

    if nd2.y == nd3.y:
        midptH.x = (nd2.x + nd3.x) / 2.0
        midptH.y = nd2.y

    ptH_in = point_in_on_poly(midptH.x,midptH.y,polygonDef) 
    ptV_in = point_in_on_poly(midptV.x,midptV.y,polygonDef)
    
    # if both points are in or out (but both of them) return true, else false
    if (ptH_in == True and ptV_in == True) or (ptH_in == False and ptV_in == False):
        return True

    return False
