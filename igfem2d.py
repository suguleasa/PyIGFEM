from scipy import  sparse

import scipy.sparse.linalg.dsolve as linsolve
#import matplotlib
#import pylab
import math
from quad2d import quad2d
from isoparam import *
# import testclosure
from probstatement import *
# from caseStudies import caseA,caseB,caseC,caseD
import numpy as np
from findIntersection import *
import itertools
import sympy
from quad2d import *
# import pylab
import scipy.io

from libFcts import *

from main import *
from scipy import interpolate

import __builtin__
sum = __builtin__.sum

#definition of rhombus corners
#A = Point(0.5, 1.0/3.0)
#B = Point(2.0/3.0, 0.5)
#C = Point(1.0/3.0, 0.5)
#D = Point(0.5, 2.0/3.0)
#rhombus = [ (A.x, A.y), (B.x,B.y), (D.x,D.y), (C.x,C.y) ]
# lower rhombus
#A    = Point(1.0/6.0, 0.0)
#B = Point(2.0/3.0, 0.0)
##C = Point(3.0/4.0+1.0/6.0, 1.0/2.0)
##D = Point(1.0/4.0+1.0/6.0, 1.0/2.0)
#C = Point(3.0/4.0, 1.0/2.0)
#D = Point(1.0/4.0, 1.0/2.0)
#A = Point(0.5, 0.0)
#B = Point(2.0/3.0, 0.5)
#C = Point(0.5, 2.0/3.0)
#D = Point(1.0/3.0, 0.5)

#rhombus = [ (A.x, A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y) ]


#definition of hexagon corners
#A = Point(1.0/4.0, 1.0/6.0)
#B = Point(3.0/4.0, 1.0/6.0)
#C = Point(5.0/6.0, 1.0/2.0)
#D = Point(3.0/4.0, 5.0/6.0)
#E = Point(1.0/4.0, 5.0/6.0)
#F = Point(1.0/6.0, 1.0/2.0)
#hexagon = [ (A.x,A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y), (E.x,E.y), (F.x,F.y) ]


EPS_FACTOR = 1e-2


#class Point:
#    def __init__(self, x, y):
#        self.x = x
#        self.y = y

# Use Simpson Rule for integration over isoparametric elements
def simpson_rule(f):
# Gauss Rule with 3 points: exact for polynomial of degree 2

  I = 1.0/6.0 * f(1.0/2.0,1.0/2.0) + 1.0/6.0 * f(1.0/2.0,0.0) + 1.0/6.0 * f(0.0,1.0/2.0)
  return I


def simpson_rule4(f):
# Gauss Rule with 4 points: exact for polynomial of degree 3
    I = -27.0/96.0 * f(1.0/3.0,1.0/3.0) + 25.0/96.0 * f(0.2,0.6) + 25.0/96.0 * f(0.2, 0.2) + 25.0/96.0 * f(0.6, 0.2)
    return I

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


def gauss_integration_HN(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem_HN,detJ_HN):
    J = 0
    L = ui.size
    for i in range(0,L):
        for j in range(0,L):
            J = J + wi[i] * wi[j] * diff_fct_triangles_HN(ui[i],ui[j],x_trans_fct,y_trans_fct,pConf,tConf,UConf,uh_elem_HN,detJ_HN)
    return J

def gauss_integration(ui,wi,UConf,pConf,tConf,x_trans_fct,y_trans_fct,uh_elem,detJ):
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

def myquad(m,n,k1,k2,ui,wi,p,t,masterNode,llist,image):
#def myquad(m,n,k1,k2,loc,ui,wi,p,t,UConf,pConf,tConf):
    
    #definition of rhombus corners
#    A = Point(0.5, 1.0/3.0)
#    B = Point(2.0/3.0, 0.5)
#    C = Point(1.0/3.0, 0.5)
#    D = Point(0.5, 2.0/3.0)
#    rhombus = [ (A.x, A.y), (B.x,B.y), (D.x,D.y), (C.x,C.y) ]
    
    #A = Point(1.0/6.0, 0.0)
    #B = Point(2.0/3.0, 0.0)
    #C = Point(3.0/4.0+1.0/6.0, 1.0/2.0)
    #D = Point(1.0/4.0+1.0/6.0, 1.0/2.0)

    #rhombus = [ (A.x, A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y) ]


    #print "sizes:",m,n
    # reading data
    #filename = 'vertical2by2.txt'
    #lines = [line.strip() for line in open(filename)]

    # numbering boundary nodes
#    lbcs = [0] + range(m,m*n,m)
#    rbcs = [m-1] + range(2*m-1,m*n,m)
#    tbcs = range(1,m-1)
#    bbcs = range(lbcs[-1]+1,rbcs[-1]) 
    
    tbcs = []
    bbcs = []
    lbcs = []
    rbcs = []
    for i in range(0, len(p)):
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
                            
#     print 'tbcs', tbcs
#     print 'lbcs', lbcs
#     print 'rbcs', rbcs
#     print 'bbcs', bbcs

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

    # vector with boundary nodes: bottom, left, right, top
#    b = range(0,m) + range(m,m*n,m) + range(2*m-1,m*n,m) + range(m*n-m+2-1,m*n-1) 

    N = len(p)  #+ 5 # number of nodes
    T = len(t) #+ 3 # number of elements

#    p = numpy.vstack([p,[ 0.8125, 0.5 ]])
#    p = numpy.vstack([p,[ 0.75, 0.5625 ]])
#    p = numpy.vstack([p,[ 0.8125, 0.5625 ]])
#    p = numpy.vstack([p,[ 0.875, 0.5625 ]])
#    p = numpy.vstack([p,[ 0.8125, 0.6125 ]])

#    for e in range(0,T):
#        nodes = t[e]
#        if nodes[0] == 42 and nodes[1] == 43:
#            t =  t[0:e] + t[e+1:len(t)]
#            t = t + [[42, 111, 113, 112]]
#            t = t + [[111, 43, 114, 113]]
#            t = t + [[112, 113, 115, 51]]
#            t = t + [[113, 114, 52, 115]]

    K = sparse.lil_matrix((N,N))
    F = sparse.lil_matrix((N,1))

    list_hanging_nodes = []
    for e in range(0,T): #800, 833
#     for e in range(820, T):
#     for e in range(833, T):

#     for e in range(0,14):
        
        nodes = t[e] # row of t =  node numbers of the 4 corners of element e
        
        root = get_node_by_id(masterNode,llist[e])
        
        
        p1,p2,p3,p4 = root.rect
    
        pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
        pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
        pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
        pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    

        # 2-column matrix containing on each row the coordinates of each of the nodes
        coords = p[nodes,:]    

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

#        if root.index == '30100':#'210333':
#            print '-----------------------------------'
#            print p1.x,p1.y,p2.x,p2.y,p3.x,p3.y, p4.x,p4.y
#            enrN1 = root.enrichNodes[0]
#            enrN2 = root.enrichNodes[1]
#            print 'nodes', nodes, 'e =',e
#            print enrN1.x,enrN1.y, enrN2.x,enrN2.y
#            print root.ishomog
#            print '-----------------------------------'
              
#         root_i = get_node_by_id(masterNode,['30100'])
#         root_i.printRect()
#         print root_i.enrichNodes[0].x, root_i.enrichNodes[0].y
#         print root_i.enrichNodes[1].x, root_i.enrichNodes[1].y
#         print root_i.ishomog, root_i.index
#     #     print root_i.index == '210333'
#     
#     
#         root_i = get_node_by_id(masterNode,['30011'])
#         root_i.printRect()
#         print root_i.enrichNodes[0].x, root_i.enrichNodes[0].y
#         print root_i.enrichNodes[1].x, root_i.enrichNodes[1].y
#         print root_i.ishomog, root_i.index

#        # multiple inclusions
#        cornerA = f_circle(x0,y0)
#        cornerB = f_circle(x1,y0)
#        cornerC = f_circle(x1,y1)
#        cornerD = f_circle(x0,y1)
#        R = 1.0/3.0
#
#        cornerA_s = f_circle_s(x0,y0)
#        cornerB_s = f_circle_s(x1,y0)
#        cornerC_s = f_circle_s(x1,y1)
#        cornerD_s = f_circle_s(x0,y1)
#        Rs = 1.0/3.0
#
#        cornerA_c1 = f_circle1(x0,y0)
#        cornerB_c1 = f_circle1(x1,y0)
#        cornerC_c1 = f_circle1(x1,y1)
#        cornerD_c1 = f_circle1(x0,y1)
#        R1 = 1.0/6.0
#
#        cornerA_c2 = f_circle2(x0,y0)
#        cornerB_c2 = f_circle2(x1,y0)
#        cornerC_c2 = f_circle2(x1,y1)
#        cornerD_c2 = f_circle2(x0,y1)
#        R2 = 1.0/6.0

        #cornerA = point_in_poly(x0,y0,rhombus)
        #cornerB = point_in_poly(x1,y0,rhombus)
        #cornerC = point_in_poly(x1,y1,rhombus) 
        #cornerD = point_in_poly(x0,y1,rhombus)

        # NE corner triangle
        #domainInclusion = [(0.9,1.0), (1.0,1.0),(1.0,0.9)]
        # SE corner triangle
        #domainInclusion = [(0.9,0), (1.0,0.0),(1.0,0.1)]
        # NW corner triangle
        #domainInclusion = [(0.0,0.9), (0.1,1.0),(0.0,1.0)]
        # SW corner triangle
        #domainInclusion = [(0.0,0.0),(0.1,0.0),(0.0,0.1)]
        # SE Triangle with SE/NW corners cut
        #domainInclusion = [(2.0/3.0,0), (1.0,0.0),(1.0,1.0/3.0)]

        #cornerA = point_in_poly(x0,y0,domainInclusion)
        #cornerB = point_in_poly(x1,y0,domainInclusion)
        #cornerC = point_in_poly(x1,y1,domainInclusion) 
        #cornerD = point_in_poly(x0,y1,domainInclusion)

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
            

            # slanted interface
            # set the coefficient of the conductivity
#            if coords[0,0] <= interface_fcn(coords[0,1]) and coords[1,0]<= interface_fcn(coords[0,1]) and k1!=k2:
        #    if point_in_on_poly(x0,y0,polygonDef) and point_in_on_poly(x0,y1,polygonDef) and point_in_on_poly(x1,y0,polygonDef) and point_in_on_poly(x1,y1,polygonDef):
        #        K_cst = k1 # y <= loc
        #    else:
        #        K_cst = k2 # y > loc
            
#            polygonDef = domainInclusion
#            # rhombus inside domain:
#            if point_in_on_poly(x0,y0,polygonDef) and point_in_on_poly(x0,y1,polygonDef) and point_in_on_poly(x1,y0,polygonDef) and point_in_on_poly(x1,y1,polygonDef):
#                K_cst = k2 
#            else:
#                K_cst = k1


#            if cornerA>R*R and cornerB>R*R and cornerC>R*R and cornerD>R*R:
#                K_cst = k1
#            else:
#                if cornerA < R*R and cornerB < R*R and cornerC < R*R and cornerD < R*R:
#                    K_cst = k2
#                else:
#                    print 'ERROR!'
#                    K_cst = k1
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
                

                

#            # multiple inclusions:
#            if ( (cornerA_s>Rs*Rs and cornerB_s>Rs*Rs and cornerC_s>Rs*Rs and cornerD_s>Rs*Rs) and 
#                (cornerA_c1>R1*R1 and cornerB_c1>R1*R1 and cornerC_c1>R1*R1 and cornerD_c1>R1*R1) and 
#                (cornerA_c2>R2*R2 and cornerB_c2>R2*R2 and cornerC_c2>R2*R2 and cornerD_c2>R2*R2) ):
#                K_cst = k1
#            else:
#                K_cst = k2

            #if cornerA == True and cornerB == True and cornerC == True and cornerD == True:
            #    K_cst = k2
            #else:    
            #    if cornerA == False and cornerB == False and cornerC == False and cornerD == False:
            #        K_cst = k1
            #    else:
            #            print 'ERROR! -  inside igfem2d.py - assembling the stiffness matrix'                

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
            
#            if nodes[0] == 41 and nodes[1] == 42 and nodes[2] == 51 and nodes[3] == 50:
#                [x_fct,y_fct] = xy_fct_HN(coords[0:4,0],coords[0:4,1],0,0,1,0)
#                Jac = jacobian_mat_HN( coords[0:4,0], coords[0:4,1],0,0,1,0)
#                det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
#                print 'east edge'
#                
#            if nodes[0] == 43 and nodes[1] == 44 and nodes[2] == 53 and nodes[3] == 52:
#                [x_fct,y_fct] = xy_fct_HN(coords[0:4,0],coords[0:4,1],0,0,0,1)
#                Jac = jacobian_mat_HN( coords[0:4,0], coords[0:4,1],0,0,0,1)
#                det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
#                print 'west edge'
#                
#            if nodes[0] == 51 and nodes[1] == 52 and nodes[2] == 61 and nodes[3] == 60:
#                [x_fct,y_fct] = xy_fct_HN(coords[0:4,0],coords[0:4,1],0,1,0,0)
#                Jac = jacobian_mat_HN( coords[0:4,0], coords[0:4,1],0,1,0,0)
#                det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
#                print 'south edge'
#                
#            if nodes[0] == 33 and nodes[1] == 34 and nodes[2] == 43 and nodes[3] == 42:
#                [x_fct,y_fct] = xy_fct_HN(coords[0:4,0],coords[0:4,1],1,0,0,0)
#                Jac = jacobian_mat_HN( coords[0:4,0], coords[0:4,1],1,0,0,0)
#                det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)
#                print 'north edge'
                
                    
# Method 2 for Hanging Nodes                
#===============================================================================
#             # construct the local matrix and local components of the load vector    
#             for i in range(0,4):
#                 for j in range(0,4):
#                     if nodes[i] >=  nodes[j]:
#                         Kefunc = lambda ee,nn: K_cst * ( Nx[i](x_fct(ee,nn),y_fct(ee,nn)) * Nx[j](x_fct(ee,nn),y_fct(ee,nn)) + Ny[i](x_fct(ee,nn),y_fct(ee,nn)) * Ny[j](x_fct(ee,nn),y_fct(ee,nn)) ) * det_Jac(ee,nn)
#                         Ke[i,j] = quad2d(Kefunc,-1,1,-1,1,ui,wi)
# 
#                 # construct the local load vector
#                 fv = lambda ee,nn: rhs(x_fct(ee,nn),y_fct(ee,nn)) * Nbasis[i](x_fct(ee,nn),y_fct(ee,nn)) * det_Jac(ee,nn)
#                 Fe[i] = quad2d(fv,-1,1,-1,1,ui,wi)
#===============================================================================
            
            # Method 1
            # construct the local matrix and local components of the load vector    
            for i in range(0,4):
                for j in range(0,4):
#                     if nodes[i] >=  nodes[j]:
                        Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                        Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)

                # construct the local load vector
                fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)


            # add the local stiffness matrix and local load vector to the global K and F
            for i in range(0,4):
                for j in range(0,4):
#                     if nodes[i] >= nodes[j]:
                        K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke[i,j]
#                         K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                F[nodes[i],0] = F[nodes[i],0] + Fe[i]

    
           
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
                    print 'False positive: Odd diagonal '
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
                                Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                                Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)

                        # construct the local load vector
                        fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                        Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)

            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,4):
                        for j in range(0,4):
                                K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke[i,j]
                        F[nodes[i],0] = F[nodes[i],0] + Fe[i]
            

             # FALSE POSITIVIES: only corner is in a different material
             # HOMOGENEOUS element
                if (corner0 == True or corner2 == True):
                    print ' False positive: Even diagonal '
                    elCorner = True
                    midInside = Point( (x0+x1)/2.0, (y0+y1)/2.0 )
                                        
                    pxValMid = image.GetPixel(int(midInside.x), int(midInside.y))
                    if pxValMid > binBnd[1]:
#                    if point_in_on_poly(midInside.x,midInside.y,polygonDef):
                        K_cst = k2
                    else:
                        K_cst = k1
                    
                    Ke = np.zeros((4,4))
                    Fe = np.zeros((4,1))
            
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,4):
                        for j in range(0,4):
                                Kefunc = lambda x,y: K_cst * ( Nx[i](x,y) * Nx[j](x,y) + Ny[i](x,y) * Ny[j](x,y) )
                                Ke[i,j] = quad2d(Kefunc,x0,x1,y0,y1,ui,wi)

                        # construct the local load vector
                        fv = lambda x,y: rhs(x,y) * Nbasis[i](x,y)
                        Fe[i] = quad2d(fv,x0,x1,y0,y1,ui,wi)

            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,4):
                        for j in range(0,4):
                                K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke[i,j]
                        F[nodes[i],0] = F[nodes[i],0] + Fe[i]


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
                    print 'Diagonal, element ', e
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

                    if len(root.enrichNodes) == 2:
                        print 'DIAGONAL WITH QUADRATIC!!!!!!!!!'
                    
                        [x_fct_1, y_fct_1] = tri_xy_fct( coords_trid1[:,0], coords_trid1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( coords_trid2[:,0], coords_trid2[:,1] )
        
                        J1 = tri_jacobian_mat( coords_trid1[:,0], coords_trid1[:,1] )
                        J2 = tri_jacobian_mat( coords_trid2[:,0], coords_trid2[:,1] )
                    
                    if len(root.enrichNodes) == 3:
                                
                        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
                
                        coord_enrich = coord_enrich_comp(root, midPoint)
                               
                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [1,2,0] # local order
                        else:
                            lOrd = [0,1,2] 
                          
                        vec1_x = [ coords_trid1[ lOrd[0],0], coords_trid1[lOrd[1],0], coords_trid1[lOrd[2],0], (coords_trid1[ lOrd[0],0] + coords_trid1[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid1[lOrd[0],0] + coords_trid1[lOrd[2],0])/2.0  ]
                        vec1_y = [ coords_trid1[ lOrd[0],1], coords_trid1[lOrd[1],1], coords_trid1[lOrd[2],1], (coords_trid1[ lOrd[0],1] + coords_trid1[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid1[lOrd[0],1] + coords_trid1[lOrd[2],1])/2.0  ]
                
                        [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
                        J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
                        
                        if(corner0 == True and corner2 == True): 
                        # even diagonal: SW - NE
                            lOrd = [2,0,1] # local order
                        else:
                            lOrd = [1,2,0] 
                            
                        vec2_x = [ coords_trid2[ lOrd[0],0], coords_trid2[lOrd[1],0], coords_trid2[lOrd[2],0], (coords_trid2[ lOrd[0],0] + coords_trid2[lOrd[1],0])/2.0, coord_enrich.x, (coords_trid2[lOrd[0],0] + coords_trid2[lOrd[2],0])/2.0  ]
                        vec2_y = [ coords_trid2[ lOrd[0],1], coords_trid2[lOrd[1],1], coords_trid2[lOrd[2],1], (coords_trid2[ lOrd[0],1] + coords_trid2[lOrd[1],1])/2.0, coord_enrich.y, (coords_trid2[lOrd[0],1] + coords_trid2[lOrd[2],1])/2.0  ]
                
                        [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
                        J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )
                        
                    det_J1 = lambda e,n: determinant(J1)(e,n)
                    det_J2 = lambda e,n: determinant(J2)(e,n)
                        
                    ## FIRST TRIANGLE
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,3):
                        for j in range(0,3):
#                             if nodes_trid1[i] >=  nodes_trid1[j]:
                                Kefunc_trid1 = lambda e,n: K_cst_trid1 * ( Nx_trid1[i](x_fct_1(e,n), y_fct_1(e,n)) * Nx_trid1[j](x_fct_1(e,n), y_fct_1(e,n)) + Ny_trid1[i](x_fct_1(e,n), y_fct_1(e,n)) * Ny_trid1[j](x_fct_1(e,n), y_fct_1(e,n)) ) * det_J1(e,n)
                                Ke_trid1[i,j] = 1.0/2.0 * quad2d(Kefunc_trid1,x0,x1,y0,y1,ui,wi)

                        # construct the local load vector
                        fv_trid1 = lambda x,y: rhs(x,y) * Nbasis_trid1[i](x,y)
                        Fe_trid1[i] = 1.0/2.0 * quad2d(fv_trid1,x0,x1,y0,y1,ui,wi)

            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,3):
                        for j in range(0,3):
#                             if nodes_trid1[i] >= nodes_trid1[j]:
                                K[nodes_trid1[i],nodes_trid1[j]] = K[nodes_trid1[i],nodes_trid1[j]] + Ke_trid1[i,j]
#                                 K[nodes_trid1[j],nodes_trid1[i]] = K[nodes_trid1[i],nodes_trid1[j]]
                        F[nodes_trid1[i],0] = F[nodes_trid1[i],0] + Fe_trid1[i]
                        
                    ## SECOND TRIANGLE
                    # construct the local matrix and local components of the load vector    
                    for i in range(0,3):
                        for j in range(0,3):
#                             if nodes_trid2[i] >=  nodes_trid2[j]:
                                Kefunc_trid2 = lambda e,n: K_cst_trid2 * ( Nx_trid2[i](x_fct_2(e,n), y_fct_2(e,n)) * Nx_trid2[j](x_fct_2(e,n), y_fct_2(e,n)) + Ny_trid2[i](x_fct_2(e,n), y_fct_2(e,n)) * Ny_trid2[j](x_fct_2(e,n), y_fct_2(e,n)) )* det_J2(e,n)
                                Ke_trid2[i,j] = 1.0/2.0 * quad2d(Kefunc_trid2,x0,x1,y0,y1,ui,wi)

                        # construct the local load vector
                        fv_trid2 = lambda x,y: rhs(x,y) * Nbasis_trid2[i](x,y)
                        Fe_trid2[i] = 1.0/2.0 * quad2d(fv_trid2,x0,x1,y0,y1,ui,wi)

            
                    # add the local stiffness matrix and local load vector to the global K and F
                    for i in range(0,3):
                        for j in range(0,3):
#                             if nodes_trid2[i] >= nodes_trid2[j]:
                                K[nodes_trid2[i],nodes_trid2[j]] = K[nodes_trid2[i],nodes_trid2[j]] + Ke_trid2[i,j]
#                                 K[nodes_trid2[j],nodes_trid2[i]] = K[nodes_trid2[i],nodes_trid2[j]]
                        F[nodes_trid2[i],0] = F[nodes_trid2[i],0] + Fe_trid2[i]
                        

                else:

                    # the North-West corner is cut, 0-4-3, 2-5-3

                    if (
                        ((enrich1[0] == x0 and enrich2[1] == y1) or
                         (enrich2[0] == x0 and enrich1[1] == y1)) and
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))                         
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        ):
                        print "NW corner"
                        [Ke_NW,Fe_NW] = NW_corner(p,ui,wi,k1,k2,nodes,root,image)

                        
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_NW[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_NW[i]


                    # the South-East corner is cut, 0-4-1, 1-5-2
                    if ( 
                        ((enrich1[1] == y0 and enrich2[0] == x1) or
                          (enrich2[1] == y0 and enrich1[0] == x1)) and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))                        
                        ):
                        print 'SE corner'
                        
                        [Ke_SE,Fe_SE] = SE_corner(p,ui,wi,k1,k2,nodes,root,image)
    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_SE[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_SE[i]
    
                    # the North East corner is cut, 1-4-2, 2-5-3
                    if ( 
                        ((enrich1[0] == x1 and enrich2[1] == y1) or
                          (enrich2[0] == x1 and enrich1[1] == y1)) and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        ):
                        print "NE corner"
                        [Ke_NE,Fe_NE] = NE_corner(p,ui,wi,k1,k2,nodes,root,image)
    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_NE[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_NE[i]

                    # the South-West corner is cut, 0-4-1, and 0-5-3
                    if ( 
                        ((enrich1[1] == y0 and enrich2[0] == x0) or
                         (enrich2[1] == y0 and enrich1[0] == x0)) and
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))                         
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        ):
                        print "SW corner"
                        [Ke_SW,Fe_SW] = SW_corner(p,ui,wi,k1,k2,nodes,root,image)
        
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_SW[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_SW[i]
                            
                    # the South edge
                    if (  ((enrich1[1] == y0 and enrich2[1] == y1) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x1) ) or
                        ( (enrich2[1] == y0 and enrich1[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1 ) ) ):

                        print 'South edge'
                        if not(on_corners(enrich2,coords)) :
#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :                            
                            south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_South,Fe_South] = South_edge(p,ui,wi,k1,k2,south_nodes,root,image)

                        print nodes
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
#                                 if south_nodes[i] >= south_nodes[j]:
                                    K[south_nodes[i],south_nodes[j]] = K[south_nodes[i],south_nodes[j]] + Ke_South[i,j]
#                                     K[south_nodes[j],south_nodes[i]] = K[south_nodes[i],south_nodes[j]]
                            F[south_nodes[i],0] = F[south_nodes[i],0] + Fe_South[i]

                    # the North edge
                    if ( ( (enrich1[1] == y0 and enrich2[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1)) or
                        ( (enrich1[1] == y1 and enrich2[1] == y0) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x0) )  ):
                        print 'North edge'

                        if not(on_corners(enrich2,coords)):
#                         if not(on_corners(enrich2,x0,y0,x1,y1)):
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5] ]
                        else:
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] ]
                        [Ke_North,Fe_North] = North_edge(p,ui,wi,k1,k2,north_nodes,root,image)

                        print north_nodes
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
#                                 if north_nodes[i] >= north_nodes[j]:
                                    K[north_nodes[i],north_nodes[j]] = K[north_nodes[i],north_nodes[j]] + Ke_North[i,j]
#                                     K[north_nodes[j],north_nodes[i]] = K[north_nodes[i],north_nodes[j]]
                            F[north_nodes[i],0] = F[north_nodes[i],0] + Fe_North[i]

                    
                    # the West edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1)) or
                        ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)  ) ):
                        print 'West edge'

                        enrich1 = np.array(p[nodes[4]])
                        enrich2 = np.array(p[nodes[5]])

#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :
                        if not(on_corners(enrich2,coords)) :
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_West,Fe_West] = West_edge(p,ui,wi,k1,k2,west_nodes,root,image)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
#                                 if west_nodes[i] >= west_nodes[j]:
                                    K[west_nodes[i],west_nodes[j]] = K[west_nodes[i],west_nodes[j]] + Ke_West[i,j]
#                                     K[west_nodes[j],west_nodes[i]] = K[west_nodes[i],west_nodes[j]]
                            F[west_nodes[i],0] = F[west_nodes[i],0] + Fe_West[i]



                    # the East edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)) or
                            ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1) )  ):
                        print 'East edge'

#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :
                        if not(on_corners(enrich2,coords)) :

                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]

                        [Ke_East,Fe_East] = East_edge(p,ui,wi,k1,k2,east_nodes,root,image)

                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,5):
                            for j in range(0,5):
#                                 if east_nodes[i] >= east_nodes[j]:
                                    K[east_nodes[i],east_nodes[j]] = K[east_nodes[i],east_nodes[j]] + Ke_East[i,j]
#                                     K[east_nodes[j],east_nodes[i]] = K[east_nodes[i],east_nodes[j]]
                            F[east_nodes[i],0] = F[east_nodes[i],0] + Fe_East[i]
    

#                     if e == 74:
#                         print p1.x, p2.x, p1.y, p4.y, enrN1.x, enrN2.x
#                         print 'element', e, x0, x1, enrich1[0], enrich2[0]
#                         print (enrich1[0] == x0  and enrich2[0] == x1)
#                         print  (enrich1[0] == x1 and enrich2[0] == x0)
#                         print not(on_corners(enrich1,x0,y0,x1,y1))
#                         print not(on_corners(enrich2,x0,y0,x1,y1))
#                         print nodes
#                         x0 = coords[0,0]
#                         x1 = coords[1,0]
#                         y0 = coords[0,1]
#                         y1 = coords[2,1] 
#                         print (p[nodes[5]])
##                    if root.index == '232110':   
##                        print ' ------=========----------blah'
##                        print not(on_corners(enrich1,coords)), not(on_corners(enrich2,coords))
##                        print 'coords', coords
##                        print 'nodes', nodes
##                        print 'enrich1', enrich1[0], enrich1[1]  
##                        print 'enrich2', enrich2[0], enrich2[1]
##                        print [x0,x1], [y0,y1]
                    # interface cuts the element horizontally into two quads, 0-4-3, 1-5-2 
                    if ( ((enrich1[0] == x0  and enrich2[0] == x1) or (enrich1[0] == x1 and enrich2[0] == x0)) and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1)) ):
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):


                        print "horizontal slide: quad-quad"
                        [Ke_Horiz,Fe_Horiz] = horizontal_cut(p,ui,wi,k1,k2,nodes,root,image)
                    
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_Horiz[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_Horiz[i]
            

                    # interface cuts the element vertically into two quads, 0-4-1, 3-5-2
                    if ( ((enrich1[1] == y0 and enrich2[1] == y1) or (enrich1[1] == y1 and enrich2[1] == y0 )) and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1)) ):
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):
                        print "vertical slide: quad-quad"
                        [Ke_Vertical,Fe_Vertical] = vertical_cut(p,ui,wi,k1,k2,nodes,root,image)
                        # add the local stiffness matrix and local load vector to the global K and F
                        for i in range(0,6):
                            for j in range(0,6):
#                                 if nodes[i] >= nodes[j]:
                                    K[nodes[i],nodes[j]] = K[nodes[i],nodes[j]] + Ke_Vertical[i,j]
#                                     K[nodes[j],nodes[i]] = K[nodes[i],nodes[j]]
                            F[nodes[i],0] = F[nodes[i],0] + Fe_Vertical[i]
            
    # end of loop
    # BCs: a * U + b * dU/dx + c * dU/dy + d = 0
    # Dirichlet: b,c = 0, homogeneous Dirichlet: b,c = 0, d = 0
    # Neumann: a = 0, and b or c may be 0, but not both

    print N
    U = sparse.lil_matrix((N,1))

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

    
    #F = F - np.dot(K,U)
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

    Kb = K[:,:]
    Fb = F[:]


#     scipy.io.savemat('Kb1.mat', mdict={'Kb': Kb})
#     scipy.io.savemat('Fb1.mat', mdict={'Fb': Fb})

    # Need to reduce the Kb matrix in order to be able to use it with SpSolve
    Kbreduced = sparse.lil_matrix((len(FreeNodes),len(FreeNodes)));
    for i in range(0,len(FreeNodes)):
        for j in range(0,len(FreeNodes)):
            Kbreduced[i,j] = Kb[FreeNodes[i],FreeNodes[j]]
    Kbreduced = Kbreduced.tocsr()

    # solve for the numerical solution
    numericalSoln = linsolve.spsolve(Kbreduced,Fb[FreeNodes,0])
    U[FreeNodes,0] = np.matrix(numericalSoln).T

#     print 'matrix, vector:',Kbreduced, Fb[FreeNodes,0]
    scipy.io.savemat('Uedge.mat', mdict={'U': U})
    scipy.io.savemat('Kedge.mat', mdict={'Kbreduced': Kbreduced})
    Fbreduced = Fb[FreeNodes,0]
    scipy.io.savemat('Fedge.mat', mdict={'Fbreduced': Fbreduced})
    
    
    # Getting the analytical solution and its vector form
    #uexact = lambda x,y: ufct(x,loc,k1,k2) 

#     HANGING NODES implementation
    for i in range(0,len(list_hanging_nodes)):
        listHN = list_hanging_nodes[i]
        U[listHN[0],0] = ( U[listHN[1],0] + U[listHN[2],0] ) / 2.0

    return  np.array(U) 

def computeNorm(p,t,pConf,tConf,ui,wi,k1,k2,U,UConf,masterNode,llist):
    print 'compute the L-2 norm ...'

    T = len(t)
    all_elems_sum = 0

    en_arr = [pp for pp in itertools.product(ui,repeat=2)]
    en_arr = np.array(en_arr)
    e_arr = en_arr[:,0] # contains the epsilon coordinates
    n_arr = en_arr[:,1] # contains the niu coordinates

    yloc = 0.1
    #xloc = yloc + 2.0/3.0

    Usolution = np.zeros((len(p),1))
    polygonList = []

    # COMPUTING THE L-2 NORM
    for e in range(0,T):
#     for e in range(800, T):
#     for e in range(833,T):
#     for e in range(1814,T):

#    for e in range(34,51):
#     for e in [35,37,38,40,43,44,47,48,49,50]:
#    for e in [38, 102, 122]:
#    for e in [96,97,100,102,103,105]: 
        
        nodes = t[e] # row of t =  node numbers of the 4 corners of element e
        nodes = np.array(nodes)
    
        root = get_node_by_id(masterNode,llist[e])

        
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

#         if len(nodes) == 4 and root.ishomog == 1:
        if (len(nodes) == 4 or  root.ishomog == 1) or thru_corner == True:
            
#             print 'element: e',e
            x_coords = coords[:,0]
            y_coords = coords[:,1]

#             [x_fct,y_fct] = xy_fct(coords[0:4,0],coords[0:4,1])
#             Jac = jacobian_mat( coords[0:4,0], coords[0:4,1])
#             det_Jac = lambda ee,nn: determinant(Jac)(ee,nn)

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

#             p1h,p2h,p3h,p4h = root.hn
#             
#             [N_hn,S_hn,E_hn,W_hn] = root.nsew
#             NbasisHN = Nbases(Nbasis,x0,x1,y0,y1,p,nodes,N_hn,S_hn,E_hn,W_hn)
#             
#             uiHN = ui
#             wiHN = wi
#             
#             p1,p2,p3,p4 = root.rect
#             
#             
#             if e!= 38:
#                 # locate the hanging nodes
#                 north_node = 0
#                 south_node = 0
#                 east_node = 0
#                 west_node = 0
#                 if S_hn == 1:
#                     south_hn = root.hn[1]
#                     shn = [south_hn.x, south_hn.y]
#                     s_ind = numpy.where(numpy.all(p==shn,axis=1))
#                     south_node = s_ind[0][0]
#                     Nbases_HN_South = lambda x,y: Nbasis57(x,y,p,nodes,x0,x1,y0,y1)[0]
#                     
#     #                [uiHN,wiHN] = lgwt(6,-1,1)
#                 else:
#                     Nbases_HN_South = lambda x,y: 0.0
#                     
#                 if N_hn == 1:
#                     north_hn = root.hn[0]
#                     nhn = [north_hn.x, north_hn.y]
#                     n_ind = numpy.where(numpy.all(p==nhn,axis=1))
#                     north_node = n_ind[0][0]           
#                     Nbases_HN_North = lambda x,y: Nbasis57(x,y,p,nodes,x0,x1,y0,y1)[1]
#                     
#     #                [uiHN,wiHN] = lgwt(6,-1,1)
#                 else:
#                     Nbases_HN_North = lambda x,y: 0.0   
#                                       
#                 if E_hn == 1:
#                     east_hn = root.hn[2]
#                     ehn = [east_hn.x, east_hn.y]
#                     e_ind = numpy.where(numpy.all(p==ehn,axis=1))
#                     east_node = e_ind[0][0]   
#                     Nbases_HN_East = lambda x,y: Nbasis68(x,y,p,nodes,x0,x1,y0,y1)[0]
#                     
#     #                [uiHN,wiHN] = lgwt(6,-1,1)
#                 else:
#                     Nbases_HN_East = lambda x,y: 0.0
#                                                  
#                 if W_hn == 1:
#                     west_hn = root.hn[3]
#                     whn = [west_hn.x, west_hn.y]
#                     w_ind = numpy.where(numpy.all(p==whn,axis=1))
#                     west_node = w_ind[0][0]            
#                     Nbases_HN_West = lambda x,y: Nbasis68(x,y,p,nodes,x0,x1,y0,y1)[1]
#                     
#     #                [uiHN,wiHN] = lgwt(6,-1,1)
#                 else:
#                     Nbases_HN_West = lambda x,y: 0.0                
#     
#                 if N_hn == 0 and S_hn == 0 and E_hn == 0 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], nodes[3] ]]
#               
#                 uh_elem = lambda x,y: ( U[t[e][0],0] * NbasisHN[0](x,y) +
#                                 U[t[e][1],0] * NbasisHN[1](x,y) +
#                                 U[t[e][2],0] * NbasisHN[2](x,y) +
#                                 U[t[e][3],0] * NbasisHN[3](x,y) +
#                                 
# #                uh_elem = lambda x,y: ( U[t[e][0],0] * Nbasis[0](x,y) +
# #                                U[t[e][1],0] * Nbasis[1](x,y) +
# #                                U[t[e][2],0] * Nbasis[2](x,y) +
# #                                U[t[e][3],0] * Nbasis[3](x,y) +
#                                 
#                                 U[north_node,0] * Nbases_HN_North(x,y) +  
#                                 U[south_node,0] * Nbases_HN_South(x,y) +
#                                 U[east_node,0] * Nbases_HN_East(x,y) +
#                                 U[west_node,0] * Nbases_HN_West(x,y)
#                         ) 
#                         
#                 Usolution[nodes[0],0] = uh_elem(p[nodes[0],0],p[nodes[0],1])
#                 Usolution[nodes[1],0] = uh_elem(p[nodes[1],0],p[nodes[1],1])
#                 Usolution[nodes[2],0] = uh_elem(p[nodes[2],0],p[nodes[2],1])
#                 Usolution[nodes[3],0] = uh_elem(p[nodes[3],0],p[nodes[3],1])
#                                 
#                 # 1 - Hanging Node:
#                 # South
#                 if N_hn == 0 and S_hn == 1 and E_hn == 0 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], nodes[2], nodes[3] ]]
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                 # East
#                 if N_hn == 0 and S_hn == 0 and E_hn == 1 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], nodes[1], east_node, nodes[2], nodes[3] ]]
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                 # North  
#                 if N_hn == 1 and S_hn == 0 and E_hn == 0 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], north_node, nodes[3] ]]
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                 # West 
#                 if N_hn == 0 and S_hn == 0 and E_hn == 0 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], nodes[3], west_node ]]
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                     
#                 # 2 Hanging Nodes:
#                 # South & East
#                 if N_hn == 0 and S_hn == 1 and E_hn == 1 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], east_node, nodes[2], nodes[3] ]]
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                 # South & North 
#                 if N_hn == 1 and S_hn == 1 and E_hn == 0 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], nodes[2], north_node, nodes[3] ]]
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                 # South & West
#                 if N_hn == 0 and S_hn == 1 and E_hn == 0 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], nodes[2], nodes[3], west_node ]]
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                 # North & East 
#                 if N_hn == 1 and S_hn == 0 and E_hn == 1 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], nodes[1], east_node, nodes[2], north_node, nodes[3] ]]
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                 # East & West
#                 if N_hn == 0 and S_hn == 0 and E_hn == 1 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], nodes[1], east_node, nodes[2], nodes[3], west_node ]]
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#     #            # North & West
#                 if N_hn == 1 and S_hn == 0 and E_hn == 0 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], nodes[1], nodes[2], north_node, nodes[3], west_node ]]
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                      
#                 # 3 Hanging Nodes:
#                 # South, East & North
#                 if N_hn == 1 and S_hn == 1 and E_hn == 1 and W_hn == 0:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], east_node, nodes[2], north_node, nodes[3] ]]
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                 # South, North & West 
#                 if N_hn == 1 and S_hn == 1 and E_hn == 0 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], nodes[2], north_node, nodes[3], west_node ]]
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#                 # East, West & North  
#                 if N_hn == 1 and S_hn == 0 and E_hn == 1 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], nodes[1], east_node, nodes[2], north_node, nodes[3], west_node ]]
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                     Usolution[north_node,0] = uh_elem(p[north_node,0],p[north_node,1])
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                 # East, West & South
#                 if N_hn == 0 and S_hn == 1 and E_hn == 1 and W_hn == 1:
#                     polygonList = polygonList + [[nodes[0], south_node, nodes[1], east_node, nodes[2], nodes[3], west_node ]]
#                     Usolution[west_node,0] = uh_elem(p[west_node,0],p[west_node,1])
#                     Usolution[east_node,0] = uh_elem(p[east_node,0],p[east_node,1])
#                     Usolution[south_node,0] = uh_elem(p[south_node,0],p[south_node,1])
#      
#     
                  
# canceling the norm computation
#             # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
#             # for transformation from the parametric element to phisycal element
#             # of the Gauss nodes ui
#             [x_transform_fct,y_transform_fct] = xy_fct(x_coords,y_coords)            
# 
#             Jac = jacobian_mat( coords[:,0], coords[:,1] )
#             detJ = lambda eps,niu: determinant(Jac)(eps,niu)
# 
#     #            el_sum =  gauss_integration_HN(ui,wi,UConf,pConf,tConf,x_fct_HN,y_fct_HN,uh_elem,det_Jac)
#                 
#                 #el_sum =  gauss_integration(uiHN,wiHN,UConf,pConf,tConf,x_transform_fct,y_transform_fct,uh_elem,detJ)
#             el_sum =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct,y_transform_fct,uh_elem,detJ)
#             all_elems_sum = all_elems_sum + el_sum;
    
    
    #             if y0 <= yloc and yloc <= y1:
    #                 tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
    #                 sfct = uh_elem(tx,yloc)
    #                 pylab.plot(tx,sfct)
            
                    
#             else:
#                 
#                     # hard coded element 38
#                     nodes = [169, 171, 187, 186, 170, 178]
# 
#                  
#                   # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
#                     tri_nodes1 = [nodes[0],nodes[4],nodes[5]]
#                     tri_nodes2 = [nodes[5],nodes[2],nodes[3]]
#                     tri_nodes3 = [nodes[4],nodes[2],nodes[5]]
#                     tri_nodes4 = [nodes[4],nodes[1],nodes[2]] # the one triangle in a diff material
#     
#                     polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
#                     polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
#                     polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
#                     polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
# 
#                     Pe1 = np.zeros((3,3))
#                     Pe1[:,0] = np.ones((3,1)).transpose()
#                     Pe1[:,1] = p[tri_nodes1[0:3],0]
#                     Pe1[:,2] = p[tri_nodes1[0:3],1]
#                     C1 = np.linalg.inv(Pe1)
#                     Nbasis_tri1 = tribasisFct(C1)
#                     Nx_tri1 = triderivX(C1)
#                     Ny_tri1 = triderivY(C1)
# 
#                     Pe2 = np.zeros((3,3))
#                     Pe2[:,0] = np.ones((3,1)).transpose()
#                     Pe2[:,1] = p[tri_nodes2[0:3],0]
#                     Pe2[:,2] = p[tri_nodes2[0:3],1]
#                     C2 = np.linalg.inv(Pe2)
#                     Nbasis_tri2 = tribasisFct(C2)
#                     Nx_tri2 = triderivX(C2)
#                     Ny_tri2 = triderivY(C2)
#     
#                     Pe3 = np.zeros((3,3))
#                     Pe3[:,0] = np.ones((3,1)).transpose()
#                     Pe3[:,1] = p[tri_nodes3[0:3],0]
#                     Pe3[:,2] = p[tri_nodes3[0:3],1]
#                     C3 = np.linalg.inv(Pe3)
#                     Nbasis_tri3 = tribasisFct(C3)
#                     Nx_tri3 = triderivX(C3)
#                     Ny_tri3 = triderivY(C3)
#     
#                     Pe4 = np.zeros((3,3))
#                     Pe4[:,0] = np.ones((3,1)).transpose()
#                     Pe4[:,1] = p[tri_nodes4[0:3],0]
#                     Pe4[:,2] = p[tri_nodes4[0:3],1]
#                     C4 = np.linalg.inv(Pe4)
#                     Nbasis_tri4 = tribasisFct(C4)
#                     Nx_tri4 = triderivX(C4)
#                     Ny_tri4 = triderivY(C4)
#     
#                     tri_coords1 = p[tri_nodes1]
#                     tri_coords2 = p[tri_nodes2]
#                     tri_coords3 = p[tri_nodes3]
#                     tri_coords4 = p[tri_nodes4]
#     
#                     [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
#                     [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
#                     [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
#                     [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
#     
#                     Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri1[1],Nbasis_tri1[2]]
#                     Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri2[0]]
#                     Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri3[0],Nbasis_tri3[2]]
#                     Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],lambda x,y: 0]
#     
#                     J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
#                     J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
#                     J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
#                     J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
#                     detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
#                     detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
#                     detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
#                     detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
#         
#     
#                     node_enr_4 = [nodes[4]]
#                     node_enr_5 = [nodes[5]]
#                     coords_enr_4 = p[node_enr_4]
#                     coords_enr_5 = p[node_enr_5]
#     
#                     x4 = coords_enr_4[0,0]
#                     y4 = coords_enr_4[0,1]
#     
#                     x5 = coords_enr_5[0,0]
#                     y5 = coords_enr_5[0,1]
#     
#                     w1 = y1 - y5
#                     w2 = y5 - y0
#                     factor_W = ( 2 * min(w1,w2) / (w1 + w2)) ** 2 
#                     if factor_W > EPS_FACTOR:
#                         factor_W = 1
#     
#                     s1 = x4 - x0
#                     s2 = x1 - x4
#                     factor_S = ( 2 * min(s1,s2) / (s1 + s2)) ** 2 
#                     if factor_S > EPS_FACTOR:
#                         factor_S = 1
#     
#                     uh_elem_1 = lambda x,y: (
#                                             U[nodes[0],0] * Nbasis[0](x,y) +
#                                             U[nodes[1],0] * Nbasis[1](x,y) + 
#                                             U[nodes[2],0] * Nbasis[2](x,y) +
#                                             U[nodes[3],0] * Nbasis[3](x,y) +
#                                             U[tri_nodes1[1],0] * Nbasis_tri1[1](x,y) * factor_S+
#                                             U[tri_nodes1[2],0] * Nbasis_tri1[2](x,y) * factor_W )
#     
#                     el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
#     
#                     uh_elem_2 = lambda x,y: (
#                                             U[nodes[0],0] * Nbasis[0](x,y) +
#                                             U[nodes[1],0] * Nbasis[1](x,y) + 
#                                             U[nodes[2],0] * Nbasis[2](x,y) +
#                                             U[nodes[3],0] * Nbasis[3](x,y) +
#                                             U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_W )
#     
#                     el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
#     
#                     uh_elem_3 = lambda x,y: (
#                                             U[nodes[0],0] * Nbasis[0](x,y) +
#                                             U[nodes[1],0] * Nbasis[1](x,y) + 
#                                             U[nodes[2],0] * Nbasis[2](x,y) +
#                                             U[nodes[3],0] * Nbasis[3](x,y) +
#                                             U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_S +
#                                             U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_W )
#     
#                     el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
#     
#                     uh_elem_4 = lambda x,y: (
#                                             U[nodes[0],0] * Nbasis[0](x,y) +
#                                             U[nodes[1],0] * Nbasis[1](x,y) + 
#                                             U[nodes[2],0] * Nbasis[2](x,y) +
#                                             U[nodes[3],0] * Nbasis[3](x,y) +
#                                             U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_S )
#     
#                     el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
#     
#                     all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
#                 
#                     Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
#                     Usolution[nodes[1],0] = uh_elem_4( p[nodes[1],0], p[nodes[1],1]  )
#                     Usolution[nodes[2],0] = uh_elem_4( p[nodes[2],0], p[nodes[2],1]  )
#                     Usolution[nodes[3],0] = uh_elem_2( p[nodes[3],0], p[nodes[3],1]  )
#                     Usolution[nodes[4],0] = uh_elem_1( p[nodes[4],0], p[nodes[4],1]  )
#                     Usolution[nodes[5],0] = uh_elem_1( p[nodes[5],0], p[nodes[5],1]  )
#    
#                 
#            if nodes[0] == 102 and nodes[1] == 104:
#                print south_node, east_node, north_node
#                print uh_elem(p[south_node,0],p[south_node,1])
#                print uh_elem(p[east_node,0],p[east_node,1])
#                print uh_elem(p[north_node,0],p[north_node,1])
#                print uh_elem(p[102,0],p[102,1]), uh_elem(p[104,0],p[104,1]), uh_elem(p[121,0], p[121,1]), uh_elem(p[120,0],p[120,1])
            
 
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
#                     # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
#                     # for transformation from the parametric element to phisycal element
#                     # of the Gauss nodes ui
#                     [x_transform_fct,y_transform_fct] = xy_fct(x_coords,y_coords)            
# 
#                     Jac = jacobian_mat( coords[:,0], coords[:,1] )
#                     detJ = lambda eps,niu: determinant(Jac)(eps,niu)
# 
#                     el_sum =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct,y_transform_fct,uh_elem,detJ)
# 
#                     all_elems_sum = all_elems_sum + el_sum;


#                     if y0 <= yloc and yloc <= y1:
#                         tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
#                         sfct = uh_elem(tx,yloc)
#                         pylab.plot(tx,sfct)

                # FALSE POSITIVIES: only corner is in a different material
                # HOMOGENEOUS element
#                if (corner0 == True or corner2 == True):
#                    print ' False positive: Even diagonal '
#                    elCorner = True


#                if enrich1[0] == x1 and elCorner == False:
                

#                if enrich1[1] == y0 and elCorner == False:

#                if enrich1[0] == x0 and elCorner == False:
#                # the North edge has the enrichment: 2-4-3
#                if enrich1[0] == x0  and elCorner == False:
            else:
                enrich1 = np.array(p[nodes[4]])
                enrich2 = np.array(p[nodes[5]])

                corner0 = ( min(abs(enrich1[0] - [x0]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner1 = ( min(abs(enrich1[0] - [x1]))<=1e-12) and (min(abs(enrich1[1] - [y0])) <= 1e-12 )
                corner2 = ( min(abs(enrich2[0] - [x1]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )
                corner3 = ( min(abs(enrich2[0] - [x0]))<=1e-12) and (min(abs(enrich2[1] - [y1])) <= 1e-12 )

#                 odd_even_diag = False
#                 # testing if interface is along a diagonal
#                 if (corner0 == True and corner2 == True) or (corner1 == True and corner3 == True):
#                     print 'Diagonal'
#                     odd_even_diag = True
#                 else:
#                     odd_even_diag = False
# 
#                 print 'Diagonal element: ', e, odd_even_diag
#                 print corner0, corner1, corner2, corner3
#                 if odd_even_diag == True:
                if (corner0 == True and corner2 == True) or (corner1 == True and corner3 == True):
                    if(corner0 == True and corner2 == True):
                        nodes_trid1 = [nodes[0], nodes[1], nodes[2]]
                        nodes_trid2 = [nodes[0], nodes[2], nodes[3]]
                        tc1 = [0,1,2]
                        tc2 = [0,2,3]
                    else:
                        nodes_trid1 = [nodes[0], nodes[1], nodes[3]]
                        nodes_trid2 = [nodes[1], nodes[2], nodes[3]]
                        tc1 = [0,1,3]
                        tc2 = [1,2,3]

                    print 'Diagonal', e
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
                    polygonList = polygonList + [[nodes_trid1[0], nodes_trid1[1], nodes_trid1[2] ]]
  
# canceling the norm computation    
#                     [x_transform_fct_trid1,y_transform_fct_trid1] = tri_xy_fct(x_coords_trid1,y_coords_trid1)            
#                     Jac_trid1 = tri_jacobian_mat( coords_trid1[:,0], coords_trid1[:,1] )
#                     detJ_trid1 = lambda eps,niu: determinant(Jac_trid1)(eps,niu)
#             
#                     el_sum_trid1 =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_trid1,y_transform_fct_trid1,uh_elem_trid1,detJ_trid1)
#                     all_elems_sum = all_elems_sum + el_sum_trid1;
    
    #                     if y0 <= yloc and yloc <= y1:
    #                         tx_trid1 = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
    #                         sfct_trid1 = uh_elem_trid1(tx_trid1,yloc)
    #                         pylab.plot(tx_trid1,sfct_trid1)
    
                    Usolution[nodes_trid2[0],0] = uh_elem_trid2(p[nodes_trid2[0],0],p[nodes_trid2[0],1])
                    Usolution[nodes_trid2[1],0] = uh_elem_trid2(p[nodes_trid2[1],0],p[nodes_trid2[1],1])
                    Usolution[nodes_trid2[2],0] = uh_elem_trid2(p[nodes_trid2[2],0],p[nodes_trid2[2],1])
                    polygonList = polygonList + [[nodes_trid2[0], nodes_trid2[1], nodes_trid2[2] ]]
    
# canceling the norm computation    
#                     [x_transform_fct_trid2,y_transform_fct_trid2] = tri_xy_fct(x_coords_trid2,y_coords_trid2)            
#                     Jac_trid2 = tri_jacobian_mat( coords_trid2[:,0], coords_trid2[:,1] )
#                     detJ_trid2 = lambda eps,niu: determinant(Jac_trid2)(eps,niu)
#         
#                     el_sum_trid2 =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_trid2,y_transform_fct_trid2,uh_elem_trid2,detJ_trid2)
#                     all_elems_sum = all_elems_sum + el_sum_trid2
    
    #                     if y0 <= yloc and yloc <= y1:
    #                         tx_trid2 = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
    #                         sfct_trid2 = uh_elem_trid2(tx_trid2,yloc)
    #                         pylab.plot(tx_trid2,sfct_trid2)

                else:
                    if ( ( (enrich1[1] == y0 and enrich2[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1)) or
                          ( (enrich1[1] == y1 and enrich2[1] == y0) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x0) )  ):
                        print 'norm computation: North edge'
    
#                         if not(on_corners(enrich2,x0,y0,x1,y1)):
                        if not(on_corners(enrich2,coords)):
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5] ]
                        else:
                            north_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4] ]
    
                        tri_nodes1 = [north_nodes[0],north_nodes[4],north_nodes[3]]
                        tri_nodes2 = [north_nodes[0],north_nodes[1],north_nodes[4]]
                        tri_nodes3 = [north_nodes[1],north_nodes[2],north_nodes[4]]
                        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
        
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
    
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
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
                                                
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_N)
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
#        
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_2( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
                    
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
        
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    
        
                    # the South edge
                    if (  ((enrich1[1] == y0 and enrich2[1] == y1) and (enrich2[0] == x0 or enrich2[0] == x1) and (enrich1[0] != x0 and enrich1[0] != x1) ) or
                          ( (enrich2[1] == y0 and enrich1[1] == y1) and (enrich1[0] == x0 or enrich1[0] == x1) and (enrich2[0] != x0 and enrich2[0] != x1 ) ) ):
                        print 'norm computation: South edge'
                        
#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :
                        if not(on_corners(enrich2,coords)) :
                          south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                          south_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
                        tri_nodes1 = [south_nodes[0],south_nodes[4],south_nodes[3]]
                        tri_nodes2 = [south_nodes[4],south_nodes[2],south_nodes[3]]
                        tri_nodes3 = [south_nodes[4],south_nodes[1],south_nodes[2]] # the one triangle in a diff material
                        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
        
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
        
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
            
        
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_S )
        
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
        
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_S )
# canceling the norm computation    
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
    
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_3( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
        
                        txP = np.arange(x0,x1+0.00001,0.001)
            
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            for idx in range(0,len(tx)):
            
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                            
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    
                    # the West edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1)) or
                          ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)  ) ):
                        print 'norm computation: West edge'
#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :
                        if not(on_corners(enrich2,coords)) :
                          west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            west_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
#                         print 'element ',e, west_nodes
                        tri_nodes1 = [west_nodes[0],west_nodes[1],west_nodes[4]]
                        tri_nodes2 = [west_nodes[1],west_nodes[2],west_nodes[4]]
                        tri_nodes3 = [west_nodes[4],west_nodes[2],west_nodes[3]] 
        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
        
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
    
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
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_W )
                                                
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_W)
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
#                         print U[tri_nodes1[2],0], U[tri_nodes2[2],0], U[tri_nodes3[0],0]    
#                         print tri_nodes2
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_3( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
                    
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
        
    
    
                    # the East edge
                    if ( ( (enrich1[0] == x0 and enrich2[0] == x1) and (enrich1[1] == y0 or enrich1[1] == y1) and (enrich2[1] != y0 and enrich2[1] != y1)) or
                          ( (enrich2[0] == x0 and enrich1[0] == x1) and (enrich2[1] == y0 or enrich2[1] == y1) and (enrich1[1] != y0 and enrich1[1] != y1) )  ):
                        print 'norm computation: East edge'
    
#                         if not(on_corners(enrich2,x0,y0,x1,y1)) :
                        if not(on_corners(enrich2,coords)) :
                          east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[5]]
                        else:
                            east_nodes = [nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]]
    
    
                        tri_nodes1 = [east_nodes[0],east_nodes[1],east_nodes[4]]
                        tri_nodes2 = [east_nodes[0],east_nodes[4],east_nodes[3]]
                        tri_nodes3 = [east_nodes[4],east_nodes[2],east_nodes[3]] # the one triangle in a diff material
        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
    
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
        
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
            
        
# canceling the norm computation
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E )
        
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
        
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_E )
# canceling the norm computation    
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_3( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_2( p[nodes[4],0], p[nodes[4],1]  )
                        txP = np.arange(x0,x1+0.00001,0.001)
            
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            for idx in range(0,len(tx)):
            
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
    
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                            
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    
    
    
                    # the North-West corner is cut, 0-4-3, 2-5-3
                    if ( 
                        ((enrich1[0] == x0 and enrich2[1] == y1) or
                        (enrich2[0] == x0 and enrich1[1] == y1)) and 
#                         odd_even_diag == False and 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        ):
                        
                        print "norm computation: NW corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (lower right corner cut)
                        tri_nodes1 = [nodes[0],nodes[1],nodes[4]]
                        tri_nodes2 = [nodes[4],nodes[1],nodes[5]]
                        tri_nodes3 = [nodes[1],nodes[2],nodes[5]]
                        tri_nodes4 = [nodes[4],nodes[5],nodes[3]] 
    
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                        polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                        [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
        
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
    
                        # scaling factor
                        #sx1 = abs(p[nodes[5]][0]- x0)
                        #sx2 = abs(x1 - p[nodes[5]][0])
                        #s_factor_x = 1#( 2 * min(sx1,sx2)/(sx1+sx2)) ** 2
        
                        #sy1 = abs(p[nodes[4]][1]- y0)
                        #sy2 = abs(y1 - p[nodes[4]][1])
                        #s_factor_y = 1#( 2 * min(sy1,sy2)/(sy1+sy2)) ** 2
    
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
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_W +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
                                                
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_N)
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_W +
                                                U[tri_nodes4[1],0] * Nbasis_tri4[1](x,y) * factor_N )
# canceling the norm computation        
#                         el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_1( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_4( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
                
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            tx4 = []
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords4):
                                    tx4.append(tx[idx])
        
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
                            tx4 = np.array(tx4)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    #                         pylab.plot( tx4, uh_elem_4(tx4, yloc))
                            
        
                    # the South-East corner is cut, 0-4-1, 1-5-2
                    if (
                        ((enrich1[1] == y0 and enrich2[0] == x1) or 
                         (enrich2[1] == y0 and enrich1[0] == x1)) and 
#                          odd_even_diag == False and 
#                          not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                          not(on_corners(enrich2,x0,y0,x1,y1))
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
        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                        polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                        [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
        
                        #Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri1[1],lambda x,y: 0]
                        #Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri2[0],Nbasis_tri2[1]]
                        #Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri3[0]]
                        #Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],Nbasis_tri4[2]]
        
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
        
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
            
        
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_S +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E ) 
        
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_E )
        
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_S +
                                                U[tri_nodes4[2],0] * Nbasis_tri4[2](x,y) * factor_E )
# canceling the norm computation
#                         el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
        
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_4( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_3( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
    
    #                    print 'SE - 1 ', uh_elem_1(p[nodes[4],0], p[nodes[4],1])
    #                    print 'SE - 2', uh_elem_2(p[nodes[4],0], p[nodes[4],1])
    #                    print 'SE - 3', uh_elem_3(p[nodes[4],0], p[nodes[4],1])
    #                    print 'SE - 4', uh_elem_4(p[nodes[4],0], p[nodes[4],1])
    #                    print 'SE - 1 - natural', uh_elem_1(p[nodes[0],0], p[nodes[0],1])
    #                    print 'SE - 2 - natural', uh_elem_2(p[nodes[0],0], p[nodes[0],1])
    #                    print 'SE - I - natural', uh_elem_1(p[nodes[3],0], p[nodes[3],1])
    #                    print 'SE - II - natural', uh_elem_2(p[nodes[3],0], p[nodes[3],1])
    #                    print 'SE - III - natural', uh_elem_3(p[nodes[3],0], p[nodes[3],1])
                        
    
                        txP = np.arange(x0,x1+0.00001,0.001)
            
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            tx4 = []
                            for idx in range(0,len(tx)):
            
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords4):
                                    tx4.append(tx[idx])
                            
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
                            tx4 = np.array(tx4)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    #                         pylab.plot( tx4, uh_elem_4(tx4, yloc))
    #     
                    # the North East corner is cut, 1-4-2, 2-5-3
                    if (
                        ((enrich1[0] == x1 and enrich2[1] == y1) or
                         (enrich2[0] == x1 and enrich1[1] == y1)) and
#                         odd_even_diag == False and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1)) 
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) 
                        ):
                        
                        print "norm computation: NE corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
                        tri_nodes1 = [nodes[0],nodes[5],nodes[3]]
                        tri_nodes2 = [nodes[0],nodes[4],nodes[5]]
                        tri_nodes3 = [nodes[0],nodes[1],nodes[4]]
                        tri_nodes4 = [nodes[4],nodes[2],nodes[5]] # the one triangle in a diff material
        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                        polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                        [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
        
        
                        Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri1[1]]
                        Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri2[1],Nbasis_tri2[2]]
                        Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri3[2],lambda x,y: 0]
                        Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],Nbasis_tri4[2]]
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
        
                
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
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
                    
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[1],0] * Nbasis_tri2[1](x,y) * factor_E +
                                                U[tri_nodes2[2],0] * Nbasis_tri2[2](x,y) * factor_N )
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
                    
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_E)
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_E +
                                                    U[tri_nodes4[2],0] * Nbasis_tri4[2](x,y) * factor_N )
# canceling the norm computation        
#                         el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
        
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_3( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_4( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_1( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_4( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_4( p[nodes[5],0], p[nodes[5],1]  )
    
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            tx4 = []
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                    
                                if point_in_on_poly( tx[idx], yloc, tri_coords4):
                                    tx4.append(tx[idx])
                            
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
                            tx4 = np.array(tx4)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    #                         pylab.plot( tx4, uh_elem_4(tx4, yloc))
        
                    # the South-West corner is cut, 0-4-1, and 0-5-3
                    if (
                        ((enrich1[1] == y0 and enrich2[0] == x0) or
                         (enrich2[1] == y0 and enrich1[0] == x0)) and 
#                         odd_even_diag == False and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1))
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords))
                        ):
                        
                        print "norm computation: SW corner"
                        # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper right corner cut)
                        tri_nodes1 = [nodes[0],nodes[4],nodes[5]]
                        tri_nodes2 = [nodes[5],nodes[2],nodes[3]]
                        tri_nodes3 = [nodes[4],nodes[2],nodes[5]]
                        tri_nodes4 = [nodes[4],nodes[1],nodes[2]] # the one triangle in a diff material
        
                        polygonList = polygonList + [[tri_nodes1[0], tri_nodes1[1], tri_nodes1[2]]]
                        polygonList = polygonList + [[tri_nodes2[0], tri_nodes2[1], tri_nodes2[2]]]
                        polygonList = polygonList + [[tri_nodes3[0], tri_nodes3[1], tri_nodes3[2]]]
                        polygonList = polygonList + [[tri_nodes4[0], tri_nodes4[1], tri_nodes4[2]]]
    
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
        
                        [x_fct_1, y_fct_1] = tri_xy_fct( tri_coords1[:,0], tri_coords1[:,1] )
                        [x_fct_2, y_fct_2] = tri_xy_fct( tri_coords2[:,0], tri_coords2[:,1] )
                        [x_fct_3, y_fct_3] = tri_xy_fct( tri_coords3[:,0], tri_coords3[:,1] )
                        [x_fct_4, y_fct_4] = tri_xy_fct( tri_coords4[:,0], tri_coords4[:,1] )
        
                        Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri1[1],Nbasis_tri1[2]]
                        Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],lambda x,y: 0,Nbasis_tri2[0]]
                        Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri3[0],Nbasis_tri3[2]]
                        Nbasis_4 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3],Nbasis_tri4[0],lambda x,y: 0]
        
                        J_tri1 = tri_jacobian_mat( tri_coords1[:,0], tri_coords1[:,1] )
                        J_tri2 = tri_jacobian_mat( tri_coords2[:,0], tri_coords2[:,1] )
                        J_tri3 = tri_jacobian_mat( tri_coords3[:,0], tri_coords3[:,1] )
                        J_tri4 = tri_jacobian_mat( tri_coords4[:,0], tri_coords4[:,1] )
                        detJ_tri1 = lambda e,n: determinant(J_tri1)(e,n)
                        detJ_tri2 = lambda e,n: determinant(J_tri2)(e,n)
                        detJ_tri3 = lambda e,n: determinant(J_tri3)(e,n)
                        detJ_tri4 = lambda e,n: determinant(J_tri4)(e,n)
            
        
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
# canceling the norm computation        
#                         el_sum_1 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_1, y_fct_1, uh_elem_1, detJ_tri1)
        
                        uh_elem_2 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes2[0],0] * Nbasis_tri2[0](x,y) * factor_W )
# canceling the norm computation        
#                         el_sum_2 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_2, y_fct_2, uh_elem_2, detJ_tri2)
        
                        uh_elem_3 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes3[0],0] * Nbasis_tri3[0](x,y) * factor_S +
                                                U[tri_nodes3[2],0] * Nbasis_tri3[2](x,y) * factor_W )
# canceling the norm computation        
#                         el_sum_3 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_3, y_fct_3, uh_elem_3, detJ_tri3)
        
                        uh_elem_4 = lambda x,y: (
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) +
                                                U[tri_nodes4[0],0] * Nbasis_tri4[0](x,y) * factor_S )
# canceling the norm computation        
#                         el_sum_4 =  gauss_integration(ui,wi,UConf,pConf,tConf, x_fct_4, y_fct_4, uh_elem_4, detJ_tri4)
#         
#                         all_elems_sum = all_elems_sum + el_sum_1 + el_sum_2 + el_sum_3 + el_sum_4
                    
                        Usolution[nodes[0],0] = uh_elem_1( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_4( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_4( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_2( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_1( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_1( p[nodes[5],0], p[nodes[5],1]  )
    
                        if y0 <= yloc and yloc <= y1:
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            tx1 = []
                            tx2 = []
                            tx3 = []
                            tx4 = []
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, tri_coords1):
                                    tx1.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords2):
                                    tx2.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords3):
                                    tx3.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, tri_coords4):
                                    tx4.append(tx[idx])
        
                            tx1 = np.array(tx1)
                            tx2 = np.array(tx2)
                            tx3 = np.array(tx3)
                            tx4 = np.array(tx4)
    #                         pylab.plot( tx1, uh_elem_1(tx1, yloc))
    #                         pylab.plot( tx2, uh_elem_2(tx2, yloc))
    #                         pylab.plot( tx3, uh_elem_3(tx3, yloc))
    #                         pylab.plot( tx4, uh_elem_4(tx4, yloc))
        
        
                    # interface cuts the element horizontally into two quads, 0-4-3, 1-5-2 
                    if ( ((enrich1[0] == x0  and enrich2[0] == x1) or 
                        (enrich1[0] == x1  and enrich2[0] == x0)) and 
#                         odd_even_diag == False and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1)) ):
                        not(on_corners(enrich1,coords)) and 
                        not(on_corners(enrich2,coords)) ):
                        print "norm computation: horizontal slide: quad-quad"
        
                        # nodes on the top and bottom side of the interface
                        top_nodes = [nodes[4], nodes[5], nodes[2],nodes[3]]
                        bottom_nodes = [nodes[0],nodes[1],nodes[5],nodes[4]]
    
                        if (enrich1[0] == x1  and enrich2[0] == x0):
                            top_nodes = [nodes[5], nodes[4], nodes[2],nodes[3]]
                            bottom_nodes = [nodes[0],nodes[1],nodes[4],nodes[5]]
    
        
                        top_coords = p[top_nodes,:]
                        bottom_coords = p[bottom_nodes,:]
    
        
                        polygonList = polygonList + [[top_nodes[0],top_nodes[1],top_nodes[2],top_nodes[3]]]
                        polygonList = polygonList + [[bottom_nodes[0],bottom_nodes[1],bottom_nodes[2],bottom_nodes[3]]]
    
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
        
                        # computing the Jacobian and the determinant of the left and right children of the parent element
                        J_top = jacobian_mat( top_coords[:,0], top_coords[:,1] )
                        J_bottom = jacobian_mat( bottom_coords[:,0], bottom_coords[:,1] )
                        detJ_top = lambda e,n: determinant(J_top)(e,n)
                        detJ_bottom = lambda e,n: determinant(J_bottom)(e,n)
        
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
# canceling the norm computation        
#                         # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
#                         # for transformation from the parametric element to phisycal element
#                         # of the Gauss nodes ui
#                         [x_transform_fct_T,y_transform_fct_T] = xy_fct(x_coords_T,y_coords_T)            
#                         el_sum_T =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_T,y_transform_fct_T,uh_elem_T,detJ_top)
#         
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
#                         [x_transform_fct_B,y_transform_fct_B] = xy_fct(x_coords_B,y_coords_B)            
#                         el_sum_B =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_B,y_transform_fct_B,uh_elem_B,detJ_bottom)
#         
#                         all_elems_sum = all_elems_sum + el_sum_T + el_sum_B;
                    
                        Usolution[nodes[0],0] = uh_elem_B( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_B( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_T( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_T( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_B( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_B( p[nodes[5],0], p[nodes[5],1]  )
    
                        if y0 <= yloc and yloc <= y1:
        
                            tx = np.arange(coords[0,0],coords[1,0]+0.01,0.001)
                            txT = []
                            txB = []
        
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, top_coords):
                                    txT.append(tx[idx])
                                if point_in_on_poly( tx[idx]-0.1, yloc, top_coords):
                                    txT.append(tx[idx]-0.1)
                                
                                if point_in_on_poly( tx[idx], yloc, bottom_coords):
                                    txB.append(tx[idx])
                                
                            txT = np.array(txT)
                            txB = np.array(txB)
    #                         pylab.plot( txT, uh_elem_T(txT, yloc))
    #                         pylab.plot( txB, uh_elem_B(txB, yloc))
        
                    # interface cuts the element vertically into two quads, 0-4-1, 3-5-2
                    if ( ((enrich1[1] == y0 and enrich2[1] == y1) or 
                        (enrich1[1] == y1 and enrich2[1] == y0)) and 
#                         odd_even_diag == False and 
#                         not(on_corners(enrich1,x0,y0,x1,y1)) and 
#                         not(on_corners(enrich2,x0,y0,x1,y1)) ):
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
        
                        polygonList = polygonList + [[left_nodes[0],left_nodes[1],left_nodes[2],left_nodes[3]]]
                        polygonList = polygonList + [[right_nodes[0],right_nodes[1],right_nodes[2],right_nodes[3]]]
    
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
         
         
                        # getting x and y coordinates transformed into the [-1,1] and [-1,1] intervals
                        [x_fct_L, y_fct_L] = xy_fct( left_coords[:,0], left_coords[:,1] )
                        [x_fct_R, y_fct_R] = xy_fct( right_coords[:,0], right_coords[:,1] )
         
                        # computing the Jacobian and the determinant of the left and right children of the parent element
                        J_left = jacobian_mat( left_coords[:,0], left_coords[:,1] )
                        J_right = jacobian_mat( right_coords[:,0], right_coords[:,1] )
                        detJ_left = lambda eps,niu: determinant(J_left)(eps,niu)
                        detJ_right = lambda eps,niu: determinant(J_right)(eps,niu)
             
                        # scaling factor
                        ##x1 = abs(loc - coords[0,0])
                        ##x2 = abs(coords[1,0] - loc)
                        #x1 = abs(left_coords[1,0]- coords[0,0])
                        #x2 = abs(coords[1,0] - left_coords[1,0])
                        s_factor = 1#( 2 * min(x1,x2)/(x1+x2)) ** 2
        
            
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
                                                #U[left_nodes[0],0] * Nbasis[0](x,y) +
                                                #U[right_nodes[1],0] * Nbasis[1](x,y) +
                                                #U[right_nodes[2],0] * Nbasis[2](x,y) +
                                                #U[left_nodes[3],0] * Nbasis[3](x,y) ) #* detJ_left
    
         
                        left_coords = p[left_nodes,:]
                        x_coords_L = left_coords[:,0]
                        y_coords_L = left_coords[:,1]
# canceling the norm computation
#                         # create the x = f(epsilon,niu) and y = g(epsilon,niu) functions
#                         # for transformation from the parametric element to phisycal element
#                         # of the Gauss nodes ui
#                         [x_transform_fct_L,y_transform_fct_L] = xy_fct(x_coords_L,y_coords_L)            
#         
#                         el_sum_L =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_L,y_transform_fct_L,uh_elem_L,detJ_left)
#         
#             
        
                        #tmp_array = np.arange(coords[0,0],left_coords[1,0],0.01)
                        #tx_L = np.append(tmp_array,left_coords[1,0])
                        #sfct_L = uh_elem_L(tx_L,0.03*yloc+0.57)
                        #sfct_L = uh_elem_L(tx_L,loc)
        
                        # on the right side of the interface
                        uh_elem_R = lambda x,y: (
                                                U[right_nodes[0],0] * psi_btm_R(x,y) * factor_S +
                                                U[right_nodes[3],0] * psi_upr_R(x,y) * factor_N +
                                                U[nodes[0],0] * Nbasis[0](x,y) +
                                                U[nodes[1],0] * Nbasis[1](x,y) + 
                                                U[nodes[2],0] * Nbasis[2](x,y) +
                                                U[nodes[3],0] * Nbasis[3](x,y) )
        
                                                #U[left_nodes[0],0] * Nbasis[0](x,y) +
                                                #U[right_nodes[1],0] * Nbasis[1](x,y) +
                                                #U[right_nodes[2],0] * Nbasis[2](x,y) +
                                                #U[left_nodes[3],0] * Nbasis[3](x,y)
                                                #)# * detJ_right
        
                        right_coords = p[right_nodes,:]
                        x_coords_R = right_coords[:,0]
                        y_coords_R = right_coords[:,1]

# canceling the norm computation            
#                         [x_transform_fct_R,y_transform_fct_R] = xy_fct(x_coords_R,y_coords_R)            
#                 
#                 
#                         el_sum_R =  gauss_integration(ui,wi,UConf,pConf,tConf,x_transform_fct_R,y_transform_fct_R,uh_elem_R,detJ_right)
#         
        
                        #tx_R = np.arange(left_coords[1,0],coords[1,0]+0.01,0.01)
                        #sfct_R = uh_elem_R(tx_R,0.03*yloc+0.57)
                        #sfct_R = uh_elem_R(tx_R,loc)
        
        
#                         all_elems_sum = all_elems_sum + el_sum_R + el_sum_L;
        
                        Usolution[nodes[0],0] = uh_elem_L( p[nodes[0],0], p[nodes[0],1]  )
                        Usolution[nodes[1],0] = uh_elem_R( p[nodes[1],0], p[nodes[1],1]  )
                        Usolution[nodes[2],0] = uh_elem_R( p[nodes[2],0], p[nodes[2],1]  )
                        Usolution[nodes[3],0] = uh_elem_L( p[nodes[3],0], p[nodes[3],1]  )
                        Usolution[nodes[4],0] = uh_elem_L( p[nodes[4],0], p[nodes[4],1]  )
                        Usolution[nodes[5],0] = uh_elem_L( p[nodes[5],0], p[nodes[5],1]  )
    
    
                        if y0 <= yloc and yloc <= y1:
        
                            tx = np.arange(coords[0,0],coords[1,0]+0.00001,0.001)
                            txL = []
                            txR = []
        
                            for idx in range(0,len(tx)):
        
                                if point_in_on_poly( tx[idx], yloc, left_coords):
                                    txL.append(tx[idx])
                                
                                if point_in_on_poly( tx[idx], yloc, right_coords):
                                    txR.append(tx[idx])
                                
                            txL = np.array(txL)
                            txR = np.array(txR)
    #                         pylab.plot( txL, uh_elem_L(txL, yloc))
    #                         pylab.plot( txR, uh_elem_R(txR, yloc))
        

#    print Usolution
    print 'Writing VTK file...' 
    print_vtk_file(p,Usolution,polygonList)
    print ' Done.'
    rtx = np.arange(0,1,0.01)
    #yidx = []
    ytriangle = []
    for idx  in range(0, len(rtx)):
        #yid = uexact(rtx[idx], xloc,k1,k2)
        #yidx.append(yid)
        ytriangle_Fct = found_in_FEM(rtx[idx],yloc,pConf,tConf,UConf)
        ytriangle.append(ytriangle_Fct)

    #pylab.plot(rtx,yidx)
#     pylab.plot(rtx,ytriangle)

    #===========================================================================
    # pylab.xlabel('x')
    # pylab.ylabel('f(x,y)')
    # pylab.title('Solution for a fixed y')
    # pylab.show()
    #===========================================================================

    return  math.sqrt(all_elems_sum)

def print_vtk_file(p,Usolution,plist):
    
    P = len(p)
    filename = 'dataset' + str(P) + 'points.vtk'
    target = open(filename,'w')
    target.write('# vtk DataFile Version 3.1 \n')
    target.write('Circle example \n')
    target.write('ASCII \n')
    target.write('DATASET POLYDATA \n')
    str1 = 'POINTS ' +  str(P) + ' FLOAT \n'
    target.write(str1)

#    print plist
#    print '-------------'
#    print p
#    print len( sum (plist, [] ) ) + len(plist)

    for i in range(0,P):
        stri = str(p[i,0]) + '  ' + str(p[i,1]) + '  ' + str(Usolution[i,0]) + ' \n'
#        print 'i = ', i, ' and ',stri
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

    str3 = '\nPOINT_DATA ' + str(P) + ' \n'
    target.write(str3)
    target.write('SCALARS Temperature FLOAT \n')
    target.write('LOOKUP_TABLE default \n')
    for z in range(0,len(Usolution)):
        strz = str(Usolution[z,0])
        target.write(strz + ' \n')

    target.close()

def coord_enrich_comp(root, midPoint):
    if root.enrichNodes[0].x != root.enrichNodes[1].x:
        xx = [root.enrichNodes[0].x, root.enrichNodes[1].x, root.enrichNodes[2].x]
        x_n = [root.enrichNodes[0].x, root.enrichNodes[1].x, root.enrichNodes[2].x]
        yy = [0,0,0]
        y_n = [root.enrichNodes[0].y, root.enrichNodes[1].y, root.enrichNodes[2].y]
        xx.sort()
        for i in range(0,3):
            indx = x_n.index(xx[i])
            yy[i] = y_n[indx]  
            
        
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
    
    return coord_enrich
        
        
def NW_corner(p,ui,wi,k1,k2,nodess,root,image):
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

    
#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)

#    R = 1.0 / 3.0 

#    if cornerD <= R * R:
#        K_cst = [k1,k1,k1,k2]
#    else:
#        K_cst = [k2,k2,k2,k1]
        
#     cornerD_s = f_circle_s(x0,y1)
#     cornerD_c1 = f_circle1(x0,y1)
#     cornerD_c2 = f_circle2(x0,y1)
#     Rs = 1.0/3.0
#     R1 = 1.0/6.0
#     R2 = 1.0/6.0
# 
#     if (cornerD_s <= Rs * Rs) or (cornerD_c1 <= R1 * R1) or (cornerD_c2 <= R2 * R2) :
#         K_cst = [k1,k1,k1,k2]
#     else:
#         K_cst = [k2,k2,k2,k1]
    
    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
    #if NW
    if ( is_in_same_bin(pxVal1,pxVal4) == False and pxVal1 > binBnd[1] and
        ( is_in_same_bin(pxVal4,pxVal2)==True and is_in_same_bin(pxVal2,pxVal3)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]
        
    #if point_in_on_poly(x0,y1,domainInclusion) == True:
    #    K_cst = [k1,k1,k1,k2]
    #else: 
    #    K_cst = [k2,k2,k2,k1]


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

        midPoint = Coordinate((root.enrichNodes[0].x +root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
        
        coord_enrich = coord_enrich_comp(root, midPoint)
        
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
            

            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)
            integral4 = simpson_rule(Kefunc4)

            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
            fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

            Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) + simpson_rule(fv4)

            #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]


def SW_corner(p,ui,wi,k1,k2,nodess, root, image):
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

#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)

#    R = 1.0 / 3.0

#    if cornerA <= R * R:
#        K_cst = [k2,k1,k1,k1]
#    else:
#        K_cst = [k1,k2,k2,k2]

#    cornerA_s = f_circle_s(x0,y0)
#    cornerA_c1 = f_circle1(x0,y0)
#    cornerA_c2 = f_circle2(x0,y0)

#    Rs = 1.0/3.0
#    R1 = 1.0/6.0
#    R2 = 1.0/6.0
#
#    if (cornerA_s <= Rs * Rs) or (cornerA_c1 <= R1 * R1) or (cornerA_c2 <= R2 * R2):
#        K_cst = [k2,k1,k1,k1]
#    else:
#        K_cst = [k1,k2,k2,k2]

    #if point_in_on_poly(x0,y0,domainInclusion) == True:
    #    K_cst = [k2,k1,k1,k1]
    #else:
    #    K_cst = [k1,k2,k2,k2]
        

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
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
        
        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)

        coord_enrich = coord_enrich_comp(root, midPoint)
                
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
            
            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)
            integral4 = simpson_rule(Kefunc4)
            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
        # construct the local load vector
        fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
        fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

        Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) + simpson_rule(fv4)

    #NEUMANN BCS are zero - code not inserted here
    return [K,Fe]

def NE_corner(p,ui,wi,k1,k2,nodess,root,image):
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

#    if x0==0.25 and y1==0.25:
#        print 'triangle 1', nodes1, coords1
#        print 'triangle 2', nodes2, coords2
#        print 'triangle 3', nodes3, coords3
#        print 'triangle 4', nodes4, coords4
#        print 'element', nodes, coords


#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)

#    R = 1.0 / 3.0

#    if cornerC <= R * R:
#        K_cst = [k1,k1,k1,k2]
#    else:
#        K_cst = [k2,k2,k2,k1]

#     cornerC_s = f_circle_s(x1,y1)
#     cornerC_c1 = f_circle1(x1,y1)
#     cornerC_c2 = f_circle2(x1,y1)
# 
#     Rs = 1.0/3.0
#     R1 = 1.0/6.0
#     R2 = 1.0/6.0
#     
#     if (cornerC_s <= Rs * Rs) or (cornerC_c1 <= R1 * R1) or (cornerC_c2 <= R2 * R2):
#         K_cst = [k1,k1,k1,k2]
#     else:
#         K_cst = [k2,k2,k2,k1]
        
    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
    #if NE
    if ( is_in_same_bin(pxVal1,pxVal2) == False and pxVal2 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal3)==True and is_in_same_bin(pxVal4,pxVal3)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]

    #if point_in_on_poly(x1,y1,domainInclusion) == True:
    #    K_cst = [k1,k1,k1,k2]
    #else:
    #    K_cst = [k2,k2,k2,k1]
       
       
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
        
        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)

        coord_enrich = coord_enrich_comp(root, midPoint)
                   
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
            
            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)
            integral4 = simpson_rule(Kefunc4)

            K[i,j] = integral1 + integral2 + integral3 + integral4

        # construct the local matrix and local components of the load vector
    for i in range(0,6):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
            fv4 = lambda e,n: rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

            Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) + simpson_rule(fv4)

    #NEUMANN BCS are zero - code not inserted here
    return [K,Fe]

def SE_corner(p,ui,wi,k1,k2,nodess,root,image):
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

    
#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)

#    R = 1.0 / 3.0

#    if cornerB <= R * R:
#        K_cst = [k1,k1,k1,k2]
#    else:
#        K_cst = [k2,k2,k2,k1]
# 
#     Rs = 1.0/3.0
#     R1 = 1.0/6.0
#     R2 = 1.0/6.0
# 
#     cornerB_s = f_circle_s(x1,y0)
#     cornerB_c1 = f_circle1(x1,y0)
#     cornerB_c2 = f_circle2(x1,y0)
# 
#     if (cornerB_s <= Rs * Rs) or (cornerB_c1 <= R1 * R1) or (cornerB_c2 <= R2 * R2):
#         K_cst = [k1,k1,k1,k2]
#     else:
#         K_cst = [k2,k2,k2,k1]

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
    #if SE
    if ( is_in_same_bin(pxVal3,pxVal4) == False and pxVal3 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)==True and is_in_same_bin(pxVal2,pxVal4)) ):
        K_cst = [k1,k1,k1,k2]
    else:
        K_cst = [k2,k2,k2,k1]

#if point_in_on_poly(x1,y0,domainInclusion) == True:
    #    K_cst = [k1,k1,k1,k2]
    #else:
    #    K_cst = [k2,k2,k2,k1]


#    if len(root.enrichNodes) >2:
#        print len(root.enrichNodes)
#        print 'in  SE SE SE SE SE SE-----------------------------------'
#        root.printRect() 
#        print coords
#        print coords4
#        print root.enrichNodes[0].x,root.enrichNodes[0].y
#        print root.enrichNodes[1].x,root.enrichNodes[1].y 
#        print root.enrichNodes[2].x,root.enrichNodes[2].y
        
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
       
        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)

        coord_enrich = coord_enrich_comp(root, midPoint)
        
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
            
            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)
            integral4 = simpson_rule(Kefunc4)

            K[i,j] = integral1 + integral2 + integral3 + integral4

    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)
        fv4 = lambda e,n:  rhs(e,n) * Nbasis_4[i](x_fct_4(e,n),y_fct_4(e,n)) * det_J4(e,n)

        Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) + simpson_rule(fv4)

            #NEUMANN BCS are zero - code not inserted here

#     print 'South-East corner', K, Fe
    return [K,Fe]


def East_edge(p,ui,wi,k1,k2,nodess,root,image):
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

#    # if SE - East
#    if ( (point_in_on_poly(x1,y0,domainInclusion) == True) and
#        ( not(point_in_on_poly(x0,y1,domainInclusion))==True and not(point_in_on_poly(x1,y1,domainInclusion))==True) ):
#        K_cst = [k2,k1,k1]
#    else:
#        if ( (point_in_on_poly(x0,y1,domainInclusion) == True and point_in_on_poly(x1,y1,domainInclusion)==True) and
#            not(point_in_on_poly(x1,y0,domainInclusion))==True ):
#            K_cst = [k1,k2,k2]
#
#    # if NE - East
#    if ( point_in_on_poly(x1,y1,domainInclusion)==True and
#        ( not(point_in_on_poly(x0,y0,domainInclusion))==True and not(point_in_on_poly(x1,y0,domainInclusion))==True) ):
#        K_cst = [k1,k1,k2]
#    else:
#        if ( (point_in_on_poly(x0,y0,domainInclusion)==True and point_in_on_poly(x1,y0,domainInclusion)==True) and
#            not(point_in_on_poly(x1,y1,domainInclusion))==True ):
#            K_cst = [k2,k2,k1]

    #print 'MANY VARIABLES',nodes,nodess[4]
    #print 'coordsinats', coords,p[nodess[4]]
    #print    point_in_on_poly(x0,y1,domainInclusion),point_in_on_poly(x1,y1,domainInclusion)
    #print point_in_on_poly(x1,y0,domainInclusion)

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    pxVal14 = image.GetPixel(int( (p1.x+p4.x)/2.0),int( (p1.y+p4.y)/2.0) )
    pxVal12 = image.GetPixel(int( (p1.x+p2.x)/2.0),int( (p1.y+p2.y)/2.0) )
    pxVal23 = image.GetPixel(int( (p2.x+p3.x)/2.0),int( (p2.y+p3.y)/2.0) )
    pxVal34 = image.GetPixel(int( (p3.x+p4.x)/2.0),int( (p3.y+p4.y)/2.0) )
    
    
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

    print 'East edge: K = ', K_cst, len(root.enrichNodes)

    if len(root.enrichNodes) == 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
        coord_enrich = coord_enrich_comp(root, midPoint)
        
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
            
            
            lOrd = [0,1,2]
            vec1_x = [ coords1[ lOrd[0],0], coords1[lOrd[1],0], coords1[lOrd[2],0], (coords1[ lOrd[0],0] + coords1[lOrd[1],0])/2.0, coord_enrich.x, (coords1[lOrd[0],0] + coords1[lOrd[2],0])/2.0  ]
            vec1_y = [ coords1[ lOrd[0],1], coords1[lOrd[1],1], coords1[lOrd[2],1], (coords1[ lOrd[0],1] + coords1[lOrd[1],1])/2.0, coord_enrich.y, (coords1[lOrd[0],1] + coords1[lOrd[2],1])/2.0  ]

            [x_fct_1, y_fct_1] = tri_xy_fct_quadratic( vec1_x, vec1_y )
            J1 = tri_jacobian_mat_quadratic( vec1_x, vec1_y )
                        
            lOrd = [2,0,1] # local order    
            vec2_x = [ coords2[ lOrd[0],0], coords2[lOrd[1],0], coords2[lOrd[2],0], (coords2[ lOrd[0],0] + coords2[lOrd[1],0])/2.0, coord_enrich.x, (coords2[lOrd[0],0] + coords2[lOrd[2],0])/2.0  ]
            vec2_y = [ coords2[ lOrd[0],1], coords2[lOrd[1],1], coords2[lOrd[2],1], (coords2[ lOrd[0],1] + coords2[lOrd[1],1])/2.0, coord_enrich.y, (coords2[lOrd[0],1] + coords2[lOrd[2],1])/2.0  ]
    
            [x_fct_2, y_fct_2] = tri_xy_fct_quadratic( vec2_x, vec2_y )
            J2 = tri_jacobian_mat_quadratic( vec2_x, vec2_y )

        
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

            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)

            K[i,j] = integral1 + integral2 + integral3 

    # construct the local matrix and local components of the load vector
    for i in range(0,5):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

        Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) 

        #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]
    

def South_edge(p,ui,wi,k1,k2,nodess,root,image):
    
    Ke = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    nodes1 = [nodess[0],nodess[4],nodess[3]]
    nodes2 = [nodess[4],nodess[2],nodess[3]]
    nodes3 = [nodess[4],nodess[1],nodess[2]]
#     nodes3 = [nodess[1],nodess[4],nodess[2]]

#     print nodes
#     print nodes1
#     print nodes2
#     print nodes3

    coords = p[nodes]
    coords1 = p[nodes1]
    coords2 = p[nodes2]
    coords3 = p[nodes3]

#     print nodes
#     print nodes1
#     print nodes2
#     print nodes3
#     print coords
#     print coords1
#     print coords2
#     print coords3

    x0 = coords[0,0]
    x1 = coords[1,0]
    y0 = coords[0,1]
    y1 = coords[2,1]

#    # if SE-South
#    if ( (point_in_on_poly(x1,y0,domainInclusion) == True) and
#        (not(point_in_on_poly(x0,y0,domainInclusion)) == True and not(point_in_on_poly(x0,y1,domainInclusion))==True) ):
#        K_cst = [k1,k1,k2]
#    else:
#        if ( (point_in_on_poly(x0,y0,domainInclusion)==True and point_in_on_poly(x0,y1,domainInclusion) == True) and
#            ( not(point_in_on_poly(x1,y0,domainInclusion)==True)) ):
#            K_cst = [k2,k2,k1]
#
#    # if SW-South
#    if ( (point_in_on_poly(x0,y0,domainInclusion)==True) and 
#        (not(point_in_on_poly(x1,y1,domainInclusion))==True and not(point_in_on_poly(x1,y0,domainInclusion)) ==True) ):
#        K_cst = [k2,k1,k1]
#    else:
#        if ( (point_in_on_poly(x1,y1,domainInclusion)==True and point_in_on_poly(x1,y0,domainInclusion)==True) and
#            ( not(point_in_on_poly(x0,y0,domainInclusion)) == True) ):
#            K_cst = [k1,k2,k2]

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    pxVal14 = image.GetPixel(int( (p1.x+p4.x)/2.0),int( (p1.y+p4.y)/2.0) )
    pxVal12 = image.GetPixel(int( (p1.x+p2.x)/2.0),int( (p1.y+p2.y)/2.0) )
    pxVal23 = image.GetPixel(int( (p2.x+p3.x)/2.0),int( (p2.y+p3.y)/2.0) )
    pxVal34 = image.GetPixel(int( (p3.x+p4.x)/2.0),int( (p3.y+p4.y)/2.0) )
    
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

#     print 'South edge: K = ', K_cst

    
    if len(root.enrichNodes) == 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
        coord_enrich = coord_enrich_comp(root, midPoint)
        
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

#     factor_S = 1
    
    Nbasis_1 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis1[1](x,y)]
    Nbasis_2 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis2[0](x,y)]
    Nbasis_3 = [Nbasis[0],Nbasis[1],Nbasis[2],Nbasis[3], lambda x,y: factor_S * Nbasis3[0](x,y)]

    Nx_1 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx1[1](x,y)]
    Nx_2 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx2[0](x,y)]
    Nx_3 = [Nx[0],Nx[1],Nx[2],Nx[3], lambda x,y: factor_S * Nx3[0](x,y)]
    
    Ny_1 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny1[1](x,y)]
    Ny_2 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny2[0](x,y)]
    Ny_3 = [Ny[0],Ny[1],Ny[2],Ny[3], lambda x,y: factor_S * Ny3[0](x,y)]

#     ggg = lambda e,n:( (-2+2*n)*(-2+2*n) + (-2+2.0/3.0*e)*(-2+2.0/3.0*e)) * ( (-0.5+0.667)/2.0)
#     fff = lambda x,y: (-2+4*y)*(-2+4*y) + (-4+4*x)*(-4+4*x)

#     fff = lambda x,y: (2-4*y)*(2-4*y) + (2-4*x)*(2-4*x)
#     fff = lambda x,y: (4*y)*(4*y) + (-2+4*x)*(-2+4*x)
#     fff = lambda x,y: (-4*y)*(-4*y) + (4-4*x)*(4-4*x)
#     fff = lambda x,y: (2-4*y)*(5.98802395) + (2-4*x)*(0)
#     fff = lambda x,y: (2-4*y)*(0) + (2-4*x)*(-2)
#     fff = lambda x,y: (2-4*y)*(-3.003003) + (2-4*x)*(0)
    
    
#     fff = lambda x,y: (Nx1[1](x,y)) *(Nx1[1](x,y)) + (Ny1[1](x,y)) *(Ny1[1](x,y))
    
#     kefunct1 = lambda e,n: fff(x_fct_1(e,n), y_fct_1(e,n)) *  det_J1(e,n)
#     kefunct2 = lambda e,n: fff(x_fct_2(e,n), y_fct_2(e,n)) *  det_J2(e,n)
#     kefunct3 = lambda e,n: fff(x_fct_3(e,n), y_fct_3(e,n)) *  det_J3(e,n) * 10

#     print 'Kefunct: simpson ',simpson_rule(kefunct1) , simpson_rule(kefunct2) , simpson_rule(kefunct3)
#     print 'simspon(ggg)',simpson_rule(ggg)
    
    
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

            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)

#             if i == 4 and j == 1:
#                 print'=============>>>>>>', integral1 , integral2 , integral3
            Ke[i,j] = integral1 + integral2 + integral3 

    # construct the local matrix and local components of the load vector
    for i in range(0,5):
        # construct the local load vector
        fv1 = lambda e,n:  rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
        fv2 = lambda e,n:  rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
        fv3 = lambda e,n:  rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

        Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) 

        #NEUMANN BCS are zero - code not inserted here

#     print Ke
#     print Fe
#     print 'South EDGE:\n', Ke, Fe
    return [Ke,Fe]

def North_edge(p,ui,wi,k1,k2,nodess,root,image):
    K = numpy.zeros((5,5))
    Fe = np.zeros((5,1))

        
    nodes = [nodess[0],nodess[1],nodess[2],nodess[3]]
    # nodes are [0,1,2,3,4,5] or SW,SE,NE,NW,SMid,EMid (upper left corner cut)
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

    #cornerA = f_circle(x0,y0)
    #cornerB = f_circle(x1,y0)
    #cornerC = f_circle(x1,y1)
    #cornerD = f_circle(x0,y1)

    #R = 1.0 / 3.0 

    #if cornerD <= R * R:
    #    K_cst = [k1,k1,k1,k2]
    #else:
    #    K_cst = [k2,k2,k2,k1]
    
#    # if NW - North
#    if ( point_in_on_poly(x0,y1,domainInclusion) == True and
#        (not(point_in_on_poly(x1,y0,domainInclusion)) == True and not(point_in_on_poly(x1,y1,domainInclusion))==True )):
#        K_cst = [k2,k1,k1]
#    else: 
#        if ( not(point_in_on_poly(x0,y1,domainInclusion)) == True and
#            (point_in_on_poly(x1,y0,domainInclusion) == True and point_in_on_poly(x1,y1,domainInclusion)==True) ):
#            K_cst = [k2,k1,k1]
#
#    # if NE - North 
#    if ( point_in_on_poly(x1,y1,domainInclusion)==True and
#        ( not(point_in_on_poly(x0,y0,domainInclusion))==True and not(point_in_on_poly(x0,y1,domainInclusion))==True ) ):
#        K_cst = [k1,k1,k2]
#    else:
#        if ( not(point_in_on_poly(x1,y1,domainInclusion))==True and
#            ( point_in_on_poly(x0,y0,domainInclusion)==True and point_in_on_poly(x0,y1,domainInclusion)==True ) ):
#            K_cst = [k2,k2,k1]

#    if point_in_on_poly(x1,y1,domainInclusion) and not(point_in_on_poly(x0,y1,domainInclusion)) and not(point_in_on_poly(x0,y0,domainInclusion)):
#        K_cst = [k1,k1,k2]
#    else:
#        if not(point_in_on_poly(x1,y1,domainInclusion)) and point_in_on_poly(x0,y1,domainInclusion) and point_in_on_poly(x0,y0,domainInclusion):
#         K_cst = [k2,k2,k1]
#
#    if point_in_on_poly(x0,y1,domainInclusion) and not(point_in_on_poly(x1,y0,domainInclusion)) and not(point_in_on_poly(x1,y1,domainInclusion)):
#        K_cst = [k2,k1,k1]
#    else:
#        if not(point_in_on_poly(x0,y1,domainInclusion)) and point_in_on_poly(x1,y0,domainInclusion) and point_in_on_poly(x1,y1,domainInclusion):
#            K_cst = [k1,k2,k2]
    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    pxVal14 = image.GetPixel(int( (p1.x+p4.x)/2.0),int( (p1.y+p4.y)/2.0) )
    pxVal12 = image.GetPixel(int( (p1.x+p2.x)/2.0),int( (p1.y+p2.y)/2.0) )
    pxVal23 = image.GetPixel(int( (p2.x+p3.x)/2.0),int( (p2.y+p3.y)/2.0) )
    pxVal34 = image.GetPixel(int( (p3.x+p4.x)/2.0),int( (p3.y+p4.y)/2.0) )
    
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
    
            
    if len(root.enrichNodes) == 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
        coord_enrich = coord_enrich_comp(root, midPoint)
        
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

            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)

            K[i,j] = integral1 + integral2 + integral3 

        # construct the local matrix and local components of the load vector
    for i in range(0,5):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

            Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) 

            #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]



def West_edge(p,ui,wi,k1,k2,nodess,root,image):
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

    #cornerA = f_circle(x0,y0)
    #cornerB = f_circle(x1,y0)
    #cornerC = f_circle(x1,y1)
    #cornerD = f_circle(x0,y1)

    #R = 1.0 / 3.0 

    #if cornerD <= R * R:
    #    K_cst = [k1,k1,k1,k2]
    #else:
    #    K_cst = [k2,k2,k2,k1]
    

#    #if NW - West
#    if ( point_in_on_poly(x0,y1,domainInclusion) == True and
#        (not(point_in_on_poly(x0,y0,domainInclusion))==True and not(point_in_on_poly(x1,y0,domainInclusion))==True) ):
#        K_cst = [k1,k1,k2]
#    else: 
#        if ( not(point_in_on_poly(x0,y1,domainInclusion)) == True and
#            (point_in_on_poly(x0,y0,domainInclusion)==True and point_in_on_poly(x1,y0,domainInclusion)==True) ):
#            K_cst = [k2,k2,k1]
#
#    # if SW - West
#    if ( point_in_on_poly(x0,y0,domainInclusion) == True and
#        (not(point_in_on_poly(x0,y1,domainInclusion)) ==True and not(point_in_on_poly(x1,y1,domainInclusion))==True )  ):
#        K_cst = [k2,k1,k1]
#    else:
#        if ( not(point_in_on_poly(x0,y0,domainInclusion)) == True and
#            ( point_in_on_poly(x0,y1,domainInclusion) ==True and point_in_on_poly(x1,y1,domainInclusion)==True )  ):
#             K_cst = [k1,k2,k2]
    
    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    pxVal14 = image.GetPixel(int( (p1.x+p4.x)/2.0),int( (p1.y+p4.y)/2.0) )
    pxVal12 = image.GetPixel(int( (p1.x+p2.x)/2.0),int( (p1.y+p2.y)/2.0) )
    pxVal23 = image.GetPixel(int( (p2.x+p3.x)/2.0),int( (p2.y+p3.y)/2.0) )
    pxVal34 = image.GetPixel(int( (p3.x+p4.x)/2.0),int( (p3.y+p4.y)/2.0) )
    
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
             
#     print 'West edge: K = ', K_cst

            
#    cornerA_s = f_circle_s(x0,y0)
#    cornerB_s = f_circle_s(x1,y0)
#    cornerD_s = f_circle_s(x0,y1)
#    cornerD_c1 = f_circle1(x0,y1)
#    cornerD_c2 = f_circle2(x0,y1)
#    Rs = 1.0/3.0
#    R1 = 1.0/6.0
#    R2 = 1.0/6.0
#
#    if (cornerD_s <= Rs * Rs) or (cornerD_c1 <= R1 * R1) or (cornerD_c2 <= R2 * R2) :
#        K_cst = [k1,k1,k1,k2]
#    else:
#        K_cst = [k2,k2,k2,k1]
        
    if len(root.enrichNodes) == 2:
        [x_fct_1, y_fct_1] = tri_xy_fct( coords1[:,0], coords1[:,1] )
        [x_fct_2, y_fct_2] = tri_xy_fct( coords2[:,0], coords2[:,1] )
        [x_fct_3, y_fct_3] = tri_xy_fct( coords3[:,0], coords3[:,1] )
    
        J1 = tri_jacobian_mat( coords1[:,0], coords1[:,1] )
        J2 = tri_jacobian_mat( coords2[:,0], coords2[:,1] )
        J3 = tri_jacobian_mat( coords3[:,0], coords3[:,1] )

    if len(root.enrichNodes) == 3:

        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)
        coord_enrich = coord_enrich_comp(root, midPoint)
        
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

            integral1 = simpson_rule(Kefunc1)
            integral2 = simpson_rule(Kefunc2)
            integral3 = simpson_rule(Kefunc3)

            K[i,j] = integral1 + integral2 + integral3 

        # construct the local matrix and local components of the load vector
    for i in range(0,5):
            # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_1[i](x_fct_1(e,n),y_fct_1(e,n)) * det_J1(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_2[i](x_fct_2(e,n),y_fct_2(e,n)) * det_J2(e,n)
            fv3 = lambda e,n: rhs(e,n) * Nbasis_3[i](x_fct_3(e,n),y_fct_3(e,n)) * det_J3(e,n)

            Fe[i] = simpson_rule(fv1) + simpson_rule(fv2) + simpson_rule(fv3) 

            #NEUMANN BCS are zero - code not inserted here

    return [K,Fe]

def horizontal_cut(p,ui,wi,k1,k2,nodes,root,image):

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


#    if x0 == 0.28125 and y1 == 0.25:
#        print 'bottom ', bottom_nodes, bottom_coords
#        print 'top', top_nodes,top_coords


    Pe = np.zeros((4,4))
    Pe[:,0] = np.ones((4,1)).transpose()
    Pe[:,1] = p[nodes[0:4],0]
    Pe[:,2] = p[nodes[0:4],1]
    Pe[:,3] = p[nodes[0:4],0]*p[nodes[0:4],1]

    C = np.linalg.inv(Pe)

    Nbasis = basisFct(C)
    Nx = derivX(C)
    Ny = derivY(C)

#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)
#    R = 1.0/3.0
#    if (cornerA <= R*R and cornerB <= R*R) and (cornerC > R*R and cornerD > R*R):
#        K_cst = [k2,k1]
#    else:
#        if (cornerA > R*R and cornerB > R*R) and (cornerC <= R*R and cornerD <= R*R):
#            K_cst = [k1,k2]
#        else:
#            print 'ERROR! inside horizontal interface of stiffness matrix computation'

#     Rs = 1.0/3.0
#     R1 = 1.0/6.0
#     R2 = 1.0/6.0
#     cornerA_s = f_circle_s(x0,y0)
#     cornerB_s = f_circle_s(x1,y0)
#     cornerC_s = f_circle_s(x1,y1)
#     cornerD_s = f_circle_s(x0,y1)
# 
#     cornerA_c1 = f_circle1(x0,y0)
#     cornerB_c1 = f_circle1(x1,y0)
#     cornerC_c1 = f_circle1(x1,y1)
#     cornerD_c1 = f_circle1(x0,y1)
# 
#     cornerA_c2 = f_circle2(x0,y0)
#     cornerB_c2 = f_circle2(x1,y0)
#     cornerC_c2 = f_circle2(x1,y1)
#     cornerD_c2 = f_circle2(x0,y1)
# 
# 
#     if ( ( (cornerA_s <= Rs*Rs and cornerB_s <= Rs*Rs) and (cornerC_s > Rs*Rs and cornerD_s > Rs*Rs) ) or
#             ( (cornerA_c1 <= R1*R1 and cornerB_c1 <= R1*R1) and (cornerC_c1 > R1*R1 and cornerD_c1 > R1*R1) ) or
#             ( (cornerA_c2 <= R2*R2 and cornerB_c2 <= R2*R2) and (cornerC_c2 > R2*R2 and cornerD_c2 > R2*R2) ) ):
#         K_cst = [k2,k1]
#     else:
#             K_cst = [k1,k2]

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
    #if SE-South
    if ( is_in_same_bin(pxVal1,pxVal2) == True and pxVal1 > binBnd[1] and
         is_in_same_bin(pxVal3,pxVal4) == True and
        (is_in_same_bin(pxVal1,pxVal4) == False and is_in_same_bin(pxVal2,pxVal3) == False) ):
        K_cst = [k2,k1]
    else:
        K_cst = [k1,k2]

#definition of rhombus corners
#    A = Point(0.5, 1.0/3.0)
#    B = Point(2.0/3.0, 0.5)
#    C = Point(1.0/3.0, 0.5)
#    D = Point(0.5, 2.0/3.0)
#    rhombus = [ (A.x, A.y), (B.x,B.y), (D.x,D.y), (C.x,C.y) ]

#    if point_in_poly(x0,y0,rhombus) and point_in_poly(x1,y0,rhombus):
#        K_cst = [k2, k1]
#    else:
#        K_cst = [k1, k2]
                
    #if point_in_on_poly(x0,y0,domainInclusion) and point_in_on_poly(x1,y0,domainInclusion):
    #    K_cst = [k2, k1]
    #else:
    #    K_cst = [k1, k2]

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
       
        midPoint = Coordinate((root.enrichNodes[0].x + root.enrichNodes[1].x)/2.0, (root.enrichNodes[0].y + root.enrichNodes[1].y)/2.0)

        coord_enrich = coord_enrich_comp(root, midPoint)
        
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
        
    detJ_top = lambda e,n: determinant(J_top)(e,n)
    detJ_bottom = lambda e,n: determinant(J_bottom)(e,n)
    
    Ke_Horiz = np.zeros((6,6))
    Fe_Horiz = np.zeros((6,1))
    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        for j in range(0,6):
            if nodes[i] >= nodes[j]:

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
    
                # construct the local load vector
            fv1 = lambda e,n: rhs(e,n) * Nbasis_top[i](x_fct_T(e,n),y_fct_T(e,n)) * detJ_top(e,n)
            fv2 = lambda e,n: rhs(e,n) * Nbasis_bottom[i](x_fct_B(e,n),y_fct_B(e,n)) * detJ_bottom(e,n)
    
            Fe_Horiz[i] = quad2d(fv1,-1,1,-1,1,ui,wi) + quad2d(fv2,-1,1,-1,1,ui,wi)


    return [Ke_Horiz,Fe_Horiz]


def vertical_cut(p,ui,wi,k1,k2,nodes,root,image):
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

#    cornerA = f_circle(x0,y0)
#    cornerB = f_circle(x1,y0)
#    cornerC = f_circle(x1,y1)
#    cornerD = f_circle(x0,y1)
#    R = 1.0/3.0

#    if (cornerA <= R*R and cornerD <= R*R) and (cornerB > R*R and cornerC > R*R):
#        K_cst = [k2,k1]
#    else:
#        if (cornerA > R*R and cornerD > R*R) and (cornerB <= R*R and cornerC <= R*R):
#            K_cst = [k1,k2]
#        else:
#            print 'ERROR inside vertical quad-quad stiffness matrix implementation'

#     Rs = 1.0/3.0
#     R1 = 1.0/6.0
#     R2 = 1.0/6.0
#     cornerA_s = f_circle_s(x0,y0)
#     cornerB_s = f_circle_s(x1,y0)
#     cornerC_s = f_circle_s(x1,y1)
#     cornerD_s = f_circle_s(x0,y1)
# 
#     cornerA_c1 = f_circle1(x0,y0)
#     cornerB_c1 = f_circle1(x1,y0)
#     cornerC_c1 = f_circle1(x1,y1)
#     cornerD_c1 = f_circle1(x0,y1)
# 
#     cornerA_c2 = f_circle2(x0,y0)
#     cornerB_c2 = f_circle2(x1,y0)
#     cornerC_c2 = f_circle2(x1,y1)
#     cornerD_c2 = f_circle2(x0,y1)
# 
#     if ( ( (cornerA_s <= Rs*Rs and cornerD_s <= Rs*Rs) and (cornerB_s > Rs*Rs and cornerC_s > Rs*Rs)) or
#         ( (cornerA_c1 <= R1*R1 and cornerD_c1 <= R1*R1) and (cornerB_c1 > R1*R1 and cornerC_c1 > R1*R1)) or
#         ( (cornerA_c2 <= R2*R2 and cornerD_c2 <= R2*R2) and (cornerB_c2 > R2*R2 and cornerC_c2 > R2*R2)) ) :
#         K_cst = [k2,k1]
#     else:
#         K_cst = [k1,k2]

    p1,p2,p3,p4 = root.rect
    
    pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
    pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
    pxVal3 = image.GetPixel(int(p3.x), int(p3.y));
    pxVal4 = image.GetPixel(int(p4.x), int(p4.y));
    
    #if SE-South
    if ( is_in_same_bin(pxVal1,pxVal4) == True and is_in_same_bin(pxVal2,pxVal3) == True and pxVal4 > binBnd[1] and
        ( is_in_same_bin(pxVal1,pxVal2)== False and is_in_same_bin(pxVal4,pxVal3) == False) ):
        K_cst = [k2,k1]
    else:
        K_cst = [k1,k2]
        

#    slanted_rectangle = polygonDef#[ (0.57,0.0), (1.0,0.0), (1.0,1.0), (0.6,1.0) ]
    #if point_in_on_poly(x0,y0,domainInclusion) and point_in_on_poly(x0,y1,domainInclusion) and point_in_on_poly(x1,y0,domainInclusion) and point_in_on_poly(x1,y1,domainInclusion):

    #if not(point_in_on_poly(x0,y0,domainInclusion)) and not(point_in_on_poly(x0,y1,domainInclusion)) and point_in_on_poly(x1,y0,domainInclusion) and point_in_on_poly(x1,y1,domainInclusion):
    #    K_cst = [k1, k2]
    #else:
    #    K_cst = [k2, k1]
    
    #if ( coords[0,0]<= left_coords[1,0] and left_coords[1,0] <= coords[1,0] ) and k1!=k2:
    #    K_cst = [k1,k2]
    #else:
    #    K_cst = [k2,k1]


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
    [x_fct_L, y_fct_L] = xy_fct( left_coords[:,0], left_coords[:,1] )
    [x_fct_R, y_fct_R] = xy_fct( right_coords[:,0], right_coords[:,1] )
    
    # computing the Jacobian and the determinant of the left and right children of the parent element
    J_left = jacobian_mat( left_coords[:,0], left_coords[:,1] )
    J_right = jacobian_mat( right_coords[:,0], right_coords[:,1] )
    detJ_left = lambda e,n: determinant(J_left)(e,n)
    detJ_right = lambda e,n: determinant(J_right)(e,n)
    
    
    Ke_Vertical = np.zeros((6,6))
    Fe_Vertical = np.zeros((6,1))
    # construct the local matrix and local components of the load vector
    for i in range(0,6):
        for j in range(0,6):
            if nodes[i] >= nodes[j]:

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
                #print i,j,left_integral, right_integral
    
        # construct the local load vector
        fv1 = lambda e,n: rhs(e,n) * Nbasis_left[i](x_fct_L(e,n),y_fct_L(e,n)) * detJ_left(e,n)
        fv2 = lambda e,n: rhs(e,n) * Nbasis_right[i](x_fct_R(e,n),y_fct_R(e,n)) * detJ_right(e,n)
    
        Fe_Vertical[i] = quad2d(fv1,-1,1,-1,1,ui,wi) + quad2d(fv2,-1,1,-1,1,ui,wi)
    
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



#definition of rhombus corners
#A = Point(0.5, 1.0/3.0)
#B = Point(2.0/3.0, 0.5)
#C = Point(1.0/3.0, 0.5)
#D = Point(0.5, 2.0/3.0)
# lower rhombus
#A = Point(0.5, 0.0)
#B = Point(2.0/3.0, 1.0/4.0)
#C = Point(1.0/3.0, 1.0/4.0)
#D = Point(1.0/2.0, 1.0/2.0)

#rhombus = [ (A.x, A.y), (B.x,B.y), (D.x,D.y), (C.x,C.y) ]

#definition of hexagon corners
#A = Point(1.0/4.0, 1.0/6.0)
#B = Point(3.0/4.0, 1.0/6.0)
#C = Point(5.0/6.0, 1.0/2.0)
#D = Point(3.0/4.0, 5.0/6.0)
#E = Point(1.0/4.0, 5.0/6.0)
#F = Point(1.0/6.0, 1.0/2.0)
#hexagon = [ (A.x,A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y), (E.x,E.y), (F.x,F.y) ]

#polygonDef = [(0.0,0.0),(0.57,0.0),(0.6,1.0),(0.0,1.0)]
# SE triangle is cut
#polygonDef = [(0.0,0.0),(0.9,0.0),(1.0,0.1),(1.0,1.0),(0.0,1.0)]
# NE triangle is cut
#polygonDef = [(0.0,0.0),(1.0,0.0),(1.0,0.9),(0.9,1.0),(0.0,1.0)]
# NW triangle is cut
#polygonDef = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.1,1.0),(0.0,0.9)]
# SW triangle is cut
#polygonDef = [(0.1,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),(0.0,0.1)]
#polygonDef = rhombus
# SE/NW Triangle with SE/NW corners cut: x=[2/3,1] to y= [0,1/3]
#polygonDef = [(0.0,0.0),(2.0/3.0,0.0),(1.0,1.0/3.0),(1.0,1.0),(0.0,1.0)]
#polygonDef = rhombus
# the top part of the domain for HORIZONTAL line splitting the domain in half
#polygonDef = domainInclusion#[(0.0,0.57),(1.0,0.6),(1.0,1.0),(0.0,1.0)]
#polygonDef = hexagon

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

def simpson_rule(f):

    I = 1.0/6.0 * f(1.0/2.0,1.0/2.0) + 1.0/6.0 * f(1.0/2.0,0.0) + 1.0/6.0 * f(0.0,1.0/2.0)
    return I

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
