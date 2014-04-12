import sys
import os
import math
import SimpleITK as sitk
from libFcts import *
from globalVars import * 
import string
from itertools import izip_longest
from itertools import groupby
import numpy
from copy import deepcopy 

from igfem2d import *
from meshgeneration import *
import scipy
from lgwt import *
from readFile import *
from findIntersection import *
# from myplot import *
import scipy.io


allnodes = 0
LIST = []
D = {}
D['0'] = {
          'R': {'Quadrant':'1', 'Direction':'H'},
          'L': {'Quadrant':'1', 'Direction':'L'},
          'D': {'Quadrant':'2', 'Direction':'H'},
          'U': {'Quadrant':'2', 'Direction':'U'},
          
          'RU': {'Quadrant':'3', 'Direction':'U'},
          'RD': {'Quadrant':'3', 'Direction':'H'},
          'LD': {'Quadrant':'3', 'Direction':'L'},
          'LU': {'Quadrant':'3', 'Direction':'LU'}
          }

D['1'] = {
          'R': {'Quadrant':'0', 'Direction':'R'},
          'L': {'Quadrant':'0', 'Direction':'H'},
          'D': {'Quadrant':'3', 'Direction':'H'},
          'U': {'Quadrant':'3', 'Direction':'U'},
                    
          'RU': {'Quadrant':'2', 'Direction':'RU'},
          'RD': {'Quadrant':'2', 'Direction':'R'},
          'LD': {'Quadrant':'2', 'Direction':'H'},
          'LU': {'Quadrant':'2', 'Direction':'U'}
          }

D['2'] = {
          'R': {'Quadrant':'3', 'Direction':'H'},
          'L': {'Quadrant':'3', 'Direction':'L'},
          'D': {'Quadrant':'0', 'Direction':'D'},
          'U': {'Quadrant':'0', 'Direction':'H'},
                    
          'RU': {'Quadrant':'1', 'Direction':'H'},
          'RD': {'Quadrant':'1', 'Direction':'D'},
          'LD': {'Quadrant':'1', 'Direction':'LD'},
          'LU': {'Quadrant':'1', 'Direction':'L'}
          }

D['3'] = {
          'R': {'Quadrant':'2', 'Direction':'R'},
          'L': {'Quadrant':'2', 'Direction':'H'},
          'D': {'Quadrant':'1', 'Direction':'D'},
          'U': {'Quadrant':'1', 'Direction':'H'},
                    
          'RU': {'Quadrant':'0', 'Direction':'R'},
          'RD': {'Quadrant':'0', 'Direction':'RD'},
          'LD': {'Quadrant':'0', 'Direction':'D'},
          'LU': {'Quadrant':'0', 'Direction':'H'}
          }

class ListNode():
    def __init__(self):
        self.next = None
        self.previous = None
        self.element = None
        
class DoubleLinkedList():
    def __init__(self):
        self.n = 0
        self.last = ListNode()
        self.first = self.last
        
    def append(self,element):
        self.last.element = element
        self.last.next = ListNode()
        tmp = self.last
        self.last = self.last.next
        self.last.previous = tmp
        self.n += 1
        
    def front(self):
        if self.n == 0:
            return None
        el = self.first.element
        self.first = self.first.next
        self.n -= 1
        return el

    def back(self):
        if self.n == 0: 
            return None
        
        el = self.last.previous.element
        self.last = self.last.previous
        self.last.next = ListNode()
        self.n -= 1
        return el
    
    def size(self):
        return self.n
    
    def elements(self):
        i = self.first
        while i.element:
            yield i.element
            i = i.next
            
class Node():
    ROOT = 0
    BRANCH = 1
    LEAF = 2
    MAX_DEPTH = 0
    
    def __init__(self, parent, rect, inImage,outImage,imageSize):

        self.imsize = imageSize
        self.parent = parent
        self.children = [None, None, None, None]
        self.has_children = False

#        self.maxdepth = 0
        if parent == None:
            self.depth = 0
            
        else:
            self.depth = parent.depth + 1
            if self.depth > Node.MAX_DEPTH:
                Node.MAX_DEPTH = self.depth
                
        
        self.ishomog = 1
        self.tlist = [] # contains the numbering of the nodes in the element
        self.tpix = []       
        self.nsew = [0,0,0,0]
        self.hn = [Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1)]
        
        self.rect = rect 
        self.index = '-1'
        [p1,p2,p3,p4] = rect
        
###        dx = abs(p1.x-p2.x)+1
###        dy = abs(p1.y-p4.y)+1
####        ind_x = round(imageSize[0]/dx)
####        ind_y = round(imageSize[1]/dy)
###        self.i = p2.x / dx
###        self.j = p4.y / dy


#        self.mat = 'Epoxy'
        self.enrichNodes = []
        
        if self.parent == None:
            self.type = Node.ROOT
#        #elif abs(p1.x - p2.x) <= MIN_SIZE:
        elif ( abs(p1.x - p2.x) <= MIN_SIZE or 
            (self.children[0]==None and self.children[1]== None and self.children[2] == None and self.children[3] == None) ):
#            print self.rect
            self.type = Node.LEAF
        else:
            self.type = Node.BRANCH
        
        self.outImage = outImage
        self.inImage = inImage
        
        Lx = set_interval(imageSize[0],self.depth)
        Ly = set_interval(imageSize[1],self.depth)
        self.i = list(Lx).index(p1.x)
        self.j = list(Ly).index(p1.y)
        
    def subdivide(self): 
    # this method subdivides a node recursively if some
    # division criterion is satisfied
    
#        if self.type == Node.LEAF:
#            return
        
        p1,p2,p3,p4 = self.rect
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid24 = find_mid_point(p2,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        rects = []
        rects.append((p1,cMid12,cMid24,cMid14)) #NW child
        rects.append((cMid12,p2,cMid23,cMid24)) #NE child
        rects.append((cMid14,cMid24,cMid34,p4)) #SW child
        rects.append((cMid24,cMid23,p3,cMid34)) #SE child

#        if p1.x == 0 and p2.x == 31 and p1.y == 467 and p3.y == 499:
#            self.printRect()
#            print self.division_criterion(rects[0], self.inImage, self.outImage),self.division_criterion(rects[1], self.inImage, self.outImage),self.division_criterion(rects[2], self.inImage, self.outImage),self.division_criterion(rects[3], self.inImage, self.outImage)
#            print [ has_inclusions(self.inImage,p1,p2,p3,p4)]
        for n in range(len(rects)):
            span = self.division_criterion(rects[n], self.inImage, self.outImage)
#            pi1,pi2,pi3,pi4 = self.children[n].rect

            if span == True:
                self.children[n] = self.getinstance(rects[n], self.inImage, self.outImage,imageSize)
                self.children[n].index = str(convert_to_base_4(tomorton(self.children[n].i, self.children[n].j)))
                diff_level = abs(len(self.children[n].index) - self.children[n].depth)
                if diff_level != 0:
                    self.children[n].index = '0'*diff_level + self.children[n].index
                
                p1r,p2r,p3r,p4r = rects[n]
                L1 = search_in(LIST,p1r,p2r,self.inImage)
                L2 = search_in(LIST,p2r,p3r,self.inImage)
                L3 = search_in(LIST,p4r,p3r,self.inImage)
                L4 = search_in(LIST,p1r,p4r,self.inImage)
        
 
                if len(L1) == 1:
                    L1 = L1[0]
                    if in_child_k(rects[n],L1) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L1]
                if len(L2) == 1:
                    L2 = L2[0]
                    if in_child_k(rects[n],L2) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L2]
                if len(L3) == 1:
                    L3 = L3[0]
                    if in_child_k(rects[n],L3) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L3]
                if len(L4) == 1:
                    L4 = L4[0]
                    if in_child_k(rects[n],L4) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L4]
        
                self.children[n].subdivide()
                
#        if self.children == [None,None,None,None]:# and root.parent != None:
#           self.has_children = False
#        elif self.children != [None,None,None,None]:
#            self.has_children = True       
#        else:
#            print 'why here?', self.children
        if ( self.children[0] != None and
             self.children[1] != None and
             self.children[2] != None and
             self.children[3] != None ):
            self.has_children = True
        else:
            self.has_children = False
            
    def divideOnce(self): 
    # this method divides once a given node
        
        p1,p2,p3,p4 = self.rect
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid24 = find_mid_point(p2,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        rects = []
        rects.append((p1,cMid12,cMid24,cMid14))
        rects.append((cMid12,p2,cMid23,cMid24))
        rects.append((cMid14,cMid24,cMid34,p4))
        rects.append((cMid24,cMid23,p3,cMid34))
        
        for n in range(len(rects)):
            span = self.division_criterionOnce(rects[n], self.inImage, self.outImage)
            if span == True:
                self.children[n] = self.getinstance(rects[n], self.inImage, self.outImage,imageSize)
                
                self.children[n].index = str(convert_to_base_4(tomorton(self.children[n].i, self.children[n].j)))
                diff_level = abs(len(self.children[n].index) - self.children[n].depth)
                if diff_level != 0:
                    self.children[n].index = '0'*diff_level + self.children[n].index 
        
                p1r,p2r,p3r,p4r = rects[n]
                
                L1 = linear_search(self.inImage,p1r,p2r);
                L2 = linear_search(self.inImage,p2r,p3r);
                L3 = linear_search(self.inImage,p4r,p3r);
                L4 = linear_search(self.inImage,p1r,p4r);
#                L1 = search_in(LIST,p1,p2,self.inImage)
#                L2 = search_in(LIST,p2,p3,self.inImage)
#                L3 = search_in(LIST,p4,p3,self.inImage)
#                L4 = search_in(LIST,p1,p4,self.inImage)
        
                l1 = ends_in_same_bin(self.inImage,p1r,p2r);
                l2 = ends_in_same_bin(self.inImage,p2r,p3r);
                l3 = ends_in_same_bin(self.inImage,p4r,p3r);
                l4 = ends_in_same_bin(self.inImage,p1r,p4r);
                
                
                if len(L1) == 1:
                    L1 = L1[0]
                    if in_child_k(rects[n],L1) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L1]
                        self.children[n].ishomog = 0
                if len(L2) == 1:
                    L2 = L2[0]
                    if in_child_k(rects[n],L2) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L2]
                        self.children[n].ishomog = 0
                if len(L3) == 1:
                    L3 = L3[0]
                    if in_child_k(rects[n],L3) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L3]
                        self.children[n].ishomog = 0
                if len(L4) == 1:
                    L4 = L4[0]
                    if in_child_k(rects[n],L4) == True:
                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L4]
                        self.children[n].ishomog = 0
        
#                 # NW
#                 if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L1, L4)
#                     
#                 # NE
#                 if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE):
#                     
#                     draw_line(self.outImage, L1, L2)
#                 # SE
#                 if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L2, L3)
#                 # SW
#                 if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L3, L4)
#                 # vertical
#                 if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L1, L3)
#                 # horizontal
#                 if (l1==1 and l2==0 and l3==1 and l4==0) and (abs(p1r.x-p2r.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L4, L2)
#                 
#                 cMid12r = find_mid_point(p1r,p2r)
#                 cMid14r = find_mid_point(p1r,p4r)
#                 cMid24r = find_mid_point(p2r,p4r)
#                 cMid23r = find_mid_point(p2r,p3r)
#                 cMid34r = find_mid_point(p3r,p4r)
#                 
#                 # case 1: interface crossing through L1 and L4
#                 if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) :
#                 #print "case 1"
#                     vecCoord1 = case_NW_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L1,L4,self);
#                     self.ishomog = 0
#                     if ( vecCoord1[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) ):
# 
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);               
#                         return True
#                     elif vecCoord1[0] != -1:
#                         self.enrichNodes = vecCoord1#[L1, L4]   
#                 
#                 # case 2: interface crossing through L1 and L2
#                 if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE):
#                     vecCoord2 = case_NE_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L1,L2);
#                     self.ishomog = 0
#                     if(vecCoord2[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) ):
#                     #print "case 2"
# 
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);
#                         return True
#                     elif vecCoord2[0] != -1:
#                         self.enrichNodes = vecCoord2
# #                        self.mat = 'Fluid'
#                     
#                         
#                 # case 3: interface crossing through L2 and L3
#                 if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) :
#                 #print "case 3"
#                     vecCoord3 = case_SE_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L2,L3);
#                     if(vecCoord3[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE)):
# 
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);
#                         return True
#                     elif vecCoord3[0].x != -1:
#                         self.enrichNodes = vecCoord3
# #                        self.mat = 'Fluid'
#                     self.ishomog = 0   
#                         
#                 # case 4: interface crossing through L4 and L3
#                 if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) :
#                 #print "case 4"
#                     vecCoord4 = case_SW_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L3,L4);
#                     if(vecCoord4[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE)):
#  
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);
#                         return True
#                     elif vecCoord4[0].x != -1:
#                         self.enrichNodes = vecCoord4
# #                        self.mat = 'Fluid'
#                     self.ishomog = 0
# 
#                 # case 5: interface crossing through L1 and L3
#                 if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) :
#                     #print "case 5"
#                     vecCoord5 = case_vertical_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L1,L3);
#                     if(vecCoord5[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE)):
# 
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);
#                         return True
#                     elif vecCoord5[0].x != -1:
#                         self.enrichNodes = vecCoord5
# #                        self.mat = 'Fluid'
#                     self.ishomog = 0   
#                     
#                 # case 6: interface crossing through L4 and L2
#                 if (l1==True and l2==False and l3==True and l4==False) and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE) :
#                 #print "case 6"
#                     vecCoord6 = case_horizontal_polynomial_test(self.inImage,self.outImage,p1r,p2r,p3r,p4r,L2,L4);
#                     if(vecCoord6[0].x == -1 and (abs(p1r.x-p2r.x) >= 2*MIN_SIZE)):
# 
#                         draw_line(self.outImage,cMid12r,cMid34r);
#                         draw_line(self.outImage,cMid14r,cMid23r);
#                         return True
#                     elif vecCoord6[0].x != -1:
#                         self.enrichNodes = vecCoord6
#                     self.ishomog = 0    
# #                        self.mat = 'Fluid'
#                                             
#                  # case 7: one line crossing through L1 and L4 and one line crossing through L2 and L3
#                  # 2-2 non adjacent corners are the same color (diagonally opposed)
#                  # the case of 3 consecutive-adjacent materials
#                 if (l1==0 and l2==0 and l3==0 and l4==0) :
# #                    print "case 7"
#                     self.ishomog = 0
#                     draw_line(self.outImage,cMid12r,cMid34r);
#                     draw_line(self.outImage,cMid14r,cMid23r);
#                     return True
#                 
#                 if (l1==1 and l2==1 and l3==0 and l4==1) or (l1==1 and l2==1 and l3==1 and l4==0) or (l1==0 and l2==1 and l3==1 and l4==1) or (l1==1 and l2==0 and l3==1 and l4==1) :
#                     self.ishomog = 0
#                     draw_line(self.outImage,cMid12r,cMid34r);
#                     draw_line(self.outImage,cMid14r,cMid23r);
#                     return True
                
                
        if ( self.children[0] != None and
             self.children[1] != None and
             self.children[2] != None and
             self.children[3] != None ):
                self.has_children = True
        else:
                self.has_children = False  
                      
    def printRect(self):
        print [self.rect[0].x, self.rect[1].x, self.rect[0].y, self.rect[3].y]
        
           
    def getinstance(self, rect, inImage, outImage):
        return Node(self, rect, inImage, outImage)
    
    def division_criterion(self, rect, inImage, outImage):
        return False
    
    def division_criterionOnce(self, rect, inImage, outImage):
        return False
         
class QuadTree(Node):
    maxdepth = 1
    leaves = []
    allnodes = []
    
    def __init__(self,rootnode):
        rootnode.subdivide() # constructs the network or nodes

    def count_nodes(self,root):
        
        allnodes = 0

        if root.has_children == False:
            return 0 
        else:
            if root.children[0] != None:
                allnodes += self.count_nodes(root.children[0]) +1
            if root.children[1] != None:
                allnodes += self.count_nodes(root.children[1]) +1
            if root.children[2] != None:
                allnodes += self.count_nodes(root.children[2]) +1
            if root.children[3] != None:
                allnodes += self.count_nodes(root.children[3]) +1

        return allnodes
    
    def leveling(self,root):
        if root.parent == 0:
            return 0
        if root.children != [None,None,None,None]:
            h1 = self.leveling(root.children[0])
            h2 = self.leveling(root.children[1])
            

    def neighbors_of_SE(self,root,masterNode):
        p1,p2,p3,p4 = root.rect
        
        # looking West:
        west_sibling = root.parent.children[2]
        if west_sibling.has_children:
             if west_sibling.children[1].has_children or west_sibling.children[3].has_children:
#                 print ' neighbors of SE: west'
                 root.divideOnce()
        
        
        # looking North:
        north_sibling = root.parent.children[1]
        if north_sibling.has_children:
            if north_sibling.children[2].has_children or north_sibling.children[3].has_children:
#                print ' neighbors of SE: north'
                root.divideOnce()       
                
        # looking East:
        # look up index of East neighbor
        if p2.x < root.imsize[0]-1: # check that we are not at the edge of the image
            east_neigh_index = str(find_neighbor_of(root.index,'R'))
            # checking to see if the east neighbor exists or is a ghost
            if it_exists(east_neigh_index,masterNode):
                east_neighbor = get_node_of_neighbor(root, root.index, east_neigh_index)
                if east_neighbor.has_children:
                    if east_neighbor.children[0].has_children or east_neighbor.children[2].has_children:
#                        print ' neighbors of SE: east'
                        root.divideOnce()

        # looking South:
        # look up index of South neighbor
        if p4.y < root.imsize[0]-1:
            south_neigh_index = str(find_neighbor_of(root.index,'D'))
            
            
            # checking to see if the south neighbor exists or is a ghost
            if it_exists(south_neigh_index,masterNode):
                south_neighbor = get_node_of_neighbor(root, root.index, south_neigh_index)
                if south_neighbor.has_children:
                    if south_neighbor.children[0].has_children or south_neighbor.children[1].has_children:
#                        print ' neighbors of SE: south'
                        root.divideOnce() 
  
        
    def neighbors_of_SW(self,root,masterNode):
        p1,p2,p3,p4 = root.rect
        
        # looking North:
        north_sibling = root.parent.children[0]
        if north_sibling.has_children:
            if north_sibling.children[2].has_children or north_sibling.children[3].has_children:
#                print ' neighbors of SW: north'
                root.divideOnce()       

        
        # looking East:
        east_sibling = root.parent.children[3]
        if east_sibling.has_children:
            if east_sibling.children[0].has_children or east_sibling.children[2].has_children:
#                print ' neighbors of SW: east'
                root.divideOnce()
         
        
        
        # looking West
        # look up index of West neighbor
        if p1.x > 0:
            west_neigh_index = str(find_neighbor_of(root.index,'L'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(west_neigh_index,masterNode):
                west_neighbor = get_node_of_neighbor(root, root.index, west_neigh_index)
                if west_neighbor.has_children:
                    if west_neighbor.children[1].has_children or west_neighbor.children[3].has_children:
#                        print ' neighbors of SW: west'
                        root.divideOnce()
        # looking South
        # look up index of South neighbor
        if p4.y < root.imsize[0]-1:
            south_neigh_index = str(find_neighbor_of(root.index,'D'))
            # checking to see if the south neighbor exists or is a ghost
            if it_exists(south_neigh_index,masterNode):
                south_neighbor = get_node_of_neighbor(root, root.index, south_neigh_index)
                if south_neighbor.has_children:
                    if south_neighbor.children[0].has_children or south_neighbor.children[1].has_children:
#                        print ' neighbors of SW: south'
                        root.divideOnce() 
  
    
    def neighbors_of_NE(self,root,masterNode):

        p1,p2,p3,p4 = root.rect
        
        # looking West
        west_sibling = root.parent.children[0]
        if west_sibling.has_children:
             if west_sibling.children[1].has_children or west_sibling.children[3].has_children:
#                 print ' neighbors of NE: west'
                 root.divideOnce()
#                 root.parent.children[1].divideOnce()
        
        # looking South
        south_sibling = root.parent.children[3]
        if south_sibling.has_children:
            if south_sibling.children[0].has_children or south_sibling.children[1].has_children:
#                print ' neighbors of NE: south'
                root.divideOnce()
          
              
        # looking East
        # look for East neighbor - on the right: R
        if p2.x < root.imsize[0]-1: # check that we are not at the edge of the image
            east_neigh_index = str(find_neighbor_of(root.index,'R'))
            # checking to see if the east neighbor exists or is a ghost
            if it_exists(east_neigh_index,masterNode):
                east_neighbor = get_node_of_neighbor(root, root.index, east_neigh_index)
                if east_neighbor.has_children:
                    if east_neighbor.children[0].has_children or east_neighbor.children[2].has_children:
#                        print ' neighbors of NE: east'
                        root.divideOnce()

        # looking North
        # look for North neighbor
        if p1.y > 0:
            north_neigh_index = str(find_neighbor_of(root.index,'U'))
            # checking to see if the north neighbor exists or is a ghost
            if it_exists(north_neigh_index,masterNode):
                north_neighbor = get_node_of_neighbor(root, root.index, north_neigh_index)
                if north_neighbor.has_children:
                    if north_neighbor.children[2].has_children or north_neighbor.children[3].has_children:
#                        print ' neighbors of NE: north'
                        root.divideOnce()

            
    def neighbors_of_NW(self,root,masterNode):
        p1,p2,p3,p4 = root.rect

        # looking East
        east_sibling = root.parent.children[1]
        if east_sibling.has_children:
            if east_sibling.children[0].has_children or east_sibling.children[2].has_children:
#                print ' neighbors of NW: east'
                root.divideOnce()

        # looking South
        south_sibling = root.parent.children[2]
        if south_sibling.has_children:
            if south_sibling.children[0].has_children or south_sibling.children[1].has_children:
#                print ' neighbors of NW: south'
                root.divideOnce()

        
        # looking West
        # look up West neighbor
        if p1.x > 0:
            west_neigh_index = str(find_neighbor_of(root.index,'L'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(west_neigh_index,masterNode):
                west_neighbor = get_node_of_neighbor(root, root.index, west_neigh_index)
                if west_neighbor.has_children:
                    if west_neighbor.children[1].has_children or west_neighbor.children[3].has_children:
#                        print ' neighbors of NW: west'
                        root.divideOnce()
            
                
        # looking North
        # look up North neighbor 
        if p1.y > 0:
            north_neigh_index = str(find_neighbor_of(root.index,'U'))
            # checking to see if the north neighbor exists or is a ghost
            if it_exists(north_neigh_index,masterNode):
                north_neighbor = get_node_of_neighbor(root, root.index, north_neigh_index)
                if north_neighbor.has_children:
                    if north_neighbor.children[2].has_children or north_neighbor.children[3].has_children:
#                        print ' neighbors of NW: north'
                        root.divideOnce()
        

    def balance_tree(self,root,masterNode):
   
        p1,p2,p3,p4 = root.rect
            
#        if root.i ==0 and root.j == 0:
#            root.printRect()
#            print root.index, root.depth
            
        if root.has_children == True:
            if root.children[0].has_children == False:
                self.neighbors_of_NW(root.children[0],masterNode)
            if root.children[1].has_children == False:
                self.neighbors_of_NE(root.children[1],masterNode)
            if root.children[2].has_children == False:
                self.neighbors_of_SW(root.children[2],masterNode)  
            if root.children[3].has_children == False:
                self.neighbors_of_SE(root.children[3],masterNode)              
        
        if root.children[0] != None:
            self.balance_tree(root.children[0],masterNode)
        if root.children[1] != None:
            self.balance_tree(root.children[1],masterNode)
        if root.children[2] != None:
            self.balance_tree(root.children[2],masterNode)
        if root.children[3] != None:
            self.balance_tree(root.children[3],masterNode)
#        print 'maximum depth: ', CNode.MAX_DEPTH   

#def height(root):
#    if root.parent == None:
#        return 0
#    h1 = height(root.children[0])
#    h2 = height(root.children[1])
#        
#    if (h1 == -1 or h2 == -1 or h1 < h2-1 or h1 > h2+1):
#        return -1
#    if h1>h2:
#        return h1 + 1
#    else:
#        return h2 + 1
        
class CNode(Node):
    
    def getinstance(self,rect,inImage,outImage,imageSize):
        return CNode(self,rect,inImage,outImage,imageSize)

    def division_criterionOnce(self, rect, inImage, outImage):
        p1,p2,p3,p4 = self.rect
        
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid24 = find_mid_point(p2,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)

        if abs(p1.x-p2.x) >= ALT_MIN_SIZE:
            draw_line(self.outImage,cMid12,cMid34)
            draw_line(self.outImage,cMid14,cMid23)
            
            return True
    
        return False
    
    def division_criterion(self, rect, inImage, outImage):

        p1,p2,p3,p4 = self.rect
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid24 = find_mid_point(p2,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        

        if abs(p1.x - p2.x) >= MAX_SIZE:               
            draw_line(self.outImage,cMid12,cMid34);
            draw_line(self.outImage,cMid14,cMid23);

            return True
        else:
            pxVal1 = self.inImage.GetPixel(p1.x,p1.y,0)
            pxVal2 = self.inImage.GetPixel(p2.x,p2.y,0)
            pxVal3 = self.inImage.GetPixel(p3.x,p3.y,0)
            pxVal4 = self.inImage.GetPixel(p4.x,p4.y,0)
            

            # are the 4 corners of the element in the same bin? i.e. homogeneous?
            isHomogeneous = four_corners_test(pxVal1,pxVal2,pxVal3,pxVal4);

            # if the four corners test fails
            if ( isHomogeneous == 0) and has_inclusions(self.inImage,p1,p2,p3,p4):
                
                
                l1 = ends_in_same_bin(self.inImage,p1,p2);
                l2 = ends_in_same_bin(self.inImage,p2,p3);
                l3 = ends_in_same_bin(self.inImage,p4,p3);
                l4 = ends_in_same_bin(self.inImage,p1,p4);
            
#                L1 = linear_search(self.inImage,p1,p2);
#                L2 = linear_search(self.inImage,p2,p3);
#                L3 = linear_search(self.inImage,p4,p3);
#                L4 = linear_search(self.inImage,p1,p4);
                L1 = search_in(LIST,p1,p2,self.inImage)
                L2 = search_in(LIST,p2,p3,self.inImage)
                L3 = search_in(LIST,p4,p3,self.inImage)
                L4 = search_in(LIST,p1,p4,self.inImage)                
                
                if ( 
                    len(L2) > 1 or len(L4) > 1 or len(L1) > 1 or len(L3) > 1 
                    # or len(L2) < 1 or len(L4) < 1 or len(L1) < 1 or len(L3) < 1
                     ):
                    # interface croses one edge multiple times
                    draw_line(self.outImage,cMid12,cMid34);
                    draw_line(self.outImage,cMid14,cMid23)
                    return True
                else:
                    if len(L1) == 1:
                        L1 = L1[0]
                    if len(L2) == 1:
                        L2 = L2[0]
                    if len(L3) == 1:
                        L3 = L3[0]
                    if len(L4) == 1:
                        L4 = L4[0]
                        
#                 print L1
#                 print L1.x, L1.y
#                 print L2.x, L2.y 
#                 print L3.x, L3.y
#                 print L4.x, L4.y

                # NW
                if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    self.ishomog = 0
#                     draw_line(self.outImage, L1, L4)
                    
                # NE
                if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE):
#                     draw_line(self.outImage, L1, L2)
                    self.ishomog = 0
                # SE
                if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L2, L3)
                    self.ishomog = 0
                # SW
                if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L3, L4)
                    self.ishomog = 0
                # vertical
                if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L1, L3)
                    self.ishomog = 0
                # horizontal
                if (l1==1 and l2==0 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                     draw_line(self.outImage, L4, L2)
                    self.ishomog = 0
                
                # case 1: interface crossing through L1 and L4
                if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                #print "case 1"
                    vecCoord1 = case_NW_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L1,L4,self);
                    self.ishomog = 0
                    if ( vecCoord1[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE) ):

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);               
                        return True
                    elif vecCoord1[0] != -1:
                        self.enrichNodes = vecCoord1#[L1, L4]
                       
#                        self.mat = 'Fluid'
#                        print self.en[0].x, self.en[0].y, self.en[1].x, self.en[1].y

                # case 2: interface crossing through L1 and L2
                if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1.x-p2.x) >= 2*MIN_SIZE):
                    vecCoord2 = case_NE_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L1,L2);
                    self.ishomog = 0
                    if(vecCoord2[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE) ):
                    #print "case 2"

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);
                        return True
                    elif vecCoord2[0] != -1:
                        self.enrichNodes = vecCoord2
#                        self.mat = 'Fluid'
                    
                        
                # case 3: interface crossing through L2 and L3
                if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                #print "case 3"
                    vecCoord3 = case_SE_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L2,L3);
                    if(vecCoord3[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE)):

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);
                        return True
                    elif vecCoord3[0].x != -1:
                        self.enrichNodes = vecCoord3
#                        self.mat = 'Fluid'
                    self.ishomog = 0   
                        
                # case 4: interface crossing through L4 and L3
                if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                #print "case 4"
                    vecCoord4 = case_SW_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L3,L4);
                    if(vecCoord4[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE)):
 
                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);
                        return True
                    elif vecCoord4[0].x != -1:
                        self.enrichNodes = vecCoord4
#                        self.mat = 'Fluid'
                    self.ishomog = 0

                # case 5: interface crossing through L1 and L3
                if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                    #print "case 5"
                    vecCoord5 = case_vertical_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L1,L3);
                    if(vecCoord5[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE)):

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);
                        return True
                    elif vecCoord5[0].x != -1:
                        self.enrichNodes = vecCoord5
#                        self.mat = 'Fluid'
                    self.ishomog = 0   
                # case 6: interface crossing through L4 and L2
                if (l1==True and l2==False and l3==True and l4==False) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                #print "case 6"
                    vecCoord6 = case_horizontal_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L2,L4);
                    if(vecCoord6[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE)):

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);
                        return True
                    elif vecCoord6[0].x != -1:
                        self.enrichNodes = vecCoord6
                    self.ishomog = 0    
#                        self.mat = 'Fluid'
                                            
                 # case 7: one line crossing through L1 and L4 and one line crossing through L2 and L3
                 # 2-2 non adjacent corners are the same color (diagonally opposed)
                 # the case of 3 consecutive-adjacent materials
                if (l1==0 and l2==0 and l3==0 and l4==0) :
#                    print "case 7"
                    self.ishomog = 0
                    draw_line(self.outImage,cMid12,cMid34);
                    draw_line(self.outImage,cMid14,cMid23);
                    return True
                
                if (l1==1 and l2==1 and l3==0 and l4==1) or (l1==1 and l2==1 and l3==1 and l4==0) or (l1==0 and l2==1 and l3==1 and l4==1) or (l1==1 and l2==0 and l3==1 and l4==1) :
                    self.ishomog = 0
                    draw_line(self.outImage,cMid12,cMid34);
                    draw_line(self.outImage,cMid14,cMid23);
                    return True


        return False
 
class CQuadTree(QuadTree):
    def __init__(self,rootnode):
        QuadTree.__init__(self, rootnode)
 
def not_a_corner(root,A):
    p1,p2,p3,p4 = root.rect

    if ( (A.x != p1.x and A.x != p2.x) and
         (A.y != p1.y and A.y != p4.y) ):
        return True
    return False
     
def set_interval(imSize,level):
    my_arr = [0,imSize-1]
    
    for i in range(0,level):
        mdpt = numpy.zeros((len(my_arr)-1,1))
        new_arr = numpy.zeros(( len(my_arr) + len(mdpt), 1  ))
        for j in range(0,len(my_arr)-1):
            mdpt[j] = (my_arr[j] + my_arr[j+1] ) // 2

        for k in range(0,len(new_arr)):
            
            if  k % 2 == 0:
                new_arr[k] = my_arr[k/2]
            else:
            
                new_arr[k] = mdpt[k/2]
        my_arr = new_arr
        
    return my_arr 

def in_child_k(rects,L):
    p1,p2,p3,p4 = rects
    if p1.x <= L.x and L.x <= p2.x and p1.y <= L.y and L.y <= p4.y:
        return True
    
    return False
  
def it_exists(index,masterNode):      
    llen = len(index)
    child = masterNode
    for i in range(0,llen):
        if child.children[int(index[i])].has_children == False:
            return False
        child = child.children[int(index[i])]
    return True
 
def node_exists(index,list):
    # look up index in list, if it exists return true, else false 
    for i in range(0,len(list)):
        if llist[i] == [str(index)]:
            return True
        
    return False
         
        
def tree_balance(tree, root,masterNode):
   
        p1,p2,p3,p4 = root.rect
            
#        if root.i ==0 and root.j == 0:
#            root.printRect()
#            print root.index, root.depth
            
        if root.has_children == True:
            
            if root.children[0].has_children == False:
                tree.neighbors_of_NW(root.children[0],masterNode)
            if root.children[1].has_children == False:
                tree.neighbors_of_NE(root.children[1],masterNode)
            if root.children[2].has_children == False:
                tree.neighbors_of_SW(root.children[2],masterNode)  
            if root.children[3].has_children == False:
                tree.neighbors_of_SE(root.children[3],masterNode)              
#
        if root.children[0] != None:
            tree_balance(tree,root.children[0],masterNode)
        if root.children[1] != None:
            tree_balance(tree,root.children[1],masterNode)
        if root.children[2] != None:
            tree_balance(tree,root.children[2],masterNode)
        if root.children[3] != None:
            tree_balance(tree,root.children[3],masterNode)


def ghost_nodes_enrichment_nodes(tree, root, masterNode):

        p1,p2,p3,p4 = root.rect
            
        west_has_children = False
        east_has_children = False
        south_has_children = False
        north_has_children = False
        
        # if root has no children look at his neighbors
        # if all of them do have children, 
        # root needs to be subdivided
        if root.has_children == False:
            
            west_neigh_index = str(find_neighbor_of(root.index,'L'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(west_neigh_index, masterNode):
                west_neighbor = get_node_of_neighbor(root, root.index, west_neigh_index)
                if west_neighbor.has_children == True:
                    west_has_children = True
                    
            east_neigh_index = str(find_neighbor_of(root.index,'R'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(east_neigh_index, masterNode):
                east_neighbor = get_node_of_neighbor(root, root.index, east_neigh_index)
                if east_neighbor.has_children == True:
                    east_has_children = True

            south_neigh_index = str(find_neighbor_of(root.index,'D'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(south_neigh_index, masterNode):
                south_neighbor = get_node_of_neighbor(root, root.index, south_neigh_index)
                if south_neighbor.has_children == True:
                    south_has_children = True

            north_neigh_index = str(find_neighbor_of(root.index,'U'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(north_neigh_index, masterNode):
                north_neighbor = get_node_of_neighbor(root, root.index, north_neigh_index)
                if north_neighbor.has_children == True:
                    north_has_children = True
                    
            if (len(root.enrichNodes) > 0 and 
                (west_has_children == True or east_has_children == True or
                 south_has_children == True or north_has_children == True)):
#                print root.index
                root.divideOnce()      

        if root.children[0] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[0],masterNode)
        if root.children[1] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[1],masterNode)
        if root.children[2] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[2],masterNode)
        if root.children[3] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[3],masterNode)

def number_of_generations(tree, root, masterNode):

#        if root.has_children == False:
#            nrchildren = 0
#        else:
#            nrchildren = nrchildren + 1
    maxx = 0
    if root.has_children == True:
        
        if root.children[0] != None:
            val = 1 + number_of_generations(tree,root.children[0],masterNode)
            if maxx < val:
                maxx = val 
        if root.children[1] != None:
            val = 1 + number_of_generations(tree,root.children[1],masterNode)
            if maxx < val:
                maxx = val
        if root.children[2] != None:
            val = 1 + number_of_generations(tree,root.children[2],masterNode)
            if maxx < val:
                maxx = val
        if root.children[3] != None:
            val = 1 + number_of_generations(tree,root.children[3],masterNode)
            if maxx < val:
                maxx = val
    return maxx


def get_list_of_nodes(tree, root, masterNode,llist):

#        if root.has_children == False:
#            nrchildren = 0
#        else:
#            nrchildren = nrchildren + 1
        if root.has_children == False:
            llist.append([root.index]) 
            
        if root.children[0] != None:
            get_list_of_nodes(tree,root.children[0],masterNode,llist)
        if root.children[1] != None:
            get_list_of_nodes(tree,root.children[1],masterNode,llist)
        if root.children[2] != None:
            get_list_of_nodes(tree,root.children[2],masterNode,llist)
        if root.children[3] != None:
            get_list_of_nodes(tree,root.children[3],masterNode,llist)

        return llist
    
def three_neighbor_rule(tree, root, masterNode):

        p1,p2,p3,p4 = root.rect
            
        west_has_children = False
        east_has_children = False
        south_has_children = False
        north_has_children = False
        
        counter = 0
        # if root has no children look at his neighbors
        # if all of them do have children, 
        # root needs to be subdivided
        if root.has_children == False:
            
            west_neigh_index = str(find_neighbor_of(root.index,'L'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(west_neigh_index, masterNode):
                west_neighbor = get_node_of_neighbor(root, root.index, west_neigh_index)
                if west_neighbor.has_children == True:
                    west_has_children = True
                    counter += 1
                    
            east_neigh_index = str(find_neighbor_of(root.index,'R'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(east_neigh_index, masterNode):
                east_neighbor = get_node_of_neighbor(root, root.index, east_neigh_index)
                if east_neighbor.has_children == True:
                    east_has_children = True
                    counter += 1

            south_neigh_index = str(find_neighbor_of(root.index,'D'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(south_neigh_index, masterNode):
                south_neighbor = get_node_of_neighbor(root, root.index, south_neigh_index)
                if south_neighbor.has_children == True:
                    south_has_children = True
                    counter += 1

            north_neigh_index = str(find_neighbor_of(root.index,'U'))    
            # checking to see if the west neighbor exists or is a ghost
            if it_exists(north_neigh_index, masterNode):
                north_neighbor = get_node_of_neighbor(root, root.index, north_neigh_index)
                if north_neighbor.has_children == True:
                    north_has_children = True
                    counter += 1
            
            #if west_has_children == True and east_has_children == True and south_has_children == True and north_has_children == True:
            if counter > 2:
#                print '3 neigh-rule - ', root.index
                root.divideOnce()      

                        
        if root.children[0] != None:
            three_neighbor_rule(tree,root.children[0],masterNode)
        if root.children[1] != None:
            three_neighbor_rule(tree,root.children[1],masterNode)
        if root.children[2] != None:
            three_neighbor_rule(tree,root.children[2],masterNode)
        if root.children[3] != None:
            three_neighbor_rule(tree,root.children[3],masterNode)

#def stress_concentration_constraint(tree, root, masterNode):
#
#    print 'graph'
##         if root.children[0] != None:
##             stress_concentration_constraint(tree,root.children[0],masterNode)
##         if root.children[1] != None:
##             stress_concentration_constraint(tree,root.children[1],masterNode)
##         if root.children[2] != None:
##             stress_concentration_constraint(tree,root.children[2],masterNode)
##         if root.children[3] != None:
##             stress_concentration_constraint(tree,root.children[3],masterNode)
#  
def set_graph(masterNode,llist):
    # each node called ROOT has 8 neighbors: 
    #  NW   NC   NE 
    #  WC  ROOT  EC
    #  SW   SC   SE
    
    mesh_graph = {}
    simple_graph = {}

    n = len(llist)
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])
        
        p1,p2,p3,p4 = root.rect
        
        # NORTH-WEST neighbor index:
#         if p1.y > 0 and p1.x > 0:
#             north_west_neigh_index = str(find_neighbor_of(root.index,'LU'))
#             # checking to see if the north-west neighbor exists or is a ghost
#             if node_exists(north_west_neigh_index,llist) == False: # node does not exist
#                 
#                 # check if one level above, the node exists:
#                 if node_exists(str(north_west_neigh_index[:-1]),llist) == True:
#                     north_west_neighbor = get_node_by_id(masterNode, [str(north_west_neigh_index)])
#                     parent_north_west_neighbor = get_node_by_id(masterNode, [str(north_west_neigh_index[:-1])])
#                     if parent_north_west_neighbor.has_children:
#                         if north_west_neighbor.has_children:
#                             north_west_neigh_index = [str(north_west_neighbor.children[3].index)]
#                     else:
#                         north_west_neigh_index = [str(parent_north_west_neighbor.index)]
# #                 else:
# #                     # check if two levels above the node exists:
# #                     if node_exists(str(north_west_neigh_index[:-2]),llist) == True:
#             else: # the node exists: done
#                 north_west_neigh_index = [str(north_west_neigh_index)]
#         else:
#             north_west_neigh_index = ['-1']

        # NORTH-WEST neighbor index:
        if p1.y > 0 and p1.x > 0:
            north_west_neigh_index = str(find_neighbor_of(root.index,'LU'))
            # checking to see if the north-west neighbor exists or is a ghost
            if node_exists(north_west_neigh_index,llist) == False: # node does not exist
                
                # check if one level above, the node exists:
                if node_exists(str(north_west_neigh_index[:-1]),llist) == False:
                    #check if two levels above, the node exists:
                    if node_exists(str(north_west_neigh_index[:-2]),llist) == False:
                        print 'ERROR: this case should NOT exist'
                    else:
                        north_west_neigh_index = [str(north_west_neigh_index[:-2])]
                else:
                    north_west_neigh_index = [str(north_west_neigh_index[:-1])]
            else: # the node exists: done
                north_west_neigh_index = [str(north_west_neigh_index)]
        else:
            north_west_neigh_index = ['-1']

#         # SOUTH-EAST neighbor index:
#         if p3.y < 999 and p3.x < 999:
#             south_east_neigh_index = str(find_neighbor_of(root.index,'RD'))
#             south_east_neighbor = get_node_by_id(masterNode, [str(south_east_neigh_index)])
# 
#             # checking to see if the south-east neighbor exists or is a ghost
#             if node_exists(str(south_east_neigh_index),llist) == False and south_east_neighbor.has_children == False:
#                 #check if one level above, the node exisits:
#                 if node_exists(str(south_east_neigh_index[:-1]),llist) == False:
#                     #check if two levels above the node existS:
#                     if node_exists(str(south_east_neigh_index[:-2]),llist) == False:
#                         print 'ERROR: this case should NOT exist'
#                         print str(south_east_neigh_index)
#                         print root.index
#                     else:
#                         south_east_neigh_index = [str(south_east_neigh_index[:-2])]
#                 else:
#                     south_east_neigh_index = [str(south_east_neigh_index[:-1])]
#             else:
#                 south_east_neigh_index = [str(south_east_neigh_index)]
#         else:
#             south_east_neigh_index = ['-1']
# 
#         print south_east_neigh_index, [str(south_east_neigh_index)]
        
        # SOUTH-WEST neighbor index:
        if p3.y < 999 and p1.x > 0:
            south_west_neigh_index = str(find_neighbor_of(root.index,'LD'))
            # checking to see if the south-west neighbor exists or is a ghost
            if node_exists(south_west_neigh_index,llist) == False:
                #check if one level above the node exists:
                if node_exists(str(south_west_neigh_index[:-1]),llist) == False:
                    # check if two levels above the node exists:
                    if node_exists(str(south_west_neigh_index[:-2]),llist) == False:
                        print 'ERROR: this case should NOT exist'
                    else:
                        south_west_neigh_index = [str(south_west_neigh_index[:-2])]
                else:
                    south_west_neigh_index = [str(south_west_neigh_index[:-1])]
            else:
                south_west_neigh_index = [str(south_west_neigh_index)]
        else:
            south_west_neigh_index = ['-1']

        # NORTH-EAST neighbor index:
        if p1.y > 0 and p2.x < 999:
            north_east_neigh_index = str(find_neighbor_of(root.index,'RU'))
            # checking to see if the north-east neighbor exists or is a ghost
            if node_exists(north_east_neigh_index,llist) == False:
                
                #check if one level above, the node exists:
                if node_exists(str(north_east_neigh_index[:-1]),llist) == False:
                    # check if two levels above, the node exists:
                    if node_exists(str(north_east_neigh_index[:-2]),llist) == False:
                        print 'ERROR: this case should NOT exist'
                    else:
                        north_east_neigh_index = [str(north_east_neigh_index[:-2])]
                else:
                    north_east_neigh_index = [str(north_east_neigh_index[:-1])]
            else:
                north_east_neigh_index = [str(north_east_neigh_index)]
        else:
            north_east_neigh_index = ['-1']    
                    
        # NORTH neighbor index:
        if p1.y > 0:
            north_neigh_index = str(find_neighbor_of(root.index,'U'))
            north_neighbor = get_node_by_id(masterNode, [str(north_neigh_index)])
            # checking to see if the north neighbor exists or is a ghost
            if node_exists(north_neigh_index,llist) == False:
                parent_north_neighbor = get_node_by_id(masterNode, [str(north_neigh_index[:-1])])
                if parent_north_neighbor.has_children:
                    if north_neighbor.has_children:
                        north_neigh_index = [str(north_neighbor.children[2].index), str(north_neighbor.children[3].index)]
                else:
                    north_neigh_index = [str(parent_north_neighbor.index)]
            else:
                north_neigh_index = [str(north_neigh_index)]
        else:
            north_neigh_index = ['-1']
            
#         # NORTH-EAST neighbor index:
#         if p1.y > 0 and p2.x < 999:
#             north_east_neigh_index = str(find_neighbor_of(root.index,'RU'))
#             north_east_neighbor = get_node_by_id(masterNode, [str(north_east_neigh_index)])
#             # checking to see if the north-east neighbor exists or is a ghost
#             if node_exists(north_east_neigh_index,llist) == False:
#                 parent_north_east_neighbor = get_node_by_id(masterNode, [str(north_east_neigh_index[:-1])])
#                 if parent_north_east_neighbor.has_children:
#                     if north_east_neighbor.has_children:
#                         north_east_neigh_index = [str(north_east_neighbor.children[2].index)]
#                 else:
#                     north_east_neigh_index = [str(parent_north_east_neighbor.index)]
#             else:
#                 north_east_neigh_index = [str(north_east_neigh_index)]
#         else:
#             north_east_neigh_index = ['-1']            


        # WEST neighbor index:
        if p1.x > 0:
            west_neigh_index = str(find_neighbor_of(root.index,'L'))
            west_neighbor = get_node_by_id(masterNode, [str(west_neigh_index)])
            # checking to see if the west neighbor exists or is a ghost
            if node_exists(west_neigh_index,llist) == False:
                parent_west_neighbor = get_node_by_id(masterNode, [str(west_neigh_index[:-1])])
                if parent_west_neighbor.has_children:
                    if west_neighbor.has_children:
                        west_neigh_index = [str(west_neighbor.children[0].index), str(west_neighbor.children[3].index)]
                else:
                   west_neigh_index = [str(parent_west_neighbor.index)]
            else:
                west_neigh_index = [str(west_neigh_index)]
        else:
            west_neigh_index = ['-1']
            
        # EAST neighbor index:
        if p2.x < 999:
            east_neigh_index = str(find_neighbor_of(root.index,'R'))
            east_neighbor = get_node_by_id(masterNode, [str(east_neigh_index)])
            # checking to see if the east neighbor exists or is a ghost
            if node_exists(east_neigh_index,llist) ==  False:
                parent_east_neighbor = get_node_by_id(masterNode, [str(east_neigh_index[:-1])])
                if parent_east_neighbor.has_children:
                    if east_neighbor.has_children:
                        east_neigh_index = [str(east_neighbor.children[0].index), str(east_neighbor.children[2].index)]
                else:
                    east_neigh_index = [str(parent_east_neighbor.index)]
            else:
                east_neigh_index = [str(east_neigh_index)]
        else:
            east_neigh_index = ['-1'] 

        
#         # SOUTH-WEST neighbor index:
#         if p3.y < 999 and p1.x > 0:
#             south_west_neigh_index = str(find_neighbor_of(root.index,'LD'))
#             south_west_neighbor = get_node_by_id(masterNode, [str(south_west_neigh_index)])
#             # checking to see if the south-west neighbor exists or is a ghost
#             if node_exists(south_west_neigh_index,llist) == False:
#                 parent_south_west_neighbor = get_node_by_id(masterNode, [str(south_west_neigh_index[:-1])])                
#                 if parent_south_west_neighbor.has_children: 
#                     if south_west_neighbor.has_children:
#                         south_west_neigh_index = [str(south_west_neighbor.children[1].index)]
#                 else:
#                     south_west_neigh_index = [str(parent_south_west_neighbor.index)]
#             else:
#                 south_west_neigh_index = [str(south_west_neigh_index)]
#         else:
#             south_west_neigh_index = ['-1']


        # SOUTH neighbor index:
        if p3.y < 999:
            south_neigh_index = str(find_neighbor_of(root.index,'D'))
            south_neighbor = get_node_by_id(masterNode,[str(south_neigh_index)])
            # checking to see if the south neighbor exists or is a ghost
            if node_exists(south_neigh_index,llist) == False:
                parent_south_neighbor = get_node_by_id(masterNode,[str(south_neigh_index[:-1])])
                if parent_south_neighbor.has_children:
                    if south_neighbor.has_children:
                        south_neigh_index = [str(south_neighbor.children[0].index), str(south_neighbor.children[1].index)]
                else:
                    south_neigh_index = [str(parent_south_neighbor.index)]
            else:
                south_neigh_index = [str(south_neigh_index)]
        else:
            south_neigh_index = ['-1']
            
 
        # SOUTH-EAST neighbor index:
        if p3.y < 999 and p3.x < 999:
             
            south_east_neigh_index = str(find_neighbor_of(root.index,'RD'))
             
            south_east_neighbor = get_node_by_id(masterNode, [str(south_east_neigh_index)])
             
            # checking to see if the south-east neighbor exists or is a ghost
            if node_exists(south_east_neigh_index,llist) == False:
                 
                parent_south_east_neighbor = get_node_by_id(masterNode, [str(south_east_neigh_index[:-1])])
             
                if parent_south_east_neighbor.has_children:
                     
                    if south_east_neighbor.has_children:
                        south_east_neigh_index = [str(south_east_neighbor.children[0].index)]
               
                else:
                    south_east_neigh_index = [str(parent_south_east_neighbor.index)]
            else:
                south_east_neigh_index = [str(south_east_neigh_index)]
        else:
            south_east_neigh_index = ['-1'] 
#         
        mesh_graph[str(root.index)] = {
                          'I': root.ishomog,
                          'NW': north_west_neigh_index,
                          'N': north_neigh_index,
                          'NE': north_east_neigh_index,
                          'W': west_neigh_index,
                          'E': east_neigh_index,
                          'SW': south_west_neigh_index,
                          'S': south_neigh_index,
                          'SE': south_east_neigh_index
                           }
#         graph_list = []
#         if north_west_neigh_index != ['-1']:
#             graph_list = graph_list + north_west_neigh_index
#         if north_neigh_index != ['-1']:
#             graph_list = graph_list + north_neigh_index
#         if north_east_neigh_index != ['-1']:
#             graph_list = graph_list +  north_east_neigh_index
#         if west_neigh_index != ['-1']:
#             graph_list = graph_list + west_neigh_index
#         if east_neigh_index != ['-1']:
#             graph_list = graph_list + east_neigh_index
#         if south_west_neigh_index != ['-1']:
#             graph_list = graph_list + south_west_neigh_index
#         if south_neigh_index != ['-1']:
#             graph_list = graph_list + south_neigh_index
#         if south_east_neigh_index != ['-1']:
#             graph_list = graph_list + south_east_neigh_index
            
             
#         simple_graph[str(root.index)] = { graph_list}
    return mesh_graph
 
def find_shortest_path(graph, start, end, path = []):
        start = '00'
        end = '11'
#        print [start, end]
        path = path + [start]
        if start == end:
            return path
        if not graph.has_key(start):
            print 'no key', start
            return None
        shortest = None
#         for node in graph[start]:
#             if node not in path:
#                 newpath = find_shortest_path(graph, node, end, path)
#                 if newpath:
#                     if not shortest or len(newpath) < len(shortest):
#                         shortest = newpath
        return shortest
           
def tomorton(x,y):
#  http://www.thejach.com/view/id/207
  x = bin(x)[2:]
  lx = len(x)
  y = bin(y)[2:]
  ly = len(y)
  L = max(lx, ly)
  m = 0
  for j in xrange(1, L+1):
    # note: ith bit of x requires x[lx - i] since our bin numbers are big endian
    xi = int(x[lx-j]) if j-1 < lx else 0
    yi = int(y[ly-j]) if j-1 < ly else 0
    m += 2**(2*j)*xi + 2**(2*j+1)*yi
  return m/4
  
def convert_to_base_4(n):
    # converting number n in base 10 to base 4

    result = []
    (n,remainder) = divmod(n,4)
    result.append(str(remainder))
    while n:
        (n,remainder) = divmod(n,4)
        result.append(str(remainder))
        
    result.reverse()
    
#    print 'reverse', result.reverse()
    
    return ''.join(result)

def morton_id(i,j):
    bi = bin(i)[2:]
    bj = bin(j)[2:]
    
    si = str(bi)
    sj = str(bj)
    sisj = [item for slist in izip_longest(si, sj) for item in slist if item is not None]
    sisj = ''.join(sisj)
    ind_base10 = int(sisj,2)

    
#    return convert_to_base_4(ind_base10)
    return tomorton(i,j)

def find_neighbor_of(index,direction):
    loc = str(index)
    llist_str = list(loc)
    for i in range(len(loc)-1,-1,-1):
        new_quadrant =  D[str(loc[i])][direction]['Quadrant']
        new_direction = D[str(loc[i])][direction]['Direction']
        if new_direction != 'H':
            direction = new_direction
            llist_str[i] = str(new_quadrant)
        else:
            llist_str[i] = str(new_quadrant)
            return str("".join(llist_str))

    return str("".join(llist_str))

def get_node_of_neighbor(root,my_ind,neigh_ind):
    
    p = root
    for i in range(0,len(my_ind)):
        p = p.parent
        
    r = p
    for j in range(0,len(neigh_ind)):
        r = r.children[int(neigh_ind[j])]
    return r

def get_node_by_id(node,id):
        # returns node, given index 
        # index could be ['00101'], thus it'a list
        
        index = id[0]
        ll = len(index)

        p = node
        for i in range(0,ll):
            p = p.children[int(index[i])]
            
        return p  
def remove_duplicates(seq): 
   # order preserving

   noDupes = []
   [noDupes.append(i) for i in seq if not noDupes.count(i)]
   return noDupes


def process_list_of_elements(llist,root):
    n = len(llist)
    pvec = numpy.zeros((0,2))
    
    coordsList1 = []
    #first, construct the pvec of just the physical coordinates of the element corners
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,llist[i])
        p1,p2,p3,p4 = root_i.rect
        coordsList1 = coordsList1 + [[p1.x, p1.y]]
        coordsList1 = coordsList1 + [[p2.x, p2.y]]
        coordsList1 = coordsList1 + [[p3.x, p3.y]]
        coordsList1 = coordsList1 + [[p4.x, p4.y]]

    cList1 = deepcopy(coordsList1)
    
    for i in range(0, len(coordsList1)):
        coordsList1[i][0] = coordsList1[i][0] / 1000.0
        coordsList1[i][1] = coordsList1[i][1] / 1000.0
        if coordsList1[i][0] != 0.0:
            coordsList1[i][0] += 0.001
        if coordsList1[i][1] != 0.0:
            coordsList1[i][1] += 0.001
             
        coordsList1[i][1] = 1 - coordsList1[i][1] 
        
        
    # sort by y, and then by x
    coordsList1 = sorted(coordsList1,key=lambda x: (x[1],x[0]))
    # remove duplicates from the list:
    coordsList1 = [ key for key,_ in groupby(coordsList1)]

    lenClist1 = len(coordsList1)
#    print 'inside process list',lenClist1
    
    coordsList2 = []
    #second, construct the pvec of the enrichment nodes
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,llist[i])
        l = len(root_i.enrichNodes)
        if l > 0:
            for j in range(0,l):
                enrN = root_i.enrichNodes[j]
                coordsList2 = coordsList2 + [[enrN.x, enrN.y]]

    cList2 = deepcopy(coordsList2)
    
    
    for i in range(0, len(coordsList2)):
        coordsList2[i][0] = coordsList2[i][0] / 1000.0
        coordsList2[i][1] = coordsList2[i][1] / 1000.0
        if coordsList2[i][0] != 0.0:
            coordsList2[i][0] += 0.001
        if coordsList2[i][1] != 0.0:
            coordsList2[i][1] += 0.001
             
        coordsList2[i][1] = 1 - coordsList2[i][1] 

    # sort by y, and then by x
    coordsList2 = sorted(coordsList2,key=lambda x: (x[1],x[0]))
    # remove duplicates from the list:
    coordsList2 = [ key for key,_ in groupby(coordsList2)]
    
    coordsList = coordsList1 + coordsList2
    coordsList = remove_duplicates(coordsList)
    pvec = numpy.vstack([pvec,coordsList])

 
    cList1 = sorted(cList1,key=lambda x: (x[1],x[0]))
    cList1 = [ key for key,_ in groupby(cList1)]
    cList2 = sorted(cList2,key=lambda x: (x[1],x[0]))
    # remove duplicates from the list
    cList2 = [ key for key,_ in groupby(cList2)]
#     cList2 = [[999, 499], [916, 624], [874, 687], [833, 749], [749, 874], [749, 876], [666, 999]]
     
    cList = cList1 + cList2
    cList = remove_duplicates(cList)
    pvecCList = numpy.zeros((0,2))
    pvecCList = numpy.vstack([pvecCList,cList])

#     pvec = pvec / 1000.0
#     for i in range(0,len(pvec)):
#         if pvec[i,0] != 0.0:
#             pvec[i,0] += 0.001
#         if pvec[i,1] != 0.0:
#             pvec[i,1] += 0.001
#              
#     for i in range(0,len(pvecv)):
#         pvec[i,1] = 1 - pvec[i,1] 

#     print pvecCList
    
    # construct the corners list t:
#     t = []
# #     for i in range(0,n):
# #         root_i = get_node_by_id(masterNode,llist[i])
# #         p1,p2,p3,p4 = root_i.rect
# #         
# # #         p1.x = p1.x/1000.0
# # #         p1.y = p1.y/1000.0
# # #         if p1.x != 0.0:
# # #             p1.x += 0.001
# # #         if p1.y != 0.0:
# # #             p1.y += 0.001
# # #         p1.y = 1 - p1.y
# #         
# # #         p2.x = p2.x/1000.0
# # #         p2.y = p2.y/1000.0
# # #         if p2.x != 0.0:
# # #             p2.x += 0.001
# # #         if p2.y != 0.0:
# # #             p2.y += 0.001
# # #         p2.y = 1 - p2.y
# #         
# #         
# # #         p3.x = p3.x/1000.0
# # #         p3.y = p3.y/1000.0
# # #         if p3.x != 0.0:
# # #             p3.x += 0.001
# # #         if p3.y != 0.0:
# # #             p3.y += 0.001
# # #         p3.y = 1 - p3.y
# #         
# #         
# #         
# # #         p4.x = p4.x/1000.0
# # #         p4.y = p4.y/1000.0
# # #         if p4.x != 0.0:
# # #             p4.x += 0.001
# # #         if p4.y != 0.0:
# # #             p4.y += 0.001
# # #         p4.y = 1 - p4.y
# #         
# #         
# #         b1 = [p1.x, p1.y]
# #         ind1 = numpy.where(numpy.all(pvecCList==b1,axis=1))
# #         c1 = ind1[0][0]
# #          
# #         b2 = [p2.x, p2.y]
# #         ind2 = numpy.where(numpy.all(pvecCList==b2,axis=1))
# #         c2 = ind2[0][0]
# #         
# #         b3 = [p3.x, p3.y]
# #         ind3 = numpy.where(numpy.all(pvecCList==b3,axis=1))
# #         c3 = ind3[0][0]
# #         
# #         b4 = [p4.x, p4.y]
# #         ind4 = numpy.where(numpy.all(pvecCList==b4,axis=1))
# #         c4 = ind4[0][0]
# #         
# #         l = len(root_i.enrichNodes)
# #         if l > 0:
# #             if l == 1:
# #                 enrN1 = root_i.enrichNodes[0]
# #                 
# # #                 enrN1.x = enrN1.x/1000.0
# # #                 enrN1.y = enrN1.y/1000.0
# # #                 if enrN1.x != 0.0:
# # #                     enrN1.x += 0.001
# # #                 if enrN1.y != 0.0:
# # #                     enrN1.y += 0.001
# # #                 enrN1.y = 1 - enrN1.y
# #                 
# #                 b5 = [enrN1.x, enrN1.y]
# #                 ind5 = numpy.where(numpy.all(pvecCLisr==b5,axis=1))
# #                 c5 = ind5[0][0]
# # #                 t = t + [[c1,c2,c3,c4,c5]]
# #                 t = t + [[c4,c3,c2,c1, c5]]
# #                 
# #             if l == 2:
# #                 enrN1 = root_i.enrichNodes[0]
# #                 
# # #                 enrN1.x = enrN1.x/1000.0
# # #                 enrN1.y = enrN1.y/1000.0
# # #                 if enrN1.x != 0.0:
# # #                     enrN1.x += 0.001
# # #                 if enrN1.y != 0.0:
# # #                     enrN1.y += 0.001
# # #                 enrN1.y = 1 - enrN1.y
# #                 
# #                 b5 = [enrN1.x, enrN1.y]
# #                 ind5 = numpy.where(numpy.all(pvecCList==b5,axis=1))
# #                 c5 = ind5[0][0]
# #                 
# #                 enrN2 = root_i.enrichNodes[1]
# #                 
# # #                 enrN2.x = enrN2.x/1000.0
# # #                 enrN2.y = enrN2.y/1000.0
# # #                 if enrN2.x != 0.0:
# # #                     enrN2.x += 0.001
# # #                 if enrN2.y != 0.0:
# # #                     enrN2.y += 0.001
# # #                 enrN2.y = 1 - enrN2.y
# #                 
# #                 b6 = [enrN2.x, enrN2.y]
# #                 ind6 = numpy.where(numpy.all(pvecCList==b6,axis=1))
# #                 c6 = ind6[0][0]   
# #                 
# # #                 t = t + [[c1,c2,c3,c4,c5,c6]]
# #                 t = t + [[c4,c3,c2,c1, c5, c6]]                             
# #         else:
# # #             t = t + [[c1,c2,c3,c4]]
# #             t = t + [[c4,c3,c2,c1]]
            

        
    return [pvec,pvecCList,lenClist1]
 
def numbering(pvec,pvecCList, llist, masterNode):
    n = len(llist)
    
#     for i in range(0,len(pvec)):
#         pvec[i,1] = 1 - pvec[i,1] 
#         
#         
#     for i in range(0,len(p_reg)):
#         if pvec[i,0] != 0.0:
#             pvec[i,0] -= 0.001
#         if pvec[i,1] != 0.0:
#             pvec[i,1] -= 0.001
#     pvec = pvec * 1000.0
       
    t = []
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,llist[i])
        p1,p2,p3,p4 = root_i.rect
        
        
        b1 = [p1.x, p1.y]
        ind1 = numpy.where(numpy.all(pvecCList==b1,axis=1))
        c1 = ind1[0][0]

        b2 = [p2.x, p2.y]
        ind2 = numpy.where(numpy.all(pvecCList==b2,axis=1))
        c2 = ind2[0][0]
         
        b3 = [p3.x, p3.y]
        ind3 = numpy.where(numpy.all(pvecCList==b3,axis=1))
        c3 = ind3[0][0]
         
        b4 = [p4.x, p4.y]
        ind4 = numpy.where(numpy.all(pvecCList==b4,axis=1))
        c4 = ind4[0][0]
                 
        l = len(root_i.enrichNodes)
        if l > 0:
            if l == 1:
                enrN1 = root_i.enrichNodes[0]
                   
                b5 = [enrN1.x, enrN1.y]
                ind5 = numpy.where(numpy.all(pvecCList==b5,axis=1))
                c5 = ind5[0][0]
                if  ( (abs(enrN1.x - p1.x) <= TOL_error and abs(enrN1.y - p1.y) <= TOL_error) or
                      ( abs(enrN1.x - p2.x) <= TOL_error and abs(enrN1.y - p2.y) <= TOL_error ) or
                      ( abs(enrN1.x - p3.x) <= TOL_error and abs(enrN1.y - p3.y) <= TOL_error ) or
                      ( abs(enrN1.x - p4.x) <= TOL_error and abs(enrN1.y - p4.y) <= TOL_error) ):
                    t = t + [[c1,c2,c3,c4]]
                else:
                    t = t + [[c1,c2,c3,c4,c5]]
#                 t = t + [[c4,c3,c2,c1, c5]]
                  
            if l == 2:

                enrN1 = root_i.enrichNodes[0]
                b5 = [enrN1.x, enrN1.y]
                ind5 = numpy.where(numpy.all(pvecCList==b5,axis=1))
                c5 = ind5[0][0]
                  
                enrN2 = root_i.enrichNodes[1]
                b6 = [enrN2.x, enrN2.y]
                ind6 = numpy.where(numpy.all(pvecCList==b6,axis=1))
                c6 = ind6[0][0]   
                  
                if abs(enrN1.x - enrN2.x)<=TOL_error and abs(enrN1.y - enrN2.y)<=TOL_error:
                    t = t + [[c1,c2,c3,c4]]
                else:
#                     if c5>c6:
#                         t = t + [[c1,c2,c3,c4,c5,c6]]
#                     else:
#                         t = t + [[c1,c2,c3,c4,c6,c5]]
                    if pvecCList[c5][1] < pvecCList[c6][1]:
                        mtl = [[c1,c2,c3,c4,c6,c5]]
                    else:
                        if pvecCList[c5][1] == pvecCList[c6][1]:
#                             if c1 == 1147 and c2 == 1148:
#                                 print pvecCList[1786], pvecCList[1785]
#                                 print c5, c6         
                    
                            if pvecCList[c5][0] < pvecCList[c6][0]:
                                mtl = [[c1,c2,c3,c4,c5,c6]]
#                                 mtl = [[c1,c2,c3,c4,c6,c5]]
                            else:
                                mtl =  [[c1,c2,c3,c4,c6,c5]]
#                                 mtl = [[c1,c2,c3,c4,c6,c5]]
                        else:       
                            mtl =  [[c1,c2,c3,c4,c5,c6]]
                    
                    t = t + mtl
                    
#                 t = t + [[c4,c3,c2,c1, c5, c6]]                             
        else:
            t = t + [[c1,c2,c3,c4]]
            
#             t = t + [[c4,c3,c2,c1]]
        
    
        
    # convert from numbering system for Image Coordinate with (0,0) in the NW corner, to 
    # Euclidean coordinate system with (0,0) in the SW corner
    tvec = []
    for j in range(0,len(t)):
        el = t[j]
        new_corners = []
        for k in range(0,len(el)):
            tind = el[k]
            old_coords = pvecCList[tind]
            
            new_coords = old_coords/1000.0
            if new_coords[1] != 0.0:
                new_coords[1] += 0.001
            if new_coords[0] != 0.0:
                new_coords[0] += 0.001
            new_coords[1] = 1 - new_coords[1]
            
    
            new_indx = numpy.where(numpy.all(pvec==new_coords,axis=1))
            indx = new_indx[0][0]
            new_corners = new_corners + [indx]
        

        if len(el) == 4:
                tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0]]
                
        if len(el) == 5:
                tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[4]]
        if len(el) == 6:
#             if new_corners[5] < new_corners[4]:
#                 tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[5], new_corners[4]]
#             else:
#                 tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[4], new_corners[5]]
            if pvec[new_corners[4]][1] < pvec[new_corners[5]][1]:
                tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[4], new_corners[5]]
            else:
                if pvec[new_corners[4]][1] == pvec[new_corners[5]][1]:
                    if pvec[new_corners[4]][0] < pvec[new_corners[5]][0]:
                        tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[4], new_corners[5]]
                    else:
                        tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[5], new_corners[4]]
                else:
                        tk = [new_corners[3], new_corners[2], new_corners[1], new_corners[0], new_corners[5], new_corners[4]]
                 
#             if tk[0] == 270 and tk[1] == 271:
#                     print tk       
#                     print new_corners
#                     print pvec[1529],pvec[1528] 
                    
#             print '----',pvec[new_corners[4]], pvec[new_corners[5]]
#             print tk
#             print pvec[1041], pvec[1042]
                       
        tvec = tvec + [tk]

    for i in range(0,n):
        root_i = get_node_by_id(masterNode,llist[i])
        root_i.tlist = tvec[i]
        root_i.tpix = t[i]
    return tvec,t

def set_nsew(llist, masterNode, full_vec):
    
    n = len(llist)
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])

        p1,p2,p3,p4 = root.rect
                        
        west_neigh_index = str(find_neighbor_of(root.index,'L'))    
        # checking to see if the west neighbor exists or is a ghost
        if it_exists(west_neigh_index, masterNode):
            west_neighbor = get_node_of_neighbor(root, root.index, west_neigh_index)
            p1w,p2w,p3w,p4w = west_neighbor.rect
            if west_neighbor.has_children == True and p1w.x < p1.x:
                root.nsew[3] = 1
                p1wchild1,p2wchild1,p3wchild1,p4wchild1 = west_neighbor.children[1].rect
                x1 = p3wchild1.x / 1000.0
                y1 = p3wchild1.y / 1000.0
                if x1 != 0.0:
                    x1 += 0.001
                if y1 != 0.0:
                    y1 += 0.001            
                y1 = 1 -  y1    
                x1 = find_nearest(full_vec,x1)
                y1 = find_nearest(full_vec,y1)  
                root.hn[3] = Coordinate(x1, y1)
                west_has_children = True
                
        east_neigh_index = str(find_neighbor_of(root.index,'R'))    
        # checking to see if the west neighbor exists or is a ghost
        if it_exists(east_neigh_index, masterNode):
            east_neighbor = get_node_of_neighbor(root, root.index, east_neigh_index)
            p1e,p2e,p3e,p4e = east_neighbor.rect
            if east_neighbor.has_children == True and p2.x < p2e.x:
                root.nsew[2] = 1
                p1echild0, p2echild0, p3echild0, p4echild0 = east_neighbor.children[0].rect
                x4 = p4echild0.x / 1000.0
                y4 = p4echild0.y / 1000.0
                if x4 != 0.0:
                    x4 += 0.001
                if y4 != 0.0:
                    y4 += 0.001            
                y4 = 1 - y4     
                x4 = find_nearest(full_vec,x4)
                y4 = find_nearest(full_vec,y4)
                root.hn[2] = Coordinate(x4, y4)
#                east_has_children = True

        south_neigh_index = str(find_neighbor_of(root.index,'D'))  
        # checking to see if the west neighbor exists or is a ghost
        if it_exists(south_neigh_index, masterNode):
            south_neighbor = get_node_of_neighbor(root, root.index, south_neigh_index)
            p1s,p2s,p3s,p4s = south_neighbor.rect
            if south_neighbor.has_children == True and p4.y < p4s.y:
                root.nsew[1] = 1
                p1schild0, p2schild0, p3schild0, p4schild0 = south_neighbor.children[0].rect
                x2 = p2schild0.x / 1000.0
                y2 = p2schild0.y / 1000.0
                if x2 != 0.0:
                    x2 += 0.001
                if y2 != 0.0:
                    y2 += 0.001            
                y2 = 1 - y2
                x2 = find_nearest(full_vec,x2)
                y2 = find_nearest(full_vec,y2)
                root.hn[1] = Coordinate(x2, y2)
#        
        north_neigh_index = str(find_neighbor_of(root.index,'U'))    
        # checking to see if the west neighbor exists or is a ghost
        if it_exists(north_neigh_index, masterNode):
            north_neighbor = get_node_of_neighbor(root, root.index, north_neigh_index)
            p1n,p2n,p3n,p4n = north_neighbor.rect
            if north_neighbor.has_children == True and p1n.y < p1.y:
                root.nsew[0] = 1
                p1nchild2, p2nchild2, p3nchild2, p4nchild2 = north_neighbor.children[2].rect
                x3 = p3nchild2.x / 1000.0
                y3 = p3nchild2.y / 1000.0
                if x3 != 0.0:
                    x3 += 0.001
                if y3 != 0.0:
                    y3 += 0.001            
                y3 = 1 - y3
                x3 = find_nearest(full_vec,x3)
                y3 = find_nearest(full_vec,y3)
                root.hn[0] = Coordinate(x3, y3 )

def correct_pvec(p,full_vec,lenClist1,llist,pvecPx):
# correct p_vec of coordinates because of lost pixels in integer division
# now pixel at location 62 is correctly set to be 0.0625 and not 0.63
    
#     print 'start p-945',p[945]
       
    # go over the coordinates of the regular grid in the p vector  
    for i in range(0, lenClist1):
        for j in [0,1]:
            val = find_nearest(full_vec,p[i,j])

#             # go over the coordinates of the intersection nodes in the pvecor
#             for k in range(lenClist1,len(p)):
#                 # if any of the coordinates of the intersection nodes corresponds
#                 # to a coordinate value on the grid, correct it
#                 if p[k,0] == p[i,j]: 
#                     p[k,0] = val
#                 if p[k,1] == p[i,j]:
#                     p[k,1] = val
            p[i,j] = val



#     print 'mid p-945', p[945]
    # go over the coordinates of the intersection nodes
    # and update them accordingly

    n = len(llist)
    vec_indices = []
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])
        p1,p2,p3,p4 = root.rect
        
        tl = root.tlist
        tp = root.tpix

        
        # look to see if the element has enrichment nodes
        if len(tl) > 4:
            
#             enrich1 = [root.enrichNodes[0].x, root.enrichNodes[0].y]
#             enrich2 = [root.enrichNodes[1].x, root.enrichNodes[1].y]

            enrich1 = pvecPx[tp[4]]
            enrich2 = pvecPx[tp[5]]
            coords = numpy.array([[p1.x, p1.y],[p2.x,p2.y],[p3.x,p3.y],[p4.x,p4.y]])
            
#             if pvecPx[tp[4]][1] < pvecPx[tp[5]][1]:
#                 enrich1 = pvecPx[tp[5]]
#                 enrich2 = pvecPx[tp[4]]
#             else:
#                 if pvecPx[tp[4]][1] == pvecPx[tp[5]][1]:
#                     if pvecPx[tp[4]][0] < pvecPx[tp[5]][0]:
#                         
#                         enrich1 = pvecPx[tp[4]]
#                         enrich2 = pvecPx[tp[5]]
#                     else:
#                         enrich1 = pvecPx[tp[5]]
#                         enrich2 = pvecPx[tp[4]]
                        
            ind1 = tl[4]
            ind2 = tl[5]
            
#            if ind2 == 1059:
#                print tl, tp
#                print enrich1
#                print enrich2
#                root.printRect()
#                print not(on_corners(enrich1,coords)) , not(on_corners(enrich2,coords))
#             if ind2 == 1059:
#                 print p[1059]
#                 print 'ce faci aici ba?', tl, root.printRect()
#                 print 'iteration 2', i
          
#             if ind1 == 1059:
#                 print 'iteration 1', i
#                 print p[1528]
#                 print enrich1[0], enrich1[1]
#                 print not(on_corners(enrich1,coords)) , not(on_corners(enrich2,coords))

#             # enrichment node 1 is one of the corners of the element
#             if on_corners(enrich1,p1.x,p1.y,p3.x,p3.y):
#                 x1 = enrich1[0] / 1000.0
#                 y1 = enrich1[1] / 1000.0
#                 if x1 != 0.0:
#                     x1 += 0.001
#                 if y1 != 0.0:
#                     y1 += 0.001            
#                 y1 = 1 - y1
#                 valx1 = find_nearest(full_vec,x1)
#                 valy1 = find_nearest(full_vec,y1)
#   
#                 # update the coordinates to the nearest "point on the grid"                
#                 p[ind1,0] = valx1
#                 p[ind1,1] = valy1
#                
#                
#             # enrichment node 2 is one of the corners of the element   
#             if on_corners(enrich2,p1.x,p1.y,p3.x,p3.y):
#                 x2 = enrich2[0] / 1000.0
#                 y2 = enrich2[1] / 1000.0
#                 if x2 != 0.0:
#                     x2 += 0.001
#                 if y2 != 0.0:
#                     y2 += 0.001            
#                 y2 = 1 - y2
#                 valx2 = find_nearest(full_vec,x2)
#                 valy2 = find_nearest(full_vec,y2)
#                   
#                 # update the coordinates to the nearest "point on the grid"
#                 p[ind2,0] = valx2
#                 p[ind2,1] = valy2
# #                  
            # enrichment node 1 is on an edge
#             if not(on_corners(enrich1,p1.x,p1.y,p3.x,p3.y)):

            if ( not(on_corners(enrich1,coords)) and (on_corners(enrich2,coords)) ):
                x1 = enrich1[0] / 1000.0
                y1 = enrich1[1] / 1000.0
                if x1 != 0.0:
                    x1 += 0.001
                if y1 != 0.0:
                    y1 += 0.001            
                y1 = 1 - y1
                
                        
                if ((enrich1[0] == p1.x or enrich1[0] == p3.x) and 
                 (enrich1[1] != p1.y and enrich1[1] != p3.y)):
                
                    valx1 = find_nearest(full_vec,x1)
                    valy1 = y1
                    
                    # update the X - coordinate to the nearest "point on the grid"
                    p[ind1,0] = valx1
                    p[ind1,1] = valy1
#                     if ind1 == 1528 or ind2 == 1528:
#                         print 'rerwrewrwerwerewrewrwrew'
#                         print valx1, find_nearest(full_vec,x1)
#                         print x1, y1
#                         print p[ind1,0], ind1, ind2
#                         print 'tl', tl
#                         print 'tp', tp
#                         print enrich1[0], enrich1[1]
#                         print enrich2[0], enrich2[1]
#                         print p1.x, p3.x
#                         print p1.y, p3.y
                                    
                if ((enrich1[0] != p1.x and enrich1[0] != p3.x) and 
                 (enrich1[1] == p1.y or enrich1[1] == p3.y)):  
                     
                    valx1 = x1
                    valy1 = find_nearest(full_vec,y1)
                                    
                    # update the Y - coordinate to the nearest "point on the grid"    
                    p[ind1,0] = valx1
                    p[ind1,1] = valy1
                    

            if ( on_corners(enrich1,coords) and not(on_corners(enrich2,coords)) ):
                x1 = enrich2[0] / 1000.0
                y1 = enrich2[1] / 1000.0
                if x1 != 0.0:
                    x1 += 0.001
                if y1 != 0.0:
                    y1 += 0.001            
                y1 = 1 - y1
                
                if  ((enrich2[0] == p1.x or enrich2[0] == p3.x) and 
                 (enrich2[1] != p1.y and enrich2[1] != p3.y)):
                
                    valx1 = find_nearest(full_vec,x1)
                    valy1 = y1

                    # update the X - coordinate to the nearest "point on the grid"
                    p[ind2,0] = valx1
                    p[ind2,1] = valy1
                    
#                     if ind1 == 1528 or ind2 == 1528:
#                         print 'inside first part'
#                         print find_nearest(full_vec,0.203125)
#                         print x1,y1
#                         print valx1
                        
                if ((enrich2[0] != p1.x and enrich2[0] != p3.x) and 
                 (enrich2[1] == p1.y or enrich2[1] == p3.y)):                
                    valx1 = x1
                    valy1 = find_nearest(full_vec,y1)

                    # update the Y - coordinate to the nearest "point on the grid"    
                    p[ind2,0] = valx1
                    p[ind2,1] = valy1                  
                    

            if (not(on_corners(enrich1,coords)) and not(on_corners(enrich2,coords))
               # and not(ind1 in vec_indices) #and not(ind2 in vec_indices)
                ):

                x1 = enrich1[0] / 1000.0
                y1 = enrich1[1] / 1000.0
                if x1 != 0.0:
                    x1 += 0.001
                if y1 != 0.0:
                    y1 += 0.001            
                y1 = 1 - y1
                
                        
                if ((enrich1[0] == p1.x or enrich1[0] == p3.x) and 
                 (enrich1[1] != p1.y and enrich1[1] != p3.y)):
                
                    valx1 = find_nearest(full_vec,x1)
                    valy1 = y1
                    
                    # update the X - coordinate to the nearest "point on the grid"
                    p[ind1,0] = valx1
                    p[ind1,1] = valy1
#                     if ind1 == 1528 or ind2 == 1528:
#                         print 'rerwrewrwerwerewrewrwrew'
#                         print valx1, find_nearest(full_vec,x1)
#                         print x1, y1
#                         print p[ind1,0], ind1, ind2
#                         print 'tl', tl
#                         print 'tp', tp
#                         print enrich1[0], enrich1[1]
#                         print enrich2[0], enrich2[1]
#                         print p1.x, p3.x
#                         print p1.y, p3.y
                                    
                if ((enrich1[0] != p1.x and enrich1[0] != p3.x) and 
                 (enrich1[1] == p1.y or enrich1[1] == p3.y)):  
                     
                    valx1 = x1
                    valy1 = find_nearest(full_vec,y1)
                                    
                    # update the Y - coordinate to the nearest "point on the grid"    
                    p[ind1,0] = valx1
                    p[ind1,1] = valy1
                
            # enrichment node 2 is on an edge
            coords = numpy.array( [[p1.x, p1.y],[p2.x,p2.y],[p3.x,p3.y],[p4.x,p4.y]])
                            
            if (not(on_corners(enrich2,coords)) and not(on_corners(enrich1,coords))
                #and not(ind2 in vec_indices) 
                ):
#                 if ind1 == 1528:
#                         print 'inside third after',p[1528]  
#                         print ((enrich2[0] == p1.x or enrich2[0] == p3.x) and 
#                  (enrich2[1] != p1.y and enrich2[1] != p3.y))
#                         print ((enrich2[0] != p1.x and enrich2[0] != p3.x) and 
#                  (enrich2[1] == p1.y or enrich2[1] == p3.y))
                        
                x1 = enrich2[0] / 1000.0
                y1 = enrich2[1] / 1000.0
                if x1 != 0.0:
                    x1 += 0.001
                if y1 != 0.0:
                    y1 += 0.001            
                y1 = 1 - y1
                
                if ((enrich2[0] == p1.x or enrich2[0] == p3.x) and 
                 (enrich2[1] != p1.y and enrich2[1] != p3.y)):
                
                    valx1 = find_nearest(full_vec,x1)
                    valy1 = y1

                    # update the X - coordinate to the nearest "point on the grid"
                    p[ind2,0] = valx1
                    p[ind2,1] = valy1
                    
#                     if ind1 == 1528 or ind2 == 1528:
#                         print 'inside first part'
#                         print find_nearest(full_vec,0.203125)
#                         print x1,y1
#                         print valx1
                        
                if ((enrich2[0] != p1.x and enrich2[0] != p3.x) and 
                 (enrich2[1] == p1.y or enrich2[1] == p3.y)):                
                    valx1 = x1
                    valy1 = find_nearest(full_vec,y1)

                    # update the Y - coordinate to the nearest "point on the grid"    
                    p[ind2,0] = valx1
                    p[ind2,1] = valy1                  
#                     if ind2 == 1528 or ind1==1528:
#                         print 'inside second part'
#                         print find_nearest(full_vec,0.203125)
#                         print x1,y1
#                         print valx1
    
            if enrich1[0] == enrich2[0] or enrich1[1] == enrich2[1]:
                root.ishomog = 1
            else:
                root.ishomog = 0
                
#             if root.index == '323210':
#                 print '--------323210----------------'
#                 print [p1.x,p3.x], [p1.y,p3.y]
#                 print 'tl', tl
#                 print 'tp', tp
#                 print enrich1, p[ind1]
#                 print enrich2, p[ind2]
#                 print root.index, root.ishomog
            
            vec_indices = vec_indices + [ind1, ind2]

    
    return p

def set_homogOLD(masterNode,llist,pvecPx):
    n = len(llist)
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])

        tp = root.tpix
    
        if len(tp)>4:
            enrich1 = pvecPx[tp[4]]
            enrich2 = pvecPx[tp[5]]
            if enrich1[0] == enrich2[0] or enrich1[1] == enrich2[1]:
                root.ishomog = 1
            else:
                root.ishomog = 0
        else:
            root.ishomog = 1

#        print root.index, root.ishomog

def find_index(array,value):
    idx = (numpy.abs(array-value)) == 0
    return idx

def find_nearest(array,value):
    idx = (numpy.abs(array-value)).argmin()
    return array[idx]
 
                
def element_normal_intersection(pt1,pt2,node,image):
    
    NoneINT = -9999
    midpt = find_mid_point(pt1,pt2)
    
    dx = pt1.x - pt2.x
    dy = pt1.y - pt2.y
    
    N1 = [-dy, dx] #normal 1
    N2 = [dy, -dx] #normal 2
    
    ptN1 = Coordinate( midpt.x + N1[0], midpt.y + N1[1])
    ptN2 = Coordinate( midpt.x + N2[0], midpt.y + N2[1])

    
    if 0 <= ptN1.x and ptN1.x <= node.imsize[0] and 0 <= ptN1.y and ptN1.y <= node.imsize[1]:
        ptN = Coordinate(ptN1.x, ptN1.y)
    elif 0 <= ptN2.x and ptN2.x <= node.imsize[0] and 0 <= ptN2.y and ptN2.y <= node.imsize[1]:
        ptN = Coordinate(ptN2.x, ptN2.y)
        
    dx_m = ptN.x - midpt.x
    dy_m = ptN.y - midpt.y

    
    p1,p2,p3,p4 = node.rect

#    print dx, dy, dx_m, dy_m
#    print ptN.x, ptN.y
#    print 'first pt',pt1.x, pt1.y
#    print pt2.x, pt2.y
#    print midpt.x, midpt.y
    
    #Compute the intersection of the normal with the 4 sides of an element
    if dx_m == 0: #vertical
        side1 = Coordinate(midpt.x, p1.y)
        side2 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
        side3 = Coordinate(midpt.x, p4.y)
        side4 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
#        draw_line_normals(image, side1,side3)
    elif dy_m == 0: #horizontal
        side1 = Coordinate(NoneINT, NoneINT) # normal runs parallel with the edge
        side2 = Coordinate(p2.x, midpt.y)
        side3 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
        side4 = Coordinate(p1.x, midpt.y)
#        draw_line_normals(image, side2,side4)
    else:
        m_slope = float(dy_m) / dx_m
        b = midpt.y - m_slope * midpt.x
#         print 'suntem aici!!!!', m_slope, b
#         node.printRect()
#         print 'side1', float(p1.y - b) / m_slope
#         print 'side4', float(p4.y - b) / m_slope 
        # SIDE 1
        if p1.x <= float(p1.y - b) / m_slope and float(p1.y - b) / m_slope <= p2.x:
            side1 = Coordinate( float(p1.y - b) / m_slope, p1.y)
            side1.x = int(side1.x)
            side1.y = int(side1.y)
#            draw_line_normals(image, midpt, side1)
        else: # intersection happens outisde the element
            side1 = Coordinate(NoneINT, NoneINT)
        
        # SIDE 2    
        if p1.y <= m_slope * p2.x + b and  m_slope * p2.x + b <= p4.y:
            side2 = Coordinate( p2.x, m_slope * p2.x + b)
            side2.x = int(side2.x)
            side2.y = int(side2.y)
#            draw_line_normals(image, midpt, (side2))
        else: # intersection happens outisde the element
            side2 = Coordinate(NoneINT, NoneINT)
        
        # SIDE 3
        if p1.x <= float(p4.y - b) / m_slope and float(p4.y - b) / m_slope <= p2.x:
            side3 = Coordinate( float(p4.y - b) / m_slope, p4.y)
            side3.x = int(side3.x)
            side3.y = int(side3.y)
#            draw_line_normals(image, midpt, side3)
        else: # intersection happens outisde the element
            side3 = Coordinate(NoneINT, NoneINT)
        # SIDE 4    
        if p1.y <= m_slope * p1.x + b and m_slope * p1.x + b <= p4.y:
            side4 = Coordinate( p1.x, m_slope * p1.x + b)
            side4.x = int(side4.x)
            side4.y = int(side4.y)
#            draw_line_normals(image, midpt, side4)
        else: # intersection happens outisde the element
            side4 = Coordinate(NoneINT, NoneINT)

    side1.x = int(side1.x)
    side1.y = int(side1.y)
    side2.x = int(side2.x)
    side2.y = int(side2.y)
    side3.x = int(side3.x)
    side3.y = int(side3.y)
    side4.x = int(side4.x)
    side4.y = int(side4.y)
    
    N_edge = 0
    S_edge = 0
    W_edge = 0
    E_edge = 0
    NW_edge = 0
    NE_edge = 0
    SE_edge = 0
    SW_edge = 0
    direction_list = []
    whichSide = []
    
#    node.printRect()
#    print side1.x, side1.y
#    print side2.x, side2.y
#    print side3.x, side3.y
#    print side4.x, side4.y
    
    
    if (side1.x != NoneINT and side1.x != p1.x and side1.x != p2.x and
        side1.y != NoneINT and (side1.y == p1.y or side1.y == p2.y) ):
#         print 'normal going through North edge' 
        N_edge = 1
        direction_list.append('U')
        whichSide.append(side1)
        
    if (side2.x != NoneINT and (side2.x == p3.x or side2.x == p2.x) and
        side2.y != NoneINT and side2.y != p3.y and side2.y != p2.y):
#         print 'normal going through East edge' 
#         print p2.x, p3.x, p3.y, p2.y, side2.x, side2.y
        E_edge = 1
        direction_list.append('R')
        whichSide.append(side2)
        
    if (side3.x != NoneINT and side3.x != p3.x and side3.x != p4.x and
        side3.y != NoneINT and (side3.y == p3.y or side3.y == p4.y) ):
#         print 'normal going through South edge' 
        S_edge = 1
        direction_list.append('D')
        whichSide.append(side3)
        
    if (side4.x != NoneINT and (side4.x == p1.x or side4.x == p4.x) and
        side4.y != NoneINT and side4.y != p1.y and side4.y != p4.y):
#         print 'normal going through West edge' 
        W_edge = 1
        direction_list.append('L')
        whichSide.append(side4)
        
    if (side1.x == p2.x and side1.y == p2.y) or (side2.x == p2.x and side2.y == p2.y):
#         print 'normal going through North East corner' 
        NE_edge = 1
        direction_list.append('RU')
        if (side1.x == p2.x and side1.y == p2.y) :
            whichSide.append(side1)
        else:
            whichSide.append(side2)
        
    if (side1.x == p1.x and side1.y == p1.y) or (side4.x == p1.x and side4.y == p1.y):
#         print 'normal going through North West corner'
        NW_edge = 1
        direction_list.append('LU')
        if (side1.x == p1.x and side1.y == p1.y):
            whichSide.append(side1)
        else:
            whichSide.append(side4)
    
    if (side2.x == p3.x and side2.y == p3.y) or (side3.x == p3.x and side3.y == p3.y):
#         print 'normal going through South East corner' 
        SE_edge = 1
        direction_list.append('RD')
        if (side2.x == p3.x and side2.y == p3.y):
            whichSide.append(side2)
        else:
            whichSide.append(side3)
        
    if (side3.x == p4.x and side3.y == p4.y) or (side4.x == p4.x and side4.y == p4.y):
#         print 'normal going through South West edge' 
        SW_edge = 1
        direction_list.append('LD')
        if (side3.x == p4.x and side3.y == p4.y) :
            whichSide.append(side3)
        else:
            whichSide.append(side4)
        
    neigh_list = [NW_edge,N_edge,NE_edge,E_edge,SE_edge,S_edge,SW_edge,W_edge]
    return [side1,side2,side3,side4,neigh_list, direction_list,whichSide]
        
def set_homog(masterNode,llist):
    n = len(llist)
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])

#        if len(root.enrichNodes)>1:
#            print root.index
            
def find_neighbor_index_of(index,direction, masterNode, llist):
    loc = str(index)
    llist_str = list(loc)

    
    
    
    for i in range(len(loc)-1,-1,-1):
#        print str(loc[i]), direction, loc
        new_quadrant =  D[str(loc[i])][direction]['Quadrant']
        new_direction = D[str(loc[i])][direction]['Direction']
        
        
        if new_direction != 'H':
        
            direction = new_direction
            llist_str[i] = str(new_quadrant)
        
        else:
            
            llist_str[i] = str(new_quadrant)
#             return str("".join(llist_str))
            break

    node = get_node_by_id(masterNode,[str(index)])
    p1,p2,p3,p4 = node.rect
    neigh_index = str("".join(llist_str))    
    
    lengthSize = node.imsize[1]-1
    widthSize = node.imsize[0]-1
    
    if node_exists(neigh_index,llist):
       
        if direction == 'U' and p1.y != 0:
            return [str(neigh_index)]
        if direction == 'D' and p4.y != lengthSize:
            return [str(neigh_index)]
        if direction == 'L' and p1.x != 0:
            return [str(neigh_index)]
        if direction == 'R' and p3.x != widthSize:
            return [str(neigh_index)]
        if direction == 'LU' and p1.x != 0 and p1.y != 0:
            return [str(neigh_index)]
        if direction == 'RU' and p2.x != widthSize and p1.y != 0:
            return [str(neigh_index)]
        if direction == 'LD' and p4.x != 0 and p4.y != lengthSize:
            return [str(neigh_index)]
        if direction == 'RD' and p3.x != widthSize and p3.y != lengthSize:
            return [str(neigh_index)]
        
    else:
        parent_of_neigh_index = neigh_index[:-1]
        uncle_node = get_node_by_id(masterNode,[str(parent_of_neigh_index)])
    
        if direction == 'U' and p1.y == 0:
            return []
        if direction == 'D' and p4.y == lengthSize:
            return []
        if direction == 'L' and p1.x == 0:
            return []
        if direction == 'R' and p3.x == widthSize:
            return []
        if direction == 'LU' and p1.x == 0 and p1.y == 0:
            return []
        if direction == 'RU' and p2.x == widthSize and p1.y == 0:
            return []
        if direction == 'LD' and p4.x == 0 and p4.y == lengthSize:
            return []
        if direction == 'RD' and p3.x == widthSize and p3.y == lengthSize:
            return []
                    
#        print str(index), neigh_index, parent_of_neigh_index
        if uncle_node.has_children == False:
   
            if direction == 'U' and p1.y != 0:
                return [str(parent_of_neigh_index)]
            if direction == 'D' and p4.y != lengthSize:
                return [str(parent_of_neigh_index)]
            if direction == 'L' and p1.x != 0:
                return [str(parent_of_neigh_index)]
            if direction == 'R' and p3.x != widthSize:
#                 print 'yayayayayayayayaydadadadad'
                return [str(parent_of_neigh_index)]
            if direction == 'LU' and p1.x != 0 and p1.y != 0:
                return [str(parent_of_neigh_index)]
            if direction == 'RU' and p2.x != widthSize and p1.y != 0:
                return [str(parent_of_neigh_index)]
            if direction == 'LD' and p4.x != 0 and p4.y != lengthSize:
                return [str(parent_of_neigh_index)]
            if direction == 'RD' and p3.x != widthSize and p3.y != lengthSize:
                return [str(parent_of_neigh_index)]
        else:
            if direction == 'L' and p1.x != 0:
                return [str(neigh_index) + '1', str(neigh_index) + '3']
            if direction == 'R' and p3.x != widthSize:
                return [str(neigh_index) + '0', str(neigh_index) + '2']
            if direction == 'U' and p1.y != 0:
                return [str(neigh_index) + '2', str(neigh_index) + '3']
            if direction == 'D' and p4.y != lengthSize:
                return [str(neigh_index) + '0', str(neigh_index) + '1']
            if direction == 'LU' and p1.x != 0 and p1.y != 0:
                return [str(neigh_index) + '3']
            if direction == 'RU' and p2.x != widthSize and p1.y != 0:
                return [str(neigh_index) + '2']
            if direction == 'RD' and p3.x != widthSize and p3.y != lengthSize:
                return [str(neigh_index) + '0']
            if direction == 'LD' and p4.x != 0 and p4.y != lengthSize:
                return [str(neigh_index) + '1']

        
    return []
  
def draw_interface(image, tree_list, masterNode):
    
    n = len(tree_list)

    # for each node in the tree:
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,tree_list[i])    
         
        
        
        if len(root_i.enrichNodes) > 1:
            if root_i.enrichNodes[0].x <= root_i.enrichNodes[1].x:
                P1 = root_i.enrichNodes[0]
                P2 = root_i.enrichNodes[1]
            else:
                P1 = root_i.enrichNodes[1]
                P2 = root_i.enrichNodes[0]
            
            if len(root_i.enrichNodes) == 2:
                draw_line(image,P1, P2)
                
            if len(root_i.enrichNodes) == 3:
                draw_line(image,P1, root_i.enrichNodes[2])
                draw_line(image,root_i.enrichNodes[2], P2)
#                 draw_line(image,root_i.enrichNodes[1], root_i.enrichNodes[2])
#                 draw_line(image,root_i.enrichNodes[2], root_i.enrichNodes[0])
                print ' quadratics!!!!'
            if len(root_i.enrichNodes) == 4:
                print ' cubics!!!!'
                if root_i.enrichNodes[2].x <= root_i.enrichNodes[3].x:
                    E = root_i.enrichNodes[2]
                    F = root_i.enrichNodes[3]
                else:
                    E = root_i.enrichNodes[3]
                    F = root_i.enrichNodes[2]
                draw_line(image,P1,E)
                draw_line(image,E,F)
                draw_line(image,F,P2)
                    
        
                 
def stress_concentration_constraint(tree_list, masterNode, image):

    n = len(tree_list)

    full_list = []
    
    # for each node in the tree:
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,tree_list[i])    
         
        if len(root_i.enrichNodes) > 1: # for each non-hom node
        # root_i.index=='312':
            [side1,side2,side3,side4,neigh_list,dir_list,whichSide] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], root_i, image)
#            print '-----------------'
#            print root_i.index
#            root_i.printRect()
#            print root_i.enrichNodes[0].x, root_i.enrichNodes[0].y
#            print root_i.enrichNodes[1].x, root_i.enrichNodes[1].y
#            print find_neighbor_index_of(root_i.index,dir_list[0], masterNode, tree_list)
#            print find_neighbor_index_of(root_i.index,dir_list[1], masterNode, tree_list)
#            print dir_list
#            print '++++++++++++++++'

#             print neigh_list
            counter1 = 0
            counter2 = 0
            
            list1 = []
            list2 = []
            
            whichEdge1 = dir_list[0]
            whichEdge2 = dir_list[1]
            
            currentIndex1 = root_i.index
            currentIndex2 = root_i.index
            list1.append(currentIndex1)
            list2.append(currentIndex2)

            whichSidePrev1 = copy_list_of_sides(whichSide)

            while counter1 <= 4:
                
                neighs = find_neighbor_index_of(currentIndex1,whichEdge1, masterNode, tree_list)
                
                if len(neighs) == 2: # there are 2 neighbors sharing an edge with me
                    neigh1 = get_node_by_id(masterNode,[str(neighs[0])]) 
                    neigh2 = get_node_by_id(masterNode,[str(neighs[1])])

                    p1n1,p2n1,p3n1,p4n1 = neigh1.rect
                    
                    if whichEdge1 == 'U':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side1.x >= p1n1.x and side1.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge1 == 'D':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side3.x >= p1n1.x and side3.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge1 == 'L':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side4.y >= p1n1.y and side4.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                            
                    if whichEdge1 == 'R':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side2.y >= p1n1.y and side2.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                    
                else:
                    if len(neighs) == 0: # there are no neighbors of mine, perhaps we are at the margins of the image
                        break
                    else:
                        neighIndex = str(neighs[0])

                neigh_node = get_node_by_id(masterNode,[str(neighIndex)])
                if not(neigh_node.ishomog):
                    #neighbor is non-homogeneous:                    
                    list1.append(neighIndex)
                    break
                else:
                    counter1 += 1
                    [side1N,side2N,side3N,side4N,neigh_listN,dir_listN,whichSideN] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], neigh_node, image)
                    currentIndex1 = neighIndex
#                     whichEdge1 = swap_edges(whichEdge1, dir_listN)
#                     print 'counter 1, whichEdge=', whichEdge1, 'neigh index = ', neighIndex
#                     print 'dirlistN=', dir_listN, 'whichSideN=', whichSideN
#                     print side1N.x, side1N.y, side2N.x,side2N.y, side3N.x, side3N.y, side4N.x,side4N.y
                    
                    # problem when channelsCircles is tested
                    if len(dir_listN) < 1 or len(whichSideN)<2:
                        list1.append(neighIndex)
                        break
#                     root_i.printRect()
#                     neigh_node.printRect()
#                    print dir_listN,whichSidePrev1,whichSideN
#                    print whichEdge1
                    whichEdge1 = swap_directions(dir_listN,whichSidePrev1,whichSideN)
#                    print whichEdge1
                    
                    # problem when channelsCircles is tested
                    if whichEdge1 == None:
                        list1.append(neighIndex)
                        break
                    
                    whichSidePrev1 = copy_list_of_sides(whichSideN)
                    list1.append(neighIndex)
                    side1.x = side1N.x
                    side1.y = side1N.y
                    side2.x = side2N.x
                    side2.y = side2N.y
                    side3.x = side3N.x
                    side3.y = side3N.y
                    side4.x = side4N.x
                    side4.y = side4N.y
                    
                 
#             print root_i.index
#             print list1   
            #===================================================================
#             process_list(list1,masterNode)
            #===================================================================
            full_list.append(list1)
            
            [side1,side2,side3,side4,neigh_list,dir_list,whichSide] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], root_i,image)

            whichSidePrev2 = copy_list_of_sides(whichSide)
            dir_listN = dir_list
            while counter2 <= 4:
                
                neighs = find_neighbor_index_of(currentIndex2,whichEdge2, masterNode, tree_list)
                
                if len(neighs) == 2: # there are 2 neighbors sharing an edge with me
                    neigh1 = get_node_by_id(masterNode,[str(neighs[0])]) 
                    neigh2 = get_node_by_id(masterNode,[str(neighs[1])])

                    p1n1,p2n1,p3n1,p4n1 = neigh1.rect
                    
                    if whichEdge2 == 'U':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side1.x >= p1n1.x and side1.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge2 == 'D':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side3.x >= p1n1.x and side3.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge2 == 'L':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side4.y >= p1n1.y and side4.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                            
                    if whichEdge2 == 'R':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side2.y >= p1n1.y and side2.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                    
                else:
                    if len(neighs) == 0: # there are no neighbors of mine, perhaps we are at the margins of the image
                        break
                    else:
                        neighIndex = str(neighs[0])

                neigh_node = get_node_by_id(masterNode,[str(neighIndex)])
                if not(neigh_node.ishomog):
                    #neighbor is non-homogeneous:                    
                    list2.append(neighIndex)
                    break
                else:
                    counter2 += 1
                    [side1N,side2N,side3N,side4N,neigh_listN,dir_listN,whichSideN] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], neigh_node, image)
                    
                    # problem when channelsCircles is tested
                    if len(dir_listN) < 1 or len(whichSideN)<2:
                        list2.append(neighIndex)
                        break
                    
                    currentIndex2 = neighIndex
                    
#                    neigh_node.printRect()
#                    print 'dirlist', dir_listN, neighIndex
#                    print side1N.x, side1N.y, side2N.x,side2N.y, side3N.x, side3N.y, side4N.x,side4N.y
#                    print 'old edge:', whichEdge2
#                    print dir_listN,whichSidePrev2,whichSideN
#                     whichEdge2 = swap_edges(whichEdge2, dir_listN)
                    whichEdge2 = swap_directions(dir_listN,whichSidePrev2,whichSideN)
                    
                    # problem when testing channelsCircles
                    if whichEdge2 == None:
                        list2.append(neighIndex)
                        break
                    
                    whichSidePrev2 = copy_list_of_sides(whichSideN)
                    list2.append(neighIndex)
                    side1.x = side1N.x
                    side1.y = side1N.y
                    side2.x = side2N.x
                    side2.y = side2N.y
                    side3.x = side3N.x
                    side3.y = side3N.y
                    side4.x = side4N.x
                    side4.y = side4N.y
                    
#             print list2
            #===================================================================
#             masterNode = process_list(list2,masterNode)
            #===================================================================
            full_list.append(list2)

    return full_list

def copy_list_of_sides(sides1):
    sides2 = [Coordinate(0,0), Coordinate(0,0)]
    
    sides2[0].x = deepcopy(sides1[0].x)
    sides2[0].y = deepcopy(sides1[0].y)
    
    sides2[1].x = deepcopy(sides1[1].x)
    sides2[1].y = deepcopy(sides1[1].y)
    
    return sides2

        
def process_list(full_list, masterNode,image):
    # processing the list with elements in between interfaces
    
    extra_list = []
    for k in range(0, len(full_list)):
        llist = full_list[k]
        
        if len(llist) < STRESS_MIN + 2:
            last_index = str(llist[-1])
            last_node = get_node_by_id(masterNode, [str(last_index)])
            if not(last_node.ishomog): 
            
                if len(llist) == 2:
                    node1 = get_node_by_id(masterNode, [str(llist[0])])
                    node1.divideOnce()
                    node2 = get_node_by_id(masterNode, [str(llist[1])])
                    node2.divideOnce()
                if len(llist) == 2 + 1:
    #                 print 'only one homog elem in between'
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    extra_list.append(node1)
                    
                        
                if len(llist) == 2 + 2:
    #                 print 'only 2 homog elems in between'
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    node2 =  get_node_by_id(masterNode, [str(llist[2])])
                    node2.divideOnce()
                    
                if len(llist) == 2 + 3:
    #                 print 'only 3 homog elems in between'
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node2 = get_node_by_id(masterNode, [str(llist[2])])
                    node3 = get_node_by_id(masterNode, [str(llist[3])])
                    
                    node1.divideOnce()
                    node2.divideOnce()
                    node3.divideOnce()

    for j in range(0, len(extra_list)):
        node1 = extra_list[j]
        node1.children[0].divideOnce()
        node1.children[1].divideOnce()
        node1.children[2].divideOnce()
        node1.children[3].divideOnce()
                    

def swap_directions(dirlistN,sides,sidesN):
#   
#    print 'inside swap_direction'
    myTOL = 2
#    print dirlistN
#    print sides[0].x, sides[0].y, sides[1].x, sides[1].y
#    print sidesN[0].x, sidesN[0].y, sidesN[1].x, sidesN[1].y
    
    if abs(sidesN[0].x - sidesN[1].x)<=1 and abs(sidesN[0].y - sidesN[1].y)<=1:
        if abs(sides[0].x - sidesN[0].x)<=1 and abs(sides[0].y - sidesN[0].y)<=1:
            return swap_edges(dirlistN[1],dirlistN)
        if abs(sides[1].x - sidesN[1].x)<=1 and abs(sides[1].y - sidesN[1].y)<=1:
            return swap_edges(dirlistN[0],dirlistN)    
        
    if abs(sides[0].x - sidesN[0].x)<=myTOL and abs(sides[0].y - sidesN[0].y)<=myTOL:
        return dirlistN[1]
    
    if abs(sides[0].x - sidesN[1].x)<=myTOL and abs(sides[0].y - sidesN[1].y)<=myTOL:
        return dirlistN[0]
    
    if abs(sides[1].x - sidesN[0].x)<=myTOL and abs(sides[1].y - sidesN[0].y)<=myTOL:
        return dirlistN[1]
    
    if abs(sides[1].x - sidesN[1].x)<=myTOL and abs(sides[1].y - sidesN[1].y)<=myTOL:
        return dirlistN[0]
    
#    return swap_edges(dirlistN[0],dirlistN)
def swap_edges(whichEdge, dirlist):
# swapping whichEdge in dirlist with the other edge left in dirlist
# dirlist always has two edges     
    if whichEdge == 'U':
        newEdge = 'D'
        
    if whichEdge == 'D':
         newEdge = 'U'
         
    if whichEdge == 'L':
         newEdge = 'R'
         
    if whichEdge == 'R':
         newEdge = 'L'
         
    if whichEdge == 'LU':
        newEdge = 'RD'
        
    if whichEdge == 'RU':
         newEdge = 'LD'
         
    if whichEdge == 'LD':
         newEdge = 'RU'
         
    if whichEdge == 'RD':
         newEdge = 'LU'
         
    if newEdge == str(dirlist[0]):
        return dirlist[1]
    else:
        return dirlist[0]
    
             
if __name__ == "__main__":
    print "Reading image in..."
#     inputImage = sitk.ReadImage("images/channelsCircles.png");
#     outputImage = sitk.ReadImage("images/channelsCircles.png");
    inputImage = sitk.ReadImage("images/circles.png");
    outputImage = sitk.ReadImage("images/circles.png");
#    inputImage = sitk.ReadImage("images/channelsSharp.png");
#    outputImage = sitk.ReadImage("images/channelsSharp.png");
#    inputImage = sitk.ReadImage((sys.argv[1]));
#    outputImage = sitk.ReadImage((sys.argv[1]));

    nameOutputImage = "images/outCircles_AllConstr.png"
    
    imageSize = inputImage.GetSize();
    print "Image size:", imageSize
 
    # setting the 4 corners coordinates
    p1 = Coordinate(0,0);
    p2 = Coordinate(imageSize[0]-1,0);
    p3 = Coordinate(imageSize[0]-1,imageSize[1]-1);
    p4 = Coordinate(0,imageSize[1]-1);
 
    rect = [p1,p2,p3,p4]
    rootNode = CNode(None,rect,inputImage,outputImage,imageSize)
    tree = CQuadTree(rootNode)
#    
    masterNode = CNode(None,rect,inputImage,outputImage,imageSize)    
#    tree = CQuadTree(masterNode)
#    masterNode = rootNode
#    tree.balance_tree(rootNode,masterNode)
 
    
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print 'No enrichment nodes and hanging nodes in the same element '
        totalNumberOfNodes = newTotalNumberOfNodes
        masterNode = rootNode
        ghost_nodes_enrichment_nodes(tree, rootNode, masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
          
    masterNode = rootNode
      
        
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
       
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print 'Rebalancing tree by multiple passes '
        masterNode = rootNode
        totalNumberOfNodes = newTotalNumberOfNodes
        tree_balance(tree,rootNode,masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
   
      
    masterNode = rootNode
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
  
      
            
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print '3 neighbor rule'
        totalNumberOfNodes = newTotalNumberOfNodes
        masterNode = rootNode
        three_neighbor_rule(tree, rootNode, masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
     
    print 'total number of element nodes', newTotalNumberOfNodes
     
    masterNode = rootNode
    
    llist = []
    tree_list_of_nodes = get_list_of_nodes(tree,masterNode,masterNode,llist)
    
    full_list = stress_concentration_constraint(tree_list_of_nodes, rootNode,outputImage)
    process_list(full_list,rootNode, outputImage)
    
    masterNode = rootNode
       
   
    llist = []
    tree_list_of_nodes = get_list_of_nodes(tree,masterNode,masterNode,llist)
  
  
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print 'No enrichment nodes and hanging nodes in the same element '
        totalNumberOfNodes = newTotalNumberOfNodes
        masterNode = rootNode
        ghost_nodes_enrichment_nodes(tree, rootNode, masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
            
    masterNode = rootNode
          
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
         
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print 'Rebalancing tree by multiple passes '
        masterNode = rootNode
        totalNumberOfNodes = newTotalNumberOfNodes
        tree_balance(tree,rootNode,masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
     
        
        
    masterNode = rootNode
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
    
        
              
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print '3 neighbor rule'
        totalNumberOfNodes = newTotalNumberOfNodes
        masterNode = rootNode
        three_neighbor_rule(tree, rootNode, masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
         
    print 'total number of element nodes', newTotalNumberOfNodes
     
    llist = []
    tree_list = get_list_of_nodes(tree,masterNode,masterNode,llist)
    draw_interface(outputImage, tree_list, masterNode)  
    
    

    
    print 'writing the image out'
 
    
    sitk.WriteImage(outputImage,nameOutputImage);

## Commenting out the solver 
#    [p_reg,p_regCList,lenClist1] = process_list_of_elements(llist,masterNode)
#      
#    [t_reg,t_px] = numbering(p_reg,p_regCList,llist, masterNode)
#      
#          
#    full_vec = numpy.linspace(0,1.0, pow(2,masterNode.MAX_DEPTH)+1)
#      
#    set_nsew(llist,masterNode,full_vec)
#      
#    p_reg = correct_pvec( p_reg, full_vec, lenClist1, llist, p_regCList)
#    # material conductivities
#    k1 = 1
#    k2 = 10
#    # generate Legendre-Gauss nodes and weights:
#    ruleOrder = 4
#    [ui,wi] = lgwt(ruleOrder,-1,1)
#          
#    # get triangular mesh data
#    f = open("multipleinclusions.res", "r")
#    f2 = open("multipleinclusions.ele", "r")
#    [pTri,UTri] = read_p_U(f)
#    tTri = read_corners(f2)
#    f.close()
#    f2.close()
#          
#           
#    UU = myquad(ndim,ndim,k1,k2,ui,wi,p_reg,t_reg,masterNode,llist,inputImage)
#    aa1 = numpy.array([UU])
#    ww1 = numpy.array(aa1[()])
#    UU = ww1[0].item()[:,:]
#      
#    print 'L-2 Norm: ',  computeNorm(p_reg,t_reg,pTri,tTri,ui,wi,k1,k2,UU,UTri,masterNode,llist)
