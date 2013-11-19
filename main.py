import sys
import os
import math
import SimpleITK as sitk
from libFcts import *
from globalVars import *
import string
from itertools import izip_longest
import numpy

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
        
        L = set_interval(imageSize[0],self.depth)
        self.i = list(L).index(p1.x)
        self.j = list(L).index(p1.y)
        
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
        
#        if self.children != [None,None,None,None]:# and root.parent != None:
#           self.has_children = True
#        else:
#            self.has_children = False               

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
       # self.prune(rootnode)
       # self.traverse(rootnode)
##    #_______________________________________________________
##    # Sets children of 'node' to None if they do not have any
##    # LEAF nodes.        
##    def prune(self, node):
###        if node.type == Node.LEAF:
###            return 1
##        leafcount = 0
##        removals = []
##        for child in node.children:
##            if child != None:
##                leafcount += self.prune(child)
##                if leafcount == 0:
##                    removals.append(child)
##        for item in removals:
##            n = node.children.index(item)
##            node.children[n] = None        
##        return leafcount
##    #_______________________________________________________
##    # Appends all nodes to a "generic" list, but only LEAF 
##    # nodes are appended to the list of leaves.
##    def traverse(self, node):
##        QuadTree.allnodes.append(node)
##        print QuadTree.maxdepth,node.depth
###        if node.type == Node.LEAF:
###            QuadTree.leaves.append(node)
##        if node.depth > QuadTree.maxdepth:
##            QuadTree.maxdepth = node.depth
##        for child in node.children:
##            if child != None:
##                self.traverse(child) # << recursion      
##    
    def count_nodes(self,root):
        
#    allnodes = allnodes + 4
        allnodes = 0

#        for child in root.children:
#            if child != None:
#                allnodes += self.count_nodes(child)
                
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
            
#    def height(self,root):
#
#        if root.parent == 0:
#            return 0
#        h1 = self.height(root.children[0])
#        h2 = self.height(root.children[1])
## 
#        if (h1 == -1 or h2 == -1 or h1 < h2-1 or h1 > h2+1):
#            return -1
#        if h1>h2:
#            return h1 + 1
#        else:
#            return h2 + 1       
#        

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

    def get_child(node,index):
        # returns child of node, given index of child
        ll = len(index)
        child = node
        for i in range(0,ll):
            child = child.children[int(index[i])]
        return child

    def division_criterionOnce(self, rect, inImage, outImage):
        p1,p2,p3,p4 = self.rect
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid24 = find_mid_point(p2,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
    
        draw_line(self.outImage,cMid12,cMid34)
        draw_line(self.outImage,cMid14,cMid23)
        return True
    
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
            
                L1 = linear_search(self.inImage,p1,p2);
                L2 = linear_search(self.inImage,p2,p3);
                L3 = linear_search(self.inImage,p4,p3);
                L4 = linear_search(self.inImage,p1,p4);
                
                if len(L2)>1 or len(L4) > 1 or len(L1) > 1 or len(L3) > 1:
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
                
                # NW
                if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    draw_line(self.outImage, L1, L4)
                # NE
                if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE):
                    draw_line(self.outImage, L1, L2)
                # SE
                if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    draw_line(self.outImage, L2, L3)
                # SW
                if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    draw_line(self.outImage, L3, L4)
                # vertical
                if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    draw_line(self.outImage, L1, L3)
                # horizontal
                if (l1==1 and l2==0 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
                    draw_line(self.outImage, L4, L2)
                    
                # case 1: interface crossing through L1 and L4
                if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1.x-p2.x) >= 2*MIN_SIZE) :
                #print "case 1"
                    vecCoord1 = case_NW_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L1,L4);

                    if ( vecCoord1[0].x == -1 and (abs(p1.x-p2.x) >= 2*MIN_SIZE) ):

                        draw_line(self.outImage,cMid12,cMid34);
                        draw_line(self.outImage,cMid14,cMid23);               
                        return True
                    elif vecCoord1[0] != -1:
                        self.enrichNodes = vecCoord1
#                        self.mat = 'Fluid'
#                        print self.en[0].x, self.en[0].y, self.en[1].x, self.en[1].y

                # case 2: interface crossing through L1 and L2
                if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1.x-p2.x) >= 2*MIN_SIZE):
                    vecCoord2 = case_NE_polynomial_test(self.inImage,self.outImage,p1,p2,p3,p4,L1,L2);
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
                        
#                        self.mat = 'Fluid'
                                            
                 # case 7: one line crossing through L1 and L4 and one line crossing through L2 and L3
                 # 2-2 non adjacent corners are the same color (diagonally opposed)
                 # the case of 3 consecutive-adjacent materials
                if (l1==0 and l2==0 and l3==0 and l4==0) :
#                    print "case 7"

                    draw_line(self.outImage,cMid12,cMid34);
                    draw_line(self.outImage,cMid14,cMid23);
                    return True
                
                if (l1==1 and l2==1 and l3==0 and l4==1) or (l1==1 and l2==1 and l3==1 and l4==0) or (l1==0 and l2==1 and l3==1 and l4==1) or (l1==1 and l2==0 and l3==1 and l4==1) :

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
  
def it_exists(index,masterNode):      
    llen = len(index)
    child = masterNode
    for i in range(0,llen):
        if child.children[int(index[i])].has_children == False:
            return False
        child = child.children[int(index[i])]
    return True
 
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
            llist.append([int(root.index)]) 
            
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
                    
            if west_has_children == True and east_has_children == True and south_has_children == True and north_has_children == True:
                root.divideOnce()      

        if root.children[0] != None:
            three_neighbor_rule(tree,root.children[0],masterNode)
        if root.children[1] != None:
            three_neighbor_rule(tree,root.children[1],masterNode)
        if root.children[2] != None:
            three_neighbor_rule(tree,root.children[2],masterNode)
        if root.children[3] != None:
            three_neighbor_rule(tree,root.children[3],masterNode)

def stress_concentration_constraint(tree, root, masterNode):

        print len(root.enrichNodes)

        if root.children[0] != None:
            stress_concentration_constraint(tree,root.children[0],masterNode)
        if root.children[1] != None:
            stress_concentration_constraint(tree,root.children[1],masterNode)
        if root.children[2] != None:
            stress_concentration_constraint(tree,root.children[2],masterNode)
        if root.children[3] != None:
            stress_concentration_constraint(tree,root.children[3],masterNode)
  

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

    
               
if __name__ == "__main__":
    print "Reading image in..."
#    inputImage = sitk.ReadImage("images/channels.png");
#    outputImage = sitk.ReadImage("images/channels.png");
    inputImage = sitk.ReadImage("images/circles.png");
    outputImage = sitk.ReadImage("images/circles.png");
#    inputImage = sitk.ReadImage((sys.argv[1]));
#    outputImage = sitk.ReadImage((sys.argv[1]));


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

    
#    node = CNode.get_child(masterNode,'0233')
#    node.printRect()
#    print 'has kids', node.has_children
#    masterNode.children[0].children[2].children[3].children[3].printRect()

    masterNode = rootNode
    totalNumberOfNodes = tree.count_nodes(rootNode)
    newTotalNumberOfNodes = -1
     
    while totalNumberOfNodes != newTotalNumberOfNodes:
        print '3 neighbor rule'
        totalNumberOfNodes = newTotalNumberOfNodes
        masterNode = rootNode
        three_neighbor_rule(tree, rootNode, masterNode)
        newTotalNumberOfNodes = tree.count_nodes(rootNode)
    
    print 'total number of nodes', newTotalNumberOfNodes
#     masterNode = rootNode
#     stress_concentration_constraint(tree,rootNode,masterNode)
    
    masterNode = rootNode
#     node = CNode.get_child(masterNode,'2')
#     print number_of_generations(tree, node, masterNode), node.depth
    llist = []
    tree_list_of_nodes = get_list_of_nodes(tree,rootNode,masterNode,llist)
    print 'length of list of nodes: ',len(tree_list_of_nodes)
#    print CNode.get_child(masterNode,'21213').mat
#    
    print 'writing the image out'

    sitk.WriteImage(outputImage,"outCircles.png");
#    sitk.WriteImage(outputImage,"outChannels.png");
