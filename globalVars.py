binBnd = [0,110,256]
AREA_INCLUSION = 16.0 # or a 4x4 small inclusion of pixels
PROB = 0.01
MAX_SIZE = 128
MIN_SIZE = 32#6
ALT_MIN_SIZE = 32#3
TOL_LINEARS = 2.0
TOL_QUADRATICS = -2.0 # cancel the cubic approximation by making a negative number
TOL_CUBICS = -2.0 # cancel the cubic approximation by making a negative number

class Coordinate(object):
    def __init__(self,x=-1,y=-1):
        self.x = x
        self.y = y

