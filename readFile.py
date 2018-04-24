import sys
import re
import numpy

def read_p_U(f):
	line = f.readline()

	nr_of_nodes = -1
	line_iter = -1 
	start_flag = 0
	u_line_iter = -1 
	u_start_flag = 0
	ii = 0
	jj = 0
	for line in iter(f):

		line = line.strip()

		if line.startswith("NUMNP  ="):
			nr_of_nodes = int(line.strip().split()[2])
			p = numpy.zeros((nr_of_nodes,2))
			U = numpy.zeros((nr_of_nodes,1))

		if line.startswith("NODE#       X           Y"):
			line_iter = 0 
			start_flag = 1

		if (line_iter <= nr_of_nodes) and nr_of_nodes>0 and line_iter>0 and start_flag==1:
			x = float(line.strip().split()[1])
			y = float(line.strip().split()[2])
			p[ii,0] = x;
			p[ii,1] = y;
			ii += 1

		if line.startswith("NODE#      TEMPERATURE"):
			u_line_iter = 0
			u_start_flag = 1

		if (u_line_iter <= nr_of_nodes) and nr_of_nodes>0 and u_line_iter>0 and u_start_flag==1:
			val = float(line.strip().split()[1])
			U[jj,0] = val
			jj += 1

		line_iter = line_iter + 1
		u_line_iter += 1
		
	return [p,U]

def read_corners(f2):
	t = []
	line = f2.readline()

	for line in iter(f2):
		if line.startswith("#"):
			break
		line = line.strip()
		c1 = int(line.strip().split()[1])-1
		c2 = int(line.strip().split()[2])-1
		c3 = int(line.strip().split()[3])-1
		t = t + [[c1,c2,c3]]	

	return t

def return_ptU():
	f = open(sys.argv[1], "r")
	f2 = open(sys.argv[2], "r")
	[p,U] = read_p_U(f)

	t = read_corners(f2)

	f.close()
	f2.close()

	for e in range(0,len(t)):
		nodes = t[e]
		nodes = numpy.array(nodes)
		coords = p[nodes,:]

	return [p,t,U]

def main():
	[p,t,U] = return_ptU()
if __name__=="__main__":
	main()
