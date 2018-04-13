
import sys
from operator import itemgetter

class Vertex:
	def __init__(self,coor):
		self.coor = coor
		self.edges = []

	def addEdge(self,v,sp):
		self.edges.append((v,sp))

	def getLast(self):
		return max(self.edges,key=itemgetter(1))[0]

'''
class Edge:
	def __init__(self,vert,sp):
		self.verts = vert
		self.sp = sp
'''

def loadVs():

	Vs = {}
	i = 0 
	for line in sys.stdin:
		print(i)
		i+=1
		l = line.split()
		if l[3]=='0':
			r1 = l[2]
			sp = int(l[6])
			r2 = l[1]
		else:
			r1 = l[1]
			sp = int(l[3])
			r2 = l[2]
		if r1 not in Vs:
			Vs[r1] = Vertex(r1)
		Vs[r1].addEdge(r2,sp)

	return Vs

print(sys.getsizeof(loadVs))