
import sys
import pickle
from operator import itemgetter
import heapq as hq

class Vertex:
	def __init__(self,name,seq):
		self.name = name
		self.seq = seq
		self.edges = [] # list of edges in the following format: ([read name],[index of the first nucleotide of the overlap])
		self.IN = 0     # only suffix=prefix edges are remembered 
		self.bc = name.split(':')[-1]
		self.chunk = None
	def addEdge(self,v,sp):
		self.edges.append((v,sp))

	def addIn(self):
		self.IN +=1

class Chunk:
	def __init__(self):
		self.A = 0
		self.B = 0
		self.Vs = []
		self.seq = ''
		self.edges = []   # list of tuples ([chunck],[index of the first nucleotide of the overlap])
		self.edgesIN = [] # list of tuples ([chunk that has an edge going to this chunk],[length of the overlap])

	def addVert(self,v):
		self.Vs.append(v)
		v.chunk = self
		

	def size(self):
		return len(self.Vs)

	def getLength(self):
		return len(self.seq)

	# evaluate() should be called after the chunk graph is modified
	# the fuction retrieves the sequence of the chunk and calculates the fold change
	def evaluate(self):
		v = self.Vs[0]
		self.seq += v.seq
		if len(v.edges)!=0:
			e = v.edges[0][1]
			i = 1
			while i<len(self.Vs):
				v = self.Vs[i]
				self.seq += v.seq[-e:]
				try:
					e = v.edges[0][1]
				except IndexError:
					break
				i+=1
		self.A = 0
		self.B = 0
		for v in self.Vs:
			if v.bc == bc[0]:
				self.A += 1
			else:
				self.B += 1
		
		i=0
		while i<len(self.edges):
			e = self.edges[i]
			e.edgesIN.append(self)
			i+=1
		
### Chunk graph construction ###

def traceChunk(v):
	vert = Vs[v]
	if vert.chunk!=None:
		return [vert],True

	if len(vert.edges) == 0:
		return [vert],False

	if len(vert.edges) == 1:
		l,s = traceChunk(vert.edges[0][0])
		return [vert]+l,s

	return [vert],False

def merge(v,C2):
	
	chunk = Chunk()
	C1 = v.chunk
	chunk.edges = C1.edges.copy()
	C1.edges = [chunk]
	C2.edges = [chunk]
	i = C1.Vs.index(v)
	chunk.Vs = C1.Vs[i:]
	for v in chunk.Vs:
		v.chunk = chunk
	C1.Vs = C1.Vs[:i]
	Cs.append(chunk)

def buildChunkFromEdge(edge):
	e = edge[0]
	sp = edge[1]
	origin = edge[1]
	chunk = Chunk()
	v = e[0]
	l,s = traceChunk(v)

	if s:
		end = l.pop()
		if end.chunk.Vs[0] == end:
			chunk.edges.append(end.chunk)

		else:
			merge(end,chunk)

	for i in l:
		chunk.addVert(i)
	if len(chunk.Vs)>0:
		for e in chunk.Vs[-1].edges:
			St.append((e,chunk))
		origin.edges.append(chunk)
		return chunk
	else:
		return False

def buildChunk(vert):
	chunk = Chunk()
	l,s = traceChunk(vert)

	if s:
		end = l.pop()
		if end.chunk.Vs[0] == end:
			chunk.edges.append(end.chunk)

		else:
			merge(end,chunk)

	for i in l:
		chunk.addVert(i)

	for e in chunk.Vs[-1].edges:
		St.append((e,chunk))
	return chunk
	
def readIN(f):
	print('readIN')
	
	Vs={}
	for line in f:
		l = line.split()
		if l[0] == 'VT':
			Vs[l[1]]=Vertex(l[1],l[2])
		elif l[0] == 'ED' and l[-2] != '1':
			if l[3]=='0':
				r1 = l[2]
				sp = int(l[6])
				r2 = l[1]
			else:
				r1 = l[1]
				sp = int(l[3])
				r2 = l[2]
			Vs[r1].addEdge(r2,sp)
			Vs[r2].addIn()

	in0 = []
	single = []
	for  v in Vs:
		if Vs[v].IN == 0:
			if len(Vs[v].edges) == 0:
				single.append(v)
			else:
				in0.append(v)
	for v in single:
		Vs.pop(v)

	print('done')
	return(Vs,in0)

def getBC(V):
	bc = set()
	for r in V:	
		if len(bc) == 2:
			break
		A = r.split(':')[-1]
		bc.add(A)

	return list(bc)

def evaluateAll():
	for c in Cs:
		c.evaluate()

def collapse():
	print('collapse')
	while len(St)>0:
		x = St.pop()
		if type(x) == str:
			Cs.append(buildChunk(x))
		else:
			c = buildChunkFromEdge(x)
			if c:
				Cs.append(c)
	evaluateAll()
	print('done')

### Retrieving sequences ###

def makeLengthHeap():
	heap = []
	i = 0
	for c in Cs:
		try:
			hq.heappush(heap,(c.getLength(),i,c))
		except TypeError:
			i+=1
			hq.heappush(heap,(c.getLength(),i,c))
	return heap

def getLongestSeqs():
	pass



'''f = open('./isolate/merged.asqg','r')
Vs,St = readIN(f)
bc = getBC(Vs)
'''
Cs = []

#collapse()
