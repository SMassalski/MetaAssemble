
import sys
import pickle
from operator import itemgetter
from functools import total_ordering
import heapq as hq



class Vertex:
	def __init__(self,name,size):
		self.name = name
		#self.seq = seq
		self.edges = [] # list of edges in the following format: ([read name],[index of the first nucleotide of the overlap])
		self.IN = 0     # only suffix=prefix edges are remembered 
		self.dups = {}
		self.chunk = None
		self.size = size
		self.merged = False
		

		s = self.name.split(':')
		self.dups[s[-3]]=int(s[-1])

	def addEdge(self,v,sp):
		self.edges.append((v,sp))

	def addIn(self):
		self.IN +=1

	def mergeDuplicates(self,dup):
		s = dup.split(':')
		self.dups[s[-3]]=int(s[-1])

@total_ordering
class Chunk:
	idCount = 0
	def __init__(self):
		self.Vs = []
		self.seq = ''
		self.edges = []   # list of tuples ([chunck],[index of the first nucleotide of the overlap])
		self.edgesIN = [] # list of tuples ([chunk that has an edge going to this chunk],[length of the overlap])
		self.dups = {'6685_04-06-2015':1,'6685_16-06-2015':1}
		self.fc = 0
		self.num = Chunk.idCount
		Chunk.idCount += 1

	def __eq__(self, other):
		return self.size() == other.size()

	def __ne__(self, other):
		return not (self==other)

	def __lt__(self, other):
		return -self.size() < -other.size()
	def addVert(self,v):
		self.Vs.append(v)
		v.chunk = self
		

	def size(self):
		return len(self.Vs)

	def getLength(self):
		return len(self.seq)

	# evaluate() should be called after the chunk graph is modified
	# the fuction retrieves the sequence of the chunk and calculates the fold change
	def evaluate(self,s):
		v = self.Vs[0]
		self.seq += s[v.name]
		if len(v.edges)!=0:
			e = v.edges[0][1]
			i = 1
			while i<len(self.Vs):
				v = self.Vs[i]
				self.seq += s[v.name][e:]
				try:
					e = v.edges[0][1]
				except IndexError:
					break
				i+=1
		for v in self.Vs:
			for d in v.dups:
				self.dups[d]+=v.dups[d]

		self.fc = self.dups['6685_16-06-2015']/self.dups['6685_04-06-2015']

		i=0
		while i<len(self.edges):
			e = self.edges[i]
			e[0].edgesIN.append((self,e[1]))
			i+=1
		
### Chunk graph construction ###

def traceChunk(v,chunk):
	l = []
	while 1:
		if len(l)>1000: ###debug loop
			log = open('trace.log','w')
			for v in l:
				log.write(v.name+'\n')
			print('check log')
			sys.exit()
		vert = Vs[v]
		l.append(vert)
		if not (vert.chunk is None):
			return l,True

		elif len(vert.edges) == 0:
			vert.chunk = chunk
			return l,False

		elif len(vert.edges) == 1:
			vert.chunk = chunk
			v = vert.edges[0][0]
				
		else:
			vert.chunk = chunk
			return l,False

def merge(v,C2,ol):
	
	chunk = Chunk()
	C1 = v.chunk
	chunk.edges = C1.edges.copy()

	i = C1.Vs.index(v)
	chunk.Vs = C1.Vs[i:]
	C1.Vs = C1.Vs[:i]

	C1.edges = [(chunk,C1.Vs[-1].edges[0][1])]
	C2.edges.append((chunk,ol))


	
	for v in chunk.Vs:
		v.chunk = chunk
	
	Cs.append(chunk)

def buildChunkFromEdge(edge):
	e = edge[0]
	origin = edge[1]
	chunk = Chunk()
	v = e[0]
	l,s = traceChunk(v,chunk)
	new_chunk = None
	

	if s: 						### s == True when the last vert already belongs to a chunk
		end = l.pop()
		if end.chunk is chunk:
			i = l.index(end)
			new = l[i:]
			l = l[:i]
			new_chunk = Chunk()
			for i in new:
				new_chunk.addVert(i)
			new_chunk.edges.append((new_chunk,new[-1].edges[0][1]))

		elif end.chunk.Vs[0] == end:		### if its the first vert of the other chunk simply add an edge to the chunk
			try:
				chunk.edges.append((end.chunk,l[-1].edges[0][1]))
			except IndexError:
				if (end.chunk,e[1]) not in origin.edges:
					origin.edges.append((end.chunk,e[1]))
					

		else:
			try:
				merge(end,chunk,l[-1].edges[0][1])
			except IndexError:
				merge(end,origin,e[1])

	for i in l:
		chunk.addVert(i)
	
	if len(chunk.Vs)>0:
		for e in chunk.Vs[-1].edges:
			stack.append((e,chunk))

		if new_chunk is not None:
			chunk.edges.append((new_chunk,l[-1].edges[0][1]))

		else:
			for e in chunk.Vs[-1].edges:
				stack.append((e,chunk))
		origin.edges.append((chunk,e[1]))


		return chunk
	elif new_chunk is not None:
		origin.edges.append((new_chunk,e[1]))
		return new_chunk
	else:
		return False

def buildChunk(vert):
	chunk = Chunk()
	l,s = traceChunk(vert,chunk)
	new_chunk = None

	if s:
		end = l.pop()
		if end.chunk is chunk:
			i = l.index(end)
			new = l[i:]
			l = l[:i]
			new_chunk = Chunk()
			for i in new:
				new_chunk.addVert(i)
			new_chunk.edges.append((new_chunk,new[-1].edges[0][1]))

		elif end.chunk.Vs[0] == end:
			chunk.edges.append((end.chunk,l[-1].edges[0][1]))

		else:
			merge(end,chunk,l[-1].edges[0][1])

	for i in l:
		chunk.addVert(i)

	if new_chunk is not None:
		chunk.edges.append((new_chunk,l[-1].edges[0][1]))
	else:
		for e in chunk.Vs[-1].edges:
			stack.append((e,chunk))
	return chunk
	

### Overlap graph representation ###

def readIN(f):
	print('readIN')
	
	Vs={}
	contained = {}
	edgesIn = {}
	
	for line in f:
		l = line.split()
		if l[0] == 'ED' and l[-2] != '1':
			if l[3]=='0':
				r1 = l[2]
				r1s = int(l[8])
				ol = int(l[8])-int(l[6]) # Overlap length
				r2 = l[1]
				r2s = int(l[5])
			else:
				r1 = l[1]
				r1s = int(l[5])
				ol = int(l[5])-int(l[3])
				r2 = l[2]
				r2s=int(l[8])

			if r1 in contained or r2 in contained:			
				pass
			elif ol == r2s: # if r2 is contained by r1
				if r1 not in Vs:
					Vs[r1] = Vertex(r1,r1s)
				Vs[r1].mergeDuplicates(r2)
				contained[r2] = r1
				try:
					edgesIn[r2].append(r1)
				except KeyError:
					edgesIn[r2] = [r1]
			
			else:
				if r2 not in Vs:
					Vs[r2] = Vertex(r2,r2s)
				if r1 not in Vs:
					Vs[r1] = Vertex(r1,r1s)

				try:
					edgesIn[r2].append(r1)
				except KeyError:
					edgesIn[r2] = [r1]

				Vs[r1].addEdge(r2,ol)
				Vs[r2].addIn()

	print('cleanup')

	for v in contained:
		for vert in edgesIn[v]:
			for e in Vs[vert].edges:
				if v == e[0]:
					Vs[vert].edges.remove(e)
					break
			if v in Vs:
				Vs.pop(v)
				
	in0 = []
	
	for  v in Vs:
		if Vs[v].IN == 0:
			in0.append(v)

	print('done')
	return Vs,in0,edgesIn

def readSeqs(f):
	print('readSeqs')
	f.seek(0)
	Vs={}
	for line in f:
		l = line.split()
		if l[0] == 'VT':
			Vs[l[1]]=l[2]
		elif l[0] == 'ED':
			break
	print('done')
	return Vs

def evaluateAll(f):
	s = readSeqs(f)
	for c in Cs:
		c.evaluate(s)

def collapse():
	print('collapse')
	while len(stack)>0:
		print(len(stack))
		x = stack.pop()
		if type(x) == str:
			Cs.append(buildChunk(x))
		else:
			c = buildChunkFromEdge(x)
			if c:
				Cs.append(c)

	print('done')

### Retrieving sequences ###
'''
def makeLengthHeap():
	heap = []
	for c in Cs:
		hq.heappush(heap,(-c.getLength(),c))
		
	return heap

def getLongestSeqs(heap=None,n=-1):
	if heap is None:
		heap = makeLengthHeap()
		print('heap created')
	if n == -1:
		n = len(heap)
	l = []
	C = []
	for i in range(n):
		try:
			c = hq.heappop(heap)[1]
			s = c.seq
			
			cB = c
			cF = c
			chunks = [c]

			
			while len(cF.edges)>0:
				eO = cF.edges.copy()
				eO.sort(key=lambda t: -t[0].getLength())
				eF = None
				for x in eO:
					if (-x[0].getLength(),x[0]) in heap:
						eF = x
						break
				if eF is None:
					break
				
				cF= eF[0]
				ol = eF[1]
				heap.remove((-cF.getLength(),cF))
				s+=cF.seq[ol:]
				chunks.append(cF)


			while len(cB.edgesIN)>0:
				eI = cB.edgesIN.copy()
				eI.sort(key=lambda t: -t[0].getLength())
				eB = None
				for x in eI:
					if (-x[0].getLength(),x[0]) in heap:
						eB = x
						break
				if eB is None:
					break
				cB = eB[0]
				ol = eB[1]
				heap.remove((-cB.getLength(),cB))
				s = cB.seq[:-ol] + s
				chunks = [cB] + chunks

			l.append((s))
			C.append(chunks)
			hq.heapify(heap)

		except IndexError:
			break
	return l,C

def makeFcHeap(inverse=False):
	heap = []
	if inverse:
		for c in Cs:
			hq.heappush(heap,(-1/c.fc,c))
	else:
		for c in Cs:
			hq.heappush(heap,(-c.fc,c))
	return heap
		
def getMaxFcSeqs(heap=None,n=-1,cut_off=1.5,verbose=False,inverse=False):
	if heap is None:
		heap = makeFcHeap(inverse)
		print('heap created')
	if n == -1:
		n = len(heap)
	l = []

	if inverse:
		for i in range(n):
			try:
				c = hq.heappop(heap)[1]
				s = c.seq
				FC = 1/c.fc
				if FC < cut_off:
					break
				control = c.dups['6685_04-06-2015']
				treated = c.dups['6685_16-06-2015']

				back = len(c.edgesIN)>0
				forward = len(c.edges)>0

				cF = c
				cB = c
				#chunks = [c]

				eO = cF.edges.copy()
				eO.sort(key=lambda t: -1/t[0].fc)
				eI = cB.edgesIN.copy()
				eI.sort(key=lambda t: -1/t[0].fc)
				eF = None
				eB = None
				for x in eO:
					if (-1/x[0].fc,x[0]) in heap:
						eF = x
						break	
				if eF is None:
					forward = False

				for x in eI:
					if (-1/x[0].fc,x[0]) in heap:
						eB = x
						break		
				if eB is None:
					back = False

				while back or forward:
					if back and forward:
						if eF[0].fc < eB[0].fc:
							if cut_off <= (control+eF[0].dups['6685_04-06-2015'])/(treated + eF[0].dups['6685_16-06-2015']):
								cF = eF[0]
								control += cF.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015']
								heap.remove((-1/cF.fc,cF))
								s+=cF.seq[eF[1]:]
								forward = len(cF.edges)>0
								#chunks.append(cF)
							else:
								break

						elif eF[0].fc > eB[0].fc:
							if cut_off <= (control+eB[0].dups['6685_04-06-2015'])/(treated + eB[0].dups['6685_16-06-2015']):
								cB = eB[0]
								control += cB.dups['6685_04-06-2015']
								treated += cB.dups['6685_16-06-2015']
								heap.remove((-1/cB.fc,cB))
								s = cB.seq[:-eB[1]] + s
								back = len(cB.edgesIN)>0
								#chunks = [cB]+chunks
							else:
								break

						else:
							if cut_off <= (control+eF[0].dups['6685_04-06-2015'])/(treated + eF[0].dups['6685_16-06-2015']):
								cF = eF[0]
								cB = eB[0]
								control += cF.dups['6685_04-06-2015'] + cB.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015'] + cB.dups['6685_16-06-2015']
								heap.remove((-1/cF.fc,cF))
								heap.remove((-1/cB.fc,cB))
								s = cB.seq[:-eB[1]] + s + cF.seq[eF[1]:]
								back = len(cB.edgesIN)>0
								forward = len(cF.edges)>0
								#chunks = [cB] + chunks + [cF]
							else:
								break
						FC = control/treated

						eO = cF.edges.copy()
						eO.sort(key=lambda t: -1/t[0].fc)
						eI = cB.edgesIN.copy()
						eI.sort(key=lambda t: -1/t[0].fc)
						eF = None
						eB = None
						for x in eO:
							if (-1/x[0].fc,x[0]) in heap:
								eF = x
								break
							
						if eF is None:
							forward = False
						for x in eI:
							if (-1/x[0].fc,x[0]) in heap:
								eB = x
								break
							
						if eB is None:
							back = False

					elif forward:
						
						c = eF[0]
						if cut_off <= (control+c.dups['6685_04-06-2015'])/(treated + c.dups['6685_16-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = control/treated
							heap.remove((-1/c.fc,c))
							s+=c.seq[eF[1]:]
							forward = len(c.edges)>0
							#chunks = chunks + [cF]
						else:
							break

						eO = c.edges.copy()
						eO.sort(key=lambda t: -1/t[0].fc)
						eF = None
						for x in eO:
							if (-1/x[0].fc,x[0]) in heap:
								eF = x
								break
							
						if eF is None:
							break

					else:
						c = eB[0]
						if cut_off <= (control+c.dups['6685_04-06-2015'])/(treated + c.dups['6685_16-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = control/treated
							heap.remove((-1/c.fc,c))
							s = c.seq[:-eB[1]] + s
							back = len(c.edgesIN)>0
							#chunks = [cB] + chunks
						else:
							break

						eI = c.edgesIN.copy()
						eI.sort(key=lambda t: -1/t[0].fc)
						eB = None
						for x in eI:
							if (-1/x[0].fc,x[0]) in heap:
								eB = x
								break
							
						if eB is None:
							break

				l.append((FC,s))
				hq.heapify(heap)

			except IndexError:
				break
	else:
		for i in range(n):
			try:
				c = hq.heappop(heap)[1]
				s = c.seq
				FC = c.fc
				if FC < cut_off:
					break
				control = c.dups['6685_04-06-2015']
				treated = c.dups['6685_16-06-2015']

				back = len(c.edgesIN)>0
				forward = len(c.edges)>0

				cF = c
				cB = c
				#chunks = [c]

				eO = cF.edges.copy()
				eO.sort(key=lambda t: -t[0].fc)
				eI = cB.edgesIN.copy()
				eI.sort(key=lambda t: -t[0].fc)
				eF = None
				eB = None
				for x in eO:
					if (-x[0].fc,x[0]) in heap:
						eF = x
						break
				if eF is None:
					forward = False

				for x in eI:
					if (-x[0].fc,x[0]) in heap:
						eB = x
						break
					
				if eB is None:
					back = False
				
				while back or forward:
					
					if back and forward:
						if eF[0].fc > eB[0].fc:
							if cut_off <= (treated + eF[0].dups['6685_16-06-2015'])/(control+eF[0].dups['6685_04-06-2015']):
								cF = eF[0]
								control += cF.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015']
								heap.remove((-cF.fc,cF))
								s+=cF.seq[eF[1]:]
								forward = len(cF.edges)>0
								#chunks = chunks + [cF]
							else:
								break

						elif eF[0].fc < (treated + eB[0].dups['6685_16-06-2015'])/(control+eB[0].dups['6685_04-06-2015']):
							if cut_off <= eB[0].fc:
								cB = eB[0]
								control += cB.dups['6685_04-06-2015']
								treated += cB.dups['6685_16-06-2015']
								heap.remove((-cB.fc,cB))
								s = cB.seq[:-eB[1]] + s
								back = len(cB.edgesIN)>0
								#chunks = [cB] + chunks
							else:
								break

						else:
							if cut_off <= (treated + eF[0].dups['6685_16-06-2015'])/(control+eF[0].dups['6685_04-06-2015']):
								cF = eF[0]
								cB = eB[0]
								control += cF.dups['6685_04-06-2015'] + cB.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015'] + cB.dups['6685_16-06-2015']
								heap.remove((-cF.fc,cF))
								heap.remove((-cB.fc,cB))
								s = cB.seq[:-eB[1]] + s + cF.seq[eF[1]:]
								back = len(cB.edgesIN)>0
								forward = len(cF.edges)>0
								#chunks = [cB] + chunks + [cF]
							else:
								break
						FC = treated/control

						eO = cF.edges.copy()
						eO.sort(key=lambda t: -t[0].fc)
						eI = cB.edgesIN.copy()
						eI.sort(key=lambda t: -t[0].fc)
						eF = None
						eB = None
						for x in eO:
							if (-x[0].fc,x[0]) in heap:
								eF = x
								break
								
						if eF is None:
							forward = False
						for x in eI:
							if (-x[0].fc,x[0]) in heap:
								eB = x
								break
							
						if eB is None:
							back = False

					elif forward:
						c = eF[0]
						if cut_off <= (treated + c.dups['6685_16-06-2015'])/(control+c.dups['6685_04-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = treated/control
							heap.remove((-c.fc,c))
							s+=c.seq[eF[1]:]
							forward = len(c.edges)>0
							#chunks = chunks + [cF]
						else:
							break

						eO = c.edges.copy()
						eO.sort(key=lambda t: -t[0].fc)
						eF = None
						for x in eO:
							if (-x[0].fc,x[0]) in heap:
								eF = x
								break
							
						if eF is None:
							break
					else:
						c = eB[0]
						if cut_off <= (treated + c.dups['6685_16-06-2015'])/(control+c.dups['6685_04-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = treated/control
							heap.remove((-c.fc,c))
							s = c.seq[:-eB[1]] + s
							back = len(c.edgesIN)>0
							#chunks = [cB] + chunks
						else:
							break

						eI = c.edgesIN.copy()
						eI.sort(key=lambda t: -t[0].fc)
						eB = None
						for x in eI:
							if (-x[0].fc,x[0]) in heap:
								eB = x
								break
							
						if eB is None:
							break
						
				l.append((FC,s))
				hq.heapify(heap)

			except IndexError:
				break
	return l
'''
def maxFc(L=None,n=-1,cut_off=1.5,inverse=False):
	if L is None:
		L = Cs.copy()
		L.sort(key=lambda c: c.fc,reverse=inverse)
		print('sorted')
	if n == -1:
		n = len(L)
	l = []
	C = []
	d = {}
	for i in range(len(L)):
		d[L[i].num] = i
	print("indexed")
	if inverse:
		for i in range(n):
			try:
				c = L.pop()
				while d[c.num] == -1:
					c = L.pop()
				s = c.seq
				chunks = [c]
				cB = c
				cF = c
				FC = 1/c.fc
				if FC < cut_off:
					break
				control = c.dups['6685_04-06-2015']
				treated = c.dups['6685_16-06-2015']

				back = len(c.edgesIN)>0
				forward = len(c.edges)>0

				if forward:
					eF=max(cF.edges,key=lambda t: d[t[0].num])
					forward = d[eF[0].num]!=-1
				if back:
					eB=max(cB.edgesIN,key=lambda t: d[t[0].num])
					back = d[eB[0].num]!=-1

				while back or forward:
					if back and forward:
						if eF[0].fc < eB[0].fc:
							if cut_off <= (control+eF[0].dups['6685_04-06-2015'])/(treated + eF[0].dups['6685_16-06-2015']):
								cF = eF[0]
								control += cF.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015']
								d[cF.num]=-1
								s+=cF.seq[eF[1]:]
								forward = len(cF.edges)>0
								#chunks.append(cF)
							else:
								break
						elif eF[0].fc > eB[0].fc:
							if cut_off <= (control+eB[0].dups['6685_04-06-2015'])/(treated + eB[0].dups['6685_16-06-2015']):
								cB = eB[0]
								control += cB.dups['6685_04-06-2015']
								treated += cB.dups['6685_16-06-2015']
								d[cB.num]=-1
								s = cB.seq[:-eB[1]] + s
								back = len(cB.edgesIN)>0
								#chunks = [cB]+chunks
							else:
								break

						else:
							if cut_off <= (control+eF[0].dups['6685_04-06-2015'])/(treated + eF[0].dups['6685_16-06-2015']):
								cF = eF[0]
								cB = eB[0]
								control += cF.dups['6685_04-06-2015'] + cB.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015'] + cB.dups['6685_16-06-2015']
								d[cF.num]=-1
								d[cB.num]=-1
								s = cB.seq[:-eB[1]] + s + cF.seq[eF[1]:]
								back = len(cB.edgesIN)>0
								forward = len(cF.edges)>0
								#chunks = [cB] + chunks + [cF]
							else:
								break
						FC = control/treated
						if forward:
							eF=max(cF.edges,key=lambda t: d[t[0].num])
							forward = d[eF[0].num]!=-1
						if back:
							eB=max(cB.edgesIN,key=lambda t: d[t[0].num])
							back = d[eB[0].num]!=-1

					elif forward:
						
						c = eF[0]
						if cut_off <= (control+c.dups['6685_04-06-2015'])/(treated + c.dups['6685_16-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = control/treated
							d[c.num]=-1
							s+=c.seq[eF[1]:]
							forward = len(c.edges)>0
							#chunks = chunks + [cF]
						else:
							break

						if forward:
							eF=max(c.edges,key=lambda t: d[t[0].num])
							forward = d[eF[0].num]!=-1

					else:
						c = eB[0]
						if cut_off <= (control+c.dups['6685_04-06-2015'])/(treated + c.dups['6685_16-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = control/treated
							d[c.num]=-1
							s = c.seq[:-eB[1]] + s
							back = len(c.edgesIN)>0
							#chunks = [cB] + chunks
						else:
							break

						if back:
							eB=max(c.edgesIN,key=lambda t: d[t[0].num])
							back = d[eB[0].num]!=-1

				l.append((FC,s))

			except IndexError:
				break

	else:
		for i in range(n):
			try:
				c = L.pop()
				while d[c.num] == -1:
					c = L.pop()
				s = c.seq
				chunks = [c]
				cB = c
				cF = c
				FC = c.fc
				if FC < cut_off:
					break
				control = c.dups['6685_04-06-2015']
				treated = c.dups['6685_16-06-2015']

				back = len(c.edgesIN)>0
				forward = len(c.edges)>0

				if forward:
					eF=max(cF.edges,key=lambda t: d[t[0].num])
					forward = d[eF[0].num]!=-1
				if back:
					eB=max(cB.edgesIN,key=lambda t: d[t[0].num])
					back = d[eB[0].num]!=-1

				while back or forward:
					if back and forward:
						if eF[0].fc < eB[0].fc:
							if cut_off <= (treated + eF[0].dups['6685_16-06-2015'])/(control+eF[0].dups['6685_04-06-2015']):
								cF = eF[0]
								control += cF.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015']
								d[cF.num]=-1
								s+=cF.seq[eF[1]:]
								forward = len(cF.edges)>0
								#chunks.append(cF)
							else:
								break
						elif eF[0].fc > eB[0].fc:
							if cut_off <= (treated + eB[0].dups['6685_16-06-2015'])/(control+eB[0].dups['6685_04-06-2015']):
								cB = eB[0]
								control += cB.dups['6685_04-06-2015']
								treated += cB.dups['6685_16-06-2015']
								d[cB.num]=-1
								s = cB.seq[:-eB[1]] + s
								back = len(cB.edgesIN)>0
								#chunks = [cB]+chunks
							else:
								break

						else:
							if cut_off <= (treated + eF[0].dups['6685_16-06-2015'])/(control+eF[0].dups['6685_04-06-2015']):
								cF = eF[0]
								cB = eB[0]
								control += cF.dups['6685_04-06-2015'] + cB.dups['6685_04-06-2015']
								treated += cF.dups['6685_16-06-2015'] + cB.dups['6685_16-06-2015']
								d[cF.num]=-1
								d[cB.num]=-1
								s = cB.seq[:-eB[1]] + s + cF.seq[eF[1]:]
								back = len(cB.edgesIN)>0
								forward = len(cF.edges)>0
								#chunks = [cB] + chunks + [cF]
							else:
								break
						FC = treated/control
						if forward:
							eF=max(cF.edges,key=lambda t: d[t[0].num])
							forward = d[eF[0].num]!=-1
						if back:
							eB=max(cB.edgesIN,key=lambda t: d[t[0].num])
							back = d[eB[0].num]!=-1

					elif forward:
						
						c = eF[0]
						if cut_off <= (treated + c.dups['6685_16-06-2015'])/(control+c.dups['6685_04-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = treated/control
							d[c.num]=-1
							s+=c.seq[eF[1]:]
							forward = len(c.edges)>0
							#chunks = chunks + [cF]
						else:
							break

						if forward:
							eF=max(c.edges,key=lambda t: d[t[0].num])
							forward = d[eF[0].num]!=-1

					else:
						c = eB[0]
						if cut_off <= (treated + c.dups['6685_16-06-2015'])/(control+c.dups['6685_04-06-2015']):
							control += c.dups['6685_04-06-2015']
							treated += c.dups['6685_16-06-2015']
							FC = treated/control
							d[c.num]=-1
							s = c.seq[:-eB[1]] + s
							back = len(c.edgesIN)>0
							#chunks = [cB] + chunks
						else:
							break

						if back:
							eB=max(c.edgesIN,key=lambda t: d[t[0].num])
							back = d[eB[0].num]!=-1

				l.append((FC,s))
				

			except IndexError:
				break
					
	return l

def longest(L=None,n=-1,min_length=100):

	if L is None:
		L = Cs.copy()
		L.sort(key=lambda c: c.getLength())
		print('sorted')
	if n == -1:
		n = len(L)
	l = []
	C = []
	d = {}
	for i in range(len(L)):
		d[L[i].num] = i
	print("indexed")
	for i in range(n):
		try:
			c = L.pop()
			while d[c.num] == -1:
				c = L.pop()
			s = c.seq
			chunks = [c]
			cB = c
			cF = c
			while len(cF.edges)>0:
				eF = max(cF.edges,key=lambda t: d[t[0].num])
				
				cF= eF[0]
				ol = eF[1]
				if d[cF.num] == -1:
					break
				d[cF.num] = -1
				s+=cF.seq[ol:]
				chunks += [cF]

			while len(cB.edgesIN)>0:
				eB = max(cB.edgesIN,key=lambda t: d[t[0].num])	
				cB = eB[0]
				ol = eB[1]
				if d[cB.num] == -1:
					break
				d[cB.num] = -1
				s = cB.seq[:-ol] + s
				chunks = [cB] + chunks

			if len(s)>min_length:
				l.append((s))
				C.append(chunks)
			

		except IndexError:
			print(len(L))
			break
	return l,C
### Util ###
def writeOut(outdir):
	print('out')
	l,c = longest()
	i = 0
	out = open(outdir+'long.fa','w')
	for s in l:
		out.write('>seq'+str(i)+' len:'+str(len(s))+'\n')
		out.write(s+'\n')
		i+=1
	out.close()

	l = maxFc()
	i = 0
	out = open(outdir+'posFC.fa','w')
	for s in l:
		out.write('>seq'+str(i)+' FC:'+str(s[0]) +' len:'+str(len(s[1]))+'\n')
		out.write(s[1]+'\n')
		i+=1
	out.close()

	l = maxFc(inverse=True)
	i = 0
	out = open(outdir+'negFC.fa','w')
	for s in l:
		out.write('>seq'+str(i)+' FC:'+str(s[0]) +' len:'+str(len(s[1]))+'\n')
		out.write(s[1]+'\n')
		i+=1
	out.close()

Cs = []
dups = {}

f = open('./nuc/separate/merged.asqg','r')
Vs,stack,ei = readIN(f)
collapse()
del Vs
evaluateAll(f)
f.close()
writeOut('./results/')
