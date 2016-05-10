#!/usr/bin/python2.7

import numpy as np
import random
import math

class Vector:
	def __init__(self, comps):
		self.comps=np.array(comps, dtype=float)
		self.n=len(self.comps)

	def norm(self):
		return np.linalg.norm(self.comps)

	def __add__(self, other):
		if isinstance(other, (int, float)):
			return Vector(other+self.comps)
		elif isinstance(other, Vector):
			return Vector(other.comps+self.comps)

	def __radd__(self, other):
		return self.__add__(other)	
	
	def __mul__(self, other):
		if type(other) in (int, float):
			return Vector(other*self.comps)
		elif isinstance(other, Vector):
			return np.inner(self.comps, other.comps)
	
	def __rmul__(self, other):
		return self.__mul__(other)

	def __repr__(self):
		return str(self.comps)

class Instance:
	def __init__(self, x, y):
		self.x=x	#x is a vector from R^n
		self.y=y	#y in {-1; +1}

	def __repr__(self):
		return '(' + str(self.y) +', ' + str(self.x) + ')'

class Set:
	def __init__(self, insts):
		self.insts=insts
		self.m=len(self.insts)

	def getrsubs(self, k):
		return random.sample(xrange(0, self.m), k)
	
	def __repr__(self):
		return '\n'.join(str(it) for it in self.insts)

def getpossubs(S, sub, w):
	ret_list=[]
	for i in sub:
		if S.insts[i].y*(w*S.insts[i].x) < 1:
			ret_list.append(i)
	return ret_list

def check_pres(eps, Rsq, t, dlt, lam):
	E=(21*4*Rsq*math.log(t/dlt))/(lam*t)
	#print E
	if E<eps:
		return True
	else:
		return False

if __name__=='__main__':
	lam=100
	k=1
	eps=50.0
	dlt=0.5
	R=0.0
	
	inst_list=[]
	with open("./indians.dat", "r") as inp:
		for line in inp:
			attr_chlist=line[:-1].split(',')
			attr_list=[float(a) for a in attr_chlist]
			x=Vector(attr_list[:-1])
			if x.norm()>R:
				R=x.norm()
			y=attr_list[len(attr_list)-1]
			inst=Instance(x, y)
			inst_list.append(inst)
	
	dim=x.n
	S=Set(inst_list)
	wt=Vector([0.0]*dim)
	loss=Vector([0.0]*dim)
	W=Vector([0.0]*dim)

	t=1
	Rsq=R**2
	while True:
		At=S.getrsubs(k)
		Atpos=getpossubs(S, At, wt)
		stept=1.0/(lam*t)
		loss*=0.0
		for i in Atpos:
			loss+=S.insts[i].x*S.insts[i].y
		wt=(1-stept*lam)*wt + (stept/k)*loss
		W+=wt
		if check_pres(eps, Rsq, t, dlt, lam):
			break
		t+=1

	print (1.0/t)*W
