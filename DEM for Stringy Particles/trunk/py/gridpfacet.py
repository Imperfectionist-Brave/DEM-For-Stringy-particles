# encoding: utf-8
#
# 2015 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
# 2015 © Anna Effeindzourou <anna.effeindzourou@gmail.com>
# 2015 © François Kneib <francois.kneib@gmail.com>
# 2015 © Klaus Thoeni <klaus.thoeni@gmail.com>
"""
Helper functions for creating cylinders, grids and membranes. For more details on this type of elements see [Effeindzourou2016]_, [Effeindzourou2015a]_, [Bourrier2013]_,.

For examples using :yref:`GridConnections<GridConnection>`, see

* :ysrc:`examples/grids/CohesiveGridConnectionSphere.py`
* :ysrc:`examples/grids/GridConnection_Spring.py`
* :ysrc:`examples/grids/Simple_Grid_Falling.py`
* :ysrc:`examples/grids/Simple_GridConnection_Falling.py`

For examples using :yref:`PFacets<PFacet>`, see

* :ysrc:`examples/pfacet/gts-pfacet.py`
* :ysrc:`examples/pfacet/mesh-pfacet.py`
* :ysrc:`examples/pfacet/pfacetcreators.py`

"""
from __future__ import print_function
import random
import math
import numpy as np
from builtins import zip
from builtins import range
import math, random, doctest, geom, numpy
from yade.wrapper import *
try:  # use psyco if available
	import psyco
	psyco.full()
except ImportError:
	pass

from yade import utils
from yade import pack

from yade._utils import createInteraction
from yade.utils import box
from numpy import linspace
from yade.minieigenHP import *
#****************************************
def node(pos,density, radius,mask=1):
	v=4.0/3*math.pi*radius**3
	mass=density*v
	b = Node()
	b.mass=mass
	b.pos=pos
	b.mask=mask
	return b
def thread_Se(head,id1, id2, radius,sphereid=None, wire=False, color=None, highlight=False, material=-1, mask=1, cellDist=None):#head==1 first 2 last
	b = Body()
	st=b.state
	node1 = O._nodes[id1]
	node2 = O._nodes[id2]
	if(head==1):
		if sphereid is None:
			print("sphere is none!!!")
			raise ValueError("The sphere must be a positive number")
		sp=O.bodies[sphereid]
		b.shape = Segment(state=st,radius=radius,node1=node1,node2=node2,sphere=sp,isFirst=True, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	if(head==2):
		if sphereid is None:
			print("sphere is none!!!")
			raise ValueError("The sphere must be a positive number")
		sp=O.bodies[sphereid]
		b.shape = Segment(state=st,radius=radius,node1=node1,node2=node2,sphere=sp,isLast=True, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	if(head!=1 and head!=2):
		b.shape = Segment(state=st,radius=radius,node1=node1,node2=node2, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	b.state.pos=node1.pos
	b.shape.vertices[0] = node1.pos
	b.shape.vertices[1] = node2.pos
	b.mask = mask
	b.mat = O.materials[material]
	return b
def thread(id1,id2,radius2,num,material1,ratio,density=300,color1=None,color2=None):#根据sphere创建连接id1 id2 节点id
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	sph1=O.bodies[id1]
	sph2=O.bodies[id2]
	segt=sph2.state.pos-sph1.state.pos
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	#print("norm",norm)
	#print("horizontal",horizontal)
	#print("initial",initial)
	
	#length=segt.norm()-2*radius1
	#numax=int(length/(2 * radius2))
	n=num
	random_numbers = [random.uniform(1, 360) for _ in range(n-2)] #加入随机数模拟柔索自然弯曲
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	pos1=sph1.state.pos+(sph1.shape.radius+radius2)*norm
	pos2=sph2.state.pos-(sph1.shape.radius+radius2)*norm
	posx=linspace(pos1[0],pos2[0],n)
	posy=linspace(pos1[1],pos2[1],n)
	posz=linspace(pos1[2],pos2[2],n)
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
	nodesIds=[]
	nodesIds2=[]
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O._nodes.append(node([i,j,k],density,radius2)))#确保节点位置每次具有一定随机分布
		#nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1)))#确保节点位置每次具有一定随机分布
	nodesIds2.append(O.bodies.append(Thread_Se(1,nodesIds[0],nodesIds[1],radius2,id1,material=material1)))
	#nodesIds2.append(O.bodies.append(Thread_Se(2,nodesIds[-2],nodesIds[-1],radius2,id2,material=material1)))
	nodesIds_=nodesIds
	nodesIds_=nodesIds_[1:-1]
	if(len(nodesIds_)>1):
		for i,j in zip(nodesIds_[:-1],nodesIds_[1:]):
			nodesIds2.append(O.bodies.append(Thread_Se(0,i,j,radius2,material=material1)))
	nodesIds2.append(O.bodies.append(Thread_Se(2,nodesIds[-2],nodesIds[-1],radius2,id2,material=material1)))
	#O.bodies[nodesIds[0]].shape.isFixed=True
	#O.bodies[nodesIds[-1]].shape.isFixed=True
	#O.bodies.appendClumped(id1,nodesIds[0])
#******************************************
def membrane(center, radius, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1, mask=1):
	
	b = Body()
	b.shape = MembraneNode(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	V = (4. / 3) * math.pi * radius**3
	geomInert = (2. / 5.) * V * radius**2
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=center, dynamic=dynamic, fixed=fixed)
	b.aspherical = False
	b.mask = mask
	return b
#******************************************
#********************************
def VSphere(center, radius, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1, mask=0):#**************************************6-19
	
	b = Body()
	b.shape = Vsphere(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	V = (4. / 3) * math.pi * radius**3
	geomInert = (2. / 5.) * V * radius**2
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=center, dynamic=dynamic, fixed=fixed)
	b.aspherical = False
	b.mask = mask
	return b
def VThread(id1, id2 ):
	i = createInteraction(id1, id2)
	for j in O.materials:
		if(isinstance(j,VsMat)):
			i.phys.k=j.k
			break	
	return i
#********************************
def chainedCylinder(
        begin=Vector3(0, 0, 0), end=Vector3(1., 0., 0.), radius=0.2, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1, mask=1
):
	"""
	Create and connect a chainedCylinder with given parameters. The shape generated by repeted calls of this function is the Minkowski sum of polyline and sphere.

	:param Real radius: radius of sphere in the Minkowski sum.
	:param Vector3 begin: first point positioning the line in the Minkowski sum
	:param Vector3 last: last point positioning the line in the Minkowski sum

	In order to build a correct chain, last point of element of rank N must correspond to first point of element of rank N+1 in the same chain (with some tolerance, since bounding boxes will be used to create connections.

	:return: Body object with the :yref:`ChainedCylinder` :yref:`shape<Body.shape>`.

	.. note:: :yref:`ChainedCylinder` is deprecated and will be removed in the future, use :yref:`GridConnection` instead. See :yref:`yade.gridpfacet.cylinder` and :yref:`yade.gridpfacet.cylinderConnection`.
	"""

	import warnings
	warnings.warn(
	        '\033[1;31mchainedCylinder is deprecated and will be removed in the future, use GridConnection instead. See examples/grids/CohesiveGridConnectionSphere.py.\033[1;0m',
	        category=UserWarning
	)

	segment = end - begin
	b = Body()
	b.shape = ChainedCylinder(radius=radius, length=segment.norm(), color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	b.shape.segment = segment
	V = 2 * (4. / 3) * math.pi * radius**3
	geomInert = (2. / 5.) * V * radius**2 + b.shape.length * b.shape.length * 2 * (4. / 3) * math.pi * radius**3
	b.state = ChainedState()
	b.state.addToChain(O.bodies.append(b))
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=begin, resetState=False, dynamic=dynamic, fixed=fixed)
	b.mask = mask
	b.bound = Aabb(color=[0, 1, 0])
	b.state.ori.setFromTwoVectors(Vector3(0., 0., 1.), segment)
	if (end == begin):
		b.state.ori = Quaternion((1, 0, 0), 0)
	return b


def gridNode(center, radius, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1):
	"""
	Create a :yref:`GridNode` which is needed to set up :yref:`GridConnections<GridConnection>`.

	See documentation of :yref:`yade.utils.sphere` for meaning of parameters.

	:return: Body object with the :yref:`gridNode` :yref:`shape<Body.shape>`.
	"""
	b = Body()
	b.shape = GridNode(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	#V=(4./3)*math.pi*radius**3	# will be overwritten by the connection
	V = 0.
	geomInert = (2. / 5.) * V * radius**2  # will be overwritten by the connection
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=center, dynamic=dynamic, fixed=fixed)
	b.aspherical = False
	b.bounded = False
	b.mask = 0  # avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
	return b


def gridConnection(id1, id2, radius, wire=False, color=None, highlight=False, material=-1, mask=1, cellDist=None):
	"""
	Create a :yref:`GridConnection` by connecting two :yref:`GridNodes<GridNode>`.

	:param id1,id2: the two :yref:`GridNodes<GridNode>` forming the cylinder. 
	:param float radius: radius of the cylinder. Note that the radius needs to be the same as the one for the :yref:`GridNodes<GridNode>`.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions for gridConnections are not yet implemented! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`GridConnection` :yref:`shape<Body.shape>`.

	.. note:: The material of the :yref:`GridNodes<GridNode>` will be used to set the constitutive behaviour of the internal connection, i.e., the constitutive behaviour of the cylinder. The material of the :yref:`GridConnection` is used for interactions with other (external) bodies.
	"""
	b = Body()
	b.shape = GridConnection(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	sph1 = O.bodies[id1]
	sph2 = O.bodies[id2]
	i = createInteraction(id1, id2)
	nodeMat = sph1.material
	b.shape.node1 = sph1
	b.shape.node2 = sph2
	sph1.shape.addConnection(b)
	sph2.shape.addConnection(b)
	if (O.periodic):
		if (cellDist != None):
			i.cellDist = cellDist
		segt = sph2.state.pos + O.cell.hSize * i.cellDist - sph1.state.pos
	else:
		segt = sph2.state.pos - sph1.state.pos
	L = segt.norm()
	V = 0.5 * L * math.pi * radius**2
	geomInert = (2. / 5.) * V * radius**2
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=sph1.state.pos, dynamic=False, fixed=True)
	sph1.state.mass = sph1.state.mass + V * nodeMat.density
	sph2.state.mass = sph2.state.mass + V * nodeMat.density
	for k in [0, 1, 2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert * nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert * nodeMat.density
	b.aspherical = False
	if O.periodic:
		i.phys.unp = -(sph2.state.pos + O.cell.hSize * i.cellDist - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
		b.shape.periodic = True
		b.shape.cellDist = i.cellDist
	else:
		i.phys.unp = -(sph2.state.pos - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
	i.geom.connectionBody = b
	print(id1,id2,i.phys.kn)
	I = math.pi * (2. * radius)**4 / 64.
	E = nodeMat.young
	i.phys.kn = E * math.pi * (radius**2) / L
	i.phys.kr = E * I / L
	i.phys.ks = 12. * E * I / (L**3)
	G = E / (2. * (1 + nodeMat.poisson))
	i.phys.ktw = 2. * I * G / L
	b.mask = mask
	return b
#*********************************12-26
def beads(id1,id2,radius1,radius2,num,material1,ratio):#根据颗粒创建连接
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	sph1=O.bodies[id1]
	sph2=O.bodies[id2]
	segt=sph2.state.pos-sph1.state.pos
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	#print("norm",norm)
	#print("horizontal",horizontal)
	#print("initial",initial)
	
	#length=segt.norm()-2*radius1
	#numax=int(length/(2 * radius2))
	n=num
	random_numbers = [random.uniform(1, 360) for _ in range(n-2)] #加入随机数模拟柔索自然弯曲
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	pos1=sph1.state.pos+(radius1+radius2)*norm
	pos2=sph2.state.pos-(radius1+radius2)*norm
	posx=linspace(pos1[0],pos2[0],n)
	posy=linspace(pos1[1],pos2[1],n)
	posz=linspace(pos1[2],pos2[2],n)
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
	nodesIds=[]
	nodesIds2=[]
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1)))#确保节点位置每次具有一定随机分布
	for i,j in zip(nodesIds[:-1],nodesIds[1:]):
		nodesIds2.append(O.bodies.append(lineConnection(i,j,radius2,material=material1)))
	#O.bodies[nodesIds[0]].shape.isFixed=True
	#O.bodies[nodesIds[-1]].shape.isFixed=True
	#O.bodies.appendClumped(id1,nodesIds[0])
	Sph1 = O.bodies[nodesIds[0]]
	Sph2 = O.bodies[nodesIds[-1]]
	Sph1.shape.sphid=id1#需要在创建交互之前设定id
	Sph2.shape.sphid=id2
	createInteraction(id1,nodesIds[0])
	createInteraction(id2,nodesIds[-1])
	
	

def lineNode(center, radius, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1,mask=1):#lineNode
	"""
	Create a :yref:`GridNode` which is needed to set up :yref:`GridConnections<GridConnection>`.

	See documentation of :yref:`yade.utils.sphere` for meaning of parameters.

	:return: Body object with the :yref:`gridNode` :yref:`shape<Body.shape>`.
	"""
	b = Body()
	b.shape = LineNode(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	#V=(4./3)*math.pi*radius**3	# will be overwritten by the connection
	V = (4./3)*math.pi*radius**3
	geomInert = (2. / 5.) * V * radius**2  # will be overwritten by the connection
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=center, dynamic=dynamic, fixed=fixed)
	b.aspherical = False
	#b.bounded = False
	#b.mask = 0  # avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
	b.mask = mask
	return b

def lineConnection(T,id1, id2, radius, head=None,wire=False, color=None, highlight=False, material=-1, mask=1, cellDist=None):
	"""
	Create a :yref:`GridConnection` by connecting two :yref:`GridNodes<GridNode>`.

	:param id1,id2: the two :yref:`GridNodes<GridNode>` forming the cylinder. 
	:param float radius: radius of the cylinder. Note that the radius needs to be the same as the one for the :yref:`GridNodes<GridNode>`.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions for gridConnections are not yet implemented! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`GridConnection` :yref:`shape<Body.shape>`.

	.. note:: The material of the :yref:`GridNodes<GridNode>` will be used to set the constitutive behaviour of the internal connection, i.e., the constitutive behaviour of the cylinder. The material of the :yref:`GridConnection` is used for interactions with other (external) bodies.
	"""
	b = Body()
	
	b.shape = LineConnection(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	sph1 = O.bodies[id1]
	sph2 = O.bodies[id2]
	
	
	sph1.shape.id2=id2#需要在创建交互之前指定id
	sph2.shape.id1=id1
	#sph1.shape.realIntrs+=1
	#sph2.shape.realIntrs+=1
	if (head != None):
		if(head=='First'):
			sph1.shape.isFirst=True
			#print(head)
		else:
			sph2.shape.isLast=True
			#print(head)
	else:
		sph1.shape.isMiddle=True
		#print(head)
	i = createInteraction(id1, id2)
	i.geom._Thread=T#geom存储thread
	nodeMat = sph1.material
	b.shape.node1 = sph1
	b.shape.node2 = sph2
	sph1.shape.addConnection(b)#添加节点
	sph2.shape.addConnection(b)
	if (O.periodic):
		if (cellDist != None):
			i.cellDist = cellDist
		segt = sph2.state.pos + O.cell.hSize * i.cellDist - sph1.state.pos
	else:
		segt = sph2.state.pos - sph1.state.pos
	L = segt.norm()
	V = 0.5 * L * math.pi * radius**2
	geomInert = (2. / 5.) * V * radius**2
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=sph1.state.pos, dynamic=False, fixed=True)
	sph1.state.mass = sph1.state.mass + V * nodeMat.density
	sph2.state.mass = sph2.state.mass + V * nodeMat.density
	for k in [0, 1, 2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert * nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert * nodeMat.density
	b.aspherical = False
	i.geom.connectionBody = b
	I = math.pi * (2. * radius)**4 / 64.
	E = nodeMat.young
	i.phys.kn = E * math.pi * (radius**2) / L
	i.phys.kr = E * I / L
	i.phys.ks = 12. * E * I / (L**3)
	G = E / (2. * (1 + nodeMat.poisson))
	i.phys.ktw = 2. * I * G / L
	b.mask = 1
	#b.mask = mask
	return b
def _lineConnection(id1, id2, radius, head=None,wire=False, color=None, highlight=False, material=-1, mask=1, cellDist=None):
	"""
	Create a :yref:`GridConnection` by connecting two :yref:`GridNodes<GridNode>`.

	:param id1,id2: the two :yref:`GridNodes<GridNode>` forming the cylinder. 
	:param float radius: radius of the cylinder. Note that the radius needs to be the same as the one for the :yref:`GridNodes<GridNode>`.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions for gridConnections are not yet implemented! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`GridConnection` :yref:`shape<Body.shape>`.

	.. note:: The material of the :yref:`GridNodes<GridNode>` will be used to set the constitutive behaviour of the internal connection, i.e., the constitutive behaviour of the cylinder. The material of the :yref:`GridConnection` is used for interactions with other (external) bodies.
	"""
	b = Body()
	
	b.shape = LineConnection(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
	sph1 = O.bodies[id1]
	sph2 = O.bodies[id2]
	
	
	sph1.shape.id2=id2#需要在创建交互之前指定id
	sph2.shape.id1=id1
	#sph1.shape.realIntrs+=1
	#sph2.shape.realIntrs+=1
	i = createInteraction(id1, id2)#创建交互
	if (head != None):
		if(head=='First'):
			sph1.shape.isFirst=True
			i.geom.isFirst=True
			#print(head)
		else:
			sph2.shape.isLast=True
			i.geom.isLast=True
			#print(head)
	else:
		sph1.shape.isMiddle=True
		#print(head)
	
	nodeMat = sph1.material
	b.shape.node1 = sph1
	b.shape.node2 = sph2
	sph1.shape.addConnection(b)#添加节点
	sph2.shape.addConnection(b)
	if (O.periodic):
		if (cellDist != None):
			i.cellDist = cellDist
		segt = sph2.state.pos + O.cell.hSize * i.cellDist - sph1.state.pos
	else:
		segt = sph2.state.pos - sph1.state.pos
	L = segt.norm()
	V = 0.5 * L * math.pi * radius**2
	geomInert = (2. / 5.) * V * radius**2
	utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=sph1.state.pos, dynamic=False, fixed=True)
	sph1.state.mass = sph1.state.mass + V * nodeMat.density
	sph2.state.mass = sph2.state.mass + V * nodeMat.density#增加圆柱质量
	for k in [0, 1, 2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert * nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert * nodeMat.density
	b.aspherical = False
	i.geom.connectionBody = b
	I = math.pi * (2. * radius)**4 / 64.
	E = nodeMat.young
	E_s = nodeMat.young_s
	i.phys.k = E_s * math.pi * (radius**2) / L#弹性系数k  young_s 15e6
	i.phys.kn = E * math.pi * (radius**2) / L
	i.phys.kr = E * I / L
	#i.phys.ks = 12. * E * I / (L**3)
	G = E / (2. * (1 + nodeMat.poisson))
	i.phys.ktw = 2. * I * G / L
	b.mask = 1
	#b.mask = mask   yong_s
	return [b,i]#6-11
def Beadsbynode(id1,id2,sphid1,sphid2,radius2,num,material1,ratio,color1=None,color2=None):#传入连接颗粒 前后端节点 根据节点创建连接id1 id2 节点id 6-6 7-1
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	
	snode1=O.bodies[id1]
	snode2=O.bodies[id2]
	
	
	
	segt=snode2.state.pos-snode1.state.pos
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	#print("norm",norm)
	#print("horizontal",horizontal)
	#print("initial",initial)
	
	#length=segt.norm()-2*radius1
	#numax=int(length/(2 * radius2))
	n=num
	random_numbers = [random.uniform(1, 360) for _ in range(n)] #加入随机数模拟柔索自然变形
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	pos1=snode1.state.pos
	pos2=snode2.state.pos
	posx=linspace(pos1[0],pos2[0],n+2)
	posy=linspace(pos1[1],pos2[1],n+2)
	posz=linspace(pos1[2],pos2[2],n+2)
	posx=posx[1:-1]
	posy=posy[1:-1]
	posz=posz[1:-1]
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
	nodesIds=[]
	nodesIds2=[]#line段id
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1,color=color1 if color1 else utils.randomColor())))#确保节点位置每次具有一定随机分布
	
	addsphids=nodesIds#加入连接节点
	for i in addsphids:
		O.bodies[i].shape.sphid1=sphid1#存储颗粒id
		O.bodies[i].shape.sphid2=sphid2
		
	nodesIds.insert(0,id1)#加入前节点
	nodesIds.append(id2)#加入后节点
	firstline=[]
	lastline=[]
	middleline=[]
	firstline=_lineConnection(nodesIds[0],nodesIds[1],radius2,head='First',material=material1,color=color2 if color2 else utils.randomColor()) #return [b i]
	nodesIds2.append(O.bodies.append(firstline[0]))#lineid1
	nodesIds_=nodesIds#节点id
	nodesIds_=nodesIds_[1:-1]
	if(len(nodesIds_)>1):
		for i,j in zip(nodesIds_[:-1],nodesIds_[1:]):
			middleline.append(_lineConnection(i,j,radius2,material=material1,color=color2 if color2 else utils.randomColor()))
			nodesIds2.append(O.bodies.append(middleline[-1][0]))
	if(len(nodesIds_)>0):
		lastline=_lineConnection(nodesIds[-2],nodesIds[-1],radius2,head='Last',material=material1,color=color2 if color2 else utils.randomColor())
		nodesIds2.append(O.bodies.append(lastline[0]))
		
	T=Thread()
	ids=[]
	ids.append(id1)
	ids.append(id2)
	T.addSphereId(ids)
	T.addSphereId2([sphid1,sphid2])
	T.addNodeId(nodesIds)
	T.addLineId(nodesIds2)
	T.geom1=firstline[1].geom
	T.geom2=lastline[1].geom
	T.phy1=firstline[1].phys
	T.phy2=lastline[1].phys
	T.radius=radius2
	#geoms
	geoms=[]
	geoms.append(firstline[1].geom)
	for i in middleline:
		geoms.append(i[1].geom)
	geoms.append(lastline[1].geom)
	T.addGeoms(geoms)
	#phys
	phys=[]
	phys.append(firstline[1].phys)
	for i in middleline:
		phys.append(i[1].phys)
	phys.append(lastline[1].phys)
	T.addPhys(phys)
	
	return T
def aabb_10Walls(strain,extrema=None,divide=0.75,**kw):#边界墙生成函数
    walls = []
    #centers=[]
    #extents=[]
    if not extrema:
    	extrema = aabbExtrema()
    #if not thickness: thickness=(extrema[1][0]-extrema[0][0])/10.
    #边界尺寸
    x=extrema[1][0]-extrema[0][0] #x
    _x=x*strain#最大移动量
    y=extrema[1][1]-extrema[0][1] #y
    z=extrema[1][2]-extrema[0][2] #z
    mi,ma = extrema 
    #mia=[ma[0],ma[1],ma[2]-divide*z] #中内点 mam
    #mai=[mi[0],mi[1],mi[2]+(1-divide)*z]#中外点 mim
    mam=[ma[0],ma[1],ma[2]-divide*z]
    mim=[mi[0],mi[1],mi[2]+(1-divide)*z]
    center_bottombox = [(mi[i] + mam[i]) / 2. for i in range(3)]#下盒中心
    center_upbox = [(mim[i] + ma[i]) / 2. for i in range(3)]#上盒中心
    center_bound = [(mam[i] + mim[i]) / 2. for i in range(3)]#分界面中心
    center_middle = [(mi[i] + ma[i]) / 2. for i in range(3)]#大盒子中心
    #添加box
    center=[center_bottombox[0]-0.5*x,center_bottombox[1],center_bottombox[2]]
    extents = [.5*(mam[i] - mi[i]) for i in range(3)]
    extents[0]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#0
    walls[-1].shape.wire = True
    
    center=[center_bottombox[0]+0.5*x,center_bottombox[1],center_bottombox[2]]
    extents = [.5*(mam[i] - mi[i]) for i in range(3)]
    extents[0]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#1
    walls[-1].shape.wire = True
    
    center=[center_upbox[0]-0.5*x,center_upbox[1],center_upbox[2]]
    extents = [.5*(ma[i] - mim[i]) for i in range(3)]
    extents[0]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#2
    walls[-1].shape.wire = True
    
    center=[center_upbox[0]+0.5*x,center_upbox[1],center_upbox[2]]
    extents = [.5*(ma[i] - mim[i]) for i in range(3)]
    extents[0]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#3
    walls[-1].shape.wire = True
    
    center=[center_bound[0]-0.5*x-0.5*_x,center_bound[1],center_bound[2]]
    extents[0]= 0.5*_x
    extents[1]= 0.5*(ma[1] - mim[1])
    extents[2]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#4
    walls[-1].shape.wire = True
    
    center=[center_bound[0]+0.5*x+0.5*_x,center_bound[1],center_bound[2]]
    extents[0]= 0.5*_x
    extents[1]= 0.5*(ma[1] - mim[1])
    extents[2]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#5
    walls[-1].shape.wire = True
    
    center=[center_bottombox[0],center_bottombox[1],center_bottombox[2]-0.5*(1-divide)*z]
    extents = [.5*(ma[i] - mim[i]) for i in range(3)]
    extents[2]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#6
    walls[-1].shape.wire = True
    
    center=[center_upbox[0],center_upbox[1],center_upbox[2]+0.5*divide*z]
    extents = [.5*(ma[i] - mim[i]) for i in range(3)]
    extents[2]=0
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#7
    walls[-1].shape.wire = True
    
    center=[center_middle[0],center_middle[1]-0.5*y,center_middle[2]]
    extents[0]=0.5*(ma[0]-mi[0])+_x
    extents[1]=0
    extents[2]=0.5*z
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#8 
    walls[-1].shape.wire = True
    
    center=[center_middle[0],center_middle[1]+0.5*y,center_middle[2]]
    extents[0]=0.5*(ma[0]-mi[0])+_x
    extents[1]=0
    extents[2]=0.5*z
    walls.append(box(center=center, extents=extents, fixed=True, **kw))#9
    walls[-1].shape.wire = True
   
    return walls
def beadsbysphere(sphid1,sphid2,radius2,num,material1,ratio,color1=None,color2=None):#根据sphere创建连接节点 6-6
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	
	
	sphere1=O.bodies[sphid1]#sphere
	sphere2=O.bodies[sphid2]
	R=sphere1.shape.radius
	
	snorm=sphere2.state.pos-sphere1.state.pos
	_snorm=snorm/snorm.norm()
	
	sph1=sphere1.state.pos+(R+radius2)*_snorm
	sph2=sphere2.state.pos-(R+radius2)*_snorm
	
	
	segt=sph2-sph1
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	#print("norm",norm)
	#print("horizontal",horizontal)
	#print("initial",initial)
	
	#length=segt.norm()-2*radius1
	#numax=int(length/(2 * radius2))
	n=num-2#n=num
	random_numbers = [random.uniform(1, 360) for _ in range(n)] #加入随机数模拟柔索自然变形
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	pos1=sph1
	pos2=sph2
	posx=linspace(pos1[0],pos2[0],n+2)
	posy=linspace(pos1[1],pos2[1],n+2)
	posz=linspace(pos1[2],pos2[2],n+2)
	posx=posx[1:-1]
	posy=posy[1:-1]
	posz=posz[1:-1]
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
		
	nodesIds=[]
	nodesIds2=[]
	nodesIds.append(O.bodies.append(lineNode(pos1,radius2,wire=False,fixed=False,material=material1,color=color1 if color1 else utils.randomColor())))
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1,color=color1 if color1 else utils.randomColor())))#确保节点位置每次具有一定随机分布
	nodesIds.append(O.bodies.append(lineNode(pos2,radius2,wire=False,fixed=False,material=material1,color=color1 if color1 else utils.randomColor())))
	
	
	T=Thread()
	ids=[]
	ids.append(sphid1)
	ids.append(sphid2)
	T.addSphereId(ids)
	sts=[]
	sts.append(sphere1.state)
	sts.append(sphere2.state)
	#print(sts[0],sts[1])
	T.addSphereState(sts)
	nsts=[]
	for i in nodesIds:
		nsts.append(O.bodies[i].state)
	T.addNodeState(nsts)
	T.initial()
	
	nodesIds2.append(O.bodies.append(lineConnection(T,nodesIds[0],nodesIds[1],radius2,head='First',material=material1,color=color2 if color1 else utils.randomColor())))
	nodesIds_=nodesIds
	nodesIds_=nodesIds_[1:-1]
	if(len(nodesIds_)>1):
		for i,j in zip(nodesIds_[:-1],nodesIds_[1:]):
			nodesIds2.append(O.bodies.append(lineConnection(T,i,j,radius2,material=material1,color=color2 if color1 else utils.randomColor())))
	if(len(nodesIds_)>0):
		nodesIds2.append(O.bodies.append(lineConnection(T,nodesIds[-2],nodesIds[-1],radius2,head='Last',material=material1,color=color2 if color1 else utils.randomColor())))
	

	return T
def beadsbynode(sphid1,sphid2,id1,id2,radius2,num,material1,ratio,color1=None,color2=None):#根据节点创建连接id1 id2 节点id 6-6
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	
	sph1=O.bodies[id1]
	sph2=O.bodies[id2]
	sphere1=O.bodies[sphid1]#sphere
	sphere2=O.bodies[sphid2]
	
	
	segt=sph2.state.pos-sph1.state.pos
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	#print("norm",norm)
	#print("horizontal",horizontal)
	#print("initial",initial)
	
	#length=segt.norm()-2*radius1
	#numax=int(length/(2 * radius2))
	n=num
	random_numbers = [random.uniform(1, 360) for _ in range(n)] #加入随机数模拟柔索自然变形
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	pos1=sph1.state.pos
	pos2=sph2.state.pos
	posx=linspace(pos1[0],pos2[0],n+2)
	posy=linspace(pos1[1],pos2[1],n+2)
	posz=linspace(pos1[2],pos2[2],n+2)
	posx=posx[1:-1]
	posy=posy[1:-1]
	posz=posz[1:-1]
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
	nodesIds=[]
	nodesIds2=[]
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1,color=color1 if color1 else utils.randomColor())))#确保节点位置每次具有一定随机分布
	nodesIds.insert(0,id1)#加入前节点
	nodesIds.append(id2)#加入后节点
	
	T=Thread()
	ids=[]
	ids.append(sphid1)
	ids.append(sphid2)
	T.addSphereId(ids)
	sts=[]
	sts.append(sphere1.state)
	sts.append(sphere2.state)
	print(sts[0],sts[1])
	T.addSphereState(sts)
	nsts=[]
	for i in nodesIds:
		nsts.append(O.bodies[i].state)
	T.addNodeState(nsts)
	T.initial()
	
	nodesIds2.append(O.bodies.append(lineConnection(nodesIds[0],nodesIds[1],radius2,head='First',material=material1,color=color2 if color1 else utils.randomColor())))
	nodesIds_=nodesIds
	nodesIds_=nodesIds_[1:-1]
	if(len(nodesIds_)>1):
		for i,j in zip(nodesIds_[:-1],nodesIds_[1:]):
			nodesIds2.append(O.bodies.append(lineConnection(i,j,radius2,material=material1,color=color2 if color1 else utils.randomColor())))
	nodesIds2.append(O.bodies.append(lineConnection(nodesIds[-2],nodesIds[-1],radius2,head='Last',material=material1,color=color2 if color1 else utils.randomColor())))

	return T
	
def beadsalone(pos1,pos2,radius2,num,material1,ratio):#不加入颗粒的柔索模型
	#*******
	def rotation(q,p,theata,scanl):
		RotaTheta=0.5*theata/180*math.pi
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(scanl*q2[0],scanl*q2[1],scanl*q2[2])
		return v
	#*******
	segt=pos2-pos1
	norm=segt/segt.norm()#方向单位向量 1----->2
	vertical=(0,0,1)#norm为z向量时
	Horizontal=(1,0,0)
	if (norm==(0,0,1)or norm==(0,0,-1)):
		V=Horizontal
	else:
		V=vertical
	horizontal=np.cross(norm,V)#水平向量
	initial=np.cross(horizontal,norm)
	norm_I = np.linalg.norm(initial)
	initial=initial/norm_I
	if initial[2]>0:
		initial=-initial#初始未旋转向向量，并指定向下为初始方向
	
	n=num#节点数量
	random_numbers = [random.uniform(1, 360) for _ in range(n-2)] #加入随机数模拟柔索自然弯曲
	#print(random_numbers)
	Offset=[]
	for i in random_numbers:
		Offset.append(rotation(norm,initial,i,ratio*radius2))  #仅适用于沿x轴生成颗粒
	'''for i in random_numbers:
		Offset.append([0,-radius2*math.sin(i/180*math.pi),-radius2*math.cos(i/180*math.pi)])  #仅适用于沿x轴生成颗粒'''
	#print(random.randint(1,360))
	Offset.insert(0,[0,0,0])
	Offset.append([0,0,0])
	#pos_1=sph1.state.pos+(radius1+radius2)*norm
	#pos_2=sph2.state.pos-(radius1+radius2)*norm
	posx=linspace(pos1[0],pos2[0],n)
	posy=linspace(pos1[1],pos2[1],n)
	posz=linspace(pos1[2],pos2[2],n)
	#position=linspace(pos1,pos2,n)
	
	for i in range(n):
		posx[i]=posx[i]+Offset[i][0]
		posy[i]=posy[i]+Offset[i][1]
		posz[i]=posz[i]+Offset[i][2]
	nodesIds=[]
	nodesIds2=[]
	for i,j,k in zip(posx,posy,posz):
		nodesIds.append(O.bodies.append(lineNode([i,j,k],radius2,wire=False,fixed=False,material=material1)))#确保节点位置每次具有一定随机分布
	for i,j in zip(nodesIds[:-1],nodesIds[1:]):
		nodesIds2.append(O.bodies.append(lineConnection(i,j,radius2,material=material1)))
	#O.bodies[nodesIds[0]].shape.isFixed=True
	#O.bodies[nodesIds[-1]].shape.isFixed=True
#*********************************

#TODO: find a better way of handling the Id lists for checking duplicated gridNodes or gridConnections with the same coordinates etc. It would be better to handle this globally, maybe implement something like O.bodies.getGridNodes
def cylinder(
        begin=Vector3(0, 0, 0),
        end=Vector3(1., 0., 0.),
        radius=0.2,
        nodesIds=[],
        cylIds=[],
        dynamic=None,
        fixed=False,
        wire=False,
        color=None,
        highlight=False,
        intMaterial=-1,
        extMaterial=-1,
        mask=1
):
	"""
	Create a cylinder with given parameters. The shape corresponds to the Minkowski sum of line-segment and sphere, hence, the cylinder has rounded vertices. The cylinder (:yref:`GridConnection<GridConnection>`) and its corresponding nodes (yref:`GridNodes<GridNode>`) are automatically added to the simulation. The lists with nodes and cylinder ids will be updated automatically.

	:param Vector3 begin: first point of the Minkowski sum in the global coordinate system.
	:param Vector3 end: last point of the Minkowski sum in the global coordinate system.
	:param Real radius: radius of sphere in the Minkowski sum.
	:param list nodesIds: list with ids of already existing :yref:`GridNodes<GridNode>`. New ids will be added.
	:param list cylIds: list with ids of already existing :yref:`GridConnections<GridConnection>`. New id will be added.
	:param intMaterial: :yref:`Body.material` used to create the interaction physics between the two GridNodes
	:param extMaterial: :yref:`Body.material` used to create the interaction physics between the Cylinder (GridConnection) and other bodies (e.g., spheres interaction with the cylinder)
	
	See :yref:`yade.utils.sphere`'s documentation for meaning of other parameters.

	"""
	id1 = O.bodies.append(gridNode(begin, radius, dynamic=dynamic, fixed=fixed, wire=wire, color=color, highlight=highlight, material=intMaterial))
	nodesIds.append(id1)
	id2 = O.bodies.append(gridNode(end, radius, dynamic=dynamic, fixed=fixed, wire=wire, color=color, highlight=highlight, material=intMaterial))
	nodesIds.append(id2)
	cylIds.append(
	        O.bodies.append(
	                gridConnection(id1, id2, radius=radius, wire=wire, color=color, highlight=highlight, material=extMaterial, mask=mask, cellDist=None)
	        )
	)


def cylinderConnection(
        vertices,
        radius=0.2,
        nodesIds=[],
        cylIds=[],
        dynamic=None,
        fixed=False,
        wire=False,
        color=None,
        highlight=False,
        intMaterial=-1,
        extMaterial=-1,
        mask=1
):
	"""
	Create a chain of cylinders with given parameters. The cylinders (:yref:`GridConnection<GridConnection>`) and its corresponding nodes (yref:`GridNodes<GridNode>`) are automatically added to the simulation. The lists with nodes and cylinder ids will be updated automatically.

	:param [Vector3] vertices: coordinates of vertices to connect in the global coordinate system.
	
	See :yref:`yade.gridpfacet.cylinder` documentation for meaning of other parameters.

	"""
	# create all gridNodes first
	nodesIdsCC = []
	for i in vertices:
		nodesIdsCC.append(
		        O.bodies.append(
		                gridNode(i, radius=radius, dynamic=dynamic, fixed=fixed, wire=wire, color=color, highlight=highlight, material=intMaterial)
		        )
		)
	nodesIds.extend(nodesIdsCC)
	# now create connection between the gridNodes
	for i, j in zip(nodesIdsCC[:-1], nodesIdsCC[1:]):
		cylIds.append(
		        O.bodies.append(
		                gridConnection(
		                        i, j, radius=radius, wire=wire, color=color, highlight=highlight, material=extMaterial, mask=mask, cellDist=None
		                )
		        )
		)


def pfacet(id1, id2, id3, wire=True, color=None, highlight=False, material=-1, mask=1, cellDist=None):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 :yref:`GridNodes<GridNode>` which are already connected via 3 :yref:`GridConnections<GridConnection>`:
	
	:param id1,id2,id3: already with :yref:`GridConnections<GridConnection>` connected :yref:`GridNodes<GridNode>`
	:param bool wire: if ``True``, top and bottom facet are shown as skeleton; otherwise facets are filled.
	:param Vector3-or-None color: color of the PFacet; random color will be assigned if ``None``.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions are not yet implemented for PFacets! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`PFacet<PFacet>` :yref:`shape<Body.shape>`.

	.. note:: :yref:`GridNodes<GridNode>` and :yref:`GridConnections<GridConnection>` need to have the same radius. This is also the radius used to create the :yref:`PFacet<PFacet>`

	"""
	b = Body()
	GridN1 = O.bodies[id1]
	GridN2 = O.bodies[id2]
	GridN3 = O.bodies[id3]
	b.shape = PFacet(color=color if color else utils.randomColor(), wire=wire, highlight=highlight, node1=GridN1, node2=GridN2, node3=GridN3)
	GridN1.bounded = False
	GridN2.bounded = False
	GridN3.bounded = False
	GridC1 = O.bodies[O.interactions[id1, id2].geom.connectionBody.id]
	GridC2 = O.bodies[O.interactions[id2, id3].geom.connectionBody.id]
	GridC3 = O.bodies[O.interactions[id1, id3].geom.connectionBody.id]
	GridC1.bounded = False
	GridC2.bounded = False
	GridC3.bounded = False

	b.shape.conn1 = GridC1
	b.shape.conn2 = GridC2
	b.shape.conn3 = GridC3

	b.shape.radius = GridN1.shape.radius
	GridC1.shape.addPFacet(b)
	GridC2.shape.addPFacet(b)
	GridC3.shape.addPFacet(b)
	GridN1.shape.addPFacet(b)
	GridN2.shape.addPFacet(b)
	GridN3.shape.addPFacet(b)

	V = 0

	utils._commonBodySetup(b, V, Vector3(0, 0, 0), material, pos=GridN1.state.pos, dynamic=False, fixed=True)
	b.aspherical = False  # mass and inertia are lumped into the GridNodes
	b.mask = mask
	return b


#TODO: find a better way of handling the Id lists for checking duplicated gridNodes or gridConnections with the same coordinates etc. It would be better to handle this globally, maybe implement something like O.bodies.getGridNodes
def pfacetCreator1(vertices, radius, nodesIds=[], cylIds=[], pfIds=[], wire=False, fixed=True, materialNodes=-1, material=-1, color=None):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 vertices and automatically append to simulation. The function uses the vertices to create :yref:`GridNodes<GridNode>` and automatically checks for existing nodes.
	
	:param [Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param list nodesIds: list with ids of already existing :yref:`GridNodes<GridNode>`. New ids will be added.
	:param list cylIds: list with ids of already existing :yref:`GridConnections<GridConnection>`. New ids will be added.
	:param list pfIds: list with ids of already existing :yref:`PFacets<PFacet>`. New ids will be added. 
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.
	
	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.
	"""
	n = len(nodesIds)
	k = [0, 0, 0]
	f = [0, 0, 0]
	u = 0
	nod = 0
	for i in vertices:
		u = 0
		for j in nodesIds:
			if (i == O.bodies[j].state.pos):
				f[nod] = j
				k[nod] = 1
				u += 1
		nod += 1
		test = True
		#if(u==0):
		for GN in nodesIds:
			if (i == O.bodies[GN].state.pos):
				u = 1
		if (u == 0):
			nodesIds.append(O.bodies.append(gridNode(i, radius, wire=wire, fixed=fixed, material=materialNodes, color=color)))

	if (k == [0, 0, 0]):
		pfacetCreator3(
		        nodesIds[n], nodesIds[n + 1], nodesIds[n + 2], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed
		)
	if (k == [1, 0, 0]):
		pfacetCreator3(f[0], nodesIds[n], nodesIds[n + 1], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [0, 1, 0]):
		pfacetCreator3(nodesIds[n], f[1], nodesIds[n + 1], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [0, 0, 1]):
		pfacetCreator3(nodesIds[n], nodesIds[n + 1], f[2], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [1, 1, 0]):
		pfacetCreator3(f[0], f[1], nodesIds[n], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [0, 1, 1]):
		pfacetCreator3(nodesIds[n], f[1], f[2], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [1, 0, 1]):
		pfacetCreator3(f[0], nodesIds[n], f[2], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)
	if (k == [1, 1, 1]):
		pfacetCreator3(f[0], f[1], f[2], cylIds=cylIds, pfIds=pfIds, wire=wire, material=material, color=color, fixed=fixed)


def pfacetCreator2(id1, id2, vertex, radius, nodesIds=[], wire=True, materialNodes=-1, material=-1, color=None, fixed=True):
	"""
	Create a :yref:`PFacet<PFacet>` element from 2 already existing and connected :yref:`GridNodes<GridNode>` and one vertex. The element is automatically appended to the simulation.
	
	:param int id1,id2: ids of already with :yref:`GridConnection` connected :yref:`GridNodes<GridNode>`.
	:param Vector3 vertex: coordinates of the vertex in the global coordinate system.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""
	n = len(nodesIds)
	nodesIds.append(O.bodies.append(gridNode(vertex, radius, wire=wire, fixed=fixed, material=materialNodes, color=color)))
	O.bodies.append(gridConnection(id1, nodesIds[n], radius=radius, material=materialNodes, color=color, wire=wire))
	O.bodies.append(gridConnection(id2, nodesIds[n], radius=radius, material=materialNodes, color=color, wire=wire))
	O.bodies.append(pfacet(id1, id2, nodesIds[n], wire=wire, material=material, color=color))


def pfacetCreator3(id1, id2, id3, cylIds=[], pfIds=[], wire=True, material=-1, color=None, fixed=True, mask=-1):
	"""
	Create a :yref:`PFacet` element from 3 already existing :yref:`GridNodes<GridNode>` which are not yet connected. The element is automatically appended to the simulation.
	
	:param int id1,id2,id3: id of the 3 :yref:`GridNodes<GridNode>` forming the :yref:`PFacet`.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""
	radius = O.bodies[id1].shape.radius
	try:
		cylIds.append(O.bodies.append(gridConnection(id1, id2, radius=radius, material=material, color=color, wire=wire, mask=mask)))
	except:
		pass
	try:
		cylIds.append(O.bodies.append(gridConnection(id2, id3, radius=radius, material=material, color=color, wire=wire, mask=mask)))
	except:
		pass
	try:
		cylIds.append(O.bodies.append(gridConnection(id3, id1, radius=radius, material=material, color=color, wire=wire, mask=mask)))
	except:
		pass
	pfIds.append(O.bodies.append(pfacet(id1, id2, id3, wire=wire, material=material, color=color, mask=mask)))


def pfacetCreator4(id1, id2, id3, pfIds=[], wire=True, material=-1, color=None, fixed=True, mask=-1):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 already existing :yref:`GridConnections<GridConnection>`. The element is automatically appended to the simulation.
	
	:param int id1,id2,id3: id of the 3 :yref:`GridConnections<GridConnection>` forming the :yref:`PFacet`.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""
	radius = O.bodies[id1].shape.radius
	GridN = []
	GridN1 = O.bodies[id1].shape.node1.id
	GridN2 = O.bodies[id1].shape.node2.id
	GridN.append(GridN1)
	GridN.append(GridN2)

	GridN1 = O.bodies[id2].shape.node1.id
	if (GridN1 not in GridN):
		GridN.append(GridN1)

	GridN2 = O.bodies[id2].shape.node2.id
	if (GridN2 not in GridN):
		GridN.append(GridN2)

	GridN1 = O.bodies[id3].shape.node1.id
	if (GridN1 not in GridN):
		GridN.append(GridN1)

	GridN2 = O.bodies[id3].shape.node2.id
	if (GridN2 not in GridN):
		GridN.append(GridN2)
	pfIds.append(O.bodies.append(pfacet(GridN[0], GridN[1], GridN[2], wire=wire, material=material, color=color, mask=mask)))


def gtsPFacet(meshfile, shift=Vector3.Zero, scale=1.0, radius=1, wire=True, fixed=True, materialNodes=-1, material=-1, color=None):
	"""
	Imports mesh geometry from .gts file and automatically creates connected :yref:`PFacet3<PFacet>` elements. For an example see :ysrc:`examples/pfacet/gts-pfacet.py`.

	:param string filename: .gts file to read.
	:param [float,float,float] shift: [X,Y,Z] parameter shifts the mesh.
	:param float scale: factor scales the mesh.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:returns: lists of :yref:`GridNode<GridNode>` ids `nodesIds`, :yref:`GridConnection<GridConnection>` ids `cylIds`, and :yref:`PFacet<PFacet>` ids `pfIds`
	"""
	import gts, yade.pack
	surf = gts.read(open(meshfile))
	surf.scale(scale, scale, scale)
	surf.translate(shift[0], shift[1], shift[2])
	nodesIds = []
	cylIds = []
	pfIds = []
	for face in surf.faces():
		a = face.vertices()[0].coords()
		b = face.vertices()[1].coords()
		c = face.vertices()[2].coords()
		pfacetCreator1(
		        [a, b, c],
		        radius=radius,
		        nodesIds=nodesIds,
		        cylIds=cylIds,
		        pfIds=pfIds,
		        wire=wire,
		        fixed=fixed,
		        materialNodes=materialNodes,
		        material=material,
		        color=color
		)
		#print a,b,c
	return nodesIds, cylIds, pfIds


def gmshPFacet(
        meshfile="file.mesh",
        shift=Vector3.Zero,
        scale=1.0,
        orientation=Quaternion.Identity,
        radius=1.0,
        wire=True,
        fixed=True,
        materialNodes=-1,
        material=-1,
        color=None
):
	"""
	Imports mesh geometry from .mesh file and automatically creates connected PFacet elements. For an example see :ysrc:`examples/pfacet/mesh-pfacet.py`.

	:param string filename: .gts file to read.
	:param [float,float,float] shift: [X,Y,Z] parameter shifts the mesh.
	:param float scale: factor scales the mesh.
	:param quaternion orientation: orientation of the imported geometry.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:returns: lists of :yref:`GridNode<GridNode>` ids `nodesIds`, :yref:`GridConnection<GridConnection>` ids `cylIds`, and :yref:`PFacet<PFacet>` ids `pfIds`
	
	mesh files can easily be created with `GMSH <http://www.geuz.org/gmsh/>`_.
	
	Additional examples of mesh-files can be downloaded from 
	http://www-roc.inria.fr/gamma/download/download.php
	"""
	infile = open(meshfile, "r")
	lines = infile.readlines()
	infile.close()

	nodelistVector3 = []
	findVerticesString = 0

	while (lines[findVerticesString].split()[0] != 'Vertices'):  # find the string with the number of Vertices
		findVerticesString += 1
	findVerticesString += 1
	numNodes = int(lines[findVerticesString].split()[0])

	for i in range(numNodes):
		nodelistVector3.append(Vector3(0.0, 0.0, 0.0))
	id = 0

	for line in lines[findVerticesString + 1:numNodes + findVerticesString + 1]:
		data = line.split()
		nodelistVector3[id] = orientation * Vector3(float(data[0]) * scale, float(data[1]) * scale, float(data[2]) * scale) + shift
		id += 1

	findTriangleString = findVerticesString + numNodes
	while (lines[findTriangleString].split()[0] != 'Triangles'):  # find the string with the number of Triangles
		findTriangleString += 1
	findTriangleString += 1
	numTriangles = int(lines[findTriangleString].split()[0])

	triList = []
	for i in range(numTriangles):
		triList.append([0, 0, 0, 0])

	tid = 0
	for line in lines[findTriangleString + 1:findTriangleString + numTriangles + 1]:
		data = line.split()
		id1 = int(data[0]) - 1
		id2 = int(data[1]) - 1
		id3 = int(data[2]) - 1
		triList[tid][0] = tid
		triList[tid][1] = id1
		triList[tid][2] = id2
		triList[tid][3] = id3
		tid += 1

	nodesIds = []
	cylIds = []
	pfIds = []
	for i in triList:
		a = nodelistVector3[i[1]]
		b = nodelistVector3[i[2]]
		c = nodelistVector3[i[3]]
		#print 'i',i
		#print 'a',a
		#print 'b',b
		#print 'c',c
		try:
			pfacetCreator1(
			        [a, b, c],
			        radius=radius,
			        nodesIds=nodesIds,
			        cylIds=cylIds,
			        pfIds=pfIds,
			        wire=wire,
			        fixed=fixed,
			        materialNodes=materialNodes,
			        material=material,
			        color=color
			)
		except:
			pass
	return nodesIds, cylIds, pfIds
