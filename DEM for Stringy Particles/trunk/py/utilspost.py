# enconding: utf-8
##########################################
#*************************************************************************
#*  Copyright (C) 2016 by Sway Zhao                                       *
#*  zhswee@gmail.com                                                      *
#*                                                                        *
#*  This program is free software; it is licensed under the terms of the  *
#*  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/
##########################################
"""
Module containing utinity functions for data post-processing.
"""
## all exported names
#__all__=['outputBoxPov','save_particleinfo','save_contactinfo','save_wallinfo','save_modelinfo','exportForceChains']

from yade import *
from yade.wrapper import *
#from yade import  _superquadrics_utils

#from yade._superquadrics_utils import *

import math,os
import numpy as np
from yade.minieigenHP import *
##########################################
#some auxiliary functions
#output *.POV of a cubic box for post-processing in Pov-ray
def save_particleinfo(filename):
    "output particle info to a file"
    f = open(filename,'w')
    print("###Info of each particle:id radius x y z spin",file=f)
    for i in O.bodies:
        if isinstance(i.shape,Sphere):#focus on superball
            p = i.state.pos
            r = i.shape.radius
            print(i.id,r,p[0],p[1],p[2],i.state.angVel.norm(),file=f)
    f.close()
    
    
def save_contactinfo(filename):
    "output contact info to a file"
    f = open(filename,'w')
    print("###c_type obj1.id obj2.id nforce(x,y,z) sforce(x,y,z)",file=f)
    print("### Info of contacts including particle-particle (named 1) and particle-wall (named 2)",file=f)
    w_contact = list()#store contacts on walls
    b_contact=list()
    for itr in O.interactions:
        if not isinstance(itr.geom,ScGeom):
            continue
        fn = itr.phys.normalForce
        if fn.norm() > 0.0:#repulsive force considered
            #this is a real contact
            fs = itr.phys.shearForce
            id1 = itr.id1#id1 is less than id2 by default
            id2 = itr.id2
            out = ''
            if isinstance(O.bodies[id1].shape,Facet):#wall?
                out = '2 '+str(id1)+' '+str(id2)+' '+str(fn[0])+' '+str(fn[1])+' '+str(fn[2])+' '+str(fs[0])+' '+str(fs[1])+' '+str(fs[2])
                w_contact.append(out)
            elif isinstance(O.bodies[id1].shape,Box):#wall?
                out = '3 '+str(id1)+' '+str(id2)+' '+str(fn[0])+' '+str(fn[1])+' '+str(fn[2])+' '+str(fs[0])+' '+str(fs[1])+' '+str(fs[2])
                b_contact.append(out)
            else:#particle-particle contact
                out = '1 '+str(id1)+' '+str(id2)+' '+str(fn[0])+' '+str(fn[1])+' '+str(fn[2])+' '+str(fs[0])+' '+str(fs[1])+' '+str(fs[2])
                print(out,file=f)
    for i in w_contact:
        print(i,file=f)
    for i in b_contact:
        print(i,file=f)
    f.close()
    
def out_fns_distribution():
    fns=[]
    for itr in O.interactions:
        if not isinstance(itr.geom,ScGeom):
            continue
        fn = itr.phys.normalForce.norm()
        fns.append(fn)
    min_fn=0
    max_fn=int(max(fns))
    import numpy as np
    import matplotlib.pyplot as plt
    bins = np.arange(min_fn, max_fn+1, 1)  # 6+1是因为np.arange不包括终点
    hist, bins = np.histogram(fns, bins=bins)
   # print("频率分布：")
    #for i in range(len(bins)-1):
       # print(f"{bins[i]} - {bins[i+1]}: {hist[i]}")
    plt.hist(fns, bins=bins, edgecolor='black',alpha=0.7)
    plt.xlabel('fn/N数值')
    plt.ylabel('frecucy频率')
    plt.title('数据频率分布直方图')
    plt.show()  

def out_threadfns_distribution():
    fns=[]
    for t in O.threads:
        fn=t.get_f()
        if fn>0.01:
            fns.append(fn)
    min_fn=0
    max_fn=int(max(fns))
    import numpy as np
    import matplotlib.pyplot as plt
    bins = np.arange(min_fn, max_fn+1, 1)  # 6+1是因为np.arange不包括终点
    hist, bins = np.histogram(fns, bins=bins)
   # print("频率分布：")
    #for i in range(len(bins)-1):
       # print(f"{bins[i]} - {bins[i+1]}: {hist[i]}")
    plt.hist(fns, bins=bins, edgecolor='black',alpha=0.7)
    plt.xlabel('fn/N数值')
    plt.ylabel('frecucy频率')
    plt.title('数据频率分布直方图')
    plt.show()
    plt.close()      
      
      
def exportForceChains(path,step=0):
    "export force chains with path of input files and step number."
    path = os.path.dirname(path+'/')
    if not os.path.exists(path):#not exists
        print ("The typed path dose not exist!")
        return False
    particlefile = os.path.join(path,'particleinfo_'+str(step)+'.txt')
    contactfile = os.path.join(path,'contactinfo_'+str(step)+'.txt')
    #wallfile = os.path.join(path,'wallinfo_'+str(step)+'.txt')
    #check if files exist
    if not (os.path.isfile(particlefile) and os.path.isfile(contactfile) ):
        print ("Input files are missing! Please check the path of particleinfo*.txt, contactinfo*.txt and wallinfo*.txt is set correctly.")
        return False
    ForceChains(particlefile,contactfile, comment="comment")
    return True
    
def ForceChains(inputfile,ct_filename,comment="comment"):
    '''
    inputfile: ballinfo
    ct_filename:contactinfo
    w_filename:wallinfo
    '''
    '''
    wall_file = open(w_filename, 'r')
    lines = wall_file.readlines()
    wall_file.close()
    w_pos = [i[:-2].split(' ') for i in lines]
    wall_planes=[
            [1,0,0,float(w_pos[0][0])],#wall 3,x
            [1,0,0,float(w_pos[1][0])],#wall 4,x
            [0,1,0,float(w_pos[2][1])],#wall 5,y
            [0,1,0,float(w_pos[3][1])], #wall 6,y
            [0,0,1,float(w_pos[4][2])],#wall 1,z
                [0,0,1,float(w_pos[5][2])]#wall 2,z
            ]
            '''
    f = open(inputfile,'r')
    # output file
    fName = ct_filename +'.vtp'
    outContactFile = open(fName, 'w')
    lines = f.readlines()[1:]########## the range may need to be changed
    f.close()
    #open the original file including the info of all balls
    nBodies = 0
    radius = dict()
    position =dict()
    for l in lines:
        l1 = l.lstrip()
        l1 = l1[:-2].split(' ')
        #print l1
        b_id = int(l1[0])
        position[b_id] = [float(l1[2]),float(l1[3]),float(l1[4])]
        radius[b_id] = float(l1[1])
        nBodies += 1
    # output file
    contact_file = open(ct_filename, 'r')
    lines = contact_file.readlines()[2:]
    contact_file.close()
    ##########################################################################
    ###contacts
    nIntrs = len(lines)
    # head
    outContactFile.write("<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>\n<PolyData>\n")
    outContactFile.write("<Piece NumberOfPoints='%s' NumberOfVerts='0' NumberOfLines='%s' NumberOfStrips='0' NumberOfPolys='0'>\n"%(str(2*nIntrs),str(nIntrs)))
    # write coords of intrs bodies (also taking into account possible periodicity
    outContactFile.write("<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n")
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        #print l1[:3],'sss',l1[3:]
        [contact_type,id1,id2] = [int(i) for i in l1[:3]]
        #[fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        #[fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        if contact_type < 2:#1: ball-ball contact
            #find positions of two touching balls
            #print("contact id1 is",id1)
            pos = position[id1]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
            pos = position[id2]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
        '''else:#2 ball_wall contact

            pos = PointonPlan(wall_planes[id1],position[id2])    ###???
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
            pos = position[id2]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))'''
        #print contact_type,id1,id2,fn_x,fn_y,fn_z,fs_x,fs_y,fs_z

    outContactFile.write("</DataArray>\n</Points>\n<Lines>\n<DataArray type='Int32' Name='connectivity' format='ascii'>\n")

    ss=''
    for con in range(2*nIntrs):
        ss+=' '+str(con)
    outContactFile.write(ss+'\n')
    outContactFile.write("</DataArray>\n<DataArray type='Int32' Name='offsets' format='ascii'>\n")
    ss=''
    for con in range(nIntrs):
        ss+=' '+str(con*2+2)
    outContactFile.write(ss)
    outContactFile.write("\n</DataArray>\n</Lines>\n")
    ##
    name = 'Fn'
    outContactFile.write("<PointData Scalars='%s'>\n<DataArray type='Float32' Name='%s' format='ascii'>\n"%(name,name))
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        #[contact_type,id1,id2] = [int(i) for i in l1[:3]]
        [fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        fn=(fn_x**2.+fn_y**2.+fn_z**2.)**0.5
        outContactFile.write("%g %g\n"%(fn,fn))
    outContactFile.write("</DataArray>\n</PointData>")
    outContactFile.write("\n</Piece>\n</PolyData>\n</VTKFile>")
    outContactFile.close()
    
def save_modelinfo(path,step=0):#保存
    "output all model info to three individual files corresponding to info of particles, contacts and walls."
    #create the path if it does not exists
    step=O.iter
    '''
    path=path+'/step_'+str(step)#判断是否存在
    if modle:
    	namemiddel="sphere"
    else:
    	namemiddel="LineNode"
    file_path=path+'/particleVtk_'+namemiddel+'.vtk'#path+
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
    	os.makedirs(directory)
    '''	
    filename=str(step)
    path = os.path.dirname(path+'/'+'step_'+filename+'/')
    if not os.path.exists(path):#not exists
        os.makedirs(path)
        print ("fiel create!")
    #print("path="+path)
    save_particleinfo(path+'/particleinfo_'+str(step)+'.txt')
    save_contactinfo(path+'/contactinfo_'+str(step)+'.txt')
    #save_wallinfo(path+'/wallinfo_'+str(step)+'.txt')
    print ("Model info output finished!")
    
def saveandout():#目前已弃用
	path="/home/hcl/桌面/1/savedata"
	iter_=O.iter
	save_modelinfo(path,iter_)
	exportForceChains(path,iter_)
	
def testoutCylinderVTK(filename,num):# for test only
	def rotation(q,p,theata):
		RotaTheta=0.5*theata
		QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
		QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
		QuaternionP=Quaternion(0,p[0],p[1],p[2])
		q1=QuaternionQ*QuaternionP
		q2=q1*QuaternionQ_conjugate
		v=(q2[0],q2[1],q2[2])
		return v
	radius=1.5e-3
	dist=5e-3
	number=num
	pos=Vector3(0,0,0)
	vectorz=Vector3(0,0,1)
	normal=Vector3(3,4,5)
	l=normal.norm()
	_normal=normal/l
	rotation_angle=math.acos(np.dot(vectorz,_normal))#jiaodu rad
	axis=np.cross(vectorz,_normal)
	axis=axis/np.linalg.norm(axis)#xuanzhuanzhou
	
	points=[]
	_points=[]
	theta=2*math.pi/number
	for i in range(number):
		pos=Vector3(math.cos(i*theta)*radius,math.sin(i*theta)*radius,0)
		points.append(pos)
	for i in range(number):
		pos=Vector3(math.cos(i*theta)*radius,math.sin(i*theta)*radius,dist)
		points.append(pos)
	for i in points:
		ri=rotation(axis,i,rotation_angle)
		_points.append(ri)
		#print('points oringinal',i)
		#print('points rotate',ri)
	outFile = open(filename, 'w')	
	print('writing sphere into a VTK file')
	outFile.write("# vtk DataFile Version 2.0\nvtk output\nASCII\nDATASET POLYDATA\n")
#	outFile.write(POINTS   ',str(len(_points)),' float','\n')
	print('POINTS   ',str(len(_points)),' float',file=outFile)
	for i in _points:
		print(i[0],' ',i[1],' ',i[2],file=outFile)
#	outFile.write('POLYGONS ',number,number*5,'\n')
	print('POLYGONS ',number,number*5,file=outFile)
	for i in range(number-1):
		print('4 ',i,' ',i+1,' ',i+number+1,' ',i+number,' ',file=outFile)
	print('4 ',number-1,' ',0,' ',number,' ',2*number-1,' ',file=outFile)
	outFile.close()
	

	
def outNodeSphereI(filename):
	f=open(filename+'vtk','w')
	f.write("interaction recorder begin \n")
	contact=list()
	for itr in O.interactions:
		if(isinstance(itr.phys,LineFrictPhys)):
			id1=itr.id1
			id2=itr.id2
			out=''
			fn = itr.phys.normalForce.norm()
			pe=itr.geom.penetrationDepth
			out=str(id1)+' '+str(id2)+' '+str(fn)+' '+str(pe)+'\n'
			contact.append(out)
	for i in contact:
		f.write(i)
	f.close()
def outNodeSphereI2(filename):
	f=open(filename,'w')
	f.write("interaction recorder begin \n")
	contact=list()
	for itr in O.interactions:
		if(isinstance(itr.phys,LineFrictPhys)):
			out=''
			pe=itr.geom.penetrationDepth
			out=str(pe)+'\n'
			f.write(out)	
	f.close()

##配位数计算
def CoordinateNum(BoxSize=[0.35,0.35,0.2],circle=False):  #配位数计算
    #para:BoxSize,xlength,ylength,height are x,y,z length of the container respectively.
    #out of date para:r, the search range of a particle
    #para:meaBox,measure box which specifies which particle should be taken into account.
    #para:all,if all particles are used to calculate.
    #caution:To get the right result, one should execute 'O.run(1)' before this function.Now I also don't know the reason.
    bodies=[]# bodies

    #get the polyhedrons, then put them into bodies
    if circle==True:
        for i in O.bodies:
            if isinstance(i.shape,Superquadrics):
                rr=math.sqrt((i.state.pos[0]-0.1)**2.0+(i.state.pos[1]-0.1)**2.0)
                if rr<=0.1:
                    bodies.append(i)
    else:
        for i in O.bodies:
            if isinstance(i.shape,Sphere):
                bodies.append(i)
    #process the boides
    CNlist=[]   #store CN
    for i in bodies:#循环body
        count=0
        #check others around the particle
        intrs=i.intrs() #intrs around this particle (*it).second interaction
        for j in intrs:
            if j.id1==i.id:
                id2=j.id2
            else:
                id2=j.id1

            #is it polyhedron?
            if isinstance(O.bodies[id2].shape,Box):#facet? Box
                count+=1
                print('sphere with Box interaction')
            else:
                obj1=i
                obj2=O.bodies[id2]
                if check(obj1.shape,obj2.shape,obj1.state,obj2.state) :
                    #if True,then count it.
                    count+=1
        CNlist.append(count)
    #process CN,maybe need to do a statistics
    FC=dict()   #store frenquency count 配位数频率分布
    for i in CNlist:
        if i in FC.keys():
            FC[i]=FC[i]+1
        else:
            FC[i]=1
    return [FC,CNlist]
    
def check(_shape1,_shape2,_state1,_state2):#检查位置
    r=_shape1.radius
    if sum([(_state1.pos[k]-_state1.pos[k])**2.0 for k in range(3)])<r**2.0:
        return True
    else:
        return False


#********************************************************************************************************
def rotation(q,p,theata): #输入旋转轴角度 计算旋转向量
	RotaTheta=0.5*theata
	QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
	QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
	QuaternionP=Quaternion(0,p[0],p[1],p[2])
	q1=QuaternionQ*QuaternionP
	q2=q1*QuaternionQ_conjugate
	v=(q2[0],q2[1],q2[2])
	return v 

def computePoints(_pos,_dist,_radius,_num):# for test only
    radius=_radius
    number=_num
    pos_local=Vector3(0,0,0)
    vectorz=Vector3(0,0,1)
    normal=Vector3(_dist[0],_dist[1],_dist[2])
    l=normal.norm()#length
    _normal=normal/l#danwiexiangliang
    rotation_angle=math.acos(np.dot(vectorz,_normal))#xuanzhuanjiaodu rad
    axis=np.cross(vectorz,_normal)
    axis=axis/np.linalg.norm(axis)#xuanzhuanzhou
    points=[]
    _points=[]#局部坐标
    theta=2*math.pi/number
    for i in range(number):
        pos=Vector3(math.cos(i*theta)*radius,math.sin(i*theta)*radius,0)
        points.append(pos)
    for i in range(number):
        pos=Vector3(math.cos(i*theta)*radius,math.sin(i*theta)*radius,l)
        points.append(pos)
    for i in points:
        ri=rotation(axis,i,rotation_angle)#旋转后坐标
        global_pos=[ri[0]+_pos[0],ri[1]+_pos[1],ri[2]+_pos[2]]
        _points.append(global_pos)
    return _points
'''def get_maxfn():
    maxfn=0.
    for i in O.interactions:
    	if isinstance(i.geom,ScGeom):
    		maxfn=max(maxfn,i.phys.normal)'''
def ForceChains_specific(inputfile,ct_filename,slices,maxfn,comment="comment"):#slices 圆柱片数
    path=os.getcwd()
    filename=str(O.iter)
    path= os.path.dirname(path+'/'+'step_'+filename+'/')#'/'自动删除
    sradius=1.5e-3
    f = open(inputfile,'r')#颗粒信息文件
    # output file
    #print("path="+path)#
    fName = path +'/Force_chain'+'.vtk'
    outContactFile = open(fName, 'w')
    lines = f.readlines()[1:]########## the range may need to be changed
    f.close()
    #open the original file including the info of all balls
    nBodies = 0
    radius = dict()
    position =dict()
    for l in lines:
        l1 = l.lstrip()#截掉空格
        l1 = l1[:-2].split(' ')
        #print l1
        b_id = int(l1[0])
        position[b_id] = [float(l1[2]),float(l1[3]),float(l1[4])]
        radius[b_id] = float(l1[1])
        nBodies += 1
    # output file
    contact_file = open(ct_filename, 'r')#接触信息文件
    lines = contact_file.readlines()[2:]
    contact_file.close()
    ##########################################################################
    ###contacts
    #nIntrs = len(lines)#接触数量
    #num_points = nIntrs*slices*2
    # head
    fns=[]
    #maxfn=8.#最大力如何定义以及获取方式
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        contact_type=int(l1[0])
        if contact_type < 2:
        	[fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        	fn=(fn_x**2.+fn_y**2.+fn_z**2.)**0.5
        	fns.append(fn)
    _maxfn=max(fns)
    print('_maxfn=',_maxfn)
    _meanfn=sum(fns)/len(fns)
    print('_meanfn=',_meanfn)
    f_radius=[] 
    for i in fns:
        #f_radius.append(i/maxfn*sradius)#最大半径的限制
        nradius=math.log((math.exp(1)-1)/maxfn*i+1)*sradius*0.7#最大半径的限制0.7
        f_radius.append(nradius) 	
    print('writing sphere force_chain into a VTK file')
    outContactFile.write("# vtk DataFile Version 2.0\nvtk output\nASCII\nDATASET POLYDATA\n")
    
    realnum=0
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        contact_type=int(l1[0])
        if contact_type < 2:
        	realnum+=1
    print('POINTS   ',str(realnum*slices*2),' float',file=outContactFile)
    count=0
    for l in lines:#suoyoujiechu
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        [contact_type,id1,id2] = [int(i) for i in l1[:3]]
        if contact_type < 2:
            pos1 = position[id1]#local coordnate
            pos2 = position[id2]
            dist=[pos2[0]-pos1[0],pos2[1]-pos1[1],pos2[2]-pos1[2]]#1-2
            #points=[]
            points=computePoints(pos1,dist,f_radius[count],slices)
            count+=1
            '''for i in range(len(points)):
                for j in range(3):
                    points[i][j]=points[i][j]+pos1[j]'''
            for i in points:
                print(i[0],' ',i[1],' ',i[2],file=outContactFile)

    print('POLYGONS ',realnum*slices,realnum*slices*5,file=outContactFile)#all
    for i in range(realnum):
        global_num=i*2*slices
        for j in range(slices-1):
            print('4 ',global_num+j,' ',global_num+j+1,' ',global_num+j+slices+1,' ',global_num+j+slices,' ',file=outContactFile)
        print('4 ',global_num+slices-1,' ',global_num+0,' ',global_num+slices,' ',global_num+2*slices-1,' ',file=outContactFile)
    print('CELL_DATA ',realnum*slices,file=outContactFile)#all
    print('SCALARS ','sample_scalars ','float ','1',file=outContactFile)#all
    print('LOOKUP_TABLE ','default',file=outContactFile)
    
    
    for i in fns:
    	#pstr=str(i/maxfn)+'\n'
    	for j in range(slices):
    		#print(str(i/maxfn*7),file=outContactFile)#cloro map
    		print(str(i),file=outContactFile)#cloro map
    		'''if i>_meanfn:
    			print(str(7),file=outContactFile)#cloro map
    		else:
    			print(str(0),file=outContactFile)#cloro map'''	
    
    outContactFile.close()

def exportForceChains_specific(path,num,maxfn,step=0):
    "export force chains with path of input files and step number."
    step=O.iter
    filename=str(step)
    path = os.path.dirname(path+'/'+'step_'+filename+'/')
    #path = os.path.dirname(path+'/')
    if not os.path.exists(path):#not exists
        os.makedirs(path)
        print ("fiel create!")
        #return False
    particlefile = os.path.join(path,'particleinfo_'+str(step)+'.txt')
    contactfile = os.path.join(path,'contactinfo_'+str(step)+'.txt')
    #wallfile = os.path.join(path,'wallinfo_'+str(step)+'.txt')
    #check if files exist
    if not (os.path.isfile(particlefile) and os.path.isfile(contactfile) ):
        print ("Input files are missing! Please check the path of particleinfo*.txt, contactinfo*.txt and wallinfo*.txt is set correctly.")
        return False
    ForceChains_specific(particlefile,contactfile,num,maxfn,comment="comment")
    return True

def saveandout_ForceChains_specific(num,maxfn=8.0):#导出颗粒力链总函数  
	path="/home/hcl/桌面/1/savedata"
	path=os.getcwd()
	iter_=O.iter
	save_modelinfo(path,iter_)#保存模型信息
	exportForceChains_specific(path,num,maxfn,iter_)#  
	print("maxfn setting is :"+str(maxfn))    
        
#******************************************************************************************************** 柔索力链导出模块
def save_linenodeinfo(filename):
    "output node info to a file"
    f = open(filename,'w')
    print("###Info of each node:id radius x y z spin",file=f)
    for i in O.bodies:
        if isinstance(i.shape,LineNode):#focus on superball
            p = i.state.pos
            r = i.shape.radius
            print(i.id,r,p[0],p[1],p[2],file=f)
    f.close()

def save_node_contactinfo(filename,_F_threshold=0.1):
    "output contact info to a file"
    f = open(filename,'w')
    print("###c_type obj1.id obj2.id Fn(normalforce)",file=f)
    print("### Info of contacts including node-node (named 1) and particle-wall (named 2)",file=f)
    w_contact = list()#store contacts on walls
    b_contact=list()
    F_threshold=0.1#1g 0.01N   10n 1000g      10g 0.1N                                  
    num=0
    under_f_t=list()
    for itr in O.interactions:
        if not isinstance(itr.phys,LinePhys):
            continue
        fn = itr.phys.Fn
        if fn > F_threshold:#Force threshold
            #this is a real contact
            id1 = itr.id1#id1 is less than id2 by default
            id2 = itr.id2
            out = ''
            out = '1 '+str(id1)+' '+str(id2)+' '+str(fn)
            print(out,file=f)
            num+=1
        else:
            id1 = itr.id1#id1 is less than id2 by default
            id2 = itr.id2
            out = ''
            out = '2 '+str(id1)+' '+str(id2)+' '+str(fn)
            under_f_t.append(out)
    for i in under_f_t:
        print(i,file=f)
    f.close()
    print('num of nodecontact which force under:',F_threshold,' is ',num)

def save_node_modelinfo(path,step=0):#保存
    "output all model info to three individual files corresponding to info of particles, contacts and walls."
    #create the path if it does not exists
    path = os.path.dirname(path+'/node_step_'+str(step)+'/')#无/
    if not os.path.exists(path):#not exists
        os.makedirs(path)
    print("path="+path)
    save_linenodeinfo(path+'/nodeinfo_'+str(step)+'.txt')
    save_node_contactinfo(path+'/node_contactinfo_'+str(step)+'.txt')
    #save_wallinfo(path+'/wallinfo_'+str(step)+'.txt')
    print ("Model info output finished!")
    
def saveAndOut_node_forcechain(_num,maxfn=8.0): #柔索力链导出总函数
	#path="/home/hcl/桌面/1/savedata/thread"
	path=os.getcwd()
	iter_=O.iter
	save_node_modelinfo(path,iter_)#保存模型信息
	exportForceChains_node_specific(path,_num,maxfn,iter_)#导出力链

def exportForceChains_node_specific(path,num,maxfn,step=0):
    "export force chains with path of input files and step number."
    step=O.iter
    path = os.path.dirname(path+'/node_step_'+str(step)+'/')
    if not os.path.exists(path):#not exists
        print ("The typed path dose not exist!")
        return False
    nodefile = os.path.join(path,'nodeinfo_'+str(step)+'.txt')
    node_contactfile = os.path.join(path,'node_contactinfo_'+str(step)+'.txt')
    #wallfile = os.path.join(path,'wallinfo_'+str(step)+'.txt')
    #check if files exist
    if not (os.path.isfile(nodefile) and os.path.isfile(node_contactfile) ):
        print ("Input files are missing! Please check the path of particleinfo*.txt, contactinfo*.txt and wallinfo*.txt is set correctly.")
        return False
    ForceChains_node_specific(nodefile,node_contactfile,num,maxfn,comment="comment")
    return True

def ForceChains_node_specific(inputfile,ct_filename,slices,maxfn,comment="comment"):#slices 圆柱片数 柔索力链信息读取后导出
    #sradius=8e-5
    path=os.getcwd()
    step=O.iter
    path= os.path.dirname(path+'/'+'node_step_'+str(step)+'/')#'/'自动删除
    sradius=1.5e-3
    f = open(inputfile,'r')#颗粒信息文件
    # output file
    #print("path="+path)#
    fName = path +'/Force_chain'+'.vtk'
    outContactFile = open(fName, 'w')
    '''
    sradius=1.5e-3
    f = open(inputfile,'r')#节点颗粒信息文件
    # output file
    #fName = ct_filename +'threads_chain.vtk'
    print("path="+path_save)
    fName = '/threads_chain'+'.vtk'
    path=path_save+fName
    outContactFile = open(fName, 'w')'''
    lines = f.readlines()[1:]########## the range may need to be changed
    f.close()
    #open the original file including the info of all balls
    nBodies = 0
    radius = dict()
    position =dict()
    for l in lines:
        l1 = l.lstrip()#截掉端空格
        l1 = l1[:-2].split(' ')#按中间空格分隔组
        #print l1
        b_id = int(l1[0])
        position[b_id] = [float(l1[2]),float(l1[3]),float(l1[4])]
        radius[b_id] = float(l1[1])
        nBodies += 1

    # output file
    contact_file = open(ct_filename, 'r')#接触信息文件
    lines = contact_file.readlines()[2:]#去掉1 2 行
    contact_file.close()
   
    fns=[]
    #maxfn=8.#最大接触力定义 8N 主函数传入
    for l in lines:
        l=l.lstrip()
        #l1 = l[:-2].split(' ')
        l1 = l.split(' ')
        contact_type=int(l1[0])
        fn=float(l1[3])
        if(contact_type<2):
            fns.append(fn)
            
    _maxfn=max(fns)
    print('max threads_fns=',_maxfn)
    print('num of nodecontact:',len(fns))

    f_radius=[] #力链半径大小  后期添加半径大小映射函数
    for i in fns:
        f_radius.append(i/maxfn*sradius*0.8)  	
    print('writing sphere into a VTK file')
    outContactFile.write("# vtk DataFile Version 2.0\nvtk output\nASCII\nDATASET POLYDATA\n")
    realnum=len(fns)
    print('POINTS   ',str(realnum*slices*2),' float',file=outContactFile)
    count=0
    for l in lines:#suoyoujiechu
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        [contact_type,id1,id2] = [int(i) for i in l1[:3]]
        fn=float(l1[3])
        if contact_type < 2 :
            pos1 = position[id1]#local coordnate
            pos2 = position[id2]
            dist=[pos2[0]-pos1[0],pos2[1]-pos1[1],pos2[2]-pos1[2]]#1-2
            #points=[]
            points=computePoints(pos1,dist,f_radius[count],slices)
            count+=1
            for i in points:
                print(i[0],' ',i[1],' ',i[2],file=outContactFile)

    print('POLYGONS ',realnum*slices,realnum*slices*5,file=outContactFile)#all
    for i in range(realnum):
        global_num=i*2*slices
        for j in range(slices-1):
            print('4 ',global_num+j,' ',global_num+j+1,' ',global_num+j+slices+1,' ',global_num+j+slices,' ',file=outContactFile)
        print('4 ',global_num+slices-1,' ',global_num+0,' ',global_num+slices,' ',global_num+2*slices-1,' ',file=outContactFile)
    print('CELL_DATA ',realnum*slices,file=outContactFile)#all
    print('SCALARS ','sample_scalars ','float ','1',file=outContactFile)#all
    print('LOOKUP_TABLE ','default',file=outContactFile)
    
    for i in fns:
    	#pstr=str(i/maxfn)+'\n'
    	for j in range(slices):
    		#print(str(i/maxfn*7),file=outContactFile)#cloro map
    		print(str(i),file=outContactFile)#cloro map
    outContactFile.close()
#******************************************************************************************************** 
#******************************************************************************************************** 导出常规颗粒 节点 柔索段的vtk文件函数
def out_particle_vtk_pure(slices,modle=1,color=2):#0 sphere pureout     1 sphere 0 node 2 lineconection
    path=os.getcwd()
    step=O.iter
    path=path+'/step_'+str(step)#判断是否存在
    if modle:
    	namemiddel="sphere"
    else:
    	namemiddel="LineNode"
    file_path=path+'/particleVtk_'+namemiddel+'.vtk'#path+
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
    	os.makedirs(directory)
    f = open(file_path,'w')
    
    print("writing sphere into a VTK file")
    print("# vtk DataFile Version 2.0",file=f)
    print("sphere vtk output",file=f)
    print("ASCII",file=f)
    print("DATASET POLYDATA",file=f)
    vertices=[]
    point_counter=0
    w=slices*2 #经度数量
    h=slices  #纬度数量
    ids=[]
    for i in O.bodies:
        if modle and isinstance(i.shape,Sphere):#focus on superball
            point_counter+=1
            ids.append(i.id)
            pos = i.state.pos
            r = i.shape.radius
            ori=i.state.se3[1]
            a=0.
            b1=0.
            hStep=math.pi/(h-1) #间隔
            wStep=2*math.pi/w
            Rot=ori.toRotationMatrix()
            #ori = i.state.se3[1].toRotationMatrix()*ore_v
            for i in range(h):
                for j in range(w):
                    Surf=getSurface(r,Vector2(b1+j*wStep,a+i*hStep-math.pi*0.5))
                    Surf=Rot*Surf+pos
                    vertices.append(Surf)
            continue
        elif not modle and isinstance(i.shape,LineNode):#focus on superball
            point_counter+=1
            ids.append(i.id)
            pos = i.state.pos
            r = i.shape.radius
            ori=i.state.se3[1]
            a=0.
            b1=0.
            hStep=math.pi/(h-1) #间隔
            wStep=2*math.pi/w
            Rot=ori.toRotationMatrix()
            #ori = i.state.se3[1].toRotationMatrix()*ore_v
            for i in range(h):
                for j in range(w):
                    Surf=getSurface(r,Vector2(b1+j*wStep,a+i*hStep-math.pi*0.5))
                    Surf=Rot*Surf+pos
                    vertices.append(Surf)
    print('POINTS   ',str(len(vertices)),' float',file=f)
    for i in vertices:
        print(i[0],' ',i[1],' ',i[2],file=f)
    print('POLYGONS ',str(point_counter*(h-1)*w),' ',str(point_counter*(h-1)*w*5),file=f)
    for i in range(point_counter):
        start_index=i*w*h
        m=0
        for mm in range(h-1):
            n=0
            for nn in range(w-1):
                print('4 ',str(start_index+m*w+n),' ',str(start_index+m*w+n+1),' ',str(start_index+(m+1)*w+n+1),' ',str(start_index+(m+1)*w+n),' ',file=f)
                n+=1
            print('4 ',str(start_index+m*w+n),' ',str(start_index+m*w),' ',str(start_index+(m+1)*w),' ',str(start_index+(m+1)*w+n),' ',file=f)
            m+=1
   
    print('CELL_DATA ',point_counter*(h-1)*w,file=f)#all polygons
    print('SCALARS ','rotation_scalars ','float ','1',file=f)#all
    print('LOOKUP_TABLE ','default',file=f)
    
    
    #maxRotation_angle=math.pi/6#最大旋转角定义
    for i in ids:#
    	for j in range((h-1)*w):
    		#print(str(i/maxRotation_angle*7),file=f)
    		print(color,file=f)
    f.close()
    print('writing sphere/node with color into a VTK file succusful!')
    return True

#******************************************************************************************************** 导出柔索段的vtk文件函数
def out_cylinder_vtk(slices,color=2):#0 sphere pureout     1 sphere 0 node 2 lineconection
    #path="/home/hcl/桌面/1/savedata"
    path=os.getcwd()
    #output particle info to a file
    step=O.iter
    path=path+'/step_'+str(step)#判断是否存在
    
    namemiddel=""
    file_path=path+'/LineConnectionVtk_'+namemiddel+'.vtk'#path+
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
    	os.makedirs(directory)
    #with open(file_path, 'w') as f
    f = open(file_path,'w')
    print("writing lineconection into a VTK file")
    print("# vtk DataFile Version 2.0",file=f)
    print("sphere vtk output",file=f)
    print("ASCII",file=f)
    print("DATASET POLYDATA",file=f)
    vertices=[]
    point_counter=0
    #w=slices*2 #经度数量
    #h=slices  #纬度数量
    refV=Vector3(0,0,1)#
    ids=[]
    realnum=0
    for i in O.bodies:
        if isinstance(i.shape,LineConnection):#focus on superball
            point_counter+=1
            ids.append(i.id)
            realnum+=1
    print('POINTS   ',str(realnum*slices*2),' float',file=f)
    for bid in ids:
        body=O.bodies[bid]
        pos1 = body.shape.node1.state.pos
        pos2 = body.shape.node2.state.pos
        director=pos2-pos1
        dist=director.norm()
        nomal_v=director/dist
        r = body.shape.radius
        points=computePoints(pos1,director,r,slices)
        for i in points:
            print(i[0],' ',i[1],' ',i[2],file=f)
            
    print('POLYGONS ',realnum*slices,realnum*slices*5,file=f)#all
    for i in range(realnum):
        global_num=i*2*slices
        for j in range(slices-1):
            print('4 ',global_num+j,' ',global_num+j+1,' ',global_num+j+slices+1,' ',global_num+j+slices,' ',file=f)
        print('4 ',global_num+slices-1,' ',global_num+0,' ',global_num+slices,' ',global_num+2*slices-1,' ',file=f)
    print('writing lineconection with color into a VTK file succusful!')
    return True
    
    
    
#******************************************************************************************************** 
#******************************************************************************************************** 导出边界墙体的vtk文件函数
def out_wall_vtk(wallids,model=1,color=2):#0 sphere pureout     1 sphere 0 node 2 lineconection  model=1为真实墙体vtk 0以试样边界为主
    #path="/home/hcl/桌面/1/savedata"
    path=os.getcwd()
    #output particle info to a file
    step=O.iter
    path=path+'/step_'+str(step)#判断是否存在
    
    namemiddel=str(len(wallids))
    file_path=path+'/wallVtk_'+namemiddel+'.vtk'#path+
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
    	os.makedirs(directory)
    #with open(file_path, 'w') as f
    f = open(file_path,'w')
    def getpoint(id):
    	body=O.bodies[id]
    	if not isinstance(body.shape,Box):
            print("id is not a box")
            return False
    	pos=body.state.pos
    	extents=body.shape.extents
    	index=0
    	for i in extents:#获取0扩展的边界方向
            if i==0:
            	break
            index+=1
    	verticles=[]
    	if index==0:#x
            verticles.append(Vector3(pos[0],pos[1]+extents[1],pos[2]-extents[2]))
            verticles.append(Vector3(pos[0],pos[1]-extents[1],pos[2]-extents[2]))
            verticles.append(Vector3(pos[0],pos[1]-extents[1],pos[2]+extents[2]))
            verticles.append(Vector3(pos[0],pos[1]+extents[1],pos[2]+extents[2]))
    	if index==1:#y
            verticles.append(Vector3(pos[0]-extents[0],pos[1],pos[2]-extents[2])) #+-  -- -+ ++
            verticles.append(Vector3(pos[0]+extents[0],pos[1],pos[2]-extents[2])) 
            verticles.append(Vector3(pos[0]+extents[0],pos[1],pos[2]+extents[2]))
            verticles.append(Vector3(pos[0]-extents[0],pos[1],pos[2]+extents[2]))
    	if index==2:#y
            verticles.append(Vector3(pos[0]+extents[0],pos[1]-extents[1],pos[2]))
            verticles.append(Vector3(pos[0]-extents[0],pos[1]-extents[1],pos[2]))
            verticles.append(Vector3(pos[0]-extents[0],pos[1]+extents[1],pos[2]))
            verticles.append(Vector3(pos[0]+extents[0],pos[1]+extents[1],pos[2]))
    	if model==0 and (id==2 or id==3):#更改2 3 墙最大z值
            verticles[2][2]=O.bodies[7].state.pos[2]#7默认墙体
            verticles[3][2]=O.bodies[7].state.pos[2]
    	return verticles   
        
    print("writing lineconection into a VTK file")
    print("# vtk DataFile Version 2.0",file=f)
    print("sphere vtk output",file=f)
    print("ASCII",file=f)
    print("DATASET POLYDATA",file=f)
    print('POINTS   ',str(len(wallids)*4),' float',file=f)
    for i in wallids:
        points=getpoint(i)
        for i in points:
            print(i[0],' ',i[1],' ',i[2],file=f)
            
    print('POLYGONS ',str(len(wallids)),len(wallids)*5,file=f)#all
    for i in range(len(wallids)):
        print('4 ',4*i,' ',4*i+1,' ',4*i+2,' ',4*i+3,' ',file=f)
        #print('4 ',global_num+slices-1,' ',global_num+0,' ',global_num+slices,' ',global_num+2*slices-1,' ',file=f)
    print('writing wall with color into a VTK file succusful!')
    return True
    
    
    
#******************************************************************************************************** 
        
#******************************************************************************************************** 导出颗粒以及颗粒旋转角度vtk文件模块
def out_particle_vtk(slices,thetas,model=0):#0 sphere pureout     1 sphere with color of Rotation out
    path="/home/hcl/桌面/1/savedata"
    #output particle info to a file
    step=O.iter
    filename=path+'/particleRotation_'+str(step)+'.vtk'
    f = open(filename,'w')
    print("writing sphere into a VTK file")
    print("# vtk DataFile Version 2.0",file=f)
    print("sphere vtk output",file=f)
    print("ASCII",file=f)
    print("DATASET POLYDATA",file=f)
    vertices=[]
    point_counter=0
    w=slices*2 #经度数量
    h=slices  #纬度数量
    ids=[]
    for i in O.bodies:
        if isinstance(i.shape,Sphere):#focus on superball
            point_counter+=1
            ids.append(i.id)
            pos = i.state.pos
            r = i.shape.radius
            ori=i.state.se3[1]
            a=0.
            b1=0.
            hStep=math.pi/(h-1) #间隔
            wStep=2*math.pi/w
            Rot=ori.toRotationMatrix()
            #ori = i.state.se3[1].toRotationMatrix()*ore_v
            for i in range(h):
                for j in range(w):
                    Surf=getSurface(r,Vector2(b1+j*wStep,a+i*hStep-math.pi*0.5))
                    Surf=Rot*Surf+pos
                    vertices.append(Surf)
    print('POINTS   ',str(len(vertices)),' float',file=f)
    for i in vertices:
        print(i[0],' ',i[1],' ',i[2],file=f)
    print('POLYGONS ',str(point_counter*(h-1)*w),' ',str(point_counter*(h-1)*w*5),file=f)
    for i in range(point_counter):
        start_index=i*w*h
        m=0
        for mm in range(h-1):
            n=0
            for nn in range(w-1):
                print('4 ',str(start_index+m*w+n),' ',str(start_index+m*w+n+1),' ',str(start_index+(m+1)*w+n+1),' ',str(start_index+(m+1)*w+n),' ',file=f)
                n+=1
            print('4 ',str(start_index+m*w+n),' ',str(start_index+m*w),' ',str(start_index+(m+1)*w),' ',str(start_index+(m+1)*w+n),' ',file=f)
            m+=1
    if(model==0):
        print('writing sphere into a VTK file succusful!')
        return True
    print('CELL_DATA ',point_counter*(h-1)*w,file=f)#all polygons
    print('SCALARS ','rotation_scalars ','float ','1',file=f)#all
    print('LOOKUP_TABLE ','default',file=f)
    
    rotations=[]
    if not len(thetas)==len(ids): 
        print("number of particles not equal")
        return False
    for i in ids:
        rotations.append(thetas[i])
    _maxRotation_angle=max(rotations)
    print('max particle rotation is ',_maxRotation_angle)
    #maxRotation_angle=math.pi/6#最大旋转角定义
    for i in rotations:#rotations  旋转角度值
    	for j in range((h-1)*w):
    		#print(str(i/maxRotation_angle*7),file=f)
    		print(str(i),file=f)
    f.close()
    print('writing sphere into a VTK file succusful!')
    return True

def getSurface(R,phi):#获取球表面点
    x=R*math.cos(phi[1])*math.cos(phi[0])
    y=R*math.cos(phi[1])*math.sin(phi[0])
    z=R*math.sin(phi[1])
    return Vector3(x,y,z)
    
def getSurface_cylinder(R,pso,vector):#获取圆柱表面点
    x=R*math.cos(phi[1])*math.cos(phi[0])
    y=R*math.cos(phi[1])*math.sin(phi[0])
    z=R*math.sin(phi[1])
    return Vector3(x,y,z)

def save_sphereRotationinfo(path,index):
    "output particle info to a file"
    step=O.iter
    middelpath=''
    if(index==0):
        middelpath='shearBefore'
    else:
        middelpath='shearAfter'
    #path="/home/hcl/桌面/1/savedata"
    #filename=path+'/particleRota_'+middelpath+str(step)+'.txt' #step 的识别? str(step)
    filename=path+'/particleRota_'+middelpath+'.txt'
    #filename=path+'/particlevtk_'+str(step)+'.vtk'
    f = open(filename,'w')
    print("###Info of each particle:id pos x y z ori w x y z",file=f)
    for i in O.bodies:
        if isinstance(i.shape,Sphere):#focus on superball
            pos = i.state.pos
            ori=i.state.se3[1]
            print(i.id,pos[0],pos[1],pos[2],ori[3],ori[0],ori[1],ori[2],file=f)#ori[3] 实部 0 1 2 虚部
    f.close()

def exportSphere_Rotation(path,num,axis,step=0):
    path = os.path.dirname(path+'/')#	返回文件路径 ==path
    if not os.path.exists(path):#not exists
        print ("The typed path dose not exist!")
        return False
    particlerotefile_Before = os.path.join(path,'particleRota_shearBefore'+'.txt')#目录和文件名合成一个路径  step 的识别? str(step)
    particlerotefile_After = os.path.join(path,'particleRota_shearAfter'+'.txt')#

    if not (os.path.isfile(particlerotefile_Before) and os.path.isfile(particlerotefile_After) ):
        print ("Input files are missing! Please check the path of particleinfo*.txt, contactinfo*.txt and wallinfo*.txt is set correctly.")
        return False
    Sphere_Rotation(path,particlerotefile_Before,particlerotefile_After,num,axis,comment="comment")#num slices 颗粒的纬度数量
    return True

def saveAndOut_Rotation(num,axis,model=0):# 0 before shear  1 after shear
    path="/home/hcl/桌面/1/savedata"
    if model==0:
        save_sphereRotationinfo(path,model);
    else:
        save_sphereRotationinfo(path,model);
        exportSphere_Rotation(path,num,axis)
    
def Sphere_Rotation(path,protefile_b,protefile_a,num,axis,comment="comment"):
    ore_v = Vector3(0,0,0)
    if axis==0:
        ore_v[0]=1
    elif axis==1:
        ore_v[1]=1
    elif axis==2:
        ore_v[2]=1
    print('axix= ',ore_v)
    f_b = open(protefile_b,'r')#颗粒剪切前旋转信息文件
    f_a = open(protefile_a,'r')#颗粒剪切前旋转信息文件
    # output file
    fName = path +'particleRotation'+'.txt'
    #outRotateFile = open(fName, 'w')
    lines_b = f_b.readlines()[1:]########## the range may need to be changed
    lines_a = f_a.readlines()[1:]
    f_b.close()
    f_a.close()
    #open the original file including the info of all balls
    #print("###Info of  particle rotation angle after shear: id theta",file=outRotateFile)
    if not len(lines_b)==len(lines_a):#not exists
        print ("The length ori of two file is not compare!")
        return False
    thetas=dict()
    nBodies=0
    for lb,la in zip(lines_b,lines_a):
        lb = lb.lstrip()#截掉端空格
        _lb = lb.split(' ')#按中间空格分隔组
        [ori3,ori0,ori1,ori2]=[float(i) for i in _lb[4:]]
        #print('ori model :',(ori3**2+ori0**2+ori1**2+ori2**2)**0.5)
        id1=int(_lb[0])
        Q_before=Quaternion(ori3,ori0,ori1,ori2)
        la = la.lstrip()
        _la = la.split(' ')
        [ori3,ori0,ori1,ori2]=[float(i) for i in _la[4:]]
        #print('ori model :',(ori3**2+ori0**2+ori1**2+ori2**2)**0.5)
        id2=int(_la[0])
        if not id1==id2:
            print("id of particle is not compared")
            return False
        Q_after=Quaternion(ori3,ori0,ori1,ori2)
        diff_theta=Q_after*Q_before.conjugate()
        #dtheta = diff_theta.toRotationMatrix()*ore_v#3*1
        #thetas[id1]=dtheta[axis]# 实部计算角度
        w_r=diff_theta[3]
        y_r=diff_theta[1]
        thetay=2 * np.arctan2(y_r, w_r)
        if thetay<0:
        	thetay=0
        thetas[id1]=thetay# 实部计算角度
        #theta_b = Q_before.toRotationMatrix()*ore_v#3*1
        #theta_a = Q_after.toRotationMatrix()*ore_v#3*1
        #thetas[id1]=theta_a[1]-theta_b[1]# 实部计算角度
        #thetas[id1]=2*math.acos(diff_theta[3])# 实部计算角度
        nBodies += 1

    out_particle_vtk(num,thetas,model=1)

def compute_rote_theta(q1,q2): #输入旋转轴以及角度 计算旋转向量
	RotaTheta=0.5*theata
	QuaternionQ=Quaternion(math.cos(RotaTheta),math.sin(RotaTheta)*q[0],math.sin(RotaTheta)*q[1],math.sin(RotaTheta)*q[2])
	QuaternionQ_conjugate=Quaternion(math.cos(RotaTheta),-math.sin(RotaTheta)*q[0],-math.sin(RotaTheta)*q[1],-math.sin(RotaTheta)*q[2])
	QuaternionP=Quaternion(0,p[0],p[1],p[2])
	q1=QuaternionQ*QuaternionP
	q2=q1*QuaternionQ_conjugate
	v=(q2[0],q2[1],q2[2])
	return v
#********************************************************************************************************

#******************************************************************************************************** 各向异性 组构张量计算
def vectordot(v):#单位向量
    temp=np.zeros([3,3])
    modv=math.sqrt(sum([i**2. for i in v]))#==1
    if modv!=0:
        v=[i/modv for i in v]
    for i in range(3):
        for j in range(3):
            temp[i,j]=v[i]*v[j]
    return temp

def anisotropyfromRAW(vectorslist):
    subtensor=list()
    count=len(vectorslist)#接触数量
    subtensor = [vectordot(v/v.norm()) for v in vectorslist]
    fabtensor=np.zeros([3,3])
    fabtensor=sum(subtensor)/count
    b,c=np.linalg.eig(fabtensor)#返回特征值和特征向量
    return math.sqrt((b[0]-b[1])**2.+(b[0]-b[2])**2.+(b[1]-b[2])**2.)/math.sqrt(2.)# 组构偏量

def output_fabrics():#颗粒接触方向各向异性
    #normal contacts
    vlist_nc = list()
    #branch vectors
    vlist_bv = list()
    for i in O.interactions:
        if not isinstance(i.geom,ScGeom):
            continue
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <10:#wall or box
            continue
        if fn.norm()>0.0:
            vlist_nc.append(fn)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos# 2->1
            vlist_bv.append(p)
    #calculate the anisotropy
    ani_nc = anisotropyfromRAW(vlist_nc)#接触力组构张量各项异性大小
    ani_bv = anisotropyfromRAW(vlist_bv)#接触方向组构张量各项异性大小
    return (ani_nc,ani_bv)
    
def output_thread_fabrics(threshold=0.1):#柔索接触力方向的各向异性
    #normal contacts
    
    vlist_nc = list()
    minthreshold=threshold
    for i in O.threads:
    	fn=i.get_f()
    	if fn>=minthreshold :
    		vlist_nc.append(i.geom1.normal)
    		vlist_nc.append(i.geom2.normal)
			
    #calculate the anisotropy
    ani_nc=0.
    if len(vlist_nc)==0:
    	return 0.
    else:
    	ani_nc = anisotropyfromRAW(vlist_nc)#接触力组构张量各项异性大小
    
    return ani_nc
 #***************************************************************
  #***************************************************************  柔索段拉伸状态评估
def computeVectorsCoplanarError(vectors,tolerance=0.01):
	if len(vectors)<3:
		print("max vector is three!!")
		return False
	v1=vectors[0]
	v2=vectors[1]
	normal_v=np.cross(v1,v2)
	#normal_v.nomalized();
	vlen=len(vectors);
	errors=[]
	for i in range(vlen-2):
		v3=vectors[i+2]
		#v3.nomalized()
		errors.append(abs(np.dot(normal_v,v3)))
		
	return [errors,np.average(errors)]
		
def getThreadsIds(minthreshold,maxthreshold=10):#设定力阈值柔索段
	#num=0
	ff=dict()
	for i in O.threads:
		fn=i.get_f()
		if fn>=minthreshold and fn<maxthreshold:
			ff[i.id]=i.get_f()
		#num+=1
	return ff
	
def Vdeformation():
	vd=dict()
	for i in O.threads:
		if i.length<i.length0 and i.get_f()<=0:
			vd[i.id-1]=[i.length0,i.length]
	return vd
def getThreadsGeoms(_id):#获取接触法向量
	Vectors=[]
	if _id<0 and _id>=len(O.threads):
		print("index out of range!!!")
		return False
	geoms=O.threads[_id].get_geoms()
	Vnum=len(geoms)
	for i in range(Vnum):
		Vectors.append(geoms[i].normal)
	return Vectors
	
#***************************
#***************************************************
def delete_linenode():#删除柔索及节点
	O.threads.clear()
	for i in O.engines:
		if(isinstance(i,ThreadForceController)):
			i.dead=True
			break
	num_del=0
	num_del2=0
	num_del3=0
	#for t in O.threads:
		#i.erase
	for b in O.bodies:
		if isinstance(b.shape,LineConnection):
			O.bodies.erase(b.id)
			num_del2+=1
	for b in O.bodies:	
		if isinstance(b.shape,LineNode) and not b.isClumpMember:
			#O.bodies.erase(b.id,False)
			O.bodies.erase(b.id)
			num_del+=1
	for i in O.bodies: 
		if isinstance(i.shape,LineNode): 
			O.bodies.erase(i.id)
			num_del3+=1
	print(num_del2,"lineconnection has been deleted!")
	print(num_del,"lineNode has been deleted!")
	print(num_del3,"after lineNode has been deleted!")
	
def quiet_system():#颗粒速度置零
    for i in O.bodies:
        i.state.vel=(0.,0.,0.)
        i.state.angVel=(0.,0.,0.)
#**************************************
#**************************************配位数
def CoordinateNum2():
	snum=0
	for i in O.bodies:
		if isinstance(i.shape,Sphere):
			snum+=1	
	intr_num=0
	for i in O.interactions:
		if isinstance(i.geom,ScGeom):
			intr_num+=1
	return intr_num*2/snum
#***********************	
#**************************************柔索数量及平均长度统计
def Threads_initial_len():
	num=0
	length0=[]
	if len(O.threads)==0:
		return 0
	for i in O.threads:
		length0.append(i.length0)
		num+=1
	average_l=sum(length0)/num
	minl=min(length0)
	maxl=max(length0)
	dx=maxl-minl
	'''FC=dict()   #store frenquency count 配位数频率分布
    	for i in CNlist:
    	if i in FC.keys():
    	 FC[i]=FC[i]+1
        	else:
            	FC[i]=1
  	  return [FC,CNlist]'''
	return -average_l*1000
	
def Threads_info(threshold=0.1):#*******************柔索作用占比 平均拉力  最大拉力  平均应变
	totalnum=len(O.threads)
	minthreshold=threshold
	filter_num=0
	fns=[]
	ids=[]
	for i in O.threads:
    		fn=i.get_f()
    		if fn>=minthreshold :
    			filter_num+=1
    			fns.append(fn)
    			ids.append(i.id-1)#threads id-1== O.threads[id']
	ratio=filter_num*1.0/totalnum
	average_fn=sum(fns)/len(fns)
	max_fn=max(fns)
	min_fn=min(fns)
	strain = [abs(O.threads[i].length-O.threads[i].length0)/abs(O.threads[i].length0) for i in ids]
	average_strain=sum(strain)/len(strain)
	return [filter_num,totalnum,ratio,average_fn,max_fn,average_strain]
#***********************************改变接触切向刚度振动固结
def change_ks():#循环固结后			
	for i in O.materials:
		if i.label=='spheremat' or i.label=='frictionless':
			i.poisson=10
	for i in O.interactions:
		if isinstance(i.geom,ScGeom):
			i.phys.ks=i.phys.ks/0.4*10
def change_ks_to_origin():#大步长振动后			
	for i in O.materials:
		if i.label=='spheremat' or i.label=='frictionless':
			i.poisson=0.4
	for i in O.interactions:
		if isinstance(i.geom,ScGeom):
			i.phys.ks=i.phys.ks/10*0.4	
	
def two_dimention_Dhistogram_initial(threshold,num=20):#threshold 接触力阈值,num 分段数量(0-180度) 导出接触力直方图统计
	intrs=[] 
	for i in O.interactions:
		if isinstance(i.geom,ScGeom) and i.isReal:
			minid=min(i.id1,i.id2)
			if minid<10:
				continue
			intrs.append(i)
	filt_inters_num=0
	normal_force=[]
	director_theta=[]
	vector_x=Vector3(1,0,0)#对应剪切面与xy平面平行
	vector_z=Vector3(0,0,1)
	for i in intrs:
		Fn=(i.phys.normalForce).norm()
		if Fn>threshold:
			filt_inters_num+=1
			normal_force.append(Fn)
			vector=i.geom.normal
			x=vector[0]
			z=vector[2]
			if z==0. and x>0.:
				director_theta.append(0.)
			if z==0. and x<0.:
				director_theta.append(math.pi)
			if z!=0:
				#print("x/z= ",x,"divide",z)
				Sign=np.sign(x)/np.sign(z)
				director_theta.append(math.acos(Sign*math.sqrt(x**2)/math.sqrt(x**2+z**2)))#0--pi
	#统计频率
	maxFn=max(normal_force)
	segment = np.array([(i+1)*180/num/180*math.pi  for i in range(num)])	#0-171	9-180
	FC=dict()   #store frenquency count 方向角度的频率分布  1-20 个单位  
	for i in director_theta:
		for j in range(num):
			if i<=segment[j]:
				if j+1 in FC.keys():
					FC[j+1]=FC[j+1]+1
				else:
					FC[j+1]=1
				break
			#break
	for j in range(num):
		if j+1 not in FC.keys():
			FC[j+1]=0
			
	path=os.getcwd()
	step=O.iter
	path= os.path.dirname(path+'/'+'contact_FC_step_'+str(step)+'/')#'/'自动删除
    	#f = open(inputfile,'r')#颗粒信息文件
    	# output file
    	#print("path="+path)#
	if not os.path.exists(path):#not exists
		os.makedirs(path)
	fName = path +'/contactForce_FC'+'.txt'
	outContactFile = open(fName, 'w')
	for i in director_theta:
		print(str(round(i,4)),file=outContactFile)
		print(str(round((i+math.pi),4)),file=outContactFile)
	outContactFile.close()	
	print("wrtting to txt file over!!!")
	return len(director_theta)
	
def two_dimention_Dhistogram(force_chain_num,num=20):#threshold 接触力阈值,num 分段数量(0-180度) 导出接触力直方图统计
	if force_chain_num==0:
		return False
	intrs=[] 
	for i in O.interactions:
		if isinstance(i.geom,ScGeom) and i.isReal:
			minid=min(i.id1,i.id2)
			if minid<10:
				continue
			intrs.append(i)
	filt_inrs=[]
	#force_chain_num
	normal_force=dict()
	countnum=0
	for i in intrs:
		Fn=(i.phys.normalForce).norm()
		normal_force[countnum]=Fn
		countnum+=1
	normal_force_sort=sorted(normal_force.items(),key=lambda x:(-x[1],x[0]))
	normal_force_sort=normal_force_sort[:force_chain_num]
	for i in range(force_chain_num):
		filt_inrs.append(intrs[normal_force_sort[i][0]])
		
	#filt_inters_num=0
	#normal_force=[]
	director_theta=[]
	vector_x=Vector3(1,0,0)#对应剪切面与xy平面平行
	vector_z=Vector3(0,0,1)
	for i in filt_inrs:
		vector=i.geom.normal
		x=vector[0]
		z=vector[2]
		if z==0. and x>0.:
			director_theta.append(0.)
		if z==0. and x<0.:
			director_theta.append(math.pi)
		if z!=0:
			#print("x/z= ",x,"divide",z)
			Sign=np.sign(x)/np.sign(z)
			director_theta.append(math.acos(Sign*math.sqrt(x**2)/math.sqrt(x**2+z**2)))#0--pi
	#统计频率
	#maxFn=max(normal_force)
	segment = np.array([(i+1)*180/num/180*math.pi  for i in range(num)])	#0-171	9-180
	FC=dict()   #store frenquency count 方向角度的频率分布  1-20 个单位  
	for i in director_theta:
		for j in range(num):
			if i<=segment[j]:
				if j+1 in FC.keys():
					FC[j+1]=FC[j+1]+1
				else:
					FC[j+1]=1
				break
			#break
	for j in range(num):
		if j+1 not in FC.keys():
			FC[j+1]=0
			
	path=os.getcwd()
	step=O.iter
	path= os.path.dirname(path+'/'+'contact_FC_step_'+str(step)+'/')#'/'自动删除
    	#f = open(inputfile,'r')#颗粒信息文件
    	# output file
    	#print("path="+path)#
	if not os.path.exists(path):#not exists
		os.makedirs(path)
	fName = path +'/contactForce_FC'+'.txt'
	outContactFile = open(fName, 'w')
	for i in director_theta:
		print(str(round(i,4)),file=outContactFile)
		print(str(round((i+math.pi),4)),file=outContactFile)
	outContactFile.close()	
	print("wrtting to txt file over!!!")
	return True	
	
	
	
	

