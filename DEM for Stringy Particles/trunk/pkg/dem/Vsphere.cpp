/*************************************************************************
*  Copyright (C) 2010 by Klaus Thoeni                                    *
*  klaus.thoeni@newcastle.edu.au                                         *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Vsphere.hpp"
#include "../../lib/base/Math.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>


namespace yade { // Cannot have #include directive inside.

	//YADE_PLUGIN((Vsphere)(VsMat)(WirePhys)(Ip2_WireMat_WireMat_WirePhys)(Law2_ScGeom_WirePhys_WirePM));

	Vsphere::~Vsphere() { }
	YADE_PLUGIN((Vsphere));

	VSphereGeom::~VSphereGeom() { }
	YADE_PLUGIN((VSphereGeom));

	VsMat::~VsMat() { }
	YADE_PLUGIN((VsMat));

	VsPhys::~VsPhys() { }
	YADE_PLUGIN((VsPhys));
	
	void VsMat::postLoad(VsMat&)
	{
		std::cout << "VsMat postload" << endl;
	}

	void VSphereGeom::updatefixedpoint(
		const State&                   rbp1,
		const State&                   rbp2) {

		const Se3r&   se31 = rbp1.se3;
		const Se3r&   se32 = rbp2.se3;

		Matrix3r Rotation = (se31.orientation).toRotationMatrix();
		Vector3r localV = Rotation * startV1;
		branch1 = -localV;
		fixedpoint1 = se31.position+ localV;

		Rotation = (se32.orientation).toRotationMatrix();
		localV = Rotation * startV2;
		branch2 = -localV;
		fixedpoint2 = se32.position + localV;
		
	}

	Quaternionr Ig2_Vsphere_Sphere_VSphereGeom::rotationQuaternion(Vector3r& axis, Real& theta) {
		axis.normalize(); // 归一化旋转轴
		Real half_theta = theta / 2.0;
		Real sin_half_theta = sin(half_theta);

		Real w = cos(half_theta);
		Real x = axis(0) * sin_half_theta;
		Real y = axis(1) * sin_half_theta;
		Real z = axis(2) * sin_half_theta;

		return Quaternionr(w, x, y, z);
	}
	//ig2 go 
	bool Ig2_Vsphere_Sphere_VSphereGeom::go(
		const shared_ptr<Shape>&       cm1,
		const shared_ptr<Shape>&       cm2,
		const State&                   state1,
		const State&                   state2,
		const Vector3r&                /*shift2*/,
		const bool&                    force,
		const shared_ptr<Interaction>& c)
	{
		Vsphere*		Vs= YADE_CAST<Vsphere*>(cm1.get());//虚拟球
		Sphere*        Sp = YADE_CAST<Sphere*>(cm2.get());//接触的球
		Real R1 = Vs->radius;
		Real R2 = Sp->radius;
		Vector3r   pos1 = state1.pos;
		Vector3r   pos2 = state2.pos;
		//Vector3r   segment = pos2 - pos1;
		shared_ptr<VSphereGeom> vsm;
		bool      isNew = !c->geom;
		if (isNew) {
			if(!force)return false;
			else {
				vsm = shared_ptr<VSphereGeom>(new VSphereGeom());
				c->geom = vsm;
				Vector3r   segment = pos2 - pos1;segment.normalize();
				//Vector3r normal = segment.normalize();
				vsm->startV1 = (R1 + Tradius) * segment;
				vsm->startV2 = -(R2 + Tradius) * segment;
				vsm->fixedpoint1 = pos1 + vsm->startV1;//how fixedpoint update?
				vsm->fixedpoint2 = pos2 + vsm->startV2;
			}
		}
		vsm = YADE_PTR_CAST<VSphereGeom>(c->geom);
		
		//update fixedpoint
		vsm->updatefixedpoint(state1, state2);
		//qiuejieqiedianweizhi
		Vector3r   segt = vsm->fixedpoint2 - vsm->fixedpoint1;
		//Vector3r   branch1 = pos1-vsm->fixedpoint1;
		//Vector3r   branch2 = pos2-vsm->fixedpoint2;
		Vector3r branch1=vsm->branch1 ;
		Vector3r branch2=vsm->branch2 ;
		//vsm->branch2 = branch2;
		Real len = segt.norm();
		Real relPos1 = branch1.dot(segt) / (len * len);
		Real relPos2 = branch2.dot(-segt) / (len * len);
		if (relPos1 <=0) {
			if (relPos2 <= 0) {
				vsm->TanPoint = Vector3r::Zero();
				vsm->normal1 = segt;vsm->normal1.normalize();
				vsm->normal2 = -segt;vsm->normal2.normalize();
				vsm->minlength = len;
				//vsm->contactPoint = len;
			}
			else {//sphere rotation
				Vector3r   segt2 = vsm->fixedpoint1 - pos2;
				Real D = segt2.norm();
				Real theta = acos((R2 + Tradius) / D);
				Vector3r axis = -segt2.cross(branch2);axis.normalize();
				Vector3r Rv = segt2;Rv.normalize();
				Quaternionr q = rotationQuaternion(axis, theta);//Matrix3r Rotation = q.toRotationMatrix();
				//Vector3r _Rv = Rotation * Rv;//球心到切点
				Vector3r _Rv = q * Rv;//球心到切点
				_Rv.normalize();
				vsm->contactPoint=vsm->TanPoint = pos2 + (R2 + Tradius)*_Rv;
				Vector3r lenline = vsm->TanPoint - vsm->fixedpoint1;
				Real length1 = lenline.norm();
				vsm->normal1 = lenline;vsm->normal1.normalize();
				vsm->normal2 = axis.cross(branch2); vsm->normal2.normalize();
				Real thetaArc = acos(_Rv.dot(-branch2) / (R2 + Tradius));
				Real length2 = (R2 + Tradius)*thetaArc;
				vsm->minlength = length1+ length2;
				
			}
		}
		else if (relPos1 > 0) {
			if (relPos2 > 0)cout << "geom has two tanPoints ,need to fix ig2" << endl;
			else {//relPos2<0
				Vector3r   segt2 = vsm->fixedpoint2 - pos1;
				Real D = segt2.norm();
				Real theta = acos((R2 + Tradius) / D);
				Vector3r axis = segt2.cross(branch1); axis.normalize();
				Vector3r Rv = segt2;Rv.normalize();
				Quaternionr q = rotationQuaternion(axis, theta);//Matrix3r Rotation = q.toRotationMatrix();
				//Vector3r _Rv = Rotation * Rv;//球心到切点
				Vector3r _Rv = q * Rv;//球心到切点//
				_Rv.normalize();
				vsm->contactPoint = vsm->TanPoint = pos1 + (R1 + Tradius)*_Rv;
				Vector3r lenline = vsm->TanPoint - vsm->fixedpoint2;
				Real length1 = lenline.norm();
				vsm->normal2 = lenline;vsm->normal2.normalize();
				vsm->normal1 = -axis.cross(branch1); vsm->normal1.normalize();
				Real thetaArc = acos(_Rv.dot(-branch1) / (R1 + Tradius));
				Real length2 = (R1 + Tradius)*thetaArc;
				vsm->minlength = length1 + length2;
			}
		}

		//Real penetrationDepth = (vsm->fixedpoint2 - vsm->fixedpoint1).norm();
		vsm->radius1 = R1;
		vsm->radius2 = R2;
		//vsm->penetrationDepth = penetrationDepth;
		return true;
	}

	bool Ig2_Vsphere_Sphere_VSphereGeom::goReverse(
		const shared_ptr<Shape>&       cm1,
		const shared_ptr<Shape>&       cm2,
		const State&                   state1,
		const State&                   state2,
		const Vector3r&                shift2,
		const bool&                    force,
		const shared_ptr<Interaction>& c)
	{
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	}
	YADE_PLUGIN((Ig2_Vsphere_Sphere_VSphereGeom));

	//ip2
	void Ip2_VsMat_FrictMat_VsPhys::go(
		const shared_ptr<Material>& b1 // LineMat
		,
		const shared_ptr<Material>& /*b2*/ // LineMat
		,
		const shared_ptr<Interaction>& interaction)
	{
		VsMat* mat1 = static_cast<VsMat*>(b1.get());
		//FrictMat* mat2 = static_cast<FrictMat*>(b2.get());
		
		if (interaction->phys) return;

		//LOG_TRACE("Ip2_LineMat_LineMat_LinePhys::go - create interaction physics");

		VSphereGeom*    geom = YADE_CAST<VSphereGeom*>(interaction->geom.get());
		assert(geom);

		shared_ptr<VsPhys> contactPhysics(new VsPhys());//新建指针
		Real                 initD = geom->minlength;//>0
		contactPhysics->normalForce = Vector3r::Zero();
		contactPhysics->initD = initD;//初始长度
		cout<<"mat1->k="<<mat1->k<<endl;
		contactPhysics->k = mat1->k;//3-28
		contactPhysics->normalLimit = normalLimit;//最大拉力限制
		interaction->phys = contactPhysics;//传址

	}
	YADE_PLUGIN((Ip2_VsMat_FrictMat_VsPhys));

	//law2
	bool Law2_VSphereGeom_VsPhys_VsPM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
	{

		VSphereGeom*    geom = static_cast<VSphereGeom*>(ig.get());
		VsPhys*  phys = static_cast<VsPhys*>(ip.get());
		const int& id1 = contact->getId1();
		const int& id2 = contact->getId2();
		//Body*      b1 = Body::byId(id1, scene).get();
		//Body*      b2 = Body::byId(id2, scene).get();

		Real displN = geom->minlength;//当前最小长度 >0
		//vector<Vector2r>& DFValues = phys->displForceValues;
		//vector<Real>&     kValues = phys->stiffnessValues;
		Real D = displN - phys->initD;//
		Real Fn=0;//
		if (D <= 0);
		else {
			Fn = phys->k*D;//D<0 Fn<0
			if (Fn > phys->normalLimit)return false;
		}
		phys->_normalForce1 = Fn * geom->normal1;//法向力1
		phys->_normalForce2 = Fn * geom->normal2;//法向力2
		scene->forces.addForce(id1, phys->_normalForce1);
		scene->forces.addForce(id2, phys->_normalForce2);
		scene->forces.addTorque(id1,-geom->branch1.cross( phys->_normalForce1));//矢径叉乘力施加力矩
		scene->forces.addTorque(id2,-geom->branch2.cross( phys->_normalForce2));//矢径叉乘力施加力矩
		return true;
	}
	YADE_PLUGIN((Law2_VSphereGeom_VsPhys_VsPM));
} // namespace yade
