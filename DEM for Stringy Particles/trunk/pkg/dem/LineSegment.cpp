/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@grenoble-inp.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LineSegment.hpp"
#include<core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

	using math::max;
	using math::min; // using inside .cpp file is ok.
	//*********************plugin
	ThreadMat::~ThreadMat() {}
	YADE_PLUGIN((ThreadMat));

	ThreadPhys::~ThreadPhys() {}
	YADE_PLUGIN((ThreadPhys));

	Segment::~Segment() {}
	YADE_PLUGIN((Segment));

	SeScGeom::~SeScGeom() {}
	YADE_PLUGIN((SeScGeom));

	//**************************
	Vector3r Segment::getSegment() const
	{
		//assert((node1!=NULL)&&(node2!=NULL))
		return node2->pos - node1->pos;
	}

	Real Segment::getLength() const { return getSegment().norm(); }
	Real Segment::getlength() const {
		if (incontact)return length;
		else return (node2->pos - node1->pos).norm();
	}
	void Segment::postLoad(Segment&) {
		std::cout<<"postload begin"<<std::endl;
		//if (vertices.size() != 2) { throw std::runtime_error("Segment must have two vertice."); }
		//if (isnan(vertices[0][0])) return;
		if (node1 == NULL || node2 == NULL) { throw std::runtime_error("Segment must have two nodes."); }
		if (radius <= 0) { throw std::runtime_error("radius can not be negative!!!."); }
		length = length0 = getLength();
		area = M_PI * pow(radius, 2);
		std::cout<<"postload middle"<<std::endl;
		if (isFirst || isLast) {
			assert(sphere&&(sphere->shape->getClassName() == "Sphere"));
			
			if (isFirst) {//thred head
				middle = false;
				node = node1;
				Vector3r initial_pos = node1->pos;//globoal
				Vector3r direct = initial_pos - sphere->state->pos;//initial vector
				local_coord = direct;
				id=sphere->getId();
			}
			else {
				middle = false;
				node = node2;
				Vector3r initial_pos = node2->pos;//globoal
				Vector3r direct = initial_pos - sphere->state->pos;//initial vector
				local_coord = direct;
				id=sphere->getId();
			}
			node->blockedDOFs = State::DOF_ALL;
		}
		std::cout<<"postload end"<<std::endl;
	}

	Vector3r Segment::getNormal() const {
		//return getSegment().normalized();
		return (node2->pos - node1->pos).normalized();
	}
	bool Segment::isConected() {
		return sphere != NULL;
	}
	bool Segment::isStretching() {
		return length > length0;
	}
	Vector3r Segment::gettan1()const {
		//if (incontact)return Tangential1;
		//else return getNormal();
		return getNormal();
	}
	Vector3r Segment::gettan2()const {
		if (incontact)return Tangential2;
		else return -getNormal();
	}
	/*Vector3r Segment::getlength()const {
		if (incontact)return Tangential1;
		else return getSegment();
	}*/
	void Segment::updatePos()const {
		assert(node);
		const shared_ptr<State>& St=sphere->state;
		//state->pos=node1->pos;
		Quaternionr& ori=St->ori;
		Matrix3r Rotation = ori.toRotationMatrix();
		Vector3r localpos = Rotation * local_coord;
		node->pos= St->pos+ localpos;
		//se->state->pos=se->node1->pos;
		
	}
	void Segment::updateBodyPos()const {
		//assert(node);
		//const shared_ptr<State>& St=sphere->state;
		state->pos=node1->pos;
		
	}
	void Segment::applyForce(Scene* scene,Vector3r f){
	  //scene = Omega::instance().getScene().get();
		scene->forces.addForce(id, f);
		scene->forces.addTorque(id, (node->pos-sphere->state->pos).cross(f));
	}
	//************************** Iphys Functors ************************//
	void Ip2_ThreadMat_FrictMat_ThreadPhys::go(
		const shared_ptr<Material>& b1 // 
		, const shared_ptr<Material>& //b2 // 
		, const shared_ptr<Interaction>& interaction)
	{
		if (interaction->phys) return;
		//const shared_ptr<FrictMat>& mat1 = YADE_PTR_CAST<FrictMat>(b1);
		const shared_ptr<ThreadMat>& mat1 = YADE_PTR_CAST<ThreadMat>(b1);
		interaction->phys = shared_ptr<ThreadPhys>(new ThreadPhys());
		const shared_ptr<ThreadPhys>& contactPhysics = YADE_PTR_CAST<ThreadPhys>(interaction->phys);
		contactPhysics->kn = mat1->Kn;
		contactPhysics->ks = 0;
	}
	YADE_PLUGIN((Ip2_ThreadMat_FrictMat_ThreadPhys));
	//************************** IGeom Functors ************************//
	//*************************************************1
	bool Ig2_Segment_Sphere_SeScGeom::go(
		const shared_ptr<Shape>& cm1,
		const shared_ptr<Shape>& cm2,
		const State&             /*state1*/,
		const State& state2,
		const Vector3r& shift2,//非周期性为0
		const bool& /*force*/,
		const shared_ptr<Interaction>& c)
	{
		Segment*				 se = YADE_CAST<Segment*>(cm1.get());//接触线的一段
		Sphere*                  sphere = YADE_CAST<Sphere*>(cm2.get());//接触的球
		const State*             sphereSt = YADE_CAST<const State*>(&state2);//球state
		Vector3r                 pos1 = se->node1->pos;
		//[[maybe_unused]]Vector3r                 pos2 = se->node2->pos;
		Vector3r spherePos = sphereSt->pos - shift2;
		bool                     isNew = !c->geom;
		shared_ptr<SeScGeom> scm;
		if (!isNew) scm = YADE_PTR_CAST<SeScGeom>(c->geom);
		else {
			scm = shared_ptr<SeScGeom>(new SeScGeom());
		}
		Real Sradius = sphere->radius;
		Real Seradius = se->radius;
		Vector3r normal = pos1-spherePos;//支向量1
		//Real     len = normal.norm();
		Real     penetrationDepth = Sradius + Seradius-normal.norm();
		
		
		if (isNew&&(penetrationDepth > (Sradius + Seradius))) {return false;}//距离过远不存在接触
		//存在接触的情况 求解虚拟接触节点 fictious_pos1 fictious_pos2
		if (isNew) {
			c->geom = scm;
			scm->id1 = se->node1->getId();//前节点gridnode
			//scm->id2 = se->node2->getId();//后节点gridnode
		}
		
		//scm->radius1 = sphere->radius;//
		//scm->radius2 = gridCo->radius;
		
		//scm->relPos = relPos;//投影位置
		//Vector3r normal = branchF / Dist;//法向 球指向中心线
		scm->penetrationDepth = penetrationDepth;//穿透深度
		scm->normal = normal;
		//scm->relPos = relPos;
		return true;
		//scm->fictiousState.pos = fictiousPos;//虚拟球位置
		//scm->contactPoint = spherePos + normal * (scm->radius1 - 0.5 * scm->penetrationDepth);//接触点
		//scm->fictiousState.vel = (1 - relPos) * gridNo1St->vel + relPos * gridNo2St->vel;
		
	}
	bool Ig2_Segment_Sphere_SeScGeom::goReverse(
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
	YADE_PLUGIN((Ig2_Segment_Sphere_SeScGeom));

	//************************** Law Functors ************************//
	bool Law2_SeScGeom_ThreadPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
	{
		//const int&    id1 = contact->getId1();//Segment id
		const int&    id2 = contact->getId2();
		//const shared_ptr<Body>& b = Body::byId(id1);//for test
		//std::cout << "id1 for b->shape.ClassName() is" << b->shape->getClassName() << std::endl;
		SeScGeom* geom = YADE_CAST<SeScGeom*>(ig.get());
		ThreadPhys* phys = YADE_CAST<ThreadPhys*>(ip.get());
		Real un = geom->penetrationDepth;//使用节点2计算存储力
		/*if (un < 0) {
			if (geom->penetrationDepth1 > 0) {
				phys->normalForce = Vector3r::Zero();
				Vector3r force = phys->normalForce;
				scene->forces.addForce(id1, force);
				return true;
			}
			else return false;
		}*/
		Real Fn = phys->kn * un;
		phys->normalForce = Fn * geom->normal;//
		Vector3r force = phys->normalForce;
		scene->forces.addForce(id2, -force);
		scene->nodeforces.addForce(geom->id1, force);
		//scene->nodeforces.addForce(geom->id2, (geom->relPos)*force);
		return true;
	}

	YADE_PLUGIN((Law2_SeScGeom_ThreadPhys_CundallStrack));

	//**************************************************
	//!##################	Bounds   #####################

	void Bo1_Segment_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& /*se3*/, const Body* /*b*/)
	{
		Segment* Se = static_cast<Segment*>(cm.get());
		if (!bv) { bv = shared_ptr<Bound>(new Aabb); }
		Aabb*    aabb = static_cast<Aabb*>(bv.get());
		Vector3r O = Se->node1->pos;
		Vector3r O2 = Se->node2->pos;
		if (!scene->isPeriodic) {
			for (int k = 0; k < 3; k++) {
				aabb->min[k] = min(O[k], O2[k]) - Se->radius;
				aabb->max[k] = max(O[k], O2[k]) + Se->radius;
			}
			return;
		}
		else {
			O = scene->cell->unshearPt(O);
			O2 = scene->cell->unshearPt(O2);
			O2 = O2 + scene->cell->hSize * Se->cellDist.cast<Real>();
			for (int k = 0; k < 3; k++) {
				aabb->min[k] = min(O[k], O2[k]) - Se->radius;
				aabb->max[k] = max(O[k], O2[k]) + Se->radius;
			}
		}
	}
	YADE_PLUGIN((Bo1_Segment_Aabb));

	//*****************
} // namespace yade
		/*Vector3r segt = se->getSegment();//线段向量由1-2
		Real     len = se->getLength();//线段长度 节点直线距离
		Vector3r branch = spherePos - pos1;//支向量1
		//Vector3r branchN = spherePos - pos2;//支向量2
		for (int i = 0; i < 3; i++) {
			if (math::abs(branch[i]) < 1e-14) branch[i] = 0.0;
			//if (math::abs(branchN[i]) < 1e-14) branchN[i] = 0.0;
		}
		Real relPos = branch.dot(segt) / (len * len);//投影点占segt的比例

		if (isNew&&(relPos < 0 || relPos > 1)) {//不存在与颗粒接触的情况，由前一段或后一段lineconnection处理 继续处理端节点是否碰撞?
			return false;
		}
		//在relPos处于 0-1情况下  根据情况求解虚拟节点的位置
		Vector3r fictiousPos = pos1 + relPos * segt;//球在直线柔索段上的投影点
		Vector3r branchF = fictiousPos - spherePos;//球投影垂直向量 球心指向-直线
		Real     Dist = branchF.norm();//球与连接体中心距离
		Real Sradius = sphere->radius;
		Real Seradius = se->radius;
		if (isNew&&(Dist > (Sradius + Seradius))) {return false;}//距离过远不存在接触
		//存在接触的情况 求解虚拟接触节点 fictious_pos1 fictious_pos2
		if (isNew) {
			c->geom = scm;
			scm->id1 = se->node1->getId();//前节点gridnode
			scm->id2 = se->node2->getId();//后节点gridnode
		}
		
		//scm->radius1 = sphere->radius;//
		//scm->radius2 = gridCo->radius;
		
		//scm->relPos = relPos;//投影位置
		Vector3r normal = branchF / Dist;//法向 球指向中心线
		scm->penetrationDepth = Sradius + Seradius - Dist;//穿透深度
		scm->normal = normal;
		scm->relPos = relPos;
		return true;
		//scm->fictiousState.pos = fictiousPos;//虚拟球位置
		//scm->contactPoint = spherePos + normal * (scm->radius1 - 0.5 * scm->penetrationDepth);//接触点
		//scm->fictiousState.vel = (1 - relPos) * gridNo1St->vel + relPos * gridNo2St->vel;*/
