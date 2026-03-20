/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@grenoble-inp.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Line.hpp"
//#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.


//!##################	SHAPES   #####################


//*********************************************************************************************************************************************************12-26
LineNode::~LineNode() { }
YADE_PLUGIN((LineNode));

LineConnection::~LineConnection() { }
YADE_PLUGIN((LineConnection));

LineNodeGeom::~LineNodeGeom() { }
YADE_PLUGIN((LineNodeGeom));

LineNodeGeom2::~LineNodeGeom2() { }
YADE_PLUGIN((LineNodeGeom2));

NodeScGeom2::~NodeScGeom2() { }
YADE_PLUGIN((NodeScGeom2));

//**************************1
void LineNode::addConnection(shared_ptr<Body> GC) { ConnList.push_back(GC); }

//**************************2
Vector3r LineConnection::getSegment() const
{
	if (!periodic) return node2->state->pos - node1->state->pos;
	//else
	const Scene* scene = Omega::instance().getScene().get();
	return node2->state->pos + scene->cell->hSize * cellDist.cast<Real>() - node1->state->pos;
}

Real LineConnection::getLength() const { return getSegment().norm(); }

Vector3r& LineNodeGeom::rotate(Vector3r& shearForce) const
{
	// approximated rotations
	shearForce -= shearForce.cross(orthonormal_axis);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	//shearForce -= normal.dot(shearForce)*normal;
	return shearForce;
}

void LineNodeGeom::precompute(
	        const State&                   rbp1,
	        const State&                   rbp2,
	        const Scene*                   scene,
	        const shared_ptr<Interaction>& /*c*/,
	        const Vector3r&                currentNormal,
	        bool                           isNew)
{
	if (!isNew) {orthonormal_axis = normal.cross(currentNormal);} else orthonormal_axis = Vector3r::Zero();
	//Update contact normal
	normal = currentNormal;
	//Precompute shear increment
	Vector3r relativeVelocity = getIncidentVel(&rbp1, &rbp2);
	//keep the shear part only
	relativeVelocity = relativeVelocity - normal.dot(relativeVelocity) * normal;
	shearInc         = relativeVelocity * scene->dt;
}
Vector3r LineNodeGeom::getIncidentVel(const State* rbp1, const State* rbp2) const
{
	Vector3r relativeVelocity=rbp2->vel - rbp1->vel;
	return relativeVelocity;
}
//**************************3  IGeom Functors ************************//

//!			o-o
bool Ig2_LineNode_LineNode_LineNodeGeom::go(
	const shared_ptr<Shape>&       /*cm1*/,
	const shared_ptr<Shape>&       cm2,
	const State&                   state1,
	const State&                   state2,
	const Vector3r&                /*shift2*/,
	const bool&                    force,
	const shared_ptr<Interaction>& c)
{
	//shared_ptr<Body> Ln1 = Body::byId(c->getId1(), scene);
	//shared_ptr<Body> Ln2 = Body::byId(c->getId2(), scene);
	//int id1 =c->getId1();
	//int id2 =c->getId2();
	//LineNode* LN2=YADE_CAST<LineNode*>(Ln2->shape.get());
	//LineNode*      LN1 = YADE_CAST<LineNode*>(cm1.get());
	LineNode*      LN2 = YADE_CAST<LineNode*>(cm2.get());//接触的球
	Vector3r   pos1 = state1.pos;
	Vector3r   pos2 = state2.pos;
	shared_ptr<LineNodeGeom> scm;
	//const Sphere *s1 = static_cast<Sphere*>(cm1.get()), *s2 = static_cast<Sphere*>(cm2.get());
	Vector3r      normal = pos2  - pos1;//方向距离向量1-2
	bool      isNew = !c->geom;
	if(isNew&&(!force)){
		//if(!c->isReal() && !force)Real penetrationDepthSq = pow((s1->radius + s2->radius), 2) - normal.squaredNorm();
		//if(penetrationDepthSq<0)return false;
		return false;
	}
	if (!isNew) scm = YADE_PTR_CAST<LineNodeGeom>(c->geom);
	else {
		scm     = shared_ptr<LineNodeGeom>(new LineNodeGeom());
		c->geom = scm;
	}
	Real norm = normal.norm();
	normal /= norm; // normal is unit vector now
	//Real penetrationDepth=2*LN2->radius-norm;
	Real penetrationDepth=-norm;
	scm->radius1          = LN2->radius;//半径
	scm->radius2          = LN2->radius;//
	scm->normal=normal;
	scm->penetrationDepth=penetrationDepth;
	scm->precompute(state1, state2, scene, c, normal, isNew);//计算剪切位移增量
	if (YADE_PTR_CAST<LineNodeGeom>(c->geom)->connectionBody)
	scm->connectionBody->state->pos = state1.pos;//pos update
	//YADE_PTR_CAST<LineNodeGeom>(c->geom)->connectionBody->state->pos = state1.pos;
	return true;

}

bool Ig2_LineNode_LineNode_LineNodeGeom::goReverse(
	const shared_ptr<Shape>&       cm1,
	const shared_ptr<Shape>&       cm2,
	const State&                   state1,
	const State&                   state2,
	const Vector3r&                shift2,
	const bool&                    force,
	const shared_ptr<Interaction>& c)
{
	return go(cm1, cm2, state2, state1, -shift2, force, c);
}
YADE_PLUGIN((Ig2_LineNode_LineNode_LineNodeGeom));


bool Ig2_LineNode_Sphere_NodeScGeom2::go(
	const shared_ptr<Shape>&       cm1,
	const shared_ptr<Shape>&       cm2,
	const State&                   state1,
	const State&                   state2,
	const Vector3r&                shift2,
	const bool&                    force,
	const shared_ptr<Interaction>& c)
{
	//shared_ptr<Body> Ln1 = Body::byId(c->getId1(), scene);
	//shared_ptr<Body> Ln2 = Body::byId(c->getId2(), scene);
	//int id =c->getId1();//linenode
	//LineNode* LN2=YADE_CAST<LineNode*>(Ln2->shape.get());
	//const Se3r&   se31 = state1.se3;
	//const Se3r&   se32 = state2.se3;//sphere
	
	LineNode*      LN = YADE_CAST<LineNode*>(cm1.get());//接触的节点
	Sphere*        SP = YADE_CAST<Sphere*>(cm2.get());//接触的球
	//int id1 =c->getId1();//linenode
	int id2 =c->getId2();//Sphere
	int sphid1=LN->sphid1;
	int sphid2=LN->sphid2;
	bool &NsIntrs=LN->NsIntrs;
	if(!(sphid1==id2||sphid2==id2)&&NsIntrs) return false;//取消重叠接触
	Vector3r      normal = state1.pos+shift2-state2.pos;//方向距离向量2->1 qiu->node
	Real norm = normal.norm();
	Real penetrationDepth=LN->radius+SP->radius-norm;
	//Real penetrationDepth=-norm;
	shared_ptr<NodeScGeom2> scm;
	bool      isNew = !c->geom;
	if(!c->isReal() && !force&&penetrationDepth<0){if((sphid1==id2||sphid2==id2)&&!isNew)NsIntrs=false;return false;}//(!force)
	if (!isNew) scm = YADE_PTR_CAST<NodeScGeom2>(c->geom);
	else {
		scm     = shared_ptr<NodeScGeom2>(new NodeScGeom2());
		c->geom = scm;
		if(sphid1==id2||sphid2==id2) NsIntrs=true;
	}
	normal /= norm; // normal is unit vector now
	scm->contactPoint     = state1.pos - (LN->radius - 0.5 * penetrationDepth) * normal;
	scm->radius1          = LN->radius;//半径
	scm->radius2          = SP->radius;//
	scm->normal=normal;
	scm->penetrationDepth=penetrationDepth;
	return true;

}

bool Ig2_LineNode_Sphere_NodeScGeom2::goReverse(
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
YADE_PLUGIN((Ig2_LineNode_Sphere_NodeScGeom2));

//Ig2_LineConnection_Sphere_LineNodeGeom2
bool Ig2_LineConnection_Sphere_LineNodeGeom2::go(
	const shared_ptr<Shape>&       cm1,
	const shared_ptr<Shape>&       cm2,
	const State&                   /*state1*/,
	const State&                   state2,
	const Vector3r&                shift2,
	const bool&                    /*force*/,
	const shared_ptr<Interaction>& c)
{
	const Se3r&   se32 = state2.se3;//sphere
	LineConnection*      LC = YADE_CAST<LineConnection*>(cm1.get());//接触的linconection
	Sphere*        SP = YADE_CAST<Sphere*>(cm2.get());//接触的球
	//LineNode*      LN1 = YADE_CAST<LineNode*>(LC->node1->shape.get);//linenode shape
	//LineNode*      LN2 = YADE_CAST<LineNode*>(LC->node2->shape.get);
	Vector3r      pos1=LC->node1->state->pos;
	//Vector3r      pos2=LC->node2->state->pos;
	Vector3r spherePos = se32.position - shift2;
	Vector3r segt      = LC->getSegment();//线段由1-2
	Real     len       = LC->getLength();//线段长度
	Vector3r branch    = spherePos - pos1;
	Real relPos = branch.dot(segt) / (len * len);
	bool                     isNew     = !c->geom;
	if((relPos<0||relPos>1))return false;//错误接触判断
	shared_ptr<LineNodeGeom2> scm;
	if (!isNew) scm = YADE_PTR_CAST<LineNodeGeom2>(c->geom);
	else {
		scm = shared_ptr<LineNodeGeom2>(new LineNodeGeom2());
	}
	Vector3r fictiousPos = pos1 + relPos * segt;//球投影点
	Vector3r branchF     = fictiousPos - spherePos;
	Real     dist        = branchF.norm();//球与连接体中心距离
	if (isNew && (dist > (SP->radius + LC->radius))) return false;
	if (isNew) c->geom = scm;
	scm->radius1              = LC->radius;//
	scm->radius2              = SP->radius;
	Vector3r normal           = branchF / dist;
	scm->normal		   = normal;
	scm->penetrationDepth     = SP->radius + LC->radius - dist;
	scm->relPos               = relPos;
	scm->id3                  = LC->node1->getId();
	scm->id4                  = LC->node2->getId();
	return true;

}

bool Ig2_LineConnection_Sphere_LineNodeGeom2::goReverse(
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
YADE_PLUGIN((Ig2_LineConnection_Sphere_LineNodeGeom2));
//!			O-o

//!##################	Bounds   #####################

void Bo1_LineConnection_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& /*se3*/, const Body* /*b*/)
{
	LineConnection* LC = static_cast<LineConnection*>(cm.get());
	if (!bv) { bv = shared_ptr<Bound>(new Aabb); }
	Aabb*    aabb = static_cast<Aabb*>(bv.get());
	Vector3r O = YADE_CAST<State*>(LC->node1->state.get())->pos;
	Vector3r O2 = YADE_CAST<State*>(LC->node2->state.get())->pos;
	if (!scene->isPeriodic) {
		for (int k = 0; k < 3; k++) {
			aabb->min[k] = min(O[k], O2[k]) - LC->radius;
			aabb->max[k] = max(O[k], O2[k]) + LC->radius;
		}
		return;
	}
	else {
		O = scene->cell->unshearPt(O);
		O2 = scene->cell->unshearPt(O2);
		O2 = O2 + scene->cell->hSize * LC->cellDist.cast<Real>();
		for (int k = 0; k < 3; k++) {
			aabb->min[k] = min(O[k], O2[k]) - LC->radius;
			aabb->max[k] = max(O[k], O2[k]) + LC->radius;
		}
	}
}

YADE_PLUGIN((Bo1_LineConnection_Aabb));
//!##################	LineNodeBounds   #####################2024-3-20
void Bo1_LineNode_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* /*b*/)
{
	LineNode* LN = static_cast<LineNode*>(cm.get());
	if (!bv) { bv = shared_ptr<Bound>(new Aabb); }
	Aabb*    aabb = static_cast<Aabb*>(bv.get());
	Vector3r halfSize = (aabbEnlargeFactor > 0 ? aabbEnlargeFactor : 1.) * Vector3r(LN->radius, LN->radius, LN->radius);
	if (!scene->isPeriodic) {
		aabb->min = se3.position - halfSize;
		aabb->max = se3.position + halfSize;
		return;
	}
}
YADE_PLUGIN((Bo1_LineNode_Aabb));
//*****************
} // namespace yade
