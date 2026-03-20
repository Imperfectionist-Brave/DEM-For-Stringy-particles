/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@imag.fr>         *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LineContactLaw.hpp"
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.


//****************************************12-26
YADE_PLUGIN((Law2_LineNodeGeom_LinePhys_LineMoment)(Law2_LineNodeGeom2_LineFrictPhys_CundallStrack)(Law2_NodeScGeom2_LineFrictPhys_CundallStrack)(LineMat)(LinePhys)(LineFrictPhys)(Ip2_LineMat_FrictMat_LineFrictPhys)(Ip2_LineMat_LineMat_LinePhys));


/********************** LineMat ****************************/
//CREATE_LOGGER(LineMat);
void LineMat::postLoad(LineMat&)
{
	std::cout << "linemat postload" << endl;//******************************************
	//BUG: ????? postLoad is called twice,
	//LOG_TRACE("LineMat::postLoad - update material parameters");

	// compute cross-section area for single wire
	as = pow(diameter * 0.5, 2) * M_PI;

	// check for stress strain curve for single wire
	if (strainStressValues.empty()) return; // uninitialized object, don't do nothing at all
	if (strainStressValues.size() < 2) throw std::invalid_argument("WireMat.strainStressValues: at least two points must be given.");
	if (strainStressValues[0](0) == 0. && strainStressValues[0](1) == 0.)
		throw std::invalid_argument("WireMat.strainStressValues: Definition must start with values greater than zero (strain>0,stress>0)");
}

/********************** Law2_ScGeom6D_LinePhys_LineMoment ****************************/
//CREATE_LOGGER(Law2_ScGeom6D_LinePhys_LineMoment);
bool Law2_LineNodeGeom2_LineFrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	//LOG_TRACE("Law2_ScGeom_WirePhys_WirePM::go - contact law");

	LineNodeGeom2*     geom = YADE_CAST<LineNodeGeom2*>(ig.get());
	LineFrictPhys* phys = YADE_CAST<LineFrictPhys*>(ip.get());

	//const int&    id1 = contact->getId1();
	const int&    id2 = contact->getId2();//sphere
	Real un = geom->penetrationDepth;//穿透深度
	if (un < 0) {
		return false;
	}
	Real Fn=phys->kn*un;
	if(Fn<0)Fn=0;
	phys->normalForce = Fn * geom->normal;//法向力
	scene->forces.addForce(id2, -phys->normalForce);
	scene->forces.addForce(geom->id3, (1 - geom->relPos)*phys->normalForce);
	scene->forces.addForce(geom->id4, geom->relPos*phys->normalForce);
	return true;
}


bool Law2_NodeScGeom2_LineFrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	//LOG_TRACE("Law2_ScGeom_WirePhys_WirePM::go - contact law");

	NodeScGeom2*     geom = YADE_CAST<NodeScGeom2*>(ig.get());
	LineFrictPhys* phys = YADE_CAST<LineFrictPhys*>(ip.get());

	
	//const int&    id2 = contact->getId2();//sphere
	Real& un = geom->penetrationDepth;//穿透深度
	if (un < 0) {
		return false;
	}
	const int&    id1 = contact->getId1();
	Real Fn=phys->kn*un;
	phys->normalForce = Fn * geom->normal;//法向力
	scene->forces.addForce(id1, phys->normalForce);
	//scene->forces.addForce(id2, -phys->normalForce);
	return true;
}


bool Law2_LineNodeGeom_LinePhys_LineMoment::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	

	LineNodeGeom*     geom = YADE_CAST<LineNodeGeom*>(ig.get());
	LinePhys* phys = YADE_CAST<LinePhys*>(ip.get());

	const int&    id1 = contact->getId1();
	const int&    id2 = contact->getId2();
	Real un = geom->penetrationDepth;//穿透深度
	Real D = un - phys->initD;//un<初始 D<0 受拉力 D>0则受压 不设定力
	Real Fn= 0;
	if(D>0){
		Fn=0;
		//Fn=phys->kn*un*geom->normal.norm();
	}
	else{
	 	Fn = phys->k * D;
		if(Fn>phys->normalLimit)return false;
	}
	phys->normalForce = Fn * geom->normal;//法向力
	Vector3r&       shearForce = geom->rotate(phys->shearForce);
	const Vector3r& dus        = geom->shearIncrement();
	shearForce -= phys->ks * dus;//切向力
	Real Fs    = shearForce.norm();
	Real maxFs = phys->shearAdhesion;;
	if(Fs>maxFs){
		maxFs               = maxFs / Fs;
		shearForce *= maxFs;//限制最大的切向力
	}
	/*if(!testStretch){
		if(geom->isFirst){
			if(D<0&&un>=1.1*phys->initD){//un<=phys->initD un>=1.1*phys->initD xiaolali 
				scene->forces.addForce(id1, -phys->normalForce- shearForce);
				scene->forces.addForce(id2, phys->normalForce+ shearForce);
			} 
			else {
				scene->forces.addForce(id1, -shearForce);
				scene->forces.addForce(id2, phys->normalForce+ shearForce);	
			}
		}
	
		else if(geom->isLast){
			if(D<0&&un>=1.1*phys->initD){
				scene->forces.addForce(id1, -phys->normalForce- shearForce);
				scene->forces.addForce(id2, phys->normalForce+ shearForce);
			} 
			else 	{
				scene->forces.addForce(id1, -phys->normalForce- shearForce);
				scene->forces.addForce(id2,  shearForce);
			}
		} 
		else{
			scene->forces.addForce(id1, -phys->normalForce- shearForce);
			scene->forces.addForce(id2, phys->normalForce+ shearForce);
		}
	}else{*/
		if(geom->isFirst){
			scene->forces.addForce(id1, -shearForce);scene->forces.addForce(id2, phys->normalForce+ shearForce);	
		}else if(geom->isLast)
		{
			scene->forces.addForce(id1, -phys->normalForce- shearForce);scene->forces.addForce(id2,  shearForce);
		} 
		else{
			scene->forces.addForce(id1, -phys->normalForce- shearForce);scene->forces.addForce(id2, phys->normalForce+ shearForce);
		}
	//}
	return true;
}





/********************** Ip2_LineMat_LineMat_LinePhys ****************************/
//CREATE_LOGGER(Ip2_LineMat_LineMat_LinePhys);
void Ip2_LineMat_LineMat_LinePhys::go(
	const shared_ptr<Material>& b1 // LineMat
	,
	const shared_ptr<Material>& b2 // LineMat
	,
	const shared_ptr<Interaction>& interaction)
{
	LineMat* sdec1 = static_cast<LineMat*>(b1.get());
	LineMat* sdec2 = static_cast<LineMat*>(b2.get());
	if (interaction->phys) return;
	
	LineNodeGeom*    geom = YADE_CAST<LineNodeGeom*>(interaction->geom.get());
	assert(geom);

	shared_ptr<LinePhys> contactPhysics(new LinePhys());//新建指针
	Real                 initD = geom->penetrationDepth;
	contactPhysics->normalForce = Vector3r::Zero();

	LineMat* mat1 = static_cast<LineMat*>(b1.get());
	LineMat* mat2 = static_cast<LineMat*>(b2.get());

	
	contactPhysics->initD = initD;//初始长度
	
	

	Real          Ea = mat1->young;
	Real          Eb = mat2->young;
	Real          Va = mat1->poisson;
	Real          Vb = mat2->poisson;
	Real          Da = geom->radius1;
	Real          Db = geom->radius2;
	Real          Kn = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses

	[[maybe_unused]]Real Ks;
	if (Va && Vb)
		Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb
		/ (Ea * Da * Va + Eb * Db * Vb); //harmonic average of two stiffnesses with ks=V*kn for each sphere
	else
		Ks = 0;
	Real shearAdhPreCalculated =  math::min(sdec1->shearCohesion, sdec2->shearCohesion);
	contactPhysics->shearAdhesion=shearAdhPreCalculated*pow(math::min(Db, Da), 2);
	contactPhysics->kn = Kn;//法向
	contactPhysics->ks = mat1->ks;//切向
	contactPhysics->k = mat1->k;//3-28
	contactPhysics->normalLimit=normalLimit;//最大拉力限制
	interaction->phys = contactPhysics;//传址
	
}


void Ip2_LineMat_FrictMat_LineFrictPhys::go(
	const shared_ptr<Material>& b1 // LineMat
	,
	const shared_ptr<Material>& /*b2*/ // LineMat
	,
	const shared_ptr<Interaction>& interaction)
{
	LineMat* sdec1 = static_cast<LineMat*>(b1.get());
	//FrictMat* sdec2 = static_cast<LineMat*>(b2.get());
	//assert(sdec2);
	if (interaction->phys) return;
	shared_ptr<LineFrictPhys> contactPhysics(new LineFrictPhys());//新建指针
	contactPhysics->kn = sdec1->kn;//法向
	interaction->phys = contactPhysics;//传址
	
}


//****************************************
} // namespace yade

