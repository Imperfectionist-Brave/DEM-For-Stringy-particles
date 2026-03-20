#include <lib/high-precision/Constants.hpp>
//#include <core/Clump.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/LineContactLaw.hpp>
//#include <pkg/dem/ScGeom.hpp>
#include <assert.h>
//#include <preprocessing/dem/Shop.hpp>
#include<fstream>
#include "ThreadForceController.hpp"
#ifdef YADE_OPENMP
	#include<omp.h>
#endif

namespace yade {
	YADE_PLUGIN((ThreadForceController));
	ThreadForceController::~ThreadForceController() { }
	//Real NodeForceController1::ComputeUnbalancedForce(bool maxUnbalanced) { return Shop::unbalancedForce(maxUnbalanced, scene); }
	void ThreadForceController::action() {//引擎运行模块
		/*if(!testStretch){
			if(!useallGeom)addNodeForce();
			else addNodeForceAllGeom();
		}
		else{*/
			addNodeForceAllGeom_test();//使用所有柔索段长参与计算
		//}
		
			
	}
	void ThreadForceController::updateNodePos(){
	if (firstRun) {std::cout << "First run, will initiallize for updateNodePos!" << std::endl; initialize(); firstRun = false;}
	const std::lock_guard<std::mutex> lockB(scene->bodies->drawloopmutex);
	const std::lock_guard<std::mutex> lock(scene->interactions->drawloopmutex);
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? math::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
		for (ThreadContainer::iterator it = _First; it != _End; it++) {
			const shared_ptr<Thread>& n = *it;
			n->updatePos();
		}
		
	}
	void ThreadForceController::addNodeForce(){
		scene->forces.sync();
		//const Real& dt = scene->dt;
		if (firstRun) {std::cout << "First run, will initiallize for addSphereForce!" << std::endl; initialize(); firstRun = false;}
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? math::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
		for (ThreadContainer::iterator it = _First; it != _End; it++) {
			
			const shared_ptr<Thread>& n = *it;
			LineNodeGeom*     _geom1 = YADE_CAST<LineNodeGeom*>(n->geom1.get());
			LineNodeGeom*     _geom2 = YADE_CAST<LineNodeGeom*>(n->geom2.get());
			LinePhys* _phys1=YADE_CAST<LinePhys*>(n->phy1.get());
			LinePhys* _phys2=YADE_CAST<LinePhys*>(n->phy2.get());
			Real diff1=_geom1->penetrationDepth - 1.1*n->length01;
			Real diff2=_geom2->penetrationDepth - 1.1*n->length02;
			if(diff1>=0||diff2>=0)continue;
			Real Fn = k * (diff1+diff2)*0.5+ksegment*0.1*0.5*(n->length01+n->length02);
			Vector3r directForce1 = -Fn * _geom1->normal;
			Vector3r directForce2 = Fn * _geom2->normal;
			_phys1->_normalForce=directForce1;
			_phys2->_normalForce=directForce2;
			addForce(scene, n->id1, directForce1);
			addForce(scene, n->id2, directForce2);
		}
	}
	
	void ThreadForceController::addNodeForceAllGeom(){
		scene->forces.sync();
		//const Real& dt = scene->dt;
		if (firstRun) {std::cout << "First run, will initiallize for addSphereForce!" << std::endl; initialize(); firstRun = false;}
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? math::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
		for (ThreadContainer::iterator it = _First; it != _End; it++) {
			const shared_ptr<Thread>& n = *it;
			LineNodeGeom*     geom=NULL;
			LinePhys* phys=NULL;
			Real length=0;
			
			for(int i=0;i<GeomsNum;i++){
				geom = YADE_CAST<LineNodeGeom*>(n->geoms[i].get());
				length+=geom->penetrationDepth;//length0<0
				
			}
			Real diff=length-1.1*n->length0;
			if(diff>=0)continue;
			Real Fn = k * diff/GeomsNum+ksegment*0.1*n->length0/GeomsNum;
			LineNodeGeom*     _geom1 = YADE_CAST<LineNodeGeom*>(n->geom1.get());
			LineNodeGeom*     _geom2 = YADE_CAST<LineNodeGeom*>(n->geom2.get());
			LinePhys* _phys1=YADE_CAST<LinePhys*>(n->phy1.get());
			LinePhys* _phys2=YADE_CAST<LinePhys*>(n->phy2.get());
			Vector3r directForce1 = -Fn * _geom1->normal;//Fn<0
			Vector3r directForce2 = Fn * _geom2->normal;
			_phys1->_normalForce=directForce1;_phys1->Fn=-Fn;
			_phys2->_normalForce=directForce2;_phys2->Fn=-Fn;
			addForce(scene, n->id1, directForce1);
			addForce(scene, n->id2, directForce2);
			for(int i=1;i<GeomsNum-1;i++){
				phys = YADE_CAST<LinePhys*>(n->phys[i].get());
				phys->Fn=-Fn;
			}
		}
		
	}
	
	
	void ThreadForceController::addNodeForceAllGeom_test(){
		scene->forces.sync();
		//const Real& dt = scene->dt;
		
		if (firstRun) {std::cout << "First run, will initiallize for addSphereForce!" << std::endl; initialize(); firstRun = false;}
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? math::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
		for (ThreadContainer::iterator it = _First; it != _End; it++) {
			std::vector<Vector3r> Gnormals;
			Gnormals.resize(GeomsNum);
			const shared_ptr<Thread>& n = *it;
			LineNodeGeom*     geom=NULL;
			LinePhys* phys=NULL;
			
			
			//bool change=false;
			Real &length=n->length;
			length=0;
			for(int i=0;i<GeomsNum;i++){
				geom = YADE_CAST<LineNodeGeom*>(n->geoms[i].get());
				length+=geom->penetrationDepth;//length0<0
				Gnormals[i]=geom->normal;
			}
			Real diff=length-n->length0;
			if(diff>=0){
				for(int i=0;i<GeomsNum;i++){
					phys = YADE_CAST<LinePhys*>(n->phys[i].get());
					if(phys->Fn==0)break;
					else {
						phys->Fn=0;
						phys->_normalForce=Vector3r::Zero();
					}
				}
				continue;//计算下一柔索
			}
			
			//**********  计算是否处于接触拉长状态
			bool stretch =true;
			for(int i=0;i<GeomsNum-1;i++){
				Real dot=Gnormals[i].dot(Gnormals[i+1]);
				if(dot<theta_cos){stretch=false;break;}//柔索段向量的最大转角如何确定 
				
			}
			if(!stretch){//
				for(int i=0;i<GeomsNum;i++){
					phys = YADE_CAST<LinePhys*>(n->phys[i].get());
					if(phys->Fn==0)break;
					else {
						phys->Fn=0;
						phys->_normalForce=Vector3r::Zero();
					}
				}
				continue;
			}
			//*************
			/*
			Vector3r normal=Gnormals[0].cross(Gnormals[1]);
			Real errors=0;//计算是否处于拉长状态
			for(int i=0;i<GeomsNum-2;i++){
				errors+=std::abs(normal.dot(Gnormals[i+2]));
				//Real dot=Gnormals[i].dot(Gnormals[i+1]);
				
			}
			
			if((errors/(GeomsNum-2))>2e-4){//
				for(int i=0;i<GeomsNum;i++){
					phys = YADE_CAST<LinePhys*>(n->phys[i].get());
					if(phys->Fn==0)break;
					else phys->Fn=0;
				}
				continue;
			}*/
			
			Real Fn = n->k * diff;
			LineNodeGeom*     _geom1 = YADE_CAST<LineNodeGeom*>(n->geom1.get());
			LineNodeGeom*     _geom2 = YADE_CAST<LineNodeGeom*>(n->geom2.get());
			LinePhys* _phys1=YADE_CAST<LinePhys*>(n->phy1.get());
			LinePhys* _phys2=YADE_CAST<LinePhys*>(n->phy2.get());
			Vector3r directForce1 = -Fn * _geom1->normal;//Fn<0
			Vector3r directForce2 = Fn * _geom2->normal;
			_phys1->_normalForce=directForce1;_phys1->Fn=-Fn;
			_phys2->_normalForce=directForce2;_phys2->Fn=-Fn;
			addForce(scene, n->id1, directForce1);
			addForce(scene, n->id2, directForce2);
			for(int i=1;i<GeomsNum-1;i++){
				phys = YADE_CAST<LinePhys*>(n->phys[i].get());
				phys->Fn=-Fn;
			}
		}
		
	}
	void ThreadForceController::initialize() {
		//store elements of the flexible wall
		std::cout << "Initializing threads" << std::endl;
		const shared_ptr<ThreadContainer>& TContainer = scene->threads;
		Num = TContainer->size();
		
		assert(Num > 0);
		ThreadContainer::iterator First = TContainer->begin();
		ThreadContainer::iterator End = TContainer->end();
		std::vector<shared_ptr<Thread>> _threads(First, End);
		//save nodeneighbor
		cout<<"threads before size = "<<threads.size()<<endl;//threads.size()=0
		
		for (std::vector<shared_ptr<Thread>>::iterator n = _threads.begin(); n != _threads.end(); ++n) {
				threads.push_back(*n);
		}
		cout<<"threads now size = "<<threads.size()<<endl;
		_First = threads.begin();
		_End = threads.end();
		
		radius=(*_First)->radius;
		if(useallGeom){
			GeomsNum=(*_First)->geoms.size();
			const shared_ptr<Thread>& nf = *_First;
			cout<<"nf->length0=="<<nf->length0<<endl;
			if(!std::isnan(nf->length0)){
				std::cout << "Initializing threads length0 and k has been completed, no need!" << std::endl; 
			}
			else{
				std::cout << "Initializing threads length0 and k !" << std::endl; 
				for (ThreadContainer::iterator it = _First; it != _End; it++) {
					const shared_ptr<Thread>& n = *it;
					n->length0=0;
					LineNodeGeom*     geom=NULL;
					for(int i=0;i<GeomsNum;i++){
						geom = YADE_CAST<LineNodeGeom*>(n->geoms[i].get());
						n->length0+=geom->penetrationDepth;//length0<0
					}
					n->k = young_s*M_PI*pow(radius,2) / (-n->length0);//n->k>0
				}
			}
			
		}
		else {
		
			for (ThreadContainer::iterator it = _First; it != _End; it++) {
			
				const shared_ptr<Thread>& n = *it;
				LineNodeGeom*     _geom1 = YADE_CAST<LineNodeGeom*>(n->geom1.get());
				LineNodeGeom*     _geom2 = YADE_CAST<LineNodeGeom*>(n->geom2.get());
				n->length01=_geom1->penetrationDepth;
				n->length02=_geom2->penetrationDepth;
			}
		
		}
		std::cout << "Initializing threads  completed" << std::endl;
	}

}



