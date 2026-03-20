#include <lib/high-precision/Constants.hpp>
//#include <core/Clump.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
//#include <pkg/dem/ScGeom.hpp>
#include <assert.h>
#include <preprocessing/dem/Shop.hpp>
#include<fstream>
#include "NodeForceController.hpp"

namespace yade {
	YADE_PLUGIN((NodeForceController));
	NodeForceController::~NodeForceController() { }
	Real NodeForceController::ComputeUnbalancedForce(bool maxUnbalanced) { return Shop::unbalancedForce(maxUnbalanced, scene); }
	void NodeForceController::action() {//引擎运行模块

		const Real& dt = scene->dt;
		if (fill) {
			if (firstRun) {
				initialize();//导入节点 固定膜
				get_gainz1();//填充伺服增益
				free_node();//中间节点不限制自由度
					//包括计算初始高度，体积等函数的模块
			}
			if (dofilling) {
				fillSample();//辅助墙填充试样
				if (scene->iter % echo_interval == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss1();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz1 / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
				}
			}
			else if(moving){
				movetop();
				if (scene->iter % echo_interval == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
				}
			}
		}
		else {
			if (z_servo) {//consolidation 固结
				//if (hydroStrain)//compaction with a constant strain rate in all directions
				//{
				//	hydroConsolidation();//恒定应变速率压实
				//}
				//else {
				//	generalConsolidation();//常规计算速率压实
				//}
				if (firstRun1) {
					std::cerr << "Consolidationing begins..." << std::endl;
					firstRun1 = false;
					get_gainz();//固结伺服增益获取 
				}
				generalConsolidation();//施加围压  一段时间获取伺服参数 get_gainz
				//输出信息
				//if (scene->iter % savedata_interval == 0) { recordData(); }//recordData输出信息函数
				if (scene->iter % echo_interval == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					//getStressStrain();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", Ss_x:" << wsxx / 1000.0 //average stress of the two walls, unit is kPa
						//<< ", Ss_y:" << wsyy / 1000.0
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
				}

			}
			else {//剪切 z_servo关闭
				if (firstRun2) {//consolidation completed normally
					std::cerr << "Shear begins..." << std::endl;
					updateSampleSize();
					height0 = height;
					firstRun2 = false;
					iterate_num = 0;
				}
				iterate_num += 1;
				if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles
					iterate_num = 0;
					//get_gain();//无需获取伺服参数 边界位移控制 
				}
				servo(dt);//集中力施加 包含更新边界 应变率控制 无需获取伺服参数
				//double udz = goalz * dt*height;
				//bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate 0.01
				//top_wall->state->pos[2] += -udz;//z
				//for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
				//	(*b)->selfstate->pos[2] -= udz;//z 
				//}
				//for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
				//	(*b)->selfstate->pos[2] += udz;//z 
				//}
				//updateSampleSize();
				if (log(height0 / height) > target_strain) {//
				//stop shear

					std::cerr << "Shear ends!" << std::endl;
					scene->stopAtIter = scene->iter + 1;//停止模拟
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", Str_z:" << log(height0 / height)//z
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
						<< "\n";
				}
				//output info
				if (scene->iter % echo_interval == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", Str_z:" << log(height0 / height)//z
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
						<< "\n";
				}
			}
		}
	}

	void NodeForceController::initialize() {
		std::cout << "First run, will initiallize!" << std::endl;

		//z轴面积
		z_area = M_PI * pow(wall_radius, 2.0);


		//store elements of the flexible wall
		std::cout << "Initializing FlexCompressionEngine" << std::endl;
		const shared_ptr<NodeNeighborContainer>& NContainer = scene->nodes;//scene->nodes  shared_ptr<NodeNeighborContainer>
		Num = NContainer->size();
		cows = NContainer->cows;
		cols = NContainer->cols;
		assert((cows*cols) == Num);
		BodyContainer::iterator bi1 = scene->bodies->begin();//获取上下边界box
		std::advance(bi1, Num);
		BodyContainer::iterator bi2 = scene->bodies->begin();
		std::advance(bi2, Num + 1);
		//for ( ; bi!=biEnd; ++bi )
	//	{
	//		//if((*bi)->isClump()) continue;//no clump for Superquadrics
	//		const shared_ptr<Body>& b = *bi;
	//		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";
	//		if (b->shape->getClassName()=="Superquadrics"){
	//			const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
	//			particlesVolume += A->getVolume();
	//		}
	//	}
		bottom_wall = *bi1;//bottom wall
		top_wall = *bi2;//top wall 获取边界墙
		boxid1 = bottom_wall->id;
		boxid2 = top_wall->id;
		Real top = top_wall->state->pos[2];
		Real bottom = bottom_wall->state->pos[2];
		height = height0 = top - bottom;
		min = bottom_wall->bound->min;//边界最小角点
		max = bottom_wall->bound->max;//边界最大角点
		max[2] = max[2] + height0;

		BodyContainer::iterator biend = scene->bodies->end();
		bi2++;
		for (; bi2 != biend; ++bi2) {
			const shared_ptr<Body>& b = *bi2;
			if (b->shape->getClassName() == "Box") {
				cout << "getclassname=" << b->shape->getClassName() << endl;
				third_wall = *bi2;
				boxid3 = third_wall->id;
				break;
			}

		}

		NodeNeighborContainer::iterator NeiFirst = NContainer->begin();
		NodeNeighborContainer::iterator NeiEnd = NContainer->end();
		std::vector<shared_ptr<NodeNeighbor>> Ntop(NeiFirst, NeiFirst + cols);
		std::vector<shared_ptr<NodeNeighbor>> Nfree(NeiFirst + cols, NeiEnd - cols);
		std::vector<shared_ptr<NodeNeighbor>> Nbottom(NeiEnd - cols, NeiEnd);
		//save nodeneighbor
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Ntop.begin(); n != Ntop.end(); ++n) {
			top_nodNei.push_back(*n);
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Nfree.begin(); n != Nfree.end(); ++n) {
			free_nodNei.push_back(*n);
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Nbottom.begin(); n != Nbottom.end(); ++n) {
			bottom_nodNei.push_back(*n);
		}
		heightgap = top-top_nodNei[0]->selfstate->pos[2];//间隙

		fix_node();//固定边界墙？
		//free_node();
		iterate_num = 0;
		solve_num = 0;
		firstRun = false;
		//get_gainz();

	}
	/*void NodeForceController::flexible_consolidate() {
		const Real& dt = scene->dt;


	}*/
	void NodeForceController::quiet_system() {
		BodyContainer::iterator bb = scene->bodies->begin();
		BodyContainer::iterator iiEnd = scene->bodies->end();

		for (; bb != iiEnd; ++bb)
		{
			State* state = (*bb)->state.get();
			state->vel = Vector3r(0.0, 0.0, 0.0);
			state->angVel = Vector3r(0.0, 0.0, 0.0);
		}
	}
	void NodeForceController::get_gainz() {//获取伺服参数

		if (!z_servo) { return; }
		avg_zstiffb = 0.0;
		avg_zstifft = 0.0;
		//avg_rstiff = 0.0;//total stiffness of particles contacting with the flexible wall
		InteractionContainer::iterator ii = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();
		for (; ii != iiEnd; ++ii) if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id == boxid1) {//caution: I assume that wall_id is less than 2. top and bottom walls
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstiffb += currentContactPhysics->kn;
				continue;
			}
			else if (id == boxid2) {
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstifft += currentContactPhysics->kn;
				continue;
			}
			//else if ((keepRigid) && (id < numWallFacet)) {//TriElements
			//	FrictPhys* currentContactPhysics =
			//		static_cast<FrictPhys*> (contact->phys.get());
			//	avg_rstiff += currentContactPhysics->kn;
			//}
			gain_zb = gain_alpha * z_area / (avg_zstiffb);//bottom wall 
			gain_zt = gain_alpha * z_area / (avg_zstifft);//top wall

		}
	}
	void NodeForceController::get_gainz1() {//获取伺服参数

		if (!z_servo) { return; }
		avg_zstiffb = 0.0;
		avg_zstifft1 = 0.0;
		//avg_rstiff = 0.0;//total stiffness of particles contacting with the flexible wall
		InteractionContainer::iterator ii = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();
		for (; ii != iiEnd; ++ii) if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id == boxid3) {//caution: I assume that wall_id is less than 2. top and bottom walls
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstifft1 += currentContactPhysics->kn;
				continue;
			}
			else if (id == boxid1) {
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstiffb += currentContactPhysics->kn;
				continue;
			}
			gain_zt1 = gain_alpha * z_area / (avg_zstifft1);//top1 wall
			gain_zb = gain_alpha * z_area / (avg_zstiffb);//top wall

		}
	}
	/*void NodeForceController::servo_cylinder(double dt) {

		consol_ff_cylinder();
	}*/
	void NodeForceController::consol_ff_cylinder() {//free_nodNei  膜节点集中力控制
		// sync thread storage of ForceContainer

		scene->forces.sync();
		//std::vector<shared_ptr<NodeNeighbor>> free_nodNei;
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator N = free_nodNei.begin(); N != free_nodNei.end(); ++N) {
			const shared_ptr<NodeNeighbor>& n = *N;
			Body::id_t id = n->selfId;
			Vector3r Force = n->computeUnitForce();
			//Vector3r Force=n->computeUnitForce();
			//if (id % 150 == 0)
				//std::cout <<id<< "computeUnitForce=" << Force.norm() <<"neighbor num"<<endl;
			addForce(scene, id, Force*wrap_pressure);
		}
	}
	void NodeForceController::fix_node() {//free_nodNei 固定膜颗粒
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nt = top_nodNei.begin(); Nt != top_nodNei.end(); ++Nt) {
			const shared_ptr<NodeNeighbor>& n = *Nt;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nf = free_nodNei.begin(); Nf != free_nodNei.end(); ++Nf) {
			const shared_ptr<NodeNeighbor>& n = *Nf;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nb = bottom_nodNei.begin(); Nb != bottom_nodNei.end(); ++Nb) {
			const shared_ptr<NodeNeighbor>& n = *Nb;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
	}
	void NodeForceController::free_node() {//free_nodNei 固定膜颗粒
		//for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nt = top_nodNei.begin(); Nt != top_nodNei.end(); ++Nt) {
		//	const shared_ptr<NodeNeighbor>& n = *Nt;
		//	//n->selfstate->blockedDOFs = State::DOF_X | State::DOF_Y;
		//	n->selfstate->blockedDOFs = 63;
		//}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nf = free_nodNei.begin(); Nf != free_nodNei.end(); ++Nf) {
			const shared_ptr<NodeNeighbor>& n = *Nf;
			n->selfstate->blockedDOFs = State::DOF_NONE;
		}
		//for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nb = bottom_nodNei.begin(); Nb != bottom_nodNei.end(); ++Nb) {
		//	const shared_ptr<NodeNeighbor>& n = *Nb;
		//	//	n->selfstate->blockedDOFs = State::DOF_X | State::DOF_Y;
		//	n->selfstate->blockedDOFs = 63;
		//}
	}
	void NodeForceController::servo(double dt) {//控制上墙及下墙,上下膜节点伺服移动
		consol_ss();//获取应力 wszz_top
		//get_gainz();//计算 gain_zb 伺服参数 需在外间隔获取
		consol_ff_cylinder();//选择固结时是否施加施加围压集中力控制？
		//#ifdef SOLE_GAIN
		//Real udz1;
		//udz1 = gain_z1 * (wszz_top - goalz);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
		//udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
		//top_wall->state->pos[2] += udz1;//z
		//#else
		//Real udz1;
		//udz1 = gain_z1 * (wszz - goalz);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
		//udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
		//top_wall->state->pos[2] += udz1;//z
		//#endif
		//std::cout<<"servo--moving wall:"<<udx<<" "<<udy<<std::endl;
		if (z_servo) {//stress control in z direction// z_servo为false 底面开始伺服
#ifdef SOLE_GAIN
			Real udz, udz1;
			udz = gain_zb * (wszz_bottom - goalz);///dt;
			udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			bottom_wall->state->pos[2] += -udz;//z
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
				(*b)->selfstate->pos[2] += -udz;//z 
			}
			udz1 = gain_zt * (wszz_top - goalz);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
			udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
			top_wall->state->pos[2] += udz1;//z
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
				(*b)->selfstate->pos[2] += udz1;//z +hao
			}
#else
			Real udz;
			udz = gain_zb * (wszz - goalz);///dt;
			udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			top_wall->state->pos[2] += udz;//z
			bottom_wall->state->pos[2] += -udz;//z
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
				(*b)->selfstate->pos[2] += -udz;//z 
			}
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
				(*b)->selfstate->pos[2] += udz;//z +hao
			}
#endif


		}
		else {//shear with rate
			double udz = goalzv * dt*height;
			bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate 0.01
			top_wall->state->pos[2] += -udz;//z
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
				(*b)->selfstate->pos[2] -= udz;//z  -hao 
			}
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
				(*b)->selfstate->pos[2] += udz;//z 
			}
		}
	}
	void NodeForceController::servotop(double dt) {//控制上墙伺服移动至规定压力
		consol_ss();//获取应力 wszz_top
		Real udz1;
		udz1 = gain_zt * (wszz_top - goalz);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
		udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
		top_wall->state->pos[2] += udz1;//z
		//else {//shear with rate
		//	double udz = goalzv * dt*height;
		//	bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate 0.01
		//	top_wall->state->pos[2] += -udz;//z
		//	for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
		//		(*b)->selfstate->pos[2] -= udz;//z  -hao 
		//	}
		//	for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
		//		(*b)->selfstate->pos[2] += udz;//z 
		//	}
		//}
	}
	void NodeForceController::servo1(double dt) {//控制辅助上1墙伺服移动
		consol_ss1();//获取顶面应力
		//consol_ff_cylinder();////填充不需要施加围压 可设定填充时施加围压
#ifdef SOLE_GAIN
			//Real udz, udz1;
			//udz = gain_zb * (wszz_bottom - goalz1);///dt;
			//udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			//bottom_wall->state->pos[2] += -udz;//z
		Real udz1;
		udz1 = gain_zt1 * (wszz_top1 - goalz1);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
		udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
		third_wall->state->pos[2] += udz1;//z
#else
		Real udz;
		udz = gain_zt1 * (wszz1 - goalz1);///dt;
		udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
		third_wall->state->pos[2] += udz;//z
	//	bottom_wall->state->pos[2] += -udz;//z
#endif

	}
	void NodeForceController::consol_ss() {
		// sync thread storage of ForceContainer
		scene->forces.sync();
		updateSampleSize();
		//calculate stress
		//wall forces
			//get box size
		/*updateBoxSize();
		force_left = getForce(scene, left_down_wall)[0];
		force_right = getForce(scene, right_down_wall)[0];
		shear = (force_right + force_left) / z_area_inShear;*/
#ifdef SOLE_GAIN
		wszz_top = getForce(scene, boxid2)[2] / z_area;
		wszz_bottom = -getForce(scene, boxid1)[2] / z_area;
		wszz = 0.5*(wszz_top + wszz_bottom);
#else
		wszz = 0.5*(getForce(scene, boxid2)[2] - getForce(scene, boxid1)[2]) / z_area;
#endif
	}
	void NodeForceController::consol_ss1() {
		// sync thread storage of ForceContainer
		scene->forces.sync();
		updateSampleSize();
#ifdef SOLE_GAIN
		wszz_top1 = getForce(scene, boxid3)[2] / z_area;
		wszz_bottom = -getForce(scene, boxid1)[2] / z_area;
		wszz1 = 0.5*(wszz_top1 + wszz_bottom);
#else
		wszz1 = 0.5*(getForce(scene, boxid3)[2] - getForce(scene, boxid1)[2]) / z_area;
#endif

	}
	void NodeForceController::generalConsolidation() {
		const Real& dt = scene->dt;
		iterate_num += 1;
		solve_num += 1;
		servo(dt);//控制位移量
		if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
			iterate_num = 0;
			get_gainz();
		}

		if (solve_num >= solve_num_max) {//check the stability of the sample
			solve_num = 0;
			if ((std::abs(wszz_bottom - goalz) < goalz*f_threshold) && (std::abs(wszz_top - goalz) < goalz*f_threshold)) {
				UnbalancedForce = ComputeUnbalancedForce();
				if (UnbalancedForce < unbf_tol) {//0.05 by default//stop running
					scene->stopAtIter = scene->iter + 1;
					std::cerr << "consolidation completed!" << std::endl;
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
				}
			}
		}
	}
	void NodeForceController::updateSampleSize() {
		Real top = top_wall->state->pos.z();
		Real top1 = third_wall->state->pos.z();
		Real bottom = bottom_wall->state->pos.z();
		height = top - bottom;
		height1 = top1 - bottom;
	}
	void NodeForceController::fillSample() {
		const Real& dt = scene->dt;
		iterate_num += 1;
		solve_num += 1;
		servo1(dt);//控制位移量
		if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
			iterate_num = 0;
			get_gainz1();
		}
		if (solve_num >= solve_num_max) {//check the stability of the sample
			solve_num = 0;
			if (std::abs(wszz_top1 - goalz1) < 0.1*goalz1*f_threshold) {
				UnbalancedForce = ComputeUnbalancedForce();
				if ((UnbalancedForce < unbf_tol)||(third_wall->state->pos.z()== (top_wall->state->pos.z()+2e-3))) {//0.05 by default//stop running
					//fill = false;
					dofilling = false;
					moving = true;
					consol_ss1();//获取应力
					scene->stopAtIter = scene->iter + 1;
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z top1:" << wszz_top1 / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
					std::cerr << "filling sample completed,move top wall for generalConsolidation!" << std::endl;
					deletbody();
					std::cerr << "deletbody out of bound is over ,paper moving top box!!" << std::endl;
					get_gainz();
				}
			}
		}
	}
	void NodeForceController::deletbody() {
		BodyContainer::iterator bb = scene->bodies->begin();
		BodyContainer::iterator iiEnd = scene->bodies->end();
		const shared_ptr<BodyContainer>& BodyCT = scene->bodies;
		
		for (; bb != iiEnd; ++bb){
			const shared_ptr<Body>& b = *bb;
			if (b->shape->getClassName() == "Sphere") {
				const shared_ptr<Sphere>& A = YADE_PTR_CAST<Sphere>(b->shape);
				//Real radius = A->radius;
				Vector3r position = b->state->pos;
				Body::id_t id = b->getId();
				if (((position[2] > max[2]) || (position[2] < min[2]))|| ((position[1] > 30e-3) || (position[1] < -30e-3))|| ((position[0] > 30e-3) || (position[0] < -30e-3))) {
					BodyCT->erase(id, false);
				}
			}
		}
	}
	void NodeForceController::movetop() {
		const Real& dt = scene->dt;
		iterate_num += 1;
		solve_num += 1;
		servotop(dt);//控制位移量
		if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
			iterate_num = 0;
			get_gainz();
		}

		if (solve_num >= solve_num_max) {//check the stability of the sample
			solve_num = 0;
			if ((std::abs(wszz_bottom - goalz) < goalz*f_threshold) && (std::abs(wszz_top - goalz) < goalz*f_threshold)) {
				UnbalancedForce = ComputeUnbalancedForce();
				if (UnbalancedForce < unbf_tol) {//0.05 by default//stop running
					scene->stopAtIter = scene->iter + 1;
					std::cerr << "moving top box completed ,prepare generalConsolidation!" << std::endl;
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						//<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)//体积计算难点
						<< "\n";
					fill = false;
				}
			}
			if (top_wall->state->pos[2] == (top_nodNei[0]->selfstate->pos[2] + 1e-3)) {
				scene->stopAtIter = scene->iter + 1;
				std::cerr << "moving top box completed forced!" << std::endl;
				std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
					<< ", Ss_z:" << wszz / 1000.0<< "\n";
				std::cerr << "the top is About to colliding! stop moving!!" << std::endl;
				fill = false;
			}
		}
	}
}
namespace yade {
	YADE_PLUGIN((NodeForceController1));
	NodeForceController1::~NodeForceController1() { }
	Real NodeForceController1::ComputeUnbalancedForce(bool maxUnbalanced) { return Shop::unbalancedForce(maxUnbalanced, scene); }
	void NodeForceController1::action() {//引擎运行模块
		const Real& dt = scene->dt;
		if (firstRun1) {
			//initialize();
			std::cout << "First run, will initiallize!" << std::endl;
			z_area = M_PI * pow(wall_radius, 2.0);//z轴面积
			BodyContainer::iterator bi1 = scene->bodies->begin();
			BodyContainer::iterator bi1end = scene->bodies->end();
			for (; bi1 != bi1end; ++bi1) {//初始化box top and bottom
				const shared_ptr<Body>& b = *bi1;
				if (b->shape->getClassName() == "Box") {
					cout << "getclassname=" << b->shape->getClassName() << endl;
					bottom_wall = *bi1;
					++bi1;
					cout << "getclassname=" << b->shape->getClassName() << endl;
					top_wall = *bi1;
					break;
				}
			}
			
			boxid1 = bottom_wall->id;
			boxid2 = top_wall->id;
			
			Real top = top_wall->state->pos[2];
			Real bottom = bottom_wall->state->pos[2];
			height = height0 = top - bottom;
			Init_SampleVolume = Current_SampleVolume = z_area * height0;

			iterate_num = 0;//迭代步数
			solve_num = 0;

			BodyContainer::iterator bi = scene->bodies->begin();
			BodyContainer::iterator biEnd = scene->bodies->end();
			particlesVolume = 0;
			for (; bi != biEnd; ++bi)
			{
				const shared_ptr<Body>& b = *bi;
				/*if (b->isClump()) {
					//std::cerr << "density = "<<b->material->density << "\n";
					const shared_ptr<Clump>& clump =YADE_PTR_CAST<Clump>(b->shape);
					if (clump->members.size() > 0) {
						shared_ptr<Body> subBody = Body::byId(clump->members.begin()->first);//use the material of the first member
						particlesVolume += b->state->mass / subBody->material->density;//this should be a general method to access a particle volume.
					}
					continue;
				}
				else if (b->isClumpMember()) {
					continue;
				}*/
				if (b->shape->getClassName() == "Sphere") {
					const shared_ptr<Sphere>& A = YADE_PTR_CAST<Sphere>(b->shape);
					particlesVolume += 4.0 / 3.0*M_PI*pow(A->radius, 3.0);

				}
			}
			get_gainz();//固结伺服增益获取 
			firstRun1 = false;
		}
		if (z_servo) {//consolidation 固结
			generalConsolidation();//施加围压  一段时间获取伺服参数 get_gainz
			//if (scene->iter % savedata_interval == 0) { recordData(); }//recordData输出信息函数//输出信息
			if (scene->iter % echo_interval == 0) {
				UnbalancedForce = ComputeUnbalancedForce();
				consol_ss();//包含 updateSampleSize();
				//getStressStrain();
				std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
					//<< ", Ss_x:" << wsxx / 1000.0 //average stress of the two walls, unit is kPa
					//<< ", Ss_y:" << wsyy / 1000.0
					<< ", Ss_z:" << wszz / 1000.0
					<< ", e:" << (Current_SampleVolume / particlesVolume - 1.0)//体积计算难点
					<< "\n";
			}

		}
		else {//剪切 z_servo关闭
			if (firstRun2) {//consolidation completed normally
				initialize();//获取膜节点设置
				std::cerr << "Shear begins..." << std::endl;
				updateSampleSize();
				height0 = height;
				get_gainz();
				BodyContainer::iterator bi1 = scene->bodies->begin();
				BodyContainer::iterator bi1end = scene->bodies->end();
				for (;bi1!=bi1end; ++bi1) {//初始化side_wall
					const shared_ptr<Body>& b = *bi1;
					if (b->shape->getClassName() == "BoxM") {
						cout << "getclassname=" << b->shape->getClassName() << endl;
						side_wall = *bi1;
						break;
					}
				}
				if(side_wall!=NULL){
					boxid3 = side_wall->id;
					posx0=side_wall->state->pos[0];
				}else cout<<"shear model no boxm !!"<<endl;
				if(!useside)posx0=bottom_wall->state->pos[0];
				firstRun2 = false;
			}
			servo(dt);//集中力施加 包含更新边界 应变率控制 无需获取伺服参数
			if (force) {
				iterate_num += 1;
				solve_num += 1;
				if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
					iterate_num = 0;
					get_gainz();
				}
				if(useside){
					if (side_wall->state->pos[0]-posx0 > 20e-3) {//
				//stop move

					std::cerr << "moving x ends!" << std::endl;
					scene->stopAtIter = scene->iter + 1;//停止模拟
					UnbalancedForce = ComputeUnbalancedForce();
					consol_f();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", displace:" << side_wall->state->pos[0]-posx0//x coordinate
						<< ", side force x:" << wfxx_side
						<< ", Ss_z:" << wszz / 1000.0
						<< "\n";
					}
				}else{
					
					if (bottom_wall->state->pos[0]-posx0 > 20e-3) {//
				//stop move

					std::cerr << "moving x ends!" << std::endl;
					scene->stopAtIter = scene->iter + 1;//停止模拟
					UnbalancedForce = ComputeUnbalancedForce();
					consol_f();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", displace:" << bottom_wall->state->pos[0]-posx0//x coordinate
						<< ", side force x:" << wfxx_bottomn
						<< ", Ss_z:" << wszz / 1000.0
						<< "\n";
					}
				}
				
			}
			else {
				if (log(height0 / height) > target_strain) {//
				//stop shear

					std::cerr << "Shear ends!" << std::endl;
					scene->stopAtIter = scene->iter + 1;//停止模拟
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", Str_z:" << log(height0 / height)//z
						<< ", Ss_z:" << wszz / 1000.0
						<< ", e:" << (Current_SampleVolume / particlesVolume - 1.0)
						<< "\n";
				}
			}

			//output info
			if (scene->iter % echo_interval == 0) {
				UnbalancedForce = ComputeUnbalancedForce();
				//wfxx_bottomn
				consol_ss();
				if (force) {
					consol_f();
					if(useside){
						std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", displace:" << side_wall->state->pos[0]-posx0//x coordinate
						<< ", side force x:" << wfxx_side
						<< ", Ss_z:" << wszz / 1000.0
						<< "\n";
					}
					else{
						std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", displace:" << bottom_wall->state->pos[0]-posx0//x coordinate
						<< ", side force x:" << wfxx_bottomn
						<< ", Ss_z:" << wszz / 1000.0
						<< "\n";
					}
					
				}
				else {

					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						//<< ", conf:" << (wsxx + wsyy) / 2000.0 //average stress of the two walls, unit is kPa
						<< ", Str_z:" << log(height0 / height)//z
						<< ", Ss_z:" << wszz / 1000.0
						<< ", e:" << (Current_SampleVolume / particlesVolume - 1.0)
						<< "\n";
				}

			}
		}
	}

	void NodeForceController1::initialize() {
		//store elements of the flexible wall
		std::cout << "Initializing MembraneNode" << std::endl;
		const shared_ptr<NodeNeighborContainer>& NContainer = scene->nodes;//scene->nodes  shared_ptr<NodeNeighborContainer>
		//proxee = scene->nodes;//scene->nodes  shared_ptr<NodeNeighborContainer>
		Num = NContainer->size();
		assert(Num > 0);
		cows = NContainer->cows;
		cols = NContainer->cols;
		assert((cows*cols) == Num);
		NodeNeighborContainer::iterator NeiFirst = NContainer->begin();
		NodeNeighborContainer::iterator NeiEnd = NContainer->end();
		std::vector<shared_ptr<NodeNeighbor>> Ntop(NeiFirst, NeiFirst + cols);
		std::vector<shared_ptr<NodeNeighbor>> Nfree(NeiFirst + cols, NeiEnd - cols);
		std::vector<shared_ptr<NodeNeighbor>> Nbottom(NeiEnd - cols, NeiEnd);
		//save nodeneighbor
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Ntop.begin(); n != Ntop.end(); ++n) {
			top_nodNei.push_back(*n);
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Nfree.begin(); n != Nfree.end(); ++n) {
			free_nodNei.push_back(*n);
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator n = Nbottom.begin(); n != Nbottom.end(); ++n) {
			bottom_nodNei.push_back(*n);
		}
		fix_node();//固定边界墙？
		free_node();

	}
	/*void NodeForceController1::flexible_consolidate() {
		const Real& dt = scene->dt;


	}*/
	void NodeForceController1::quiet_system() {
		BodyContainer::iterator bb = scene->bodies->begin();
		BodyContainer::iterator iiEnd = scene->bodies->end();

		for (; bb != iiEnd; ++bb)
		{
			State* state = (*bb)->state.get();
			state->vel = Vector3r(0.0, 0.0, 0.0);
			state->angVel = Vector3r(0.0, 0.0, 0.0);
		}
	}
	void NodeForceController1::get_gainz() {//获取上下墙的伺服参数G

		//if (!z_servo) { return; }
		avg_zstiffb = 0.0;
		avg_zstifft = 0.0;
		InteractionContainer::iterator ii = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();
		for (; ii != iiEnd; ++ii) if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id == boxid1) {//caution: I assume that wall_id is less than 2. top and bottom walls
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstiffb += currentContactPhysics->kn;
				continue;
			}
			else if (id == boxid2) {
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				avg_zstifft += currentContactPhysics->kn;
				continue;
			}
			gain_zb = gain_alpha * z_area / (avg_zstiffb);//bottom wall 
			gain_zt = gain_alpha * z_area / (avg_zstifft);//top wall

		}
	}
	/*void NodeForceController1::servo_cylinder(double dt) {

		consol_ff_cylinder();
	}*/
	void NodeForceController1::consol_ff_cylinder() {//free_nodNei  膜节点集中力控制
		// sync thread storage of ForceContainer

		scene->forces.sync();
		//std::vector<shared_ptr<NodeNeighbor>> free_nodNei;
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator N = free_nodNei.begin(); N != free_nodNei.end(); ++N) {
			const shared_ptr<NodeNeighbor>& n = *N;
			Body::id_t id = n->selfId;
			Vector3r Force = n->computeUnitForce();
			//Vector3r Force=n->computeUnitForce();
			//if (id % 150 == 0)
				//std::cout <<id<< "computeUnitForce=" << Force.norm() <<"neighbor num"<<endl;
			addForce(scene, id, Force*wrap_pressure);
		}
	}
	void NodeForceController1::fix_node() {//free_nodNei 固定膜颗粒
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nt = top_nodNei.begin(); Nt != top_nodNei.end(); ++Nt) {
			const shared_ptr<NodeNeighbor>& n = *Nt;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nf = free_nodNei.begin(); Nf != free_nodNei.end(); ++Nf) {
			const shared_ptr<NodeNeighbor>& n = *Nf;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nb = bottom_nodNei.begin(); Nb != bottom_nodNei.end(); ++Nb) {
			const shared_ptr<NodeNeighbor>& n = *Nb;
			n->selfstate->blockedDOFs = State::DOF_ALL;
		}
	}
	void NodeForceController1::free_node() {//free_nodNei 固定膜颗粒
		//for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nt = top_nodNei.begin(); Nt != top_nodNei.end(); ++Nt) {
		//	const shared_ptr<NodeNeighbor>& n = *Nt;
		//	//n->selfstate->blockedDOFs = State::DOF_X | State::DOF_Y;
		//	n->selfstate->blockedDOFs = 63;
		//}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nf = free_nodNei.begin(); Nf != free_nodNei.end(); ++Nf) {
			const shared_ptr<NodeNeighbor>& n = *Nf;
			n->selfstate->blockedDOFs = State::DOF_NONE;
		}
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nb = bottom_nodNei.begin(); Nb != bottom_nodNei.end(); ++Nb) {
			const shared_ptr<NodeNeighbor>& n = *Nb;
				//n->selfstate->blockedDOFs = State::DOF_Z;
				n->selfstate->blockedDOFs = State::DOF_ALL;
			//n->selfstate->blockedDOFs = 63;
		}
		if(freez){
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nt = top_nodNei.begin(); Nt != top_nodNei.end(); ++Nt) {
				const shared_ptr<NodeNeighbor>& n = *Nt;
				n->selfstate->blockedDOFs = State::DOF_X | State::DOF_Y;
			}
		
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator Nb = bottom_nodNei.begin(); Nb != bottom_nodNei.end(); ++Nb) {
				const shared_ptr<NodeNeighbor>& n = *Nb;
				n->selfstate->blockedDOFs = State::DOF_X | State::DOF_Y;
			}
			
		}
		
	}
	void NodeForceController1::servo(double dt) {//控制上墙及下墙,上下膜节点伺服移动
		consol_ss();//获取应力 wszz_top wszz_bottom 更新墙应力
		if (z_servo) {//stress control in z direction// z_servo为false 底面开始伺服
#ifdef SOLE_GAIN
			Real udz, udz1;
			udz = gain_zb * (wszz_bottom - goalz);///dt;
			udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			bottom_wall->state->pos[2] += -udz;//z
			udz1 = gain_zt * (wszz_top - goalz);///dt;//dt;伺服参数G=gain_x/dt v=G*(wsxx_left - goalx) 位移s=v*dt= gain_x * (wsxx_left - goalx)=udx
			udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
			top_wall->state->pos[2] += udz1;//z
#else
			Real udz;
			udz = gain_zb * (wszz - goalz);///dt;
			udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			top_wall->state->pos[2] += udz;//z
			bottom_wall->state->pos[2] += -udz;//z
#endif
		}
		else if (force) {
			Real udz;
			udz = gain_zb * (wszz_bottom - goalz);///dt; 获取gain_zb 
			udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
			Real udx = goalxv * dt*height;
			if(useside)side_wall->state->pos[0] += udx;//x轴正向移动
			bottom_wall->state->pos[0] += udx;
			for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
				(*b)->selfstate->pos[0] += udx;//x轴正向移动
			}
			if (movewallOnly){
				bottom_wall->state->pos[2] += -udz;//z
			}
			else{
				bottom_wall->state->pos[2] += -udz;//z
				for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
					(*b)->selfstate->pos[2] += -udz;//z 
				}
			}
			consol_ff_cylinder();//集中力施加

		}
		else {//shear with constant rate
			Real udz = goalzv * dt*height;
			if (movewallOnly) {
				bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate 0.01
				top_wall->state->pos[2] += -udz;//z
			}
			else {
				bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate 0.01
				top_wall->state->pos[2] += -udz;//z
				for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = top_nodNei.begin(); b != top_nodNei.end(); ++b) { //上节点位移控制
					(*b)->selfstate->pos[2] -= udz;//z  -hao 
				}
				for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
					(*b)->selfstate->pos[2] += udz;//z 
				}
			}
			consol_ff_cylinder();//集中力施加
		}
	}
	void NodeForceController1::consol_ss() {
		// sync thread storage of ForceContainer
		scene->forces.sync();
		updateSampleSize();
#ifdef SOLE_GAIN
		wszz_top = getForce(scene, boxid2)[2] / z_area;
		wszz_bottom = -getForce(scene, boxid1)[2] / z_area;
		wszz = 0.5*(wszz_top + wszz_bottom);
#else
		wszz = 0.5*(getForce(scene, boxid2)[2] - getForce(scene, boxid1)[2]) / z_area;
#endif
	}
	void NodeForceController1::consol_f() {
		// sync thread storage of ForceContainer
		scene->forces.sync();
		//updateSampleSize();wfxx_bottom, wfxx_bottomn
		wfxx_bottomn = 0;
		wfxx_side=0;
		wfxx_bottom = getForce(scene, boxid1)[0];
		wfxx_side=getForce(scene, boxid3)[0];
		for (std::vector<shared_ptr<NodeNeighbor>>::iterator b = bottom_nodNei.begin(); b != bottom_nodNei.end(); ++b) {
			wfxx_bottomn += getForce(scene, (*b)->selfId)[0];
		}
	}
	void NodeForceController1::generalConsolidation() {
		const Real& dt = scene->dt;
		iterate_num += 1;
		solve_num += 1;
		servo(dt);//控制位移量
		if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
			iterate_num = 0;
			get_gainz();
		}

		if (solve_num >= solve_num_max) {//check the stability of the sample
			solve_num = 0;
			if ((std::abs(wszz_bottom - goalz) < goalz*f_threshold) && (std::abs(wszz_top - goalz) < goalz*f_threshold)) {
				UnbalancedForce = ComputeUnbalancedForce();
				if (UnbalancedForce < unbf_tol) {//0.05 by default//stop running
					scene->stopAtIter = scene->iter + 1;
					consol_ss();
					std::cerr << "consolidation completed! plesea  generate Membrane now" << std::endl;
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						<< ", e:" << (Current_SampleVolume / particlesVolume - 1.0)//体积计算难点!!!
						<< "\n";
					if(firstRun){
						//deletbody();
						firstRun=false;
					}
					//this->dead = true;
				}
			}
		}
	}
	void NodeForceController1::updateSampleSize() {
		Real top = top_wall->state->pos.z();
		Real bottom = bottom_wall->state->pos.z();
		height = top - bottom;
		Current_SampleVolume = z_area * height;
	}
	void NodeForceController1::deletbody() {//delet Facet
		BodyContainer::iterator bb = scene->bodies->begin();
		BodyContainer::iterator iiEnd = scene->bodies->end();
		const shared_ptr<BodyContainer>& BodyCT = scene->bodies;
		for (; bb != iiEnd; ++bb) {
			const shared_ptr<Body>& b = *bb;//squaredNorm()  (pow(b->state->pos.x(),2)+pow(b->state->pos.y(),2))>pow(wall_radius,2)
			if ((b->shape->getClassName() == "Facet") || (b->state->pos.z() > top_wall->state->pos.z()) || (b->state->pos.z() < bottom_wall->state->pos.z()) || ((pow(b->state->pos.x(), 2) + pow(b->state->pos.y(), 2)) > pow(wall_radius, 2))) {
				Body::id_t id = b->getId();
				BodyCT->erase(id, false);
			}
		}
	}
}

