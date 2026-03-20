#include "DirectShearController.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Clump.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <assert.h>
#include <preprocessing/dem/Shop.hpp>
#include<fstream>

namespace yade { // Cannot have #include directive inside.

	CREATE_LOGGER(ShearEngine);
	YADE_PLUGIN((ShearEngine));

	ShearEngine::~ShearEngine() { }


	Real ShearEngine::ComputeUnbalancedForce(bool maxUnbalanced)
	{
		return Shop::unbalancedForce(maxUnbalanced, scene);
	}
	void ShearEngine::updateBoxSize() //更新盒子大小
	{
		Real left;
		left = box[left_down_wall]->state->pos.x();
		Real left_inShear;
		left_inShear = box[left_up_wall]->state->pos.x();
		Real right;
		right = box[right_down_wall]->state->pos.x();
		Real front;
		front = box[front_wall]->state->pos.y();
		Real back;
		back = box[back_wall]->state->pos.y();
		Real top;
		top = box[top_wall]->state->pos.z();
		Real bottom;
		bottom = box[bottom_wall]->state->pos.z();

		width = right - left;
		depth = back - front;
		height = top - bottom;
		width_inShear = right - left_inShear;
		boxVolume = width * depth*height;
		x_area = depth * height;
		y_area = width * height;
		z_area = depth * width;//z向面积
		z_area_inShear = depth * width_inShear;//剪切面面积

	}
	void ShearEngine::getBox()//box墙体
	{
		for (int i = 0; i < 14; i++) {
			box[i] = Body::byId(i);
		}
	}

	//2026-3-14
	bool ShearEngine::if_set_big_Box()//box墙体
	{
		return  big_front_wall_id != -1;
	}

	//
	void ShearEngine::servo(double dt) {
		//if (!consol_on){return;}//if the flag consol_on is false, which means not consolidating

		consol_ss();//更新边界大小
		Real udz1;
		udz1 = gain_z1 * (wszz_top - goalz);///dt;
		udz1 = math::sign(udz1)*math::min(fabs(udz1), max_vel*dt);//constrain the maximum velocity of the wall
		box[top_wall]->state->pos[2] += udz1;//z
		//std::cout<<"servo--moving wall:"<<udx<<" "<<udy<<std::endl;
		/*if (z_servo) {//stress control in z direction
			if(two_way){
				Real udz;
				udz = gain_z * (wszz_bottom - goalz);///dt;
				udz = math::sign(udz)*math::min(fabs(udz), max_vel*dt);//constrain the maximum velocity of the wall
				box[bottom_wall]->state->pos[2] += -udz;
			}

		}*/
		if(!z_servo) {//strain control during shear
			box[2]->state->pos[0] += goalx * dt*width;//z,remain a constant loading strain rate
			box[3]->state->pos[0] += goalx * dt*width;//z
			box[4]->state->pos[0] += goalx * dt*width;
			box[7]->state->pos[0] += goalx * dt*width;
		}
	}
	void ShearEngine::consol_ss() {
		// sync thread storage of ForceContainer
		scene->forces.sync();
		//calculate stress
		//wall forces
			//get box size
		updateBoxSize();
		force_left = getForce(scene, left_down_wall)[0];
		force_right = getForce(scene, right_down_wall)[0];
		shear = (force_right + force_left) / z_area_inShear;
		wszz_top = getForce(scene, top_wall)[2] / z_area;
		wszz_bottom = -getForce(scene, bottom_wall)[2] / z_area;
		wszz = 0.5*(wszz_top + wszz_bottom);

	}									//recorder
	void ShearEngine::recordData() {
		if (!out.is_open()) {
			assert(!out.is_open());

			std::string fileTemp = file;


			if (fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
			//out.open(fileTemp.c_str(), truncate(fileTemp) ? std::fstream::trunc : std::fstream::app);
			out.open(fileTemp.c_str(), std::fstream::app);
			if (!out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `" + fileTemp + "'.");
		}
		Porosity_e=(width*depth*height) / particlesVolume - 1.0;
		out << "Iter " << scene->iter << " Ubf " << UnbalancedForce
			<< ", Strain:" << log(width / width_inShear)//z
			<< ", Stress_Z:" << wszz / 1000.0
			<< ", Stress_shear:" << shear / 1000.0
			<< ", Porosity_e:" << Porosity_e
			<< ", Volume_strain:" <<log((Porosity_e+1)/(Porosity_e0+1))
			<< endl;

	}
	void ShearEngine::get_gain() {//determine servo gain parameters


		avg_zstiff = 0.0;
		avg_zstiff1 = 0.0;

		InteractionContainer::iterator ii = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();

		for (; ii != iiEnd; ++ii) if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			//Real fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce.norm();
			//if (fn!=0)
			//{
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id < 8) {//caution: I assume that wall_id is less than 6.
				Real kn;
				FrictPhys* currentContactPhysics = static_cast<FrictPhys*> (contact->phys.get());
				kn = currentContactPhysics->kn;
				if (id == 6)
					avg_zstiff += kn;
				else if (id == 7)
					avg_zstiff1 += kn;
			}
		}

		gain_z = gain_alpha * z_area / (avg_zstiff);
		gain_z1 = gain_alpha * z_area / (avg_zstiff1);

		//std::cout<<"avg_xstiffness"<<avg_xstiff<<"avg_ystiffness"<<avg_ystiff <<"avg_zstiffness"<<avg_zstiff<<std::endl;
	}
	
	void ShearEngine::move_topwall(){//向上移动上墙循环固结
		if(firstMove){
			Porosity_first_GC=(width*depth*height) / particlesVolume - 1.0;
			firstMove=false;
			cout<<"Porosity_first_GC ="<<(width*depth*height) / particlesVolume - 1.0<<endl;
		}
		if(std::abs(Poro_ratio[0]-Porosity_first_GC)<0.08)
		box[top_wall]->state->pos[2] += scalars/1000;//z
		else{
			Real compute_scalars=0.;
			if(Poro_ratio[0]-Poro_ratio_stop<0)compute_scalars=1.0;
			else compute_scalars=(scalars-1.0)/(Porosity_first_GC-0.08-Poro_ratio_stop)*(Poro_ratio[0]-Poro_ratio_stop)+1.0;
			box[top_wall]->state->pos[2] += compute_scalars/1000;//z
		}
		/*Real d_height=box[top_wall]->state->pos.z()-box[bottom_wall]->state->pos.z();
		Real Poro_ratio_now=Poro_ratio[0];
		Real Poro_ratio_front=Poro_ratio[1];
		Real udz1=0;
		if(Poro_ratio_front==0){
			udz1=scalars*d_height*0.01*Poro_ratio_now;
		}else{
			udz1=(Poro_ratio_front-Poro_ratio_now)*d_height*scalars;
		}
		box[top_wall]->state->pos[2] += udz1;//z
		if((std::abs(Poro_ratio_now-Poro_ratio_front)<0.0005)&&Poro_ratio_now>0)loop_Conso=false;*/
	}
	void ShearEngine::generalConsolidation() {//固结
		const Real& dt = scene->dt;
		iterate_num += 1;
		solve_num += 1;
		servo(dt);//控制位移量
		if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles by default
			iterate_num = 0;
			get_gain();
		}

		if (!vibrate_gc&&solve_num >= solve_num_max) {//check the stability of the sample
			solve_num = 0;
			Poro_ratio[1]=Poro_ratio[0];
			Poro_ratio[0]=((width*depth*height) / particlesVolume - 1.0);//now e
			if(Poro_ratio[0]<=Poro_ratio_stop&&loop_Conso){//到到设定堆积密度
				scene->stopAtIter = scene->iter + 1;
				std::cerr << "consolidation completed!" << std::endl;	
				UnbalancedForce = ComputeUnbalancedForce();
				std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
							<< ", Ss_z:" << wszz / 1000.0
							<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
							<< "\n";
				loop_Conso=false;
			}
			if ((std::abs(wszz_bottom - goalz) < goalz*f_threshold) && (std::abs(wszz_top - goalz) < goalz*f_threshold)) {
				UnbalancedForce = ComputeUnbalancedForce();
				if (UnbalancedForce < unbf_tol) {//o.5 by default
							 //stop running
					if(loop_Conso){
						//Conso_Active=false;
						std::cerr << "loop consolidation! move awll begin" << std::endl;
						move_topwall();
						std::cerr << "now Poro_ratio is:" << ", Ss_z:" << wszz / 1000.0
							<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
							<< "\n";
					}else{
						scene->stopAtIter = scene->iter + 1;
						std::cerr << "consolidation completed!" << std::endl;	
						UnbalancedForce = ComputeUnbalancedForce();
						//consol_ss();
						std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
							<< ", Ss_z:" << wszz / 1000.0
							<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
							<< "\n";
					}
					
					
				}
			}
			Real Poro_ratio_now=Poro_ratio[0];
			Real Poro_ratio_front=Poro_ratio[1];
			if((std::abs(Poro_ratio_now-Poro_ratio_front)<0.0001)&&Poro_ratio_now>0&&loop_Conso)
			{
				std::cerr << "loop consolidation! move awll begin" << std::endl;
						move_topwall();
						std::cerr << "now Poro_ratio is:" << ", Ss_z:" << wszz / 1000.0
							<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
							<< "\n";
			}
		}
		else if(vibrate_gc){//大步长振动固结
			vibrate_count+=1;
			if(vibrate_count>=int(real_time/dt)){
				scene->stopAtIter = scene->iter + 1;
				std::cerr << "vibrate consolidation has completed! change ks again" << std::endl;
				UnbalancedForce = ComputeUnbalancedForce();
						//consol_ss();
						std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
							<< ", Ss_z:" << wszz / 1000.0
							<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
							<< "\n";
							vibrate_gc=false;
			}
			
		}
	}
	
	void ShearEngine::computeParticlesVolume()
	{
		BodyContainer::iterator bi = scene->bodies->begin();
		BodyContainer::iterator biEnd = scene->bodies->end();
		particlesVolume = 0;
		for (; bi != biEnd; ++bi)
			{
				//if((*bi)->isClump()) continue;//no clump for Superquadrics
				const shared_ptr<Body>& b = *bi;
				//std::cerr << "watch point at particle volume calc " << "\n";
				if (b->shape->getClassName() == "Sphere") {
					const shared_ptr<Sphere>& A = YADE_PTR_CAST<Sphere>(b->shape);
					particlesVolume += 4.0 / 3.0*M_PI*pow(A->radius, 3.0);
				}
			}
	}
	void ShearEngine::action()
	{
		// sync thread storage of ForceContainer同步线程存储
		scene->forces.sync();
		if (firstRun) { // sync boundaries ids in the table

			getBox();
			Real left = box[left_down_wall]->state->pos.x();
			Real right = box[right_down_wall]->state->pos.x();
			Real front = box[front_wall]->state->pos.y();
			Real back = box[back_wall]->state->pos.y();
			Real top = box[top_wall]->state->pos.z();
			Real bottom = box[bottom_wall]->state->pos.z();

			width = width0 = right - left;
			depth = depth0 = back - front;
			height = height0 = top - bottom;

			x_area = depth * height;
			y_area = width * height;
			z_area = depth * width;
			Init_boxVolume = height0 * width0*depth0;
			//calculate particlesVolume
			computeParticlesVolume();
			firstRun = false;
			iterate_num = 0;
			solve_num = 0;
			get_gain();
			//2026-3-14*
			if (if_set_big_Box()) {
				big_length0 = box[big_right_up_wall_id]->state->pos.x()- box[big_left_up_wall_id]->state->pos.x();//初始长
				big_width0 = box[big_back_wall_id]->state->pos.y() - box[big_front_wall_id]->state->pos.y();
			}
		}
		const Real& dt = scene->dt;
		if (z_servo) {
			if (stop) {
				if (scene->iter % 1000 == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					consol_ss();
					//getStressStrain();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
						<< "\n";
				}

			}
			else if (move_side_wall) {//2026-3-14
				Real now_length= box[big_right_up_wall_id]->state->pos.x() - box[big_left_up_wall_id]->state->pos.x();//当前长
				Real now_width= box[big_back_wall_id]->state->pos.y() - box[big_front_wall_id]->state->pos.y();
				Real now_dis_x = big_length0 - now_length;
				Real now_dis_y = big_width0 - now_width;
				if (scene->iter % 1000 == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					//consol_ss();
					//getStressStrain();
					std::cerr << "Iter " << scene->iter << " now_dis_x  " << now_dis_x
						<< " now_dis_y  " << now_dis_y
						<< "\n";
				}
				if (std::abs(now_dis_x) < move_dis_threshold/1000) {
					move_side_wall = false;
					std::cerr << "Iter " << scene->iter << " now_dis_x  " << now_dis_x
						<< " now_dis_y  " << now_dis_y
						<< "\n";
					scene->stopAtIter = scene->iter + 1;
					std::cerr << "move big wall over!!!" << std::endl;

				}
				else {
					box[big_right_up_wall_id]->state->pos[0] -= m_goalx * dt*big_length0;
					box[big_left_up_wall_id]->state->pos[0] += m_goalx * dt*big_length0;
					box[big_back_wall_id]->state->pos[1] -= m_goalx * dt*big_width0;
					box[big_left_up_wall_id]->state->pos[1] += m_goalx * dt*big_width0;
				}
				/*
				if (big_length0 - now_dis_x >= 0) {
					box[big_right_up_wall_id]->state->pos[0] -= m_goalx * dt*big_length0;
					box[big_left_up_wall_id]->state->pos[0] += m_goalx * dt*big_length0;
				}

				if (big_width0 - now_dis_y >= 0) {
					box[big_back_wall_id]->state->pos[1] -= m_goalx * dt*big_width0;
					box[big_left_up_wall_id]->state->pos[1] += m_goalx * dt*big_width0;

				}
				*/
				
				
				

			}
			else {
				generalConsolidation();//侧限固结
				//if (scene->iter % savedata_interval == 0) { recordData(); }
				if (scene->iter % 1000 == 0) {
					UnbalancedForce = ComputeUnbalancedForce();
					//consol_ss();
					//getStressStrain();
					std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
						<< ", Ss_z:" << wszz / 1000.0
						<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
						<< "\n";
				}
			}

		}
		else {//shear
			if (firstRun2) {//consolidation completed normally
				std::cerr << "Shear begins..." << std::endl;
				updateBoxSize();
				width0 = width;
				height0 = height;
				depth0 = depth;
				if(Porosity_e0==0.){
					Porosity_e0=(width*depth*height) / particlesVolume - 1.0;
				
					BodyContainer::iterator bi = scene->bodies->begin();
					BodyContainer::iterator biEnd = scene->bodies->end();
					for ( ; bi!=biEnd; ++bi ){	
						const shared_ptr<Body>& b = *bi;
						//if(b->isClump())continue;
						b->state->vel=Vector3r::Zero();
						b->state->angVel=Vector3r::Zero();
					}
				}
				
				firstRun2 = false;
			}

			iterate_num += 1;
			if (iterate_num >= iterate_num_max) {//update the servo gains every 100 cycles
				iterate_num = 1;
				get_gain();

			}

			//box[bottom_wall]->state->pos[2] += goalz * dt*height;//z,remain a constant loading strain rate
			//box[top_wall]->state->pos[2] += -goalz * dt*height;//z
			servo(dt);
			updateBoxSize();
			if (log(width / width_inShear) > target_strain) {//
				//stop shear
				std::cerr << "Shear ends!" << std::endl;
				scene->stopAtIter = scene->iter + 1;//停止模拟
			}
			//output info
			if (scene->iter % echo_interval == 0) {
				UnbalancedForce = ComputeUnbalancedForce();
				recordData();
				std::cerr << "Iter " << scene->iter << " Ubf " << UnbalancedForce
					<< ", Str_z:" << log(width / width_inShear)//z log(width/(width - width_inShear))
					<< ", Ss_z:" << wszz / 1000.0
					<< ", shear:" << shear / 1000.0
					<< ", e:" << ((width*depth*height) / particlesVolume - 1.0)
					<< "\n";
			}

		}


	}
	Vector2r ShearEngine::getStress()
	{

		consol_ss();
		/*wszz_top = getForce(scene, top_wall)[2] / z_area;
		wszz_bottom = -getForce(scene, bottom_wall)[2] / z_area;*/
		return Vector2r(wszz_top, wszz_bottom);

	}
	//@@@

}
