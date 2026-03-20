#pragma once

#include <lib/base/Math.hpp>
#include <core/Scene.hpp>
#include <core/BodyContainer.hpp>
#include <pkg/common/BoundaryController.hpp>
#include<fstream>
#include<string>
#include <boost/array.hpp>
namespace yade { // Cannot have #include directive inside.

	class Scene;
	class State;


	class ShearEngine : public BoundaryController {
	private:
		bool                  firstRun;
		bool                  firstRun2;
		bool                  firstMove;
		inline const Vector3r getForce(Scene* rb, Body::id_t id) const
		{
			return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */
		}

		shared_ptr<Body> box[14];//pointers of facets of up box//@@@@
		std::ofstream out;

	public:
		//! internal index values for retrieving walls
		//enum { wall_left = 0, wall_right, wall_bottom, wall_top, wall_back, wall_front };
		enum { left_down_wall = 0, right_down_wall, left_up_wall, right_up_wall, left_cover_wall, right_cover_wall, bottom_wall, top_wall, front_wall, back_wall };//@@@@
		void updateBoxSize();
		void getBox();//box->body
		Real x_area; //area in the x axis
		Real y_area;
		Real z_area;
		Real z_area_inShear;
		void servo(double dt);
		double  wszz_top, wszz_bottom;//上下面应力
		double  force_left, force_right;
		void consol_ss();
		Real  gain_z;
		Real  gain_z1;
		void get_gain();
		Real  avg_zstiff;
		Real  avg_zstiff1;
		void generalConsolidation();
		unsigned int iterate_num, solve_num;
		Real Init_boxVolume;
		double wsxx, wsyy, wszz, shear;//stress in the three directions
		void recordData();//输出信息至文件
		//@@@
		//2026-3-14
		bool if_set_big_Box();
		//

		Real particlesVolume;
		//! Value of box volume
		Real boxVolume;
		Vector2r getStress();

		virtual ~ShearEngine();

		void action() override;
		void move_topwall();
		void computeParticlesVolume();
		Real ComputeUnbalancedForce(bool maxUnbalanced = false);

		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ShearEngine, BoundaryController, "ShearEngine",
			((int, bottom_wall_id, 6, , "id of boundary ; coordinate 1- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, top_wall_id, 7, , "id of boundary ; coordinate 1+ (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, left_down_wall_id, 0, , "id of boundary ; coordinate 0- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, right_down_wall_id, 1, , "id of boundary ; coordinate 0+ (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, front_wall_id, 8, , "id of boundary ; coordinate 2+ (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, back_wall_id, 9, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, left_up_wall_id, 2, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, right_up_wall_id, 3, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, left_cover_wall_id, 4, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, right_cover_wall_id, 5, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((Real, height, 0, Attr::readonly, "size of the box (1-axis) |yupdate|"))
			((Real, width, 0, Attr::readonly, "size of the box (0-axis) |yupdate|"))
			((Real, depth, 0, Attr::readonly, "size of the box (2-axis) |yupdate|"))
			((Real, width_inShear, 0, Attr::readonly, "size of the shear plane |yupdate|"))
			((Real, height0, 0, , "Reference size for strain definition. See :yref:`TriaxialStressController::height`"))
			((Real, width0, 0, , "Reference size for strain definition. See :yref:`TriaxialStressController::width`"))
			((Real, depth0, 0, , "Reference size for strain definition. See :yref:`TriaxialStressController::depth`"))
			//2026-3-14
			((bool, move_side_wall, false , , ""))
			((Real, x_dis, 0, , "size of the big box dis to small size"))
			((Real, y_dis, 0, , "size of the big box dis to small size"))
			((int, big_front_wall_id, -1, , "id of boundary ; coordinate 2+ (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, big_back_wall_id, -1, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, big_left_up_wall_id, -1, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, big_right_up_wall_id, -1, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((int, big_bottom_wall_id, -1, , "id of boundary ; coordinate 2- (default value is ok if aabbWalls are appended BEFORE spheres.)"))
			((Real, big_length0, 0, , "Reference size for strain definition. See :yref:`TriaxialStressController::height`"))
			((Real, big_width0, 0, , "Reference size for strain definition. See :yref:`TriaxialStressController::width`"))
			((Real, move_dis_threshold, 0.01, , "threshold for approaching the target fast."))
			((Real, m_goalx, 0.01, , "prescribed strain rate on axis 1, as defined by :yref:`TriaxialStressController::stressMask`"))
			//
			//@@@
			((bool, z_servo, true, , "stress control or strain control in the z direction"))
			((bool, stop, false, , "stop generalConsolidation for "))
			((bool, loop_Conso, true, , "stop generalConsolidation for "))
			((Real, scalars, 15, , "for loop_Conso use "))
			((bool, vibrate_gc, false, , "for Gc real time vibrate "))
			((Real, real_time, 6e-3, , "for Gc real time vibrate "))
			((int, vibrate_count, 0, , "for Gc real time vibrate "))
			((Real, Poro_ratio_stop, 0.7, , "for loop_Conso stop use "))
			((Real, Porosity_e, 0., , "for loop_Conso Porosity_e "))
			((Real, Porosity_e0, 0., , "for loop_Conso Porosity_e0 "))
			((Real, Porosity_first_GC, 0., , "for loop_Conso Porosity of first  generalConsolidation over"))
			((vector<Real>,Poro_ratio, , , "for loop_Conso use "))
			((bool, Conso_Active, true, , "stop generalConsolidation for "))
			((Real, unbf_tol, 0.01, , "tolerance of UnbalancedForce"))
			((Real, UnbalancedForce, 1, , "mean resultant forces divided by mean contact force"))
			((unsigned int, iterate_num_max, 100, , "interval of restart counting iterate_num"))
			((unsigned int, solve_num_max, 2000, , "interval of restart counting solve_num"))
			((Real, f_threshold, 0.01, , "threshold for approaching the target fast."))
			((double, gain_alpha, 0.5, , "alpha during geting gains"))
			((double, target_strain, 0.15, , "the target axial strain during shear"))
			((unsigned int, echo_interval, 1000, , "interval of iteration for printing run info."))
			((std::string, file, "recordData.txt", , "Name of file to save to; must not be empty."))
			//@@@


			((Real, goalx, 0.01, , "prescribed strain rate on axis 1, as defined by :yref:`TriaxialStressController::stressMask`"))
			((Real, goalz, 0, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))
			((unsigned int, savedata_interval, 1000, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))
			((Real, max_vel, 1, , "Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher. max_vel can be set to infinity in many cases, but sometimes helps stabilizing packings. Based on this value, different maxima are computed for each axis based on the dimensions of the sample, so that if each boundary moves at its maximum velocity, the strain rate will be isotropic (see e.g. :yref:`TriaxialStressController::max_vel1`)."))
			,
			/* extra initializers */
			,
			/* constructor */
			x_area = y_area = z_area = z_area_inShear = 0;//@@@@
			firstRun = true;
			firstRun2 = true;
			firstMove =true;
			Poro_ratio.resize(2);
			Poro_ratio.push_back(0.);
			Poro_ratio.push_back(0.);



			,
			.def_readonly("boxVolume", &ShearEngine::boxVolume, "Total packing volume.")
			.def_readonly("Init_boxVolume", &ShearEngine::Init_boxVolume, "Total packing volume.")
			.def_readonly("particlesVolume", &ShearEngine::particlesVolume, "Total volume of particles (clumps and :yref:`dynamic<Body::dynamic>` spheres). |ycomp|")
			.def_readonly("spheresVolume", &ShearEngine::particlesVolume, "Shorthand for :yref:`TriaxialStressController::particlesVolume`")
			.def("getStress", &ShearEngine::getStress, "get stress within the assembly.")
			.def("computePV", &ShearEngine::computeParticlesVolume, "computePV ParticlesVolume .")
			)
			// clang-format on
			DECLARE_LOGGER;
	};
	REGISTER_SERIALIZABLE(ShearEngine);

} // namespace yade
