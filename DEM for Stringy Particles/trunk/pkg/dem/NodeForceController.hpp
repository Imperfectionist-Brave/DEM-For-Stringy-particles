#pragma once

#include <lib/base/Math.hpp>
#include <core/Scene.hpp>
#include <core/BodyContainer.hpp>
//#include <core/NodeNeighbor.hpp>
#include <pkg/common/BoundaryController.hpp>
#include <preprocessing/dem/SpherePack.hpp>
#include<fstream>
#include<string>
#include <boost/array.hpp>
						//bool erase(Body::id_t id, bool eraseClumpMembers);
#define SOLE_GAIN
namespace yade { // Cannot have #include directive inside.
		class Scene;
		class State;
		class NodeNeighborContainer;
		SpherePack Sp;
		class NodeForceController : public BoundaryController {
		private:

		bool                  firstRun;//填充
		bool                  firstRun1;//固结
		bool                  firstRun2;//剪切
		int Num;//膜节点数量
		//shared_ptr<NodeNeighborContainer>
		//shared_ptr<Body>*  Nodelist;//pointers of paticals of Node
		shared_ptr<Body>  top_wall;//上边界墙
		shared_ptr<Body> bottom_wall;//下边界墙
		shared_ptr<Body> third_wall;//辅助压紧墙边界
		std::ofstream out;
		std::vector<shared_ptr<NodeNeighbor>> top_nodNei;//上边界节点存储
		std::vector<shared_ptr<NodeNeighbor>> bottom_nodNei;//下边界节点存储
		std::vector<shared_ptr<NodeNeighbor>> free_nodNei;
		/*	bool                  firstRun;
			bool                  firstRun2;*/

		inline const Vector3r getForce(Scene* rb, Body::id_t id) const
		{
			return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */
		}
		inline void addForce(Scene* rb, Body::id_t id, const Vector3r& f) const
		{
			// sync thread storage of ForceContainer
			//scene->forces.sync();
			return rb->forces.addForce(id, f); /* needs sync, which is done at the beginning of action */
		}

		public:
		//enum { left_down_wall = 0, right_down_wall, left_up_wall, right_up_wall, left_cover_wall, right_cover_wall, bottom_wall, top_wall, front_wall, back_wall };//@@@@
		//void updateBoxSize();
#ifdef SOLE_GAIN
		double  wszz_top, wszz_bottom;
#endif
		double  wszz1, wszz_top1;
		double  wszz;//stress in the three directions
		Real /*avg_rstiff,*/ avg_zstiffb, avg_zstifft, avg_zstifft1;
		Real z_area;
		Vector3r min,max;

		unsigned int iterate_num, solve_num;
		unsigned int rampNum;//chunck number for accelerating walls

		virtual ~NodeForceController();
		//void getNode();//box->body
		//void getBox();//box->body

		void initialize();//初始化函数
		void updateSampleSize();
		void deletbody();//删除边界外body

		void action() override;//执行接口

		Real ComputeUnbalancedForce(bool maxUnbalanced = false);

		void get_gainz();
		void get_gainz1();
		void servo(double dt);
		void servotop(double dt);
		void servo1(double dt);
		//void servo_cylinder(double dt);
		void consol_ss();
		void consol_ss1();

		void consol_ff_cylinder();
		void fix_node();
		void free_node();

		//void flexible_consolidate();//固结
		//void flexible_shear();//剪切
		void fillSample();//填充试样
		void movetop();//移动上边界
		void generalConsolidation();//固结

		void quiet_system();

		Real gain_zb, gain_zt, gain_zt1;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(NodeForceController, BoundaryController,
			"An engine maintaining constant stresses or constant strain rates on some boundaries of a parallepipedic packing. The stress/strain control is defined for each axis using :yref:`TriaxialStressController::stressMask` (a bitMask) and target values are defined by goal1,goal2, and goal3. The sign conventions of continuum mechanics are used for strains and stresses (positive traction)."
			"\n\n.. note::\n\t The algorithms used have been developed initialy for simulations reported in [Chareyre2002a]_ and [Chareyre2005]_. They have been ported to Yade in a second step and used in e.g. [Kozicki2008]_,[Scholtes2009b]_,[Jerier2010b]."
			,
			//@@@
			/*((bool, z_servo, true, , "stress control or strain control in the z direction"))
			((Real, unbf_tol, 0.01, , "tolerance of UnbalancedForce"))
			((Real, UnbalancedForce, 1, , "mean resultant forces divided by mean contact force"))
			((unsigned int, iterate_num_max, 100, , "interval of restart counting iterate_num"))
			((unsigned int, solve_num_max, 2000, , "interval of restart counting solve_num"))
			((Real, f_threshold, 0.01, , "threshold for approaching the target fast."))
			((double, gain_alpha, 0.5, , "alpha during geting gains"))
			((double, target_strain, 0.15, , "the target axial strain during shear"))
			((unsigned int, echo_interval, 1000, , "interval of iteration for printing run info."))*/
			//@@@
			((Real, unbf_tol, 0.01, , "tolerance of UnbalancedForce"))


			((unsigned int, echo_interval, 2500, , "interval of iteration for printing run info."))
			//((unsigned int, ramp_interval, 10, , "chunck bins for accelerating walls gradually."))//deprecated
			//((unsigned int, ramp_chunks, 2, , "chuncks for accelerating walls gradually."))//deprecated

			((double, target_strain, 0.15, , "the target axial strain during shear"))
			((double, gain_alpha, 0.5, , "alpha during geting gains"))
			((Real, f_threshold, 0.05, , "threshold for approaching the target fast."))//阈值

			((double, wall_radius, 20e-3, , "the radius of the cylindrical wall"))//直径
			((unsigned int, iterate_num_max, 100, , "interval of restart counting iterate_num"))
			((unsigned int, solve_num_max, 2000, , "interval of restart counting solve_num"))
			((Real, UnbalancedForce, 1, , "mean resultant forces divided by mean contact force"))
			((bool, dofilling, false, , "wrap pressure force control"))//填充试样开关
			((bool, movewallOnly, true, , "wrap pressure force control"))//仅控制上下边界位移
			((bool, moving, false, , "wrap pressure force control"))//移动上边界开关
			//((bool, free_fill, false, , "wrap pressure force control"))//不限制膜自由度填充
			((bool, fill, true, , "wrap pressure force control"))//不限制膜自由度填充
			((bool, z_servo, true, , "Force control"))
			((bool, keepRigid, true, , "keep the flexible wall rigid? Usually used during consolidation."))//固结时保持柔性墙的刚性
			((Real, height, 0, , "height"))
			((Real, height1, 0, , "height"))//辅助面至底面高度
			((Real, height0, 0, , "initial height"))
			((int, cows, -1, , "Membrane height num"))
			((int, cols, -1, , "Membrane weight num"))
			((Real, wrap_pressure, 1e5, , "wrap pressure with 10kpa"))//围压设定100000
			((Real, goalz, 1e5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//竖向应力
			((Real, goalzv, 0.5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//竖向加载速度系数
			((Real, goalz1, 1e5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//填充试样设定应力
			((Real, max_vel, 10, , "Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher."))
			((Real, heightgap, -1, , "."))


			((int, boxid1, -1, , "bottom id."))
			((int, boxid2, -1, , "top id."))
			((int, boxid3, -1, , "third id."))

			//((unsigned int, savedata_interval, 5000, , "the interval steps between two savedata operations"))
			,
			/* extra initializers */
			,
			/* constructor 初始化模块*/
			//x_area = y_area = z_area = z_area_inShear = 0;//@@@@
			z_area = 0;
		firstRun = true;
		firstRun1 = true;
		firstRun2 = true;
		rampNum = 0;
		//firstRun2 = true;
		,
			/*.def_readonly("boxVolume", &TriaxialStressController2::boxVolume, "Total packing volume.")
			.def_readonly("Init_boxVolume", &TriaxialStressController2::Init_boxVolume, "Total packing volume.")
			.def_readonly("particlesVolume", &TriaxialStressController2::particlesVolume, "Total volume of particles (clumps and :yref:`dynamic<Body::dynamic>` spheres). |ycomp|")
			.def_readonly("spheresVolume", &TriaxialStressController2::particlesVolume, "Shorthand for :yref:`TriaxialStressController::particlesVolume`")*/
			)
			// clang-format on
			DECLARE_LOGGER;
		};
		REGISTER_SERIALIZABLE(NodeForceController);

} // namespace yade
namespace yade { // Cannot have #include directive inside.
	class Scene;
		class State;
		class NodeNeighborContainer;
		//const shared_ptr<NodeNeighborContainer> proxee;//根据上下边界构造膜节点
		//SpherePack Sp;
		class NodeForceController1 : public BoundaryController {
		private:

			bool                  firstRun;//删除body
			bool                  firstRun1;//固结
			bool                  firstRun2;//剪切
			int Num;//膜节点数量
			//shared_ptr<NodeNeighborContainer>
			//shared_ptr<Body>*  Nodelist;//pointers of paticals of Node
			shared_ptr<Body>  top_wall;//上边界墙  序号1
			shared_ptr<Body> bottom_wall;//下边界墙 序号0
			shared_ptr<Body> side_wall;//下边界墙 序号0
			std::ofstream out;
			std::vector<shared_ptr<NodeNeighbor>> top_nodNei;//上边界节点存储
			std::vector<shared_ptr<NodeNeighbor>> bottom_nodNei;//下边界节点存储
			std::vector<shared_ptr<NodeNeighbor>> free_nodNei;
			/*	bool                  firstRun;
				bool                  firstRun2;*/

			inline const Vector3r getForce(Scene* rb, Body::id_t id) const
			{
				return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */
			}
			inline void addForce(Scene* rb, Body::id_t id, const Vector3r& f) const
			{
				// sync thread storage of ForceContainer
				//scene->forces.sync();
				return rb->forces.addForce(id, f); /* needs sync, which is done at the beginning of action */
			}

		public:
			//enum { left_down_wall = 0, right_down_wall, left_up_wall, right_up_wall, left_cover_wall, right_cover_wall, bottom_wall, top_wall, front_wall, back_wall };//@@@@
			//void updateBoxSize();
#ifdef SOLE_GAIN
		double  wszz_top, wszz_bottom;
#endif
		double  wfxx_bottom, wfxx_bottomn,wfxx_side;
		double  wszz1, wszz_top1;
		double  wszz;//stress in the three directions
		Real /*avg_rstiff,*/ avg_zstiffb, avg_zstifft;
		Real z_area;
		//Vector3r min, max;

		Real particlesVolume;//体积
		Real Init_SampleVolume;
		Real Current_SampleVolume;

		unsigned int iterate_num, solve_num;
		unsigned int rampNum;//chunck number for accelerating walls

		virtual ~NodeForceController1();
		//void getNode();//box->body
		//void getBox();//box->body

		void initialize();//初始化函数
		void updateSampleSize();
		void deletbody();//删除圆柱facet body

		void action() override;//执行接口

		Real ComputeUnbalancedForce(bool maxUnbalanced = false);

		void get_gainz();
		void servo(double dt);
		//void servo_cylinder(double dt);
		void consol_ss();//计算应力
		void consol_f();//计算力
		//void generateMembrane();//固结后生成膜颗粒 python生成

		void consol_ff_cylinder();
		void fix_node();
		void free_node();

		//void flexible_consolidate();//固结
		//void flexible_shear();//剪切
		//void fillSample();//填充试样
		//void movetop();//移动上边界
		void generalConsolidation();//固结

		void quiet_system();

		Real gain_zb, gain_zt, gain_zt1;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(NodeForceController1, BoundaryController,
			"An engine maintaining constant stresses or constant strain rates on some boundaries of a parallepipedic packing. The stress/strain control is defined for each axis using :yref:`TriaxialStressController::stressMask` (a bitMask) and target values are defined by goal1,goal2, and goal3. The sign conventions of continuum mechanics are used for strains and stresses (positive traction)."
			"\n\n.. note::\n\t The algorithms used have been developed initialy for simulations reported in [Chareyre2002a]_ and [Chareyre2005]_. They have been ported to Yade in a second step and used in e.g. [Kozicki2008]_,[Scholtes2009b]_,[Jerier2010b]."
			,
			//@@@
			/*((bool, z_servo, true, , "stress control or strain control in the z direction"))
			((Real, unbf_tol, 0.01, , "tolerance of UnbalancedForce"))
			((Real, UnbalancedForce, 1, , "mean resultant forces divided by mean contact force"))
			((unsigned int, iterate_num_max, 100, , "interval of restart counting iterate_num"))
			((unsigned int, solve_num_max, 2000, , "interval of restart counting solve_num"))
			((Real, f_threshold, 0.01, , "threshold for approaching the target fast."))
			((double, gain_alpha, 0.5, , "alpha during geting gains"))
			((double, target_strain, 0.15, , "the target axial strain during shear"))
			((unsigned int, echo_interval, 1000, , "interval of iteration for printing run info."))*/
			//@@@
			((Real, unbf_tol, 0.01, , "tolerance of UnbalancedForce"))


			((unsigned int, echo_interval, 2500, , "interval of iteration for printing run info."))
			//((unsigned int, ramp_interval, 10, , "chunck bins for accelerating walls gradually."))//deprecated
			//((unsigned int, ramp_chunks, 2, , "chuncks for accelerating walls gradually."))//deprecated

			((double, target_strain, 0.15, , "the target axial strain during shear"))
			((double, gain_alpha, 0.5, , "alpha during geting gains"))
			((Real, f_threshold, 0.05, , "threshold for approaching the target fast."))//阈值

			((Real, NodeRadius, 1e-3, , "top id."))//膜节点半径大小
			((Real, wall_radius, 20e-3, , "the radius of the cylindrical wall"))//直径
			((unsigned int, iterate_num_max, 100, , "interval of restart counting iterate_num"))
			((unsigned int, solve_num_max, 2000, , "interval of restart counting solve_num"))
			((Real, UnbalancedForce, 1, , "mean resultant forces divided by mean contact force"))
			((bool, dofilling, false, , "wrap pressure force control"))//填充试样开关
			((bool, movewallOnly, true, , "wrap pressure force control"))//仅控制上下边界压缩
			//((bool, free_fill, false, , "wrap pressure force control"))//不限制膜自由度填充
			((bool, freez, false, , "wrap pressure force control"))//不限制上下膜z轴自由度
			//((bool, fill, False, , "wrap pressure force control"))//不限制膜自由度填充
			((bool, z_servo, true, , "Force control"))
			((bool, force, false, , "Force control"))//横移实验
			((bool, useside, false, , "Force control"))//横移实验
			((bool, keepRigid, true, , "keep the flexible wall rigid? Usually used during consolidation."))//固结时保持柔性墙的刚性
			((Real, height, 0, , "height"))//两面高度
			((Real, height0, 0, , "initial height"))//初始高度
			((Real, posx0, 0, , "initial height"))//边墙初始x位置
			((int, cows, -1, , "Membrane height num"))
			((int, cols, -1, , "Membrane weight num"))
			((Real, wrap_pressure, 1e5, , "wrap pressure with 10kpa"))//围压设定100000
			((Real, goalz, 1e5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//竖向应力
			((Real, goalzv, 0.5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//竖向加载速度系数
			((Real, goalxv, 0.05, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//横向移动边界速度系数
			((Real, goalz1, 1e5, , "prescribed stress rate on axis 3, as defined by :yref:`TriaxialStressController::stressMask`"))//填充试样设定应力
			((Real, max_vel, 10, , "Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher."))


			((int, boxid1, -1, , "bottom id."))
			((int, boxid2, -1, , "top id."))
			((int, boxid3, -1, , "side id."))

			((Real, cylinderRadius, 0, , "third id."))
			//((double, wall_radius, 20e-3, , "the radius of the cylindrical wall"))//直径

			//((unsigned int, savedata_interval, 5000, , "the interval steps between two savedata operations"))
			,
			/* extra initializers */
			,
			/* constructor 初始化模块*/
			//x_area = y_area = z_area = z_area_inShear = 0;//@@@@
			z_area = 0;
			firstRun = true;
			firstRun1 = true;
			firstRun2 = true;
			rampNum = 0;
			wfxx_bottomn = 0;
			wfxx_bottom = 0;
			wfxx_side=0;
			//firstRun2 = true;
			,
			.def("deletbody",&NodeForceController1::deletbody,"Add a GridConnection to the GridNode.")
			/*.def_readonly("boxVolume", &TriaxialStressController2::boxVolume, "Total packing volume.")
			.def_readonly("Init_boxVolume", &TriaxialStressController2::Init_boxVolume, "Total packing volume.")
			.def_readonly("particlesVolume", &TriaxialStressController2::particlesVolume, "Total volume of particles (clumps and :yref:`dynamic<Body::dynamic>` spheres). |ycomp|")
			.def_readonly("spheresVolume", &TriaxialStressController2::particlesVolume, "Shorthand for :yref:`TriaxialStressController::particlesVolume`")*/
			)
			// clang-format on
			DECLARE_LOGGER;
	};
	REGISTER_SERIALIZABLE(NodeForceController1);
} // namespace yade
