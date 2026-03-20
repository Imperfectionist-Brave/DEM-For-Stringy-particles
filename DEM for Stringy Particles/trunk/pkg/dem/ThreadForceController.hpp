#pragma once
#include <lib/base/Math.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Line.hpp>
#include <core/ThreadContainer.hpp>
//#include <core/NodeNeighbor.hpp>
#include <pkg/common/BoundaryController.hpp>
//#include <preprocessing/dem/SpherePack.hpp>
#include<fstream>
#include<string>
#include <boost/array.hpp>
namespace yade { // Cannot have #include directive inside.
	class Scene;
	class State;
	class Thread;
	class ThreadContainer;
	class ThreadForceController : public BoundaryController {
	private:

		bool                  firstRun;//填充
		std::ofstream out;
		std::vector<shared_ptr<Thread>> threads;//线类
	
		ThreadContainer::iterator _First ;
		ThreadContainer::iterator _End ;
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
		inline void addTorque(Scene* rb, Body::id_t id, const Vector3r& f) const
		{
			// sync thread storage of ForceContainer
			//scene->forces.sync();
			return rb->forces.addTorque(id, f); /* needs sync, which is done at the beginning of action */
		}

	public:
		virtual ~ThreadForceController();
		void initialize();//初始化函数
		void addNodeForce();
		void addNodeForceAllGeom();//6-17
		void addNodeForceAllGeom_test();//6-25
		void updateNodePos();
		void action() override;//执行接口
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThreadForceController, BoundaryController,
			"An engine maintaining constant stres.",
			((Real, young_s,0.6e9 , , "Young's modulus"))
			((Real, k, 1e5, , "tolerance of UnbalancedForce"))//rousuozhengti
			((Real, ksegment, 1000, , "tolerance of UnbalancedForce"))//rousuoduan
			((Real, radius, 8e-5, , "tolerance of UnbalancedForce"))//
			((Real, theta_cos, 0.9, , "tolerance of UnbalancedForce"))//拉伸时角度最大夹角限制
			((int, Num, -1, , "Threads num"))
			((int, GeomsNum, -1, , "Threads geom num"))
			((bool, updateNode, false, , "tolerance of UnbalancedForce"))
			((bool, useallGeom, false, , "self id"))
			((bool, testStretch, false, , "use k to accumulate force"))
			,
			/* extra initializers */
			,
			/* constructor 初始化模块*/
			firstRun = true;
			,
			)
			// clang-format on
			DECLARE_LOGGER;
	};
	REGISTER_SERIALIZABLE(ThreadForceController);

} // namespace yade


