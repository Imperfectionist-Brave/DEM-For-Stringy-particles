// 2004 © Olivier Galizzi <olivier.galizzi@imag.fr>
// 2004,2019 © Janek Kozicki <cosurgi@berlios.de>
// 2010 © Václav Šmilauer <eudoxos@arcig.cz>
// 2019 © Anton Gladky <gladk@debian.org>
// 2018 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>

#pragma once

#include <lib/serialization/Serializable.hpp>
#include <core/Body.hpp>
#include <boost/iterator/filter_iterator.hpp>

namespace yade { // Cannot have #include directive inside.

class Body;
class InteractionContainer;

#ifdef YADE_OPENMP
#define YADE_PARALLEL_FOREACH_BODY_BEGIN(b_, bodies)                                                                                                           \
	bodies->updateRealBodies();                                                                                                                            \
	const vector<Body::id_t>& realBodies = bodies->realBodies;                                                                                             \
	const bool                redirect   = bodies->useRedirection;                                                                                         \
	const Body::id_t          _sz(redirect ? realBodies.size() : bodies->size());                                                                          \
	_Pragma("omp parallel for") for (int kParallelForeachIndexCounter = 0; kParallelForeachIndexCounter < _sz; kParallelForeachIndexCounter++)             \
	{                                                                                                                                                      \
		if (not redirect and not(*bodies)[kParallelForeachIndexCounter]) continue;                                                                     \
		b_((*bodies)[redirect ? realBodies[kParallelForeachIndexCounter] : kParallelForeachIndexCounter]);
#else
#define YADE_PARALLEL_FOREACH_BODY_BEGIN(b_, bodies)                                                                                                           \
	bodies->updateRealBodies();                                                                                                                            \
	const vector<Body::id_t>& realBodies = bodies->realBodies;                                                                                             \
	const bool                redirect   = bodies->useRedirection;                                                                                         \
	const Body::id_t          _sz(redirect ? realBodies.size() : bodies->size());                                                                          \
	for (int kParallelForeachIndexCounter = 0; kParallelForeachIndexCounter < _sz; kParallelForeachIndexCounter++) {                                       \
		if (not redirect and not(*bodies)[kParallelForeachIndexCounter]) continue;                                                                     \
		b_((*bodies)[redirect ? realBodies[kParallelForeachIndexCounter] : kParallelForeachIndexCounter]);
#endif
#define YADE_PARALLEL_FOREACH_BODY_END() }

/*
Container of bodies implemented as flat std::vector.
The iterator will silently skip null body pointers which may exist after removal. The null pointers can still be accessed via the [] operator.
When the container is sparse (after erasing bodies or mpi decomposition) the container supports redirection, i.e. looping on the indices of the non-void elements of the container (realBodies) rather than on elements of the container. This is integrated in the above YADE_PARALLEL_FOREACH macros.

Any alternative implementation should use the same API.
*/
class BodyContainer : public Serializable {
private:
	using ContainerT = std::vector<shared_ptr<Body>>;

public:
	friend class InteractionContainer; // accesses the body vector directly

	//An iterator that will automatically jump slots with null bodies
	struct isNonEmptySharedPtr {
		bool operator()(const shared_ptr<Body>& b) const { return b.operator bool(); }
	};
	using iterator = boost::filter_iterator<isNonEmptySharedPtr, ContainerT::iterator>;

	Body::id_t insert(shared_ptr<Body>);
	Body::id_t insertAtId(shared_ptr<Body> b, Body::id_t candidate);

	// Container operations
	void                    clear();
	iterator                begin();
	iterator                end();
	size_t                  size() const;
	shared_ptr<Body>&       operator[](unsigned int id);
	const shared_ptr<Body>& operator[](unsigned int id) const;

	bool exists(Body::id_t id) const;
	bool erase(Body::id_t id, bool eraseClumpMembers);

	void updateRealBodies();

	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BodyContainer,Serializable,"Standard body container for a scene",
		((ContainerT,body,,,"The underlying vector<shared_ptr<Body> >"))
		((bool,dirty,true,(Attr::noSave|Attr::readonly|Attr::hidden),"true if after insertion/removal of bodies, used only if collider::keepListsShort=true"))
		((bool,checkedByCollider,false,(Attr::noSave|Attr::readonly|Attr::hidden),""))
		((vector<Body::id_t>,insertedBodies,vector<Body::id_t>(),Attr::readonly,"The list of newly bodies inserted, to be used and purged by collider"))
		((vector<Body::id_t>,erasedBodies,vector<Body::id_t>(),Attr::readonly,"The list of erased bodies, to be used and purged by collider"))
		((vector<Body::id_t>,realBodies,vector<Body::id_t>(),Attr::readonly,"Redirection vector to non-null bodies, used to optimize loops after numerous insertion/erase. In MPI runs the list is restricted to bodies and neighbors present in current subdomain."))
		((bool,useRedirection,false,,"true if the scene uses up-to-date lists for boundedBodies and realBodies; turned true automatically 1/ after removal of bodies if :yref:`enableRedirection=True <BodyContainer.enableRedirection>`, and 2/ in MPI execution. |yupdate|"))
		((bool,enableRedirection,true,,"let collider switch to optimized algorithm with body redirection when bodies are erased - true by default"))
		#ifdef YADE_MPI
		((vector<Body::id_t>,subdomainBodies,vector<Body::id_t>(),,"The list of bounded bodies in the subdomain"))
		#endif
		,/*ctor*/,
		.def("updateRealBodies",&BodyContainer::updateRealBodies,"update lists realBodies and subdomainBodies. This function is called automatically by e.g. ForceContainer::reset(), it is safe to call multiple times from many places since if the lists are up-to-date he function will just return.")
		)
	// clang-format on

	DECLARE_LOGGER;

	// mutual exclusion to avoid crashes in the rendering loop
	std::mutex drawloopmutex;

private:
	bool eraseAlreadyLocked(Body::id_t id, bool eraseClumpMembers);
};
REGISTER_SERIALIZABLE(BodyContainer);

} // namespace yade
namespace yade {
	class NodeNeighbor : public Serializable {//保存膜节点邻居的类
		public:
			
			int num;//具有邻居的个数
			bool mask;//计算方向
			Vector3r UnitForce;//单位力
			Body::id_t selfId;//自身的body id 号
			shared_ptr<State> selfstate;//自身位置
			std::vector<Body::id_t> bodyids;//临近ids
			std::vector<shared_ptr<State>> bodyStates;//临近state
			

			/*vector<Vector3r> L;//位置向量序列
			vector<Real> Areas;//面积序列
			vector<Vector3r> Norm;//法向量*/
			//动态数组分配
			Vector3r *L=NULL;
			Real * Areas=NULL;
			Vector3r*Norm=NULL;
			
			void initial();
			void computeL();
			void computeAreaAndNorm();
			 Vector3r computeUnitForce();

			//std::map<int, string &> bodyState;//4或6个
		public:
			//NodeNeighbor(){};
			//virtual ~NodeNeighbor(){};
			NodeNeighbor(Body::id_t _id, shared_ptr<State> _self, bool _mask,std::vector<Body::id_t> _ids,std::vector<shared_ptr<State>> _States) {//构造函数

				selfId = _id;
				selfstate = _self;
				mask = _mask;
				bodyids = _ids;
				bodyStates = _States;
				initial();
			}
			//virtual ~NodeNeighbor() {};
			// only NodeContainer can set the id of a Node
			friend class NodeNeighborContainer;
			YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(NodeNeighbor, Serializable, "Storage node information.",
				((int, id, -1, , "Current timestep for integration."))
				//((int, mask, -1, , "Current timestep for integration."))//0 逆时针计算 1 顺时针计算法向量
				,/*ctor*/
				,/* py */
				.def_readwrite("mat",&NodeNeighbor::id,"Shorthand for :yref:`Body::material`")
				);
	};
	REGISTER_SERIALIZABLE(NodeNeighbor);//注册导出类
}
namespace yade {
	class NodeNeighborContainer : public Serializable { //// scene中 new NodeNeighbor 堆区创建
		//private:
			

		public:
			typedef std::vector<shared_ptr<NodeNeighbor>> NeiContainerT;//大小根据参数传递确定
			typedef std::vector<Body::id_t> NeiIdT;//NodeNeighbor id 容器
			typedef std::vector<std::vector<Body::id_t>> Martix2D;//二维 id 容器
			int cows;//颗粒膜行数
			int cols;//列数
			typedef NeiContainerT::iterator iterator;
			void clear();
			iterator                begin();
			iterator                end();
			int                  size() const;
			bool exists(Body::id_t id) const;
			std::vector<Body::id_t> insert(int _cows, int _cols, std::vector<std::vector<Body::id_t>> Matrix);

			//void updateNode();
			void extend();//扩充ids
			void allocat();//分配NodeNeighbor

			YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(NodeNeighborContainer, Serializable, "Storage node information.",
				((NeiContainerT, Neighbor, , , "Current timestep for integration."))
				((NeiIdT, Self, , , "Current timestep for integration."))
				((Martix2D, Ids, , , "Current timestep for integration."))
				((long, iterBorn, -1, , "Step number at which the Node was added to simulation."))
				((Real, timeBorn, -1, , "Step number at which the Node was added to simulation."))
				((vector<Body::id_t>, topNode, vector<Body::id_t>(), Attr::readonly, "The list of newly bodies inserted, to be used and purged by collider"))
				((vector<Body::id_t>, freeNode, vector<Body::id_t>(), Attr::readonly, "The list of newly bodies inserted, to be used and purged by collider"))
				((vector<Body::id_t>, bottomNode, vector<Body::id_t>(), Attr::readonly, "The list of newly bodies inserted, to be used and purged by collider"))
				//((Real, timeBorn, -1, , "Time at which the Node was added to simulation."))
				/*ctor*/
				,
				,
				.def_readwrite("mat",&NodeNeighborContainer::timeBorn,"Shorthand for :yref:`Body::material`")
				/* py */
			);
		};
REGISTER_SERIALIZABLE(NodeNeighborContainer);
	
}
