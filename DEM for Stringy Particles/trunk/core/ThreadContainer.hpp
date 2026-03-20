#pragma once

#include <lib/serialization/Serializable.hpp>
//#include <core/Body.hpp>
#include <core/State.hpp>
#include <pkg/common/Thread.hpp>
//#include <boost/iterator/filter_iterator.hpp>
//#include <pkg/common/Line.hpp>
namespace yade {

class ThreadContainer : public Serializable { //// scene中 new NodeNeighbor 堆区创建
		//private:
public:
	typedef std::vector<shared_ptr<Thread>> ContainerT;//

	typedef ContainerT::iterator iterator;
	void clear() { Threads.clear();}
	iterator                begin() { return Threads.begin(); }
	iterator                end() { return Threads.end(); }
	int                  size() const { return Threads.size(); }
	//bool exists(int id) const { return ((id >= 0) && (id < Threads.size())); }
	int insert(shared_ptr<Thread> T) {
		Threads.push_back(T);
		T->id = Threads.size();
		return T->id;
	}
	shared_ptr<Thread>&       operator[](unsigned int id);
	const shared_ptr<Thread>& operator[](unsigned int id) const;
	bool exists(unsigned int id) const;
	bool erase(unsigned int id);
	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(ThreadContainer, Serializable, "Storage node information.",
		((ContainerT, Threads, , , "Current timestep for integration."))
		((long, iterBorn, -1, , "Step number at which the Node was added to simulation."))
		((Real, timeBorn, -1, , "Step number at which the Node was added to simulation."))
		//((Real, timeBorn, -1, , "Time at which the Node was added to simulation."))
		/*ctor*/
		,
		,
		/* py */
	);
};
REGISTER_SERIALIZABLE(ThreadContainer);
	
}
