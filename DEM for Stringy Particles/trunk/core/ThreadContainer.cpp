#include "ThreadContainer.hpp"
//#include "Body.hpp"
//#include "Clump.hpp"
//#include "Scene.hpp"
//#include <lib/base/LoggingUtils.hpp>
//#ifdef YADE_OPENMP
//#include <omp.h>
//#endif
namespace yade {
	YADE_PLUGIN((ThreadContainer));
	
	shared_ptr<Thread>& ThreadContainer::operator[](unsigned int id) { return Threads[id]; }
	const shared_ptr<Thread>& ThreadContainer::operator[](unsigned int id) const { return Threads[id]; }
	//bool ThreadContainer::exists(unsigned int id) const { return ((id >= 0) && ((size_t)id < Threads.size()) && ((bool)Threads[id])); }
}


