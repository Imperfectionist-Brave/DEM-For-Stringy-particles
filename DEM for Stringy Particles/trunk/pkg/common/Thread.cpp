#include "Thread.hpp"
#include <pkg/common/Line.hpp>
#include <pkg/dem/LineContactLaw.hpp>
namespace yade {
	YADE_PLUGIN((Thread));
	void Thread::initial()
		{
			Node_num = nodesStates.size();
			LineConnect_num = nodesStates.size() - 1;
			assert(Node_num > 1);
			local_coord1 = nodesStates[0]->pos - sphereStates[0]->pos;
			local_coord2 = nodesStates[Node_num-1]->pos - sphereStates[1]->pos;
			length1=length01=computeL();
			directF=nodesStates[1]->pos - nodesStates[0]->pos;
			directL=nodesStates[Node_num - 2]->pos - nodesStates[Node_num-1]->pos;
			
		}
	void Thread::printlength(){
		length_avg=(length1+length2)/2;
		cout<<"length01"<<length01<<"  "<<"length02"<<length02<<"  "<<"length1"<<length1<<"  "<<"length2"<<length2<<endl;
		cout<<"diff length1"<<length1-length_avg<<"  "<<"diff length2"<<length2-length_avg<<endl;
		//cout<<"diff length is"<<length_avg-length0<<endl;
	}
	 void Thread::addSphereState(const py::list & pylist)
	{
		 	[[maybe_unused]]int Length = boost::python::len(pylist);
		 	assert(Length == 2);
			shared_ptr<State> st1=py::extract<shared_ptr<State>>(pylist[0]);
			shared_ptr<State> st2=py::extract<shared_ptr<State>>(pylist[1]);
			sphereStates.push_back(st1);
			sphereStates.push_back(st2);
	}
	void   Thread::addNodeState(const py::list & pylist)
		{	
		 	int Length = boost::python::len(pylist);
		 	assert(Length >=2);
		 	//vector<shared_ptr<State>> St;
		 	for(int i=0;i<Length;i++){
		 		shared_ptr<State> st=py::extract<shared_ptr<State>>(pylist[i]);
		 		nodesStates.push_back(st);
		 	}
		}
		
	void   Thread::addSphereId(vector<int> ID)
		{
			assert(ID.size()== 2);id1 = ID[0];id2 = ID[1];
		}
		
		
	void   Thread::addSphereId2(vector<int> ID)
		{
			assert(ID.size()== 2);
			sphid.clear();
			int _id1 = ID[0];
			int _id2 = ID[1];
			sphid.push_back(_id1);
			sphid.push_back(_id2);
			
		}
		
	void   Thread::addNodeId(vector<int> ID)//9-8
		{
			int num = ID.size();
			assert(num>=2);
			nodeid.clear();
			for(int i =0;i<num;i++){
			
				nodeid.push_back(ID[i]);
			}
		}
		
	void   Thread::addLineId(vector<int> ID)//9-8
		{
			int num = ID.size();
			assert(num>=2);
			lineid.clear();
			for(int i =0;i<num;i++){
				lineid.push_back(ID[i]);
			}
			
		}
		
	void Thread::addGeoms(const py::list & pylist)
	{
	
		 int Length = boost::python::len(pylist);
		 assert(Length >=2);
		 //vector<shared_ptr<IGeom>> ig;
		 for(int i=0;i<Length;i++){
		 	shared_ptr<IGeom> ig=py::extract<shared_ptr<IGeom>>(pylist[i]);
		 	geoms.push_back(ig);
		 }
	}
	
	void Thread::addPhys(const py::list & pylist)
	{
	
		 int Length = boost::python::len(pylist);
		 assert(Length >=2);
		// vector<shared_ptr<State>> St;
		 for(int i=0;i<Length;i++){
		 	shared_ptr<IPhys> ip=py::extract<shared_ptr<IPhys>>(pylist[i]);
		 	phys.push_back(ip);
		 }
	}
	Real Thread::outlength0(){
	 	return length0;
	}
	Real Thread::outlengthAll(){
		Real l=0;
		int segNum=geoms.size();
		if(segNum>0){
			for (int i=0;i<segNum;i++){
				LineNodeGeom*     geom = YADE_CAST<LineNodeGeom*>(geoms[i].get());
				l-=geom->penetrationDepth;
			}
		}
		return l;
	}
	
	py::list Thread::get_geoms()
	{
		py::list ret;
		int num=geoms.size();
		for(int i=0;i<num;i++){
			ret.append(geoms[i]);
		}
		return ret;
		
	}
	
	Real Thread::get_f(){
		assert(phys.size()>0);
		LinePhys* _phys=NULL;
		_phys = YADE_CAST<LinePhys*>(phys[0].get());
		return _phys->Fn;
	}
	
}



