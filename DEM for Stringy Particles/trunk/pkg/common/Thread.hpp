#pragma once

#include <lib/serialization/Serializable.hpp>
//#include <core/Body.hpp>
#include <core/State.hpp>
#include <core/IGeom.hpp>
#include <pkg/common/NormShearPhys.hpp>
//#include <boost/iterator/filter_iterator.hpp>
//#include <pkg/common/Line.hpp>


namespace yade {
	class Thread : public Serializable {//线类 
	public:

		int Node_num;//节点个数
		int LineConnect_num;//段个数

		//void initial();
		inline Vector3r computeDirect(shared_ptr<State>& st1,shared_ptr<State>& st2)  {return st1->pos - st2->pos;}
		inline Real computenormal(Vector3r v)  {return v.norm();}
		inline Real computeL() 
		{
			Real _length = 0;
			for (int i = 0; i < Node_num-1; ++i) {
				_length+=computenormal(computeDirect(nodesStates[i+1],nodesStates[i]));
			}
			length=_length;
			return _length;
		}
		inline Vector3r computeDirectF() {return nodesStates[1]->pos - nodesStates[0]->pos;}
		inline Vector3r computeDirectL() {return nodesStates[Node_num-2]->pos - nodesStates[Node_num-1]->pos;}
		Vector3r computeNormF() {return nodesStates[0]->pos - sphereStates[0]->pos;}
		Vector3r computeNormL() {return nodesStates[Node_num - 1]->pos - sphereStates[1]->pos;}
		inline void update(shared_ptr<State>& st ,shared_ptr<State>& nst ,Vector3r & v,Vector3r & vv) {
			Quaternionr ori = st->ori;
			Matrix3r Rotation = ori.toRotationMatrix();
			vv= Rotation * v;
			nst->pos=st->pos+vv;
		}
		void updatePos() {
			update(sphereStates[0],nodesStates[0],local_coord1,local_pos1);
			update(sphereStates[1],nodesStates[Node_num-1],local_coord2,local_pos2);
		}
	public:
		virtual ~Thread() {};
		// only NodeContainer can set the id of a Node
		void   addSphereId(vector<int> ID);
		void   addSphereId2(vector<int> ID);
		void   addNodeId(vector<int> ID);
		void   addLineId(vector<int> ID);//9-8
		
		void   addNodeState(const py::list & pylist);
		
		void   addSphereState(const py::list & pylist);
		void   addGeoms(const py::list & pylist);//6-17
		void   addPhys(const py::list & pylist);//6-17
		void initial();
		void printlength();
		Real outlength0();
		Real outlengthAll();
		py::list get_geoms();
		Real get_f();
		friend class ThreadContainer;
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Thread, Serializable, "Storage node information.",
			((Vector3r, local_coord1, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Vector3r, local_coord2, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Vector3r, local_pos1, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Vector3r, local_pos2, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Vector3r, directF, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Vector3r, directL, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Real, radius, 8e-5, , "Thread's radius"))
			((Real, k, 1e5, , "tolerance of UnbalancedForce"))//rousuozhengti k
			((Real, length0, NaN, , "initial length0 [m]"))
			((Real, length01, NaN, , "initial length0 [m]"))
			((Real, length02, NaN, , "initial length0 [m]"))
			((Real, length, NaN, , "initial length0 [m]"))
			((Real, length1, NaN, , "initial length [m]"))
			((Real, length2, NaN, , "initial length [m]"))
			((Real, length_avg, NaN, , "initial length [m]"))
			((shared_ptr<IGeom>, geom1, , , "initial length [m]"))
			((shared_ptr<IGeom>, geom2, , , "initial length [m]"))
			((bool, useallGeom, false, , "self id"))
			((vector<shared_ptr<IGeom>>, geoms, , , "First :yref:`Body` the GridConnection is connected to."))
			((vector<shared_ptr<IPhys>>, phys, , , "First :yref:`Body` the GridConnection is connected to."))
			((shared_ptr<IPhys>, phy1, , , "initial length [m]"))
			((shared_ptr<IPhys>, phy2, , , "initial length [m]"))
			((int, id, -1, , "self id"))
			((int, id1, -1, , "sphere id1"))//node id1
			((int, id2, -1, , "sphere id2"))//node id2
			((vector<int>, sphid, , , "sphere id2"))//sphere id2
			((vector<int>, nodeid, , , "sphere id2"))//node id2
			((vector<int>, lineid, , , "sphere id2"))//node id2
			//((shared_ptr<Body>, First_node1, , , "First :yref:`Body` the GridConnection is connected to."))
			//((shared_ptr<Body>, Last_node1, , , "First :yref:`Body` the GridConnection is connected to."))
			((vector<shared_ptr<State>>, nodesStates, , , "First :yref:`Body` the GridConnection is connected to."))
			((vector<shared_ptr<State>>, sphereStates, , , "First :yref:`Body` the GridConnection is connected to."))
			,/*ctor*/
			,/* py */
			.def("addNodeState", &Thread::addNodeState, (boost::python::arg("State")), "Add a GridConnection to the GridNode.")
			.def("addSphereState", &Thread::addSphereState, (boost::python::arg("State")), "Add a GridConnection to the GridNode.")
			.def("addSphereId", &Thread::addSphereId, (boost::python::arg("int")), "Add a GridConnection to the GridNode.")
			.def("addSphereId2", &Thread::addSphereId2, (boost::python::arg("int")), "Add a GridConnection to the GridNode.")//9-8
			.def("addNodeId", &Thread::addNodeId, (boost::python::arg("int")), "Add a GridConnection to the GridNode.")//9-8
			.def("addLineId", &Thread::addLineId, (boost::python::arg("int")), "Add a GridConnection to the GridNode.")//9-8
			.def("initial", &Thread::initial,(boost::python::arg("void")) , "Add a GridConnection to the GridNode.")
			.def("printlength", &Thread::printlength,(boost::python::arg("void")) , "Add a GridConnection to the GridNode.")
			.def("outlength0", &Thread::outlength0,(boost::python::arg("Real")) , "Add a GridConnection to the GridNode.")
			.def("outlengthAll", &Thread::outlengthAll,(boost::python::arg("Real")) , "Add a GridConnection to the GridNode.")
			.def("addGeoms", &Thread::addGeoms, (boost::python::arg("IGeom")), "Add a GridConnection to the GridNode.")
			.def("addPhys", &Thread::addPhys, (boost::python::arg("IPhys")), "Add a GridConnection to the GridNode.")
			.def("get_geoms", &Thread::get_geoms, (boost::python::arg("none")), "Add a GridConnection to the GridNode.")
			.def("get_f", &Thread::get_f, (boost::python::arg("none")), "Add a GridConnection to the GridNode.")
		);
	};
	REGISTER_SERIALIZABLE(Thread);
}


