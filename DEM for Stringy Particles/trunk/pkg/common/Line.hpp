/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@grenoble-inp.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/* TABLE OF CONTENT, and the minimum you need to understand.
- 2 new shapes for grids : GridNode (vertices) and GridConnection (edges)
- 2 new contact geometries :
    * GridNodeGeom6D to handle GridNode-GridNode contacts (the internal behaviour of the grid)
    * ScGridCoGeom to handle Sphere-GridConnection contacts (the interaction between the grid and an external sphere)
        Note : the Sphere-Grid contacts are always handled by the GridConnections, but the forces are applied on the related GridNodes.
        Note : there is no contact between two GridConnections, they must be linked with GridNodes.
- The 2 related Ig2 :
    * Ig2_GridNode_GridNode_GridNodeGeom6D (doing almost the same than Ig2_Sphere_Sphere_ScGeom6D)
    * Ig2_Sphere_GridConnection_ScGridCoGeom (the biggest part of the code, it handles contact detection and history when a sphere is sliding on the grid over consecutive GridConnections)
- The Law2_ScGridCoGeom_FrictPhys_CundallStrack who handles the elastic frictional Sphere-GridConnection contact. The GridNode-GridNode law is Law2_ScGeom6D_CohFrictPhys_CohesionMoment by inheritance.
*/

#pragma once
#include "Sphere.hpp"
#include <core/Body.hpp>
#include <core/Dispatching.hpp>
#include <pkg/dem/CohesiveFrictionalContactLaw.hpp>
#include <pkg/dem/ElasticContactLaw.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/common/Thread.hpp>//******6-7
#ifdef YADE_OPENGL
#include <pkg/common/GLDrawFunctors.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

//!##################	SHAPES   #####################

//**************************************************************************************************12-26
class LineConnection : public Shape {
public:
	virtual ~LineConnection();
	Real                     getLength() const;
	Vector3r                 getSegment() const;
	//void                     addPFacet(shared_ptr<Body> PF);
	//vector<shared_ptr<Body>> getPFacets() const { return pfacetList; }
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(LineConnection,Shape,"GridConnection shape (see [Effeindzourou2016]_, [Bourrier2013]_). Component of a grid designed to link two :yref:`GridNodes<GridNode>`. It is highly recommended to use :yref:`yade.gridpfacet.gridConnection` to generate correct :yref:`GridConnections<GridConnection>`.",
		((Real,length0,NaN,,"initial length0 [m]"))
		((Real,radius,NaN,,"Radius [m]"))
		((shared_ptr<Body> , node1 , ,,"First :yref:`Body` the GridConnection is connected to."))
		((shared_ptr<Body> , node2 , ,,"Second :yref:`Body` the GridConnection is connected to."))
		((shared_ptr<Body> , sphere , ,,"connected sphere id."))
		((bool, periodic, false,,"true if two nodes from different periods are connected."))
		//((vector<shared_ptr<Body> >,pfacetList,,Attr::hidden,"List of :yref:`PFacet<PFacet>` the GridConnection is connected to."))
		((Vector3i , cellDist , Vector3i(0,0,0),,"Distance of bodies in cell size units, if using periodic boundary conditions. Note that periodic boundary conditions for GridConnections have not yet been fully implemented.")),
		createIndex();, /*ctor*/
		/*py*/  
		//.def("addPFacet",&LineConnection::addPFacet,(boost::python::arg("Body")),"Add a PFacet to the GridConnection.")
		//.def("getPFacets",&LineConnection::getPFacets,"get list of linked PFacets.")
	);
	// clang-format on
	REGISTER_CLASS_INDEX(LineConnection, Shape);
};
REGISTER_SERIALIZABLE(LineConnection);

//!   0
class LineNode : public Shape {
public:
	virtual ~LineNode();
	void                     addConnection(shared_ptr<Body> GC);
	vector<shared_ptr<Body>> getConnections() const { return ConnList; }//获取节点id
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(LineNode, Shape, "GridNode shape, component of a grid.\nTo create a Grid, place the nodes first, they will define the spacial discretisation of it. It is highly recommended to use :yref:`yade.gridpfacet.gridNode` to generate correct :yref:`GridNodes<GridNode>`. Note that the GridNodes should only be in an Interaction with other GridNodes. The Sphere-Grid contact is only handled by the :yref:`GridConnections<GridConnection>`.",
		//((vector<shared_ptr<Body> >, pfacetList, , Attr::hidden, "List of :yref:`PFacets<PFacet>` the GridConnection is connected to."))
		((Real,radius,NaN,,"Radius [m]"))
		((vector<shared_ptr<Body>>, ConnList, , Attr::hidden, "List of :yref:`GridConnections<GridConnection>` the GridNode is connected to."))//节点包含的连接体（最多两个）
		((bool, NsIntrs, false, , "num of real interactions."))//与连接颗粒接触时为真
		((int, id1, -1, , "true if the node connected with sphere."))
		((int, id2,-1, , "true if the node connected with sphere."))
		((int, sphid1, -1, , "fixed sphere id1."))
		((int, sphid2, -1, , "fixed sphere id2."))
		((bool, isMiddle, false, , "true if the node connected with sphere."))
		((bool, isFirst, false, , "true if the node connected with sphere."))
		((bool, isLast, false, , "true if the node connected with sphere.")),
		/*ctor*/
		createIndex(); ,
		/*py*/
		.def("addConnection", &LineNode::addConnection, (boost::python::arg("Body")), "Add a GridConnection to the GridNode.")
		.def("getConnections", &LineNode::getConnections, "get list of linked :yref:`GridConnection`'s.")
	);
	// clang-format on
	REGISTER_CLASS_INDEX(LineNode, Shape);
};
REGISTER_SERIALIZABLE(LineNode);

//!##################	Contact Geometry   #####################

//!			o-o
class LineNodeGeom : public GenericSpheresContact {  //连接点之间的接触形状
private:
	Vector3r orthonormal_axis; //rotation vector in contact plane
public:
	virtual ~LineNodeGeom();
	Real &radius1, &radius2;
	
	void precompute(
	        const State&                   rbp1,
	        const State&                   rbp2,
	        const Scene*                   scene,
	        const shared_ptr<Interaction>& c,
	        const Vector3r&                currentNormal,
	        bool                           isNew);
	Vector3r& rotate(Vector3r& tangentVector) const;
	const Vector3r& shearIncrement() const { return shearInc; }
	Vector3r getIncidentVel(const State* rbp1, const State* rbp2 ) const;
	
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(LineNodeGeom, GenericSpheresContact, "Geometry of a :yref:`GridNode`-:yref:`GridNode` contact. Inherits almost everything from :yref:`ScGeom6D`.",
		((shared_ptr<Body>, connectionBody, , , "Reference to the :yref:`GridNode` :yref:`Body` who is linking the two :yref:`GridNodes<GridNode>`."))//连接两点的线体,用于更新cylinder位置渲染
		((bool, isFirst, false, , "true if the node connected with sphere."))
		((bool, isLast, false, , "true if the node connected with sphere."))
		((Real, penetrationDepth, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
		((Vector3r,shearInc,Vector3r::Zero(),(Attr::noSave|Attr::readonly),"Shear displacement increment in the last step")),
		/* extra initializers */((radius1,GenericSpheresContact::refR1)) ((radius2,GenericSpheresContact::refR2)),
		/* ctor */ createIndex();orthonormal_axis=Vector3r::Zero();,
		/* py */
		);
	// clang-format on
	REGISTER_CLASS_INDEX(LineNodeGeom, GenericSpheresContact);
};
REGISTER_SERIALIZABLE(LineNodeGeom);

class LineNodeGeom2 : public GenericSpheresContact {  //连接点之间的接触形状
public:
	virtual ~LineNodeGeom2();
	Real &radius1, &radius2;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(LineNodeGeom2, GenericSpheresContact, "Geometry of a :yref:`GridNode`-:yref:`GridNode` contact. Inherits almost everything from :yref:`ScGeom6D`.",
		((int,id3,0,,"id of the first :yref:`GridNode`. |yupdate|"))
		((int,id4,0,,"id of the second :yref:`GridNode`. |yupdate|"))
		((Real, relPos, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
		((Real, penetrationDepth, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)")),
		/* extra initializers */((radius1,GenericSpheresContact::refR1)) ((radius2,GenericSpheresContact::refR2)),
		/* ctor */ createIndex(); ,
		/* py */
		);
	// clang-format on
	REGISTER_CLASS_INDEX(LineNodeGeom2, GenericSpheresContact);
};
REGISTER_SERIALIZABLE(LineNodeGeom2);


class NodeScGeom2 : public GenericSpheresContact {  //连接点 颗粒之间的接触形状
public:
	virtual ~NodeScGeom2();
	Real &radius1, &radius2;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(NodeScGeom2, GenericSpheresContact, "Geometry of a :yref:`GridNode`-:yref:`GridNode` contact. Inherits almost everything from :yref:`ScGeom6D`.",
		((Real, relPos, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
		((Real, penetrationDepth, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)")),
		/* extra initializers */((radius1,GenericSpheresContact::refR1)) ((radius2,GenericSpheresContact::refR2)),
		/* ctor */ createIndex(); ,
		/* py */
		);
	// clang-format on
	REGISTER_CLASS_INDEX(NodeScGeom2, GenericSpheresContact);
};
REGISTER_SERIALIZABLE(NodeScGeom2);

//!			o-o
class Ig2_LineNode_LineNode_LineNodeGeom : public IGeomFunctor{
public:
	virtual bool
		go(const shared_ptr<Shape>&       cm1,
			const shared_ptr<Shape>&       cm2,
			const State&                   state1,
			const State&                   state2,
			const Vector3r&                /*shift2*/,
			const bool&                    force,
			const shared_ptr<Interaction>& c) override;
	virtual bool goReverse(
		const shared_ptr<Shape>&       cm1,
		const shared_ptr<Shape>&       cm2,
		const State&                   state1,
		const State&                   state2,
		const Vector3r&                shift2,
		const bool&                    force,
		const shared_ptr<Interaction>& c) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_LineNode_LineNode_LineNodeGeom, IGeomFunctor, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
		((bool, updateRotations, true, , "Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
		((bool, creep, false, , "Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."))
	);
	// clang-format on
	FUNCTOR2D(LineNode, LineNode);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(LineNode, LineNode);
};
REGISTER_SERIALIZABLE(Ig2_LineNode_LineNode_LineNodeGeom);

class Ig2_LineConnection_Sphere_LineNodeGeom2 : public IGeomFunctor{
public:
	virtual bool
		go(const shared_ptr<Shape>&       cm1,
			const shared_ptr<Shape>&       cm2,
			const State&                   state1,
			const State&                   state2,
			const Vector3r&                /*shift2*/,
			const bool&                    force,
			const shared_ptr<Interaction>& c) override;
	virtual bool goReverse(
		const shared_ptr<Shape>&       cm1,
		const shared_ptr<Shape>&       cm2,
		const State&                   state1,
		const State&                   state2,
		const Vector3r&                shift2,
		const bool&                    force,
		const shared_ptr<Interaction>& c) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_LineConnection_Sphere_LineNodeGeom2, IGeomFunctor, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
		((bool, updateRotations, true, , "Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
		((bool, creep, false, , "Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."))
	);
	// clang-format on
	FUNCTOR2D(LineConnection, Sphere);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(LineConnection, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_LineConnection_Sphere_LineNodeGeom2);

class Ig2_LineNode_Sphere_NodeScGeom2 : public IGeomFunctor{
public:
	shared_ptr<InteractionContainer> interactions;
	virtual bool
		go(const shared_ptr<Shape>&       cm1,
			const shared_ptr<Shape>&       cm2,
			const State&                   state1,
			const State&                   state2,
			const Vector3r&                /*shift2*/,
			const bool&                    force,
			const shared_ptr<Interaction>& c) override;
	virtual bool goReverse(
		const shared_ptr<Shape>&       cm1,
		const shared_ptr<Shape>&       cm2,
		const State&                   state1,
		const State&                   state2,
		const Vector3r&                shift2,
		const bool&                    force,
		const shared_ptr<Interaction>& c) override;
	
	/*bool found(int &id1,int&id2){
		if (id1 > id2) swap(id1, id2);
		shared_ptr<Interaction> in;
		const shared_ptr<Body>& b1((*(scene->bodies))[id1]);
		Body::MapId2IntrT::iterator I(b1->intrs.find(id2));
		if (I != b1->intrs.end()) in=I->second;//shared_ptr<Interaction>
		return in->isReal();
	}*/
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_LineNode_Sphere_NodeScGeom2, IGeomFunctor, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
		((bool, updateRotations, true, , "Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
		//((shared_ptr<InteractionContainer>,interactions,NULL,Attr::hidden,"All interactions between bodies."))
		((bool, creep, false, , "Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."))
		
	);
	// clang-format on
	FUNCTOR2D(LineNode, Sphere);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(LineNode, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_LineNode_Sphere_NodeScGeom2);

//!##################	Bounds   #####################

class Bo1_LineConnection_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*) override;
	FUNCTOR1D(LineConnection);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Bo1_LineConnection_Aabb, BoundFunctor, "Functor creating :yref:`Aabb` from a :yref:`GridConnection`.",
		((Real, aabbEnlargeFactor, ((void)"deactivated", -1), , "Relative enlargement of the bounding box; deactivated if negative."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Bo1_LineConnection_Aabb);
//!##################	Linenode Bounds   #####################2024-3-20
class Bo1_LineNode_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*) override;
	FUNCTOR1D(LineNode);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Bo1_LineNode_Aabb, BoundFunctor, "Functor creating :yref:`Aabb` from a :yref:`GridConnection`.",
		((Real, aabbEnlargeFactor, ((void)"deactivated", -1), , "Relative enlargement of the bounding box; deactivated if negative."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Bo1_LineNode_Aabb);

//!##################	Rendering   #####################
#ifdef YADE_OPENGL
class Gl1_LineConnection : public GlShapeFunctor {
private:
	//static int glCylinderList;
	//void subdivideTriangle(Vector3r& v1,Vector3r& v2,Vector3r& v3, int depth);
	void drawCylinder(bool wire, Real radius, Real length, const Quaternionr& shift = Quaternionr::Identity());
	//void initGlLists(void);
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	void out(Quaternionr q);
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_LineConnection, GlShapeFunctor, "Renders :yref:`Cylinder` object",
		((bool, wire, false, , "Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((bool, glutNormalize, true, , "Fix normals for non-wire rendering"))
		((int, glutSlices, 8, , "Number of cylinder slices."))
		((int, glutStacks, 4, , "Number of cylinder stacks."))
	);
	// clang-format on
	RENDERS(LineConnection);
};
REGISTER_SERIALIZABLE(Gl1_LineConnection);

class Gl1_LineNode : public GlShapeFunctor {
private:
	// for stripes
	static vector<Vector3r> vertices, faces;
	static int              glStripedSphereList;
	static int              glGlutSphereList;
	void                    subdivideTriangle(Vector3r& v1, Vector3r& v2, Vector3r& v3, int depth);
	//Generate GlList for GLUT sphere
	void initGlutGlList();
	//Generate GlList for sliced spheres
	void initStripedGlList();
	//for regenerating glutSphere or glutTorus list if needed
	static Real prevQuality;
	//for regenerating glutSphere or glutTorus list if needed
	static string prevDisplayMode;
	//for regenerating glutTorus list if needed
	static char prevCircleAllowedRotationAxis;

public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_LineNode,GlShapeFunctor,"Renders :yref:`Sphere` object",
		((Real,quality,1.0,,"Change discretization level of spheres. quality>1  for better image quality, at the price of more cpu/gpu usage, 0<quality<1 for faster rendering. If mono-color spheres are displayed (:yref:`Gl1_Sphere::stripes` = False), quality mutiplies :yref:`Gl1_Sphere::glutSlices` and :yref:`Gl1_Sphere::glutStacks`. If striped spheres are displayed (:yref:`Gl1_Sphere::stripes` = True), only integer increments are meaningfull : quality=1 and quality=1.9 will give the same result, quality=2 will give finer result."))
		((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((bool,stripes,false,,"In non-wire rendering, show stripes clearly showing particle rotation."))
		((bool,localSpecView,true,,"Compute specular light in local eye coordinate system."))
		((int,glutSlices,12,(Attr::noSave | Attr::readonly),"Base number of sphere slices, multiplied by :yref:`Gl1_Sphere::quality` before use); not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
		((int,glutStacks,6,(Attr::noSave | Attr::readonly),"Base number of sphere stacks, multiplied by :yref:`Gl1_Sphere::quality` before use; not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
		((bool,circleView,false,,"For 2D simulations : display tori instead of spheres, so they will appear like circles if the viewer is looking in the right direction. In this case, remember to disable perspective by pressing \"t\"-key in the viewer."))
		((Real,circleRelThickness,0.2,,"If :yref:`Gl1_Sphere::circleView` is enabled, this is the torus diameter relative to the sphere radius (i.e. the circle relative thickness)."))
		((char,circleAllowedRotationAxis,'z',,"If :yref:`Gl1_Sphere::circleView` is enabled, this is the only axis ('x', 'y' or 'z') along which rotation is allowed for the 2D simulation. It allows right orientation of the tori to appear like circles in the viewer. For example, if circleAllowedRotationAxis='x' is set, blockedDOFs=\"YZ\" should also be set for all your particles."))
	);
	// clang-format on
	RENDERS(LineNode);
};

REGISTER_SERIALIZABLE(Gl1_LineNode);

#endif
//******************
} // namespace yade
