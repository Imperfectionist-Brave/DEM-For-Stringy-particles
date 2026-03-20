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
#include <pkg/common/Sphere.hpp>
#include <core/Body.hpp>
#include <core/Dispatching.hpp>
#include <core/Material.hpp>
#include <pkg/dem/ElasticContactLaw.hpp>
#include <pkg/dem/FrictPhys.hpp>

#include <pkg/dem/DemXDofGeom.hpp>
#include <core/IGeom.hpp>
#include<pkg/dem/Node.hpp>//5-29
#ifdef YADE_OPENGL
#include <pkg/common/GLDrawFunctors.hpp>
#endif

namespace yade { // Cannot have #include directive inside.
//!##################	MATERIAL   #####################
	class ThreadMat : public Material {
	public:
		virtual ~ThreadMat();
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(ThreadMat, Material, "Elastic material with Coulomb friction.",
			((Real, Kn, 1e3, , "Normal stiffness (N/m)."))
			,
			/*ctor*/ createIndex();
		);
		REGISTER_CLASS_INDEX(ThreadMat, Material);
	};
	REGISTER_SERIALIZABLE(ThreadMat);
	//!##################	PHYS   #####################
	class ThreadPhys : public FrictPhys {
	public:
		virtual ~ThreadPhys();
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(ThreadPhys, FrictPhys, "Simple elastic material with friction for volumetric constitutive laws",
			,
			/*ctor*/ createIndex();
		);
		REGISTER_CLASS_INDEX(ThreadPhys, FrictPhys);
	};
	REGISTER_SERIALIZABLE(ThreadPhys);
	//!************************************** Ip2**********************************
	class Ip2_ThreadMat_FrictMat_ThreadPhys : public IPhysFunctor {
	public:
		void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;

		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_ThreadMat_FrictMat_ThreadPhys, IPhysFunctor, "",
			,
			);
		// clang-format on
		FUNCTOR2D(ThreadMat,FrictMat);
	};
	REGISTER_SERIALIZABLE(Ip2_ThreadMat_FrictMat_ThreadPhys);

	//!##################	SHAPES   #####################
	class Segment : public Shape {//lineconnection combine the thread
	
	public:
		virtual ~Segment();
		/// Normals of edges
		Real                     getLength() const;
		Real                     getlength() const;
		Vector3r                 getSegment() const;
		void postLoad(Segment&);
		Vector3r getNormal() const;
		Vector3r gettan1() const;
		Vector3r gettan2() const;
		void updatePos() const;
		void updateBodyPos() const;
		void applyForce(Scene* scene,Vector3r f);
		inline void addForce(Scene* rb, Body::id_t idt, const Vector3r& f) const
		{
			// sync thread storage of ForceContainer
			//rb->forces.sync();
			//return 
			rb->forces.addForce(idt, f);
			rb->forces.addTorque(idt, (node->pos-sphere->state->pos).cross(f)); /* needs sync, which is done at the beginning of action */
		}
		bool isConected();//是否为连接段
		bool isStretching();//是否处于伸长变形状态
		// clang-format off
		//shared_ptr<Node> node1, node2;//柔索段节点
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(Segment, Shape, "PFacet (particle facet) geometry (see [Effeindzourou2016]_, [Effeindzourou2015a]_). It is highly recommended to use the helper functions in :yref:`yade.gridpfacet` (e.g., gridpfacet.pfacetCreator1-4) to generate correct :yref:`PFacet<PFacet>` elements.",
			((vector<shared_ptr<Node>>, nodes, , , "First :yref:`Body` the Pfacet is connected to."))
			((shared_ptr<Body>, sphere, , , "third :yref:`Body` the sphere is connected to."))//可能连接的球体
			((Body::id_t,id,-1,Attr::readonly,"Unique id of this body."))//sphere id
			((shared_ptr<State>, state, , , "third :yref:`Body` the sphere is connected to."))//主体body state指针
			((shared_ptr<Node>, node, , , "third :yref:`Body` the sphere is connected to."))//可能连接的球体
			((shared_ptr<Node>, node1, , Attr::triggerPostLoad, "third :yref:`Body` the sphere is connected to."))
			((shared_ptr<Node>, node2, , Attr::triggerPostLoad, "third :yref:`Body` the sphere is connected to."))
			((vector<Vector3r>, vertices, vector<Vector3r>(2, Vector3r(NaN, NaN, NaN)), , "Vertex positions in local coordinates."))
			((Vector3r, Tangential1, Vector3r(NaN, NaN, NaN), , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))//内部力施加切线方向
			((Vector3r, Tangential2, Vector3r(NaN, NaN, NaN), , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
			((Vector3r, normal, Vector3r(NaN, NaN, NaN), (Attr::readonly | Attr::noSave), "PFacet's normal (in local coordinate system)"))
			//((Vector3r, initial_pos, Vector3r(NaN, NaN, NaN), (Attr::triggerPostLoad), "PFacet's normal (in local coordinate system)"))
			((Vector3r, local_coord, Vector3r(NaN, NaN, NaN), , "PFacet's normal (in local coordinate system)"))
			((Real, radius, 8e-5, , "Thread's radius"))
			((Real, length, -1, , "Thread's current length"))
			((Real, length0, -1, , "Thread's initial length"))
			((bool, incontact, false, , "default min number of node"))
			((bool, stretch, false, , "default min number of node"))
			((bool, broken, false, , "segment breaken"))
			((bool, isFirst, false, , "First in Thread"))
			((bool, isLast, false, , "Last in Thread"))
			((bool, middle, true, , "Last in Thread"))
			((bool, momentRotation, true, , "default min number of node"))
			((Real, area, NaN, (Attr::readonly | Attr::noSave), "PFacet's area"))
			((Vector3i, cellDist, Vector3i(0, 0, 0), , "Distance of bodies in cell size units, if using periodic boundary conditions. Note that periodic boundary conditions for PFacets have not yet been fully implemented."))
			,
			/* ctor */ createIndex();
		);
		// clang-format on
		DECLARE_LOGGER;

		REGISTER_CLASS_INDEX(Segment, Shape);
	};
	REGISTER_SERIALIZABLE(Segment);
	//****************************************GEOM*******************************
	//class SeScGeom : public GenericSpheresContact {//?
	class SeScGeom : public IGeom {//Segment with sphere geom
	public:
		virtual ~SeScGeom();
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(SeScGeom, IGeom,
			"Class representing :",
			((Vector3r, normal,Vector3r(NaN, NaN, NaN) , , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
			((Vector3r, normal1, , , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))//外部力方向
			((Vector3r, normal2, , , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
			((Vector3r, Tangential1, Vector3r(NaN, NaN, NaN), , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))//内部力方向
			((Vector3r, Tangential2, Vector3r(NaN, NaN, NaN), , "Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
			((Vector3r, contactPoint1, Vector3r(NaN, NaN, NaN), , "some reference point for the interaction (usually in the middle). |ycomp|"))
			((Vector3r, contactPoint2, Vector3r(NaN, NaN, NaN), , "some reference point for the interaction (usually in the middle). |ycomp|"))
			((Real, penetrationDepth, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
			((Real, penetrationDepth1, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
			((Real, penetrationDepth2, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
			((Real, Deformation, NaN, , "Reference radius of particle #1. |ycomp|"))//变形长度
			((bool, contact1, false, , "Reference radius of particle #1. |ycomp|"))//接触状态1
			((bool, contact2, false, , "Reference radius of particle #1. |ycomp|"))
			((int,id1,0,,"id of the first :yref:`GridNode`. |yupdate|"))
			((int,id2,0,,"id of the second :yref:`GridNode`. |yupdate|"))
			((Real, relPos, , , "Reference radius of particle #1. |ycomp|"))	
			((Real, refR1, , , "Reference radius of particle #1. |ycomp|"))
			((Real, refR2, , , "Reference radius of particle #2. |ycomp|")),
			/* extra initializers */
			,
			/* ctor */
			createIndex();,
			/* py */
			);
		// clang-format on
		REGISTER_CLASS_INDEX(SeScGeom, IGeom);
	};
	REGISTER_SERIALIZABLE(SeScGeom);
	//************************************** IG2**********************************
	class Ig2_Segment_Sphere_SeScGeom : public IGeomFunctor {
	public:
		virtual bool
			go(const shared_ptr<Shape>&       cm1,
				const shared_ptr<Shape>&       cm2,
				const State&                   state1,
				const State&                   state2,
				const Vector3r&                shift2,
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
		void penetrationDepth2();
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Ig2_Segment_Sphere_SeScGeom, IGeomFunctor, "Create/update a :yref:`ScGridCoGeom6D` instance representing the geometry of a contact point between a :yref:`GricConnection` and a :yref:`Sphere` including relative rotations.",
			((Real, interactionDetectionFactor, 1, , "Enlarge both radii by this factor (if >1), to permit creation of distant interactions."))
		);
		// clang-format on
		FUNCTOR2D(Segment, Sphere);
		DEFINE_FUNCTOR_ORDER_2D(Segment, Sphere);
	};
	REGISTER_SERIALIZABLE(Ig2_Segment_Sphere_SeScGeom);
	//!************************************** Law2**********************************
	class Law2_SeScGeom_ThreadPhys_CundallStrack : public LawFunctor {//连接体之间、连接体与球体之间
	public:
		bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_SeScGeom_ThreadPhys_CundallStrack, LawFunctor, "Law between a cohesive frictional :yref:`GridConnection` and a cohesive frictional :yref:`Sphere`. Almost the same than :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`, but THE ROTATIONAL MOMENTS ARE NOT COMPUTED.",
			((bool, neverErase, false, , "Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
			, ,
			);
		// clang-format on
		FUNCTOR2D(SeScGeom, ThreadPhys);
	};
	REGISTER_SERIALIZABLE(Law2_SeScGeom_ThreadPhys_CundallStrack);


	//class Thread : public Shape {//lineconnection combine the thread
	//public:
	//	virtual ~Thread();
	//	/// Normals of edges
	//	Vector3r ne[3];
	//	/// Inscribing cirle radius
	//	Real icr;
	//	/// Length of the vertice vectors
	//	Real vl[3];
	//	/// Unit vertice vectors
	//	Vector3r vu[3];
	//	void initial() {
	//		assert(lineconnctions.size()>0);
	//		First_lineconnction = lineconnctions[0];
	//		Last_lineconnction = lineconnctions[lineconnctions.size()-1];
	//	}
	//	// clang-format off
	//	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Thread, Shape, "PFacet (particle facet) geometry (see [Effeindzourou2016]_, [Effeindzourou2015a]_). It is highly recommended to use the helper functions in :yref:`yade.gridpfacet` (e.g., gridpfacet.pfacetCreator1-4) to generate correct :yref:`PFacet<PFacet>` elements.",
	//		((shared_ptr<Body>, node1, , , "First :yref:`Body` the Pfacet is connected to."))
	//		((shared_ptr<Body>, node2, , , "Second :yref:`Body` the Pfacet is connected to."))
	//		//((shared_ptr<Body>, node3, , , "third :yref:`Body` the Pfacet is connected to."))//柔索平面计算??
	//		((shared_ptr<Body>, conn1, , , "First :yref:`Body` the Pfacet is connected to."))
	//		((shared_ptr<Body>, conn2, , , "Second :yref:`Body` the Pfacet is connected to."))
	//		((shared_ptr<Body>, conn3, , , "third :yref:`Body` the Pfacet is connected to."))
	//		((shared_ptr<Body>, sphere1, , , "third :yref:`Body` the sphere1 is connected to."))
	//		((shared_ptr<Body>, sphere2, , , "third :yref:`Body` the sphere2 is connected to."))
	//		/*((shared_ptr<Body>, First_lineconnction, , , "third :yref:`Body` the Pfacet is connected to."))
	//		((shared_ptr<Body>, Last_lineconnction, , , "third :yref:`Body` the Pfacet is connected to."))*/
	//		((Vector3r, FirstConnctionPos, , , "third :yref:`Body` the Pfacet is connected to."))//第一球连接位置
	//		((Vector3r, LastConnctionPos, , , "third :yref:`Body` the Pfacet is connected to."))
	//		((vector<shared_ptr<Body>>, lineconnctions, , , "third :yref:`Body` the Pfacet is connected to."))
	//		((vector<shared_ptr<Body>>, linenode, , , "third :yref:`Body` the Pfacet is connected to."))
	//		((Vector3r, normal, Vector3r(NaN, NaN, NaN), (Attr::readonly | Attr::noSave), "PFacet's normal (in local coordinate system)"))
	//		((Real, radius, -1, , "Thread's radius"))
	//		((Real, length, -1, , "Thread's current length"))
	//		((Real, length0, -1, , "Thread's initial length"))
	//		((int, line_min_num, 1, , "default min number of line"))
	//		((int, node_min_num, 2, , "default min number of node"))
	//		((bool, stretch, false, , "default min number of node"))
	//		((bool, momentRotation, true, , "default min number of node"))
	//		((Real, area, NaN, (Attr::readonly | Attr::noSave), "PFacet's area"))
	//		((Vector3i, cellDist, Vector3i(0, 0, 0), , "Distance of bodies in cell size units, if using periodic boundary conditions. Note that periodic boundary conditions for PFacets have not yet been fully implemented."))
	//		,
	//		/* ctor */ createIndex();
	//	);
	//	// clang-format on
	//	DECLARE_LOGGER;
	//
	//	REGISTER_CLASS_INDEX(Thread, Shape);
	//};
	//REGISTER_SERIALIZABLE(Thread);

	//!##################	Contact Geometry   #####################

	//!##################	Bounds   #####################

	class Bo1_Segment_Aabb : public BoundFunctor {
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*) override;
		FUNCTOR1D(Segment);
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Bo1_Segment_Aabb, BoundFunctor, "Functor creating :yref:`Aabb` from a :yref:`GridConnection`.",
			((Real, aabbEnlargeFactor, ((void)"deactivated", -1), , "Relative enlargement of the bounding box; deactivated if negative."))
		);
		// clang-format on
	};
	REGISTER_SERIALIZABLE(Bo1_Segment_Aabb);

	//!##################	Linenode Bounds   #####################2024-3-20
	//!##################	Rendering   #####################

#ifdef YADE_OPENGL
	class Gl1_Segment : public GlShapeFunctor {
	private:
		//static int glCylinderList;
		//void subdivideTriangle(Vector3r& v1,Vector3r& v2,Vector3r& v3, int depth);
		void drawCylinder(bool wire, Real radius, Real length, const Quaternionr& shift = Quaternionr::Identity());
		//void initGlLists(void);
	public:
		void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
		void out(Quaternionr q);
		// clang-format off
		YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_Segment, GlShapeFunctor, "Renders :yref:`Cylinder` object",
			((bool, wire, false, , "Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
			((bool, glutNormalize, true, , "Fix normals for non-wire rendering"))
			((int, glutSlices, 8, , "Number of cylinder slices."))
			((int, glutStacks, 4, , "Number of cylinder stacks."))
		);
		// clang-format on
		RENDERS(Segment);
	};
	REGISTER_SERIALIZABLE(Gl1_Segment);
#endif
	//******************
} // namespace yade
