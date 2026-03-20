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
#include <pkg/dem/ElasticContactLaw.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/common/ElastMat.hpp>
#include <core/InteractionLoop.hpp>

//#ifdef YADE_OPENGL
//#include <pkg/common/GLDrawFunctors.hpp>
//#endif

namespace yade { // Cannot have #include directive inside.

//!##################	SHAPES   #####################

	class MembraneNode : public Shape {
	public:
		MembraneNode(Real _radius)
	        : radius(_radius)
		{
		 createIndex();
		}
		virtual ~MembraneNode() {};
		//virtual ~MembraneNode();
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(MembraneNode, Shape, "GridConnection shape (see [Effeindzourou2016]_, [Bourrier2013]_). Component of a grid designed to link two :yref:`GridNodes<GridNode>`. It is highly recommended to use :yref:`yade.gridpfacet.gridConnection` to generate correct :yref:`GridConnections<GridConnection>`.",
			((Real,radius,NaN,,"Radius [m]")),
			createIndex();, /*ctor*/
			/*py*/
		);
		// clang-format on
		REGISTER_CLASS_INDEX(MembraneNode, Shape);
	};
	REGISTER_SERIALIZABLE(MembraneNode);
	//############################bound
	class Bo1_MembraneNode_Aabb : public BoundFunctor {
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*) override;
		FUNCTOR1D(MembraneNode);
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Bo1_MembraneNode_Aabb, BoundFunctor, "Functor creating :yref:`Aabb` from a :yref:`GridConnection`.",
			((Real, aabbEnlargeFactor, ((void)"deactivated", -1), , "Relative enlargement of the bounding box; deactivated if negative."))
		);
		// clang-format on
	};
	REGISTER_SERIALIZABLE(Bo1_MembraneNode_Aabb);

	//!##################	Contact Geometry   #####################
	class MembraneMat : public FrictMat {
	public:
		virtual ~MembraneMat() {};
		/// Serialization
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(MembraneMat, FrictMat, "Material description extending :yref:`FrictMat` with cohesive properties and rotational stiffnesses. For use e.g. with :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
			//((bool, isMembrane, true, , "Whether this body can form possibly cohesive interactions (if true and depending on other parameters such as :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys.setCohesionNow`)."))
			//((Real, knfactor, 1, , "Dimensionless rolling stiffness."))
			//((Real, ksfactor, 1, , "Dimensionless rolling stiffness."))
			((Real, Kn, 5e5, , "node normal stiffness "))
			((Real, Ks, 5e5, , "node shearstiffness ")),
			createIndex();
			
			
		);
		// clang-format on
		/// Indexable
		REGISTER_CLASS_INDEX(MembraneMat, FrictMat);
	};
	REGISTER_SERIALIZABLE(MembraneMat);
	//!			O-O
	class NodeFrictPhys : public FrictPhys {
	public:
		virtual ~NodeFrictPhys() = default;
		//Vector3r getRotStiffness() const;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(NodeFrictPhys, FrictPhys, "Version of :yref:`FrictPhys` with a rotational stiffness",
			//((bool, useNodestiffness, true, , "node normal stiffness "))
			((Real, initD, 0., , "initial 初始穿透距离."))
			((Real, normalLimit, 1e5, , "法向强度."))
			((Real, Kn, 0, , "node normal stiffness "))
			((Real, Ks, 0, , "node shearstiffness ")),
			createIndex()
		);
		// clang-format on
		REGISTER_CLASS_INDEX(NodeFrictPhys, FrictPhys);
	};
	REGISTER_SERIALIZABLE(NodeFrictPhys);

	class Ip2_MembraneMat_MembraneMat_NodeFrictPhys : public IPhysFunctor {
	public:
		void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Ip2_MembraneMat_MembraneMat_NodeFrictPhys, IPhysFunctor, "Create a :yref:`FrictPhys` from two :yref:`FrictMats<FrictMat>`. The compliance of one sphere under point load is defined here as $1/(E.D)$, with $E$ the stiffness of the sphere and $D$ its diameter. The compliance of the contact itself is taken as the sum of compliances from each sphere, i.e. $1/(E_1.D_1)+1/(E_2.D_2)$ in the general case, or $2/(E.D)$ in the special case of equal sizes and equal stiffness. Note that summing compliances is equivalent to summing the harmonic average of stiffnesses. This reasoning is applied in both the normal and the tangential directions (as in e.g. [Scholtes2009a]_), hence the general form of the contact stiffness:\n\n $k = \\frac{E_1D_1*E_2D_2}{E_1D_1+E_2D_2}=\\frac{k_1*k_2}{k_1+k_2}$, with $k_i=E_iD_i$.\n\n In the above equation $E_i$ is taken equal to :yref:`FrictMat::young` of sphere $i$ for the normal stiffness, and :yref:`FrictMat::young` $\\times$ :yref:`ElastMat::poisson` for the shear stiffness. In the case of a contact between a :yref:`ViscElMat` and a :yref:`FrictMat`, be sure to set :yref:`FrictMat::young` and :yref:`FrictMat::poisson`, otherwise the default value will be used.\n\n The contact friction is defined according to :yref:`Ip2_FrictMat_FrictMat_FrictPhys::frictAngle` (minimum of the two materials by default).",
			((Real, NormalStrength,1e5 , , "Ind."))
			((Real, ShearStrength,1e5 , , "Ind."))
		);
		FUNCTOR2D(MembraneMat, MembraneMat);
		// clang-format on
	};
	REGISTER_SERIALIZABLE(Ip2_MembraneMat_MembraneMat_NodeFrictPhys);

	

	//!##################	IGeom Functors   #####################

	//!			O-O 接触几何
	class Ig2_MembraneNode_MembraneNode_ScGeom6D : public Ig2_Sphere_Sphere_ScGeom {
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
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Ig2_MembraneNode_MembraneNode_ScGeom6D, Ig2_Sphere_Sphere_ScGeom, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
			((bool, updateRotations, true, , "Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
			((bool, creep, false, , "Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."))
		);
		// clang-format on
		FUNCTOR2D(MembraneNode, MembraneNode);
		// needed for the dispatcher, even if it is symmetric
		DEFINE_FUNCTOR_ORDER_2D(MembraneNode, MembraneNode);
	};
	REGISTER_SERIALIZABLE(Ig2_MembraneNode_MembraneNode_ScGeom6D);
	
	class Ig2_MembraneNode_Sphere_ScGeom : public IGeomFunctor {
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
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Ig2_MembraneNode_Sphere_ScGeom, IGeomFunctor, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between 				two :yref:`GridNode<GridNode>`, including relative rotations.",
		((Real,interactionDetectionFactor,1,,"Enlarge both radii by this factor (if >1), to permit creation of distant interactions.\n\nInteractionGeometry will be computed when 			interactionDetectionFactor*(rad1+rad2) > distance.\n\n.. note::\n\t This parameter is functionally coupled with :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor`, which will create 				larger bounding boxes and should be of the same value."))
		((bool,avoidGranularRatcheting,true,,"Define relative velocity so that ratcheting is avoided. It applies for sphere-sphere contacts. It eventualy also apply for sphere-emulating 			interactions (i.e. convertible into the ScGeom type), if the virtual sphere's motion is defined correctly (see e.g. :yref:`Ig2_Sphere_ChainedCylinder_CylScGeom`).\n\n"
		"Short explanation of what we want to avoid :\n\n"
		"Numerical ratcheting is best understood considering a small elastic cycle at a contact between two grains : assuming b1 is fixed, impose this displacement to b2 :\n\n"
  		"#. translation *dx* in the normal direction\n"))
		);
		// clang-format on
		FUNCTOR2D(MembraneNode, Sphere);
		// needed for the dispatcher, even if it is symmetric
		DEFINE_FUNCTOR_ORDER_2D(MembraneNode, Sphere);
	};
	REGISTER_SERIALIZABLE(Ig2_MembraneNode_Sphere_ScGeom);


	//!##################	Laws   #####################
	class Law2_ScGeom6D_NodeFrictPhys_MembraneLaw : public LawFunctor {
	public:
		bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom6D_NodeFrictPhys_MembraneLaw, LawFunctor, "Law between a frictional :yref:`GridConnection` and a frictional :yref:`Sphere`. Almost the same than :yref:`Law2_ScGeom_FrictPhys_CundallStrack`, but the force is divided and applied on the two :yref:`GridNodes<GridNode>` only.",
			((bool, useK, true, , "Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
			, ,
			);
		// clang-format on
		FUNCTOR2D(ScGeom6D, NodeFrictPhys);
	};
	REGISTER_SERIALIZABLE(Law2_ScGeom6D_NodeFrictPhys_MembraneLaw);
}
