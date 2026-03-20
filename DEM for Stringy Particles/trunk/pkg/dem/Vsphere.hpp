#pragma once
#include <core/Dispatching.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.
	//shape
	class Vsphere : public Shape {
	public:
		Vsphere(Real _radius): radius(_radius)
			        { createIndex();}
		virtual ~Vsphere();
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Vsphere, Shape, "GridNode shape, component of a grid.\nTo create a Grid, place the nodes first, they will define the spacial discretisation of it. It is highly recommended to use :yref:`yade.gridpfacet.gridNode` to generate correct :yref:`GridNodes<GridNode>`. Note that the GridNodes should only be in an Interaction with other GridNodes. The Sphere-Grid contact is only handled by the :yref:`GridConnections<GridConnection>`.",
			//((vector<shared_ptr<Body> >, pfacetList, , Attr::hidden, "List of :yref:`PFacets<PFacet>` the GridConnection is connected to."))
			((Real, radius, NaN, , "Radius [m]"))
			((int, sphid, -1, , "clump's sphere id.")),
			/*ctor*/
			createIndex();,
			/*py*/
		);
		// clang-format on
		REGISTER_CLASS_INDEX(Vsphere, Shape);
	};
	REGISTER_SERIALIZABLE(Vsphere);

	//geom
	class VSphereGeom : public GenericSpheresContact {  //虚球和颗粒接触
	public:
		virtual ~VSphereGeom();
		void updatefixedpoint(
			const State&                   rbp1,
			const State&                   rbp2);
		Real &radius1, &radius2;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(VSphereGeom, GenericSpheresContact, "Geometry of a :yref:`GridNode`-:yref:`GridNode` contact. Inherits almost everything from :yref:`ScGeom6D`.",
			((shared_ptr<Body>, connectionBody, , , "Reference to the :yref:`GridNode` :yref:`Body` who is linking the two :yref:`GridNodes<GridNode>`."))//连接两点的线体,用于更新cylinder位置渲染
			((int, sphid, -1, , "clump's sphere id."))//sphere id
			((Real, penetrationDepth, NaN, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
			((Real, minlength, 0, (Attr::noSave | Attr::readonly), "Penetration distance of spheres (positive if overlapping)"))
			((Vector3r, startV1, Vector3r::Zero(), , "Force application point 1"))//起始点向量1（球为中心的局部坐标系）***
			((Vector3r, startV2, Vector3r::Zero(), , "Force application point 2"))//起始点向量2（球为中心的局部坐标系）***
			((Vector3r, normal1, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))//施加对虚球力方向
			((Vector3r, normal2, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))
			((Vector3r, branch1, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))
			((Vector3r, branch2, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))
			((Vector3r, fixedpoint1, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))//连接点1
			((Vector3r, fixedpoint2, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step"))//连接点2
			((Vector3r, TanPoint, Vector3r::Zero(), (Attr::noSave | Attr::readonly), "Shear displacement increment in the last step")),
			/* extra initializers */((radius1, GenericSpheresContact::refR1)) ((radius2, GenericSpheresContact::refR2)),
			/* ctor */ createIndex(); ,
			/* py */
			);
		// clang-format on
		REGISTER_CLASS_INDEX(VSphereGeom, GenericSpheresContact);
	};
	REGISTER_SERIALIZABLE(VSphereGeom);

	//phys
	class VsPhys : public NormShearPhys {
	public:
		virtual ~VsPhys();
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(VsPhys, NormShearPhys, "An :yref:`interaction physics<IPhys>` that extends :yref:`RotStiffFrictPhys` adding a breakable cohesive nature. Used e.g. by :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
			((Vector3r, _normalForce1, Vector3r::Zero(), , "Normal force after previous step (in global coordinates), as sustained by particle #2 (from particle #1)."))
			((Vector3r, _normalForce2, Vector3r::Zero(), , "Normal force after previous step (in global coordinates), as sustained by particle #2 (from particle #1)."))
			((Real, k, 0, , "normal stiffness k"))//柔索拉伸刚度
			((Real, initD, 0., , "initial length."))
			((Real, normalLimit, 0., , "normal strength."))
			((vector<Vector2r>, displForceValues, , Attr::readonly, "Defines the values for force-displacement curve."))
			((vector<Real>, stiffnessValues, , Attr::readonly, "Defines the values for the various stiffnesses (the elastic stiffness is stored as kn)."))
			,
			createIndex();
		);
		// clang-format on
		/// Indexable
		REGISTER_CLASS_INDEX(VsPhys, NormShearPhys);
	};
	REGISTER_SERIALIZABLE(VsPhys);


	//material
	class VsMat : public ElastMat {
	public:
		virtual ~VsMat();
		void              postLoad(VsMat&);
		/// Serialization
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(VsMat, ElastMat, "Material description extending :yref:`FrictMat` with cohesive properties and rotational stiffnesses. For use e.g. with :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
			//((bool,isCohesive,true,,"Whether this body can form possibly cohesive interactions (if true and depending on other parameters such as :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys.setCohesionNow`)."))
			((Real, diameter, 0.0005, , "Diameter of the single line in [m] (the diameter is used to compute the cross-section area of the wire)."))
			((vector<Vector2r>, strainStressValues, , , "Piecewis"))
			((Real, k, 0, Attr::triggerPostLoad, "normal stiffness k"))//线材料拉伸刚度
			,
			createIndex();
		);
		// clang-format on
		/// Indexable
		REGISTER_CLASS_INDEX(VsMat, ElastMat);
	};
	REGISTER_SERIALIZABLE(VsMat);

	//ip2
	class Ip2_VsMat_FrictMat_VsPhys : public IPhysFunctor {
	public:
		void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_VsMat_FrictMat_VsPhys, IPhysFunctor,
			"Generates cohesive-frictional interactions with moments, used in the contact law :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`. The normal/shear stiffness and friction definitions are the same as in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`, check the documentation there for details.\n\nAdhesions related to the normal and the shear components are calculated from :yref:`CohFrictMat::normalCohesion` ($C_n$) and :yref:`CohFrictMat::shearCohesion` ($C_s$). For particles of size $R_1$,$R_2$ the adhesion will be $a_i=C_i min(R_1,R_2)^2$, $i=n,s$.\n\nTwist and rolling stiffnesses are proportional to the shear stiffness through dimensionless factors alphaKtw and alphaKr, such that the rotational stiffnesses are defined by $k_s \\alpha_i R_1 R_2$, $i=tw\\,r$",
			((Real, normalLimit, 15, , "normal strength."))
			,
			//cohesionDefinitionIteration = -1;//?
			);
		// clang-format on
		FUNCTOR2D(VsMat, FrictMat);
	};
	REGISTER_SERIALIZABLE(Ip2_VsMat_FrictMat_VsPhys);

	//ig2
	class Ig2_Vsphere_Sphere_VSphereGeom : public IGeomFunctor {
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
		Quaternionr rotationQuaternion(Vector3r& axis,Real& theta);
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Ig2_Vsphere_Sphere_VSphereGeom, IGeomFunctor, "Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
			((Real, Tradius, 8e-5, , "Penetration distance of spheres (positive if overlapping)"))//
		);
		// clang-format on
		FUNCTOR2D(Vsphere, Sphere);
		// needed for the dispatcher, even if it is symmetric
		DEFINE_FUNCTOR_ORDER_2D(Vsphere, Sphere);
	};
	REGISTER_SERIALIZABLE(Ig2_Vsphere_Sphere_VSphereGeom);

	//law2
	class Law2_VSphereGeom_VsPhys_VsPM : public LawFunctor {
	public:
		bool                    go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
		// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_VSphereGeom_VsPhys_VsPM, LawFunctor, "Law for linear.",
			((bool, useK, false, , "use k to accumulate force"))
			, ,
			);
		// clang-format on
		FUNCTOR2D(VSphereGeom, VsPhys);
		DECLARE_LOGGER;
	};
	REGISTER_SERIALIZABLE(Law2_VSphereGeom_VsPhys_VsPM);

} // namespace yade
