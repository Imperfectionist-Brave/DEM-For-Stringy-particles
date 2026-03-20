#pragma once
#include <core/Dispatching.hpp>
//#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/MembraneNode.hpp>
#include <core/Aabb.hpp>
#include <pkg/common/GLDrawFunctors.hpp>
//#include <core/Dispatching.hpp>

namespace yade { // Cannot have #include directive inside.

class BoxM : public Shape {
public:
	BoxM(const Vector3r& _extents)
		: extents(_extents)
	{
	}
	virtual ~BoxM() {};
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(BoxM, Shape, "Box (cuboid) particle geometry. (Avoid using in new code, prefer :yref:`Facet` instead.)",
		((Vector3r, extents, , , "Half-size of the cuboid")),
		/* ctor */ createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(BoxM, Shape);
};
REGISTER_SERIALIZABLE(BoxM);


//class BoxM;
class Bo1_BoxM_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*) override;
	FUNCTOR1D(BoxM);
	// clang-format off
	YADE_CLASS_BASE_DOC(Bo1_BoxM_Aabb, BoundFunctor, "Create/update an :yref:`Aabb` of a :yref:`Box`.");
	// clang-format on
};

REGISTER_SERIALIZABLE(Bo1_BoxM_Aabb);

class Gl1_BoxM : public GlShapeFunctor {
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	RENDERS(BoxM);
	// clang-format off
	YADE_CLASS_BASE_DOC(Gl1_BoxM, GlShapeFunctor, "Renders :yref:`Box` object");
	// clang-format on
};
REGISTER_SERIALIZABLE(Gl1_BoxM);

class Ig2_BoxM_Sphere_ScGeom : public IGeomFunctor {
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
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_BoxM_Sphere_ScGeom,IGeomFunctor,"Create an interaction geometry :yref:`ScGeom` from :yref:`Box` and :yref:`Sphere`, representing the box with a projected virtual sphere of same radius.",
    ((Real,interactionDetectionFactor,1,,"Enlarge sphere radii by this factor (if >1), to permit creation of distant interactions.\n\nInteractionGeometry will be computed when interactionDetectionFactor*(rad) > distance.\n\n.. note::\n\t This parameter is functionally coupled with :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor`, which will create larger bounding boxes and should be of the same value.")))
	// clang-format on
	FUNCTOR2D(BoxM, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(BoxM, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_BoxM_Sphere_ScGeom);

class Ig2_BoxM_MembraneNode_ScGeom : public IGeomFunctor {
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
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_BoxM_MembraneNode_ScGeom,IGeomFunctor,"Create an interaction geometry :yref:`ScGeom` from :yref:`Box` and :yref:`Sphere`, representing the box with a projected virtual sphere of same radius.",
    ((Real,interactionDetectionFactor,1,,"Enlarge sphere radii by this factor (if >1), to permit creation of distant interactions.\n\nInteractionGeometry will be computed when interactionDetectionFactor*(rad) > distance.\n\n.. note::\n\t This parameter is functionally coupled with :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor`, which will create larger bounding boxes and should be of the same value.")))
	// clang-format on
	FUNCTOR2D(BoxM, MembraneNode);
	DEFINE_FUNCTOR_ORDER_2D(BoxM, MembraneNode);
};
REGISTER_SERIALIZABLE(Ig2_BoxM_MembraneNode_ScGeom);

} // namespace yade
