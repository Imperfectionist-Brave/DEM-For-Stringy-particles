#pragma once
#include <core/Scene.hpp>
#include <pkg/common/BoundaryController.hpp>
//#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.

class UniaxialForcer : public BoundaryController {
private:
	bool  needsInit;
	inline void addForce(Scene* rb, Body::id_t id, const Vector3r& f) const
	{
		return rb->forces.addForce(id, f); /* needs sync, which is done at the beginning of action */
	}
	inline Real  axisCoord(Body::id_t id) { return Body::byId(id, scene)->state->pos[axis]; };
	void  init();

public:
	//bool isActivated() override { return active; }
	Real sumPosForces, sumNegForces;
	//Real initAccelTime_s /* value always in s, computed from initAccelTime */;
	/** negCoords*/
	Real posCoords, negCoords;
	unsigned int solve_num;
	void action() override;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(UniaxialForcer,BoundaryController,"Axial displacing two groups of bodies in the opposite direction with given strain rate.",
			((Real, loadForce, 10.0, , "Force load on the segment"))
			((Vector3r, unitForce, Vector3r(0,0,0), , "Force load on the segment"))
			((Real, loadTime_s, 1e-3, , "time of the force arrived"))
			((long, loadStepNum, 0, , "1 for straight,2 for 1/4straight"))
			((unsigned int, force_num_max, 2000, , "interval of restart counting solve_num"))
			//((long, initialIter, NaN, , "1 for straight,2 for 1/4straight"))
			((long, currentIter, 0, , "1 for straight,2 for 1/4straight"))
			//((Real, loadRate,NaN,,"Rate of strain, starting at 0, linearly raising to strainRate. [-]"))
			//((Real,stopStrain,NaN,,"Strain at which we will pause simulation; inactive (nan) by default; must be reached from below (in absolute value)"))
			//((bool,active,true,,"Whether this engine is activated"))
			((Real,currentStrain,NaN,,"Current strain rate (update automatically). |yupdate|"))
			((int,axis,0,,"The axis which is strained (0,1,2 for x,y,z)"))
			((Body::id_t,posId,-1,,"Bodies on which strain will be applied (on the positive end along the axis)"))
			((Body::id_t,negId,-1,,"Bodies on which strain will be applied (on the negative end along the axis)"))
			((Real,originalLength,NaN,,"Distance of reference bodies in the direction of axis before straining started (computed automatically) [m]"))
			((Real,currentLength,NaN,,"Distance of reference bodies in the direction of axis before straining started (computed automatically) [m]"))
			((bool,blockDisplacements,true,,"Whether displacement of boundary bodies perpendicular to the strained axis are blocked or are free"))
			((int, stretchModel, 1, , "1 for straight,2 for 1/4straight")),
			/*ctor*/ needsInit=true;solve_num=0;
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(UniaxialForcer);

} // namespace yade
