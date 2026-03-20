#include "FrictPhys2.hpp"
#include <pkg/dem/FrictPhys.hpp>
//#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((FrictMat2)(Ip2_FrictMat2_FrictMat2_FrictPhys));

// The following code was moved from Ip2_FrictMat_FrictMat_FrictPhys.hpp

void Ip2_FrictMat2_FrictMat2_FrictPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;

	const shared_ptr<FrictMat2>& mat1 = YADE_PTR_CAST<FrictMat2>(b1);
	const shared_ptr<FrictMat2>& mat2 = YADE_PTR_CAST<FrictMat2>(b2);

	//Real Ra, Rb;                                                           //Vector3r normal;
	//assert(dynamic_cast<GenericSpheresContact*>(interaction->geom.get())); //only in debug mode
	//GenericSpheresContact* sphCont = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	//Ra                             = sphCont->refR1 > 0 ? sphCont->refR1 : sphCont->refR2;
	//Rb                             = sphCont->refR2 > 0 ? sphCont->refR2 : sphCont->refR1;

	//interaction->phys                           = shared_ptr<FrictPhys>(new FrictPhys());
	//const shared_ptr<FrictPhys>& contactPhysics = YADE_PTR_CAST<FrictPhys>(interaction->phys);
	//Real                         Ea             = mat1->young;
	//Real                         Eb             = mat2->young;
	//Real                         Va             = mat1->poisson;
	//Real                         Vb             = mat2->poisson;
	//Real                         kna            = Ea * Ra;
	//Real                         knb            = Eb * Rb;
	//Real                         ksa            = kna * Va;
	//Real                         ksb            = knb * Vb;

	////match maker or half the harmonic average of the two stiffnesses, when (2*Ri*Ei=2*kni) is the stiffness of a contact point on sphere "i"
	//Real Kn = (!kn) ? 2 * kna * knb / (kna + knb) : (*kn)(mat1->id, mat2->id, kna, knb);
	////same for shear stiffness
	//Real Ks = (!ks) ? 2 * ksa * ksb / (ksa + ksb) : (*ks)(mat1->id, mat2->id, ksa, ksb);

	//Real frictionAngle                     = (!frictAngle) ? math::min(mat1->frictionAngle, mat2->frictionAngle)
	//                                                       : (*frictAngle)(mat1->id, mat2->id, mat1->frictionAngle, mat2->frictionAngle);
	//contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	//contactPhysics->kn                     = Kn;
	//contactPhysics->ks                     = Ks;
	interaction->phys = shared_ptr<FrictPhys>(new FrictPhys());
	const shared_ptr<FrictPhys>& contactPhysics = YADE_PTR_CAST<FrictPhys>(interaction->phys);
	Real frictionAngle = std::min(mat1->frictionAngle, mat2->frictionAngle);
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	contactPhysics->kn = 2.0*mat1->Kn*mat2->Kn / (mat1->Kn + mat2->Kn);
	contactPhysics->ks = 2.0*mat1->Ks*mat2->Ks / (mat1->Ks + mat2->Ks);
};


} // namespace yade
