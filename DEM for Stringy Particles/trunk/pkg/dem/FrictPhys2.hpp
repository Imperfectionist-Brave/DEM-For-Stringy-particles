/*************************************************************************
*  Copyright (C) 2007 by Bruno CHAREYRE                                  *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include <core/Dispatching.hpp>
//#include <pkg/common/ElastMat.hpp>
#include <core/Material.hpp>
//#include <pkg/common/MatchMaker.hpp>
//#include <pkg/common/NormShearPhys.hpp>

namespace yade { // Cannot have #include directive inside.

class FrictMat2 : public Material {
public:
	virtual ~FrictMat2() {};
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(FrictMat2, Material, "Purely elastic material. The material parameters may have different meanings depending on the :yref:`IPhysFunctor` used : true Young and Poisson in :yref:`Ip2_FrictMat_FrictMat_MindlinPhys`, or contact stiffnesses in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`.",
		((Real, young, 1e9, , "elastic modulus [Pa]. It has different meanings depending on the Ip functor."))
		((Real, poisson, .25, , "Poisson's ratio or the ratio between shear and normal stiffness [-]. It has different meanings depending on the Ip functor.  "))
		((Real, Kn, 1e8, , "Normal stiffness (N/m)."))
		((Real, Ks, 1e5, , "Shear stiffness (N/m)."))
		((Real, frictionAngle, .5, , "Contact friction angle (in radians).")),
		/*ctor*/ createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(FrictMat2, Material);
};
REGISTER_SERIALIZABLE(FrictMat2);

class Ip2_FrictMat2_FrictMat2_FrictPhys : public IPhysFunctor {
public:
	 void go(const shared_ptr<Material>& b1,
		const shared_ptr<Material>& b2,
		const shared_ptr<Interaction>& interaction)override;
	FUNCTOR2D(FrictMat2, FrictMat2);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_FrictMat2_FrictMat2_FrictPhys, IPhysFunctor, "",
		//((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))//will crash:Abort (core dumped)
		//((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
		,
		);
};
REGISTER_SERIALIZABLE(Ip2_FrictMat2_FrictMat2_FrictPhys);


} // namespace yade
