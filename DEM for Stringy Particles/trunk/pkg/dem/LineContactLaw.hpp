/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>     *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <core/Dispatching.hpp>
#include <core/GlobalEngine.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/Line.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <boost/tuple/tuple.hpp>

namespace yade { // Cannot have #include directive inside.


//**************************************************************************************12-26
class LineMat : public ElastMat {
public:
	virtual ~LineMat() {};
	void              postLoad(LineMat&);
	/// Serialization
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(LineMat,ElastMat,"Material description extending :yref:`FrictMat` with cohesive properties and rotational stiffnesses. For use e.g. with :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
		//((bool,isCohesive,true,,"Whether this body can form possibly cohesive interactions (if true and depending on other parameters such as :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys.setCohesionNow`)."))
		((Real,alphaKr,2.0,,"Dimensionless rolling stiffness."))
		((Real,alphaKtw,2.0,,"Dimensionless twist stiffness."))
		((Real, diameter, 0.0005, , "Diameter of the single line in [m] (the diameter is used to compute the cross-section area of the wire)."))
		((vector<Vector2r>, strainStressValues, , Attr::triggerPostLoad, "Piecewis"))
		((Real,shearCohesion,-1,,"Shear strength, homogeneous to a pressure. If negative the shear force is purely elastic."))
		((Real, k,0 , , "线材料拉伸刚度k"))//线材料拉伸刚度
		((Real, young_s,15e6 , , "拉伸弹性模量"))//
		//((Real, _k,0 , , "线材料拉伸刚度k"))//柔索拉伸刚度
		((Real, kn,0 , , "线材料拉伸刚度k"))//
		((Real, ks,0 , , "线材料qiexiang刚度k"))//
		((Real,frictionAngle,.5,,"Contact friction angle (in radians). Hint : use 'radians(degreesValue)' in python scripts."))
		((Real, as, 0., Attr::readonly, "Cross-section area of a single wire used to transform stress into force. [m²]"))
		((bool,momentRotationLaw,false,,"Use bending/twisting moment at contact. The contact may have moments only if both bodies have this flag true. See :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment.always_use_moment_law` for details."))
		,
		createIndex();
		);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(LineMat, ElastMat);
};

REGISTER_SERIALIZABLE(LineMat);

//!**************contactphys***************************//
class LinePhys : public NormShearPhys  {
public:
	virtual ~LinePhys() = default;

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(LinePhys, NormShearPhys, "An :yref:`interaction physics<IPhys>` that extends :yref:`RotStiffFrictPhys` adding a breakable cohesive nature. Used e.g. by :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
		((bool, momentRotationLaw, false, , "set from :yref:`CohFrictMat::momentRotationLaw` in order to possibly use bending/twisting moment at contacts (if true). See :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::always_use_moment_law` for details."))
		((Vector3r,_normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates), as sustained by particle #2 (from particle #1)."))
		((Real,Fn,0,,"cohesive part of the shear strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`)"))
		((Real,shearAdhesion,0,,"cohesive part of the shear strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`)"))
		((Real, k,0, , "材料刚度"))
		//((Real, _k,0 , , "线材料拉伸刚度k"))//柔索拉伸刚度
		((Real, initD, 0., , "初始距离."))
		((Real, normalLimit, 0., , "最大法向力."))
		((vector<Vector2r>, displForceValues, , Attr::readonly, "Defines the values for force-displacement curve."))
		((vector<Real>, stiffnessValues, , Attr::readonly, "Defines the values for the various stiffnesses (the elastic stiffness is stored as kn)."))
		//((Vector3r, WireForce, Vector3r::Zero(), , "Wire force ."))//***2023-12-2
		,
		createIndex();
	);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(LinePhys, NormShearPhys);
};

REGISTER_SERIALIZABLE(LinePhys);


class LineFrictPhys : public NormShearPhys  {
public:
	virtual ~LineFrictPhys() = default;

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(LineFrictPhys, NormShearPhys, "An :yref:`interaction physics<IPhys>` that extends :yref:`RotStiffFrictPhys` adding a breakable cohesive nature. Used e.g. by :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
		((Real, kn,0 , , "线材料拉伸刚度k"))//
		//((Vector3r, WireForce, Vector3r::Zero(), , "Wire force ."))//***2023-12-2
		,
		createIndex();
	);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(LineFrictPhys, NormShearPhys);
};

REGISTER_SERIALIZABLE(LineFrictPhys);

//!****************Ip2*******************//
class Ip2_LineMat_LineMat_LinePhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	//int  cohesionDefinitionIteration;//?

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_LineMat_LineMat_LinePhys, IPhysFunctor,
		"Generates cohesive-frictional interactions with moments, used in the contact law :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`. The normal/shear stiffness and friction definitions are the same as in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`, check the documentation there for details.\n\nAdhesions related to the normal and the shear components are calculated from :yref:`CohFrictMat::normalCohesion` ($C_n$) and :yref:`CohFrictMat::shearCohesion` ($C_s$). For particles of size $R_1$,$R_2$ the adhesion will be $a_i=C_i min(R_1,R_2)^2$, $i=n,s$.\n\nTwist and rolling stiffnesses are proportional to the shear stiffness through dimensionless factors alphaKtw and alphaKr, such that the rotational stiffnesses are defined by $k_s \\alpha_i R_1 R_2$, $i=tw\\,r$",
		((Real, normalLimit, 0., , "设定最大法向力."))
		,
		//cohesionDefinitionIteration = -1;//?
	);
	// clang-format on
	FUNCTOR2D(LineMat, LineMat);
};

REGISTER_SERIALIZABLE(Ip2_LineMat_LineMat_LinePhys);



class Ip2_LineMat_FrictMat_LineFrictPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	//int  cohesionDefinitionIteration;//?

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_LineMat_FrictMat_LineFrictPhys, IPhysFunctor,
		"Generates cohesive-frictional interactions with moments,",
		,
		//cohesionDefinitionIteration = -1;//?
	);
	// clang-format on
	FUNCTOR2D(LineMat, FrictMat);
};

REGISTER_SERIALIZABLE(Ip2_LineMat_FrictMat_LineFrictPhys);


//!************law2**************//LineNodeGeom
class Law2_LineNodeGeom_LinePhys_LineMoment : public LawFunctor {
public:
	bool                    go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_LineNodeGeom_LinePhys_LineMoment, LawFunctor, "Law for linear.",
	
	((bool, testStretch, false, , "use k to accumulate force"))
		, ,
	);
	// clang-format on
	FUNCTOR2D(LineNodeGeom, LinePhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_LineNodeGeom_LinePhys_LineMoment);

class Law2_LineNodeGeom2_LineFrictPhys_CundallStrack : public LawFunctor {
public:
	bool                    go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_LineNodeGeom2_LineFrictPhys_CundallStrack, LawFunctor, "Law for linear.",
	((bool, useK, false, , "use k to accumulate force"))
		, ,
	);
	// clang-format on
	FUNCTOR2D(LineNodeGeom2, LineFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_LineNodeGeom2_LineFrictPhys_CundallStrack);


class Law2_NodeScGeom2_LineFrictPhys_CundallStrack : public LawFunctor {
public:
	bool                    go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_NodeScGeom2_LineFrictPhys_CundallStrack, LawFunctor, "Law for linear.",
	((bool, useK, false, , "use k to accumulate force"))
		, ,
	);
	// clang-format on
	FUNCTOR2D(NodeScGeom2, LineFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_NodeScGeom2_LineFrictPhys_CundallStrack);

//***************************************************
} // namespace yadeshearAdhesion
