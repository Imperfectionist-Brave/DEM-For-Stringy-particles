/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@grenoble-inp.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_OPENGL

#include "LineSegment.hpp"
#include <lib/high-precision/Constants.hpp>
#include <lib/opengl/OpenGLWrapper.hpp>


namespace yade { // Cannot have #include directive inside.

//!##################	Rendering   #####################

bool Gl1_Segment::wire;
bool Gl1_Segment::glutNormalize;
int  Gl1_Segment::glutSlices;
int  Gl1_Segment::glutStacks;

void Gl1_Segment::out(Quaternionr q)
{
	AngleAxisr aa(q);
	std::cout << " axis: " << aa.axis()[0] << " " << aa.axis()[1] << " " << aa.axis()[2] << ", angle: " << aa.angle() << " | ";
}

void Gl1_Segment::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& /*st*/, bool wire2, const GLViewInfo&)
{
	Segment*               Se     = static_cast<Segment*>(cm.get());
	Real                          r      = Se->radius;
	
	//const shared_ptr<Interaction> intr   = scene->interactions->find((int)LC->node1->getId(), (int)LC->node2->getId());//需要判断是否断裂
	//long iter = scene->iter;
	//std::cout<<"node2 pos"<<Se->node2->pos<<"  "<<"node1 pos"<<Se->node1->pos<<"length"<<length<<std::endl;
	//std::cout<<"node2 pos"<<Se->node2->pos<<"  "<<"node2 pos"<<Se->node1->pos<<std::endl;
	Vector3r                      segt   = Se->node2->pos - Se->node1->pos;
	//Real                          length = Se->getLength();
	Real                          length = segt.norm();
	//if (scene->isPeriodic && intr) segt += scene->cell->intrShiftPos(intr->cellDist);
	//glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(cm->color[0],cm->color[1],cm->color[2]));

	glColor3v(cm->color);
	if (glutNormalize) glPushAttrib(GL_NORMALIZE);
	// 	glPushMatrix();
	Quaternionr shift;
	shift.setFromTwoVectors(Vector3r::UnitZ(), segt);

	//st->ori = Quaternionr::Identity(); // Otherwise clumped connexions get rotated by the clump motion and the view is messed up (note that orientation is never used in mechanical calculations in the case of connexions and pfacets).
	
	if (!Se->broken) { drawCylinder(wire || wire2, r, length, shift); }//********5-29
	// 	if (intr && scene->isPeriodic) { glTranslatef(-segt[0],-segt[1],-segt[2]); drawCylinder(wire || wire2, r,length,-shift);}
	if (glutNormalize) glPopAttrib();
	// 	glPopMatrix();
	return;
}

void Gl1_Segment::drawCylinder(bool wireNonMember, Real radius, Real length, const Quaternionr& shift)
{
	glPushMatrix();
	GLUquadricObj* quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, (GLenum)(wireNonMember ? GLU_SILHOUETTE : GLU_FILL));
	gluQuadricNormals(quadObj, (GLenum)GLU_SMOOTH);
	gluQuadricOrientation(quadObj, (GLenum)GLU_OUTSIDE);
	AngleAxisr aa(shift);
	glRotate(aa.angle() * 180.0 / Mathr::PI, aa.axis()[0], aa.axis()[1], aa.axis()[2]);
	gluCylinder(quadObj, radius, radius, length, glutSlices, glutStacks);
	gluQuadricOrientation(quadObj, (GLenum)GLU_INSIDE);
	//glutSolidSphere(radius,glutSlices,glutStacks);
	glTranslate(0.0, 0.0, length);

	//glutSolidSphere(radius,glutSlices,glutStacks);
	//    gluDisk(quadObj,0.0,radius,glutSlices,_loops);
	gluDeleteQuadric(quadObj);
	glPopMatrix();
}
YADE_PLUGIN((Gl1_Segment));

} // namespace yade
//******************************
#endif
