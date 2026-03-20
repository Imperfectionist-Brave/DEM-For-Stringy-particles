/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  Copyright (C) 2008 by Sergei Dorofeenko                               *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  © 2008 Václav Šmilauer <eudoxos@arcig.cz>                             *
*                                                                        *
*  © 2008 Bruno Chareyre <bruno.chareyre@grenoble-inp.fr                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/*
#ifdef YADE_OPENGL

#include "Gl1_MembraneNode.hpp"
#include <lib/opengl/OpenGLWrapper.hpp>
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

//YADE_PLUGIN((GlBoundFunctor)(GlShapeFunctor)(GlIGeomFunctor)(GlIPhysFunctor)(GlStateFunctor)(GlBoundDispatcher)(GlShapeDispatcher)(GlIGeomDispatcher)(
//        GlIPhysDispatcher)(GlStateDispatcher));

// Spheres==========================================

bool             Gl1_MembraneNode::wire;
bool             Gl1_MembraneNode::stripes;
int              Gl1_MembraneNode::glutSlices;
int              Gl1_MembraneNode::glutStacks;
Real             Gl1_MembraneNode::quality;
bool             Gl1_MembraneNode::localSpecView;
bool             Gl1_MembraneNode::circleView;
Real             Gl1_MembraneNode::circleRelThickness;
vector<Vector3r> Gl1_MembraneNode::vertices, Gl1_MembraneNode::faces;
int              Gl1_MembraneNode::glStripedSphereList = -1;
int              Gl1_MembraneNode::glGlutSphereList    = -1;
Real             Gl1_MembraneNode::prevQuality         = 0;
string           Gl1_MembraneNode::prevDisplayMode     = "";
char             Gl1_MembraneNode::circleAllowedRotationAxis;
char             Gl1_MembraneNode::prevCircleAllowedRotationAxis = 'z';

void Gl1_MembraneNode::go(const shared_ptr<Shape>& cm, const shared_ptr<State>&, bool wire2, const GLViewInfo&)
{
	glClearDepth(1.0f);
	glEnable(GL_NORMALIZE);

	Real r = (static_cast<MembraneNode*>(cm.get()))->radius;
	glColor3v(cm->color);
	if (circleView) {
		bool somethingChanged
		        = (math::abs(quality - prevQuality) > 0.001 || prevDisplayMode != "torus"
		           || prevCircleAllowedRotationAxis != circleAllowedRotationAxis);
		if (somethingChanged) {
			prevCircleAllowedRotationAxis = circleAllowedRotationAxis;
			prevDisplayMode               = "torus";
			glDeleteLists(glGlutSphereList, 1);
			glGlutSphereList = glGenLists(1);
			glNewList(glGlutSphereList, GL_COMPILE)
				;
				glEnable(GL_LIGHTING);
				glShadeModel(GL_SMOOTH);
				switch (tolower(circleAllowedRotationAxis)) { //rotate the torus according to the axis from which we want to look at it.
					case 'z': break;                      //Initial torus axis is z, nothing to do
					case 'x': glRotatef(90, 0, 1, 0); break;
					case 'y': glRotatef(90, 1, 0, 0); break;
					default: cerr << "Error in Gl1_MembraneNode::go, circleAllowedRotationAxis should be \"x\", \"y\" or \"z\"." << endl;
				}
				glutSolidTorus(
				        0.5 * circleRelThickness * r,
				        r * (1.0 - circleRelThickness / 2.),
				        int(math::round(quality * glutStacks)),
				        int(math::round(quality * glutSlices))); //generate torus
			glEndList();
		}
		glCallList(glGlutSphereList);
	} else {
		if (wire || wire2) {
			glutWireSphere(r, int(math::round(quality * glutSlices)), int(math::round(quality * glutStacks)));
		} else {
			//Check if quality has been modified or if previous lists are invalidated (e.g. by creating a new qt view), then regenerate lists
			bool somethingChanged
			        = (math::abs(quality - prevQuality) > 0.001 || glIsList(glStripedSphereList) != GL_TRUE || prevDisplayMode != "sphere");
			if (somethingChanged) {
				initStripedGlList();
				initGlutGlList();
				prevQuality     = quality;
				prevDisplayMode = "sphere";
			}
			glScale(r, r, r);
			if (stripes) {
				glCallList(glStripedSphereList);
			} else {
				glCallList(glGlutSphereList);
			}
		}
	}
	return;
}
YADE_PLUGIN((Gl1_MembraneNode));

void Gl1_MembraneNode::subdivideTriangle(Vector3r& v1, Vector3r& v2, Vector3r& v3, int depth)
{
	Vector3r v;
	//Change color only at the appropriate level, i.e. 8 times in total, since we draw 8 mono-color sectors one after another
	if (depth == int(quality) || quality <= 0) {
		v = (v1 + v2 + v3) / 3.0;
		GLfloat matEmit[4];
		if (v[1] * v[0] * v[2] > 0) {
			matEmit[0] = 0.3f;
			matEmit[1] = 0.3f;
			matEmit[2] = 0.3f;
			matEmit[3] = 1.f;
		} else {
			matEmit[0] = 0.15f;
			matEmit[1] = 0.15f;
			matEmit[2] = 0.15f;
			matEmit[3] = 0.2f;
		}
		glMaterialfv(GL_FRONT, GL_EMISSION, matEmit);
	}
	if (depth == 1) { //Then display 4 triangles
		Vector3r v12 = v1 + v2;
		Vector3r v23 = v2 + v3;
		Vector3r v31 = v3 + v1;
		v12.normalize();
		v23.normalize();
		v31.normalize();
		//Use TRIANGLE_STRIP for faster display of adjacent facets
		glBegin(GL_TRIANGLE_STRIP)
			;
			glNormal3v(v1);
			glVertex3v(v1);
			glNormal3v(v31);
			glVertex3v(v31);
			glNormal3v(v12);
			glVertex3v(v12);
			glNormal3v(v23);
			glVertex3v(v23);
			glNormal3v(v2);
			glVertex3v(v2);
		glEnd();
		//terminate with this triangle left behind
		glBegin(GL_TRIANGLES)
			;
			glNormal3v(v3);
			glVertex3v(v3);
			glNormal3v(v23);
			glVertex3v(v23);
			glNormal3v(v31);
			glVertex3v(v31);
		glEnd();
		return;
	}
	Vector3r v12 = v1 + v2;
	Vector3r v23 = v2 + v3;
	Vector3r v31 = v3 + v1;
	v12.normalize();
	v23.normalize();
	v31.normalize();
	subdivideTriangle(v1, v12, v31, depth - 1);
	subdivideTriangle(v2, v23, v12, depth - 1);
	subdivideTriangle(v3, v31, v23, depth - 1);
	subdivideTriangle(v12, v23, v31, depth - 1);
}

void Gl1_MembraneNode::initStripedGlList()
{
	if (!vertices.size()) { //Fill vectors with vertices and facets
		//Define 6 points for +/- coordinates
		vertices.push_back(Vector3r(-1, 0, 0)); //0
		vertices.push_back(Vector3r(1, 0, 0));  //1
		vertices.push_back(Vector3r(0, -1, 0)); //2
		vertices.push_back(Vector3r(0, 1, 0));  //3
		vertices.push_back(Vector3r(0, 0, -1)); //4
		vertices.push_back(Vector3r(0, 0, 1));  //5
		//Define 8 sectors of the sphere
		faces.push_back(Vector3r(3, 4, 1));
		faces.push_back(Vector3r(3, 0, 4));
		faces.push_back(Vector3r(3, 5, 0));
		faces.push_back(Vector3r(3, 1, 5));
		faces.push_back(Vector3r(2, 1, 4));
		faces.push_back(Vector3r(2, 4, 0));
		faces.push_back(Vector3r(2, 0, 5));
		faces.push_back(Vector3r(2, 5, 1));
	}
	//Generate the list. Only once for each qtView, or more if quality is modified.
	glDeleteLists(glStripedSphereList, 1);
	glStripedSphereList = glGenLists(1);
	glNewList(glStripedSphereList, GL_COMPILE)
		;
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		// render the sphere now
		for (int i = 0; i < 8; i++)
			subdivideTriangle(
			        vertices[(unsigned int)faces[i][0]],
			        vertices[(unsigned int)faces[i][1]],
			        vertices[(unsigned int)faces[i][2]],
			        1 + (int)quality);
	glEndList();
}

void Gl1_MembraneNode::initGlutGlList()
{
	//Generate the "no-stripes" display list, each time quality is modified
	glDeleteLists(glGlutSphereList, 1);
	glGlutSphereList = glGenLists(1);
	glNewList(glGlutSphereList, GL_COMPILE)
		;
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		glutSolidSphere(1.0, int(math::round(math::max(quality * glutSlices, (Real)2.))), int(math::round(math::max(quality * glutStacks, (Real)3.))));
	glEndList();
}

} // namespace yade

#endif /* YADE_OPENGL */
