#define GLEW_STATIC //Using the static lib, so we need to enable it
#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include <Ant/AntTweakBar.h>
#include <memory>
#include <algorithm>
#include "helper.h"
#include "Shader.h"
#include "Skybox.h"
#include "times.h"
#include "InstancedMesh.h"
#include "FBO.h"
#include "Texture.h"
#include "glm/gtx/string_cast.hpp"
#include "marchingCubesVolume.h"

#include "HalfEdgeMesh.h"
#include "MeshResampler.h"

//Time
Timer timer;
int frame;
long timeCounter, timebase;
char timeString[50];

//Resolution (has to be changed in helper.h too)
glm::vec2 resolution = glm::vec2(1024, 768);

//Externals
cameraSystem cam(1.0f, 1.0f, glm::vec3(20.95f, 20.95f, -0.6f));
glm::mat4 projMatrix;
glm::mat4 viewMatrix;

bool leftMouseClick;
int leftMouseClickX;
int leftMouseClickY;

//Skybox
Skybox skybox;
char* negz = "C:/Dev/Assets/SkyboxTextures/Yokohama2/negz.jpg";
char* posz = "C:/Dev/Assets/SkyboxTextures/Yokohama2/posz.jpg";
char* posy = "C:/Dev/Assets/SkyboxTextures/Yokohama2/posy.jpg";
char* negy = "C:/Dev/Assets/SkyboxTextures/Yokohama2/negy.jpg";
char* negx = "C:/Dev/Assets/SkyboxTextures/Yokohama2/negx.jpg";
char* posx = "C:/Dev/Assets/SkyboxTextures/Yokohama2/posx.jpg";

//Shaders
Shader basicShader;
Shader modelLoaderShader;

//Skybox
Shader skyboxShader;

//Models
simpleModel *teaPot = 0;
HalfedgeMesh heMesh;

// tweak bar
TwBar *tweakBar;
bool wireFrameTeapot = true;
/* *********************************************************************************************************
TweakBar
********************************************************************************************************* */
void setupTweakBar() {
	TwInit(TW_OPENGL_CORE, NULL);
	tweakBar = TwNewBar("Settings");
	TwAddSeparator(tweakBar, "Wireframe", nullptr);
	TwAddVarRW(tweakBar, "Wireframe Teapot", TW_TYPE_BOOLCPP, &wireFrameTeapot, " label='Wireframe Teapot' ");
}

/* *********************************************************************************************************
Initiation
********************************************************************************************************* */

GLuint vboHalfEdgeMesh[2];
vector<glm::vec3> halfEdgeMeshVertices;
vector<glm::vec3> halfEdgeMeshColors;

void init() {
	
	teaPot = new simpleModel("C:/Dev/Assets/Teapot/teapot.obj", false);
	teaPot->upload();

	/*****************************************************************
	HalfEdgeMesh
	*****************************************************************/
	//Create
	vector< vector<Index> > polygons;
	for (int i = 0; i < teaPot->indices.size() / 3; i++) {
		vector<unsigned __int64> polygon;
		polygon.push_back(teaPot->indices[3 * i + 0]);
		polygon.push_back(teaPot->indices[3 * i + 1]);
		polygon.push_back(teaPot->indices[3 * i + 2]);

		polygons.push_back(polygon);
	}
	
	heMesh.build(polygons, teaPot->vertices);

	int currentCounter = 0;
	int faceNr = 875;

	FaceIter f_flip = heMesh.facesBegin();
	FaceIter f_flip_neighbour;
	for(int i = 0; i < faceNr; i++){
		f_flip++;
	}

	//heMesh.flipEdge(f_flip->halfedge()->edge());

	heMesh.splitEdge(f_flip->halfedge()->edge());

	//heMesh.collapseEdge(f_flip->halfedge()->edge());

	//f_flip_neighbour = f_flip->halfedge()->twin()->face();

	//MeshResampler resampler;
	//resampler.upsample(heMesh);

	for (FaceIter f = heMesh.facesBegin(); f != heMesh.facesEnd(); f++) {
		currentCounter++;
		//1. Half Edge of current face
		HalfedgeIter h = f->halfedge();
		//std::cout << "Got Halfedge " << std::endl;
		halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
		//std::cout << "Got Vertex " << std::endl;
		halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

		//2. Half Edge of current face
		h = h->next();
		//std::cout << "Got Halfedge " << std::endl;
		halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
		//std::cout << "Got Vertex " << std::endl;
		halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

		//3. Half Edge of current face
		h = h->next();		
		//std::cout << "Got Halfedge " << std::endl;
		halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
		//std::cout << "Got Vertex " << std::endl;
		halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	}

	glGenBuffers(2, vboHalfEdgeMesh);
	glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[0]);
	glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshVertices.size() * sizeof(float) * 3, halfEdgeMeshVertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[1]);
	glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshColors.size() * sizeof(float) * 3, halfEdgeMeshColors.data(), GL_STATIC_DRAW);

	/*****************************************************************
	Skybox (Only for aesthetic reasons, can be deleted)
	*****************************************************************/
	skybox.createSkybox(negz, posz, posy, negy, negx, posx);
}


void loadShader(bool init) {
	basicShader = Shader("./shader/standard.vs.glsl", "./shader/standard.fs.glsl");
	modelLoaderShader = Shader("./shader/modelLoader.vs.glsl", "./shader/modelLoader.fs.glsl");
	skyboxShader = Shader("./shader/skybox.vs.glsl", "./shader/skybox.fs.glsl");
}

/* *********************************************************************************************************
Scenes: Unit cube + Pointcloud, Results of marching cubes
********************************************************************************************************* */
void sponzaStandardScene(){
	skyboxShader.enable();
	skyboxShader.uniform("projMatrix", projMatrix);
	skyboxShader.uniform("viewMatrix", cam.cameraRotation);
	skybox.Draw(skyboxShader);
	skyboxShader.disable();

	glm::mat4 modelMatrix;

	/* ********************************************
	Half-Edge-Mesh
	**********************************************/

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	if (wireFrameTeapot) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glPolygonMode(GL_BACK, GL_LINE);
	}
	basicShader.enable();
	basicShader.uniform("projMatrix", projMatrix);
	basicShader.uniform("viewMatrix", viewMatrix);
	modelMatrix = glm::mat4(1.0f);
	modelMatrix = glm::translate(modelMatrix, glm::vec3(-10.0f, 0.0f, 0.0f));
	modelMatrix = glm::scale(modelMatrix, glm::vec3(3.0f));
	basicShader.uniform("modelMatrix", modelMatrix);
	basicShader.uniform("col", glm::vec3(1.0f, 1.0f, 0.0f));
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);

	glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[0]);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[1]);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glDrawArrays(GL_TRIANGLES, 0, halfEdgeMeshVertices.size());

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	basicShader.disable();
	if (wireFrameTeapot) {
		glPolygonMode(GL_FRONT, GL_FILL);
		glPolygonMode(GL_BACK, GL_FILL);
	}
}	
/* *********************************************************************************************************
Vector/Triangle - Intersection
********************************************************************************************************* */
float rayIntersectsTriangle(glm::vec3 origin, glm::vec3 direction, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, float& u_out, float& v_out) {
	//Intersection - Test:
	//http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
	glm::vec3 h, s, q;
	float a, f, u, v;

	glm::vec3 e1 = v1 - v0;
	glm::vec3 e2 = v2 - v0;

	h = glm::cross(direction, e2);
	a = glm::dot(e1, h);

	if (a > -0.00001f && a < 0.00001f)
		return(-1.0f);

	f = 1.0f / a;
	s = origin - v0;
	u = f * (glm::dot(s, h));

	if (u < 0.0f || u > 1.0f)
		return(-1.0f);

	q = glm::cross( s, e1);
	v = f * glm::dot(direction, q);

	if (v < 0.0f || u + v > 1.0f)
		return(-1.0f);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	float t = f * glm::dot(e2, q);

	if (t > 0.00001f) { // ray intersection
		u_out = u;
		v_out = v;
		return(t);
	}

	else // this means that there is a line intersection
		 // but not a ray intersection
		return (-1.0f);
}

/*********************************************
Mouse selection on simpleModel types (vector of vertices/indices)
*********************************************/
void mouseTriangleSelecction() {
	if (leftMouseClick) {

		leftMouseClick = false;

		float x = 2.0f * (float(leftMouseClickX) / float(WIDTH)) - 1.0f;
		float y = 1.0f - 2.0f * (float(leftMouseClickY) / float(HEIGHT));

		glm::vec4 ray_clip = glm::vec4(x, y, -1.0, 1.0);

		glm::vec4 ray_eye = glm::inverse(projMatrix) * ray_clip;
		ray_eye = glm::vec4(ray_eye.x, ray_eye.y, -1.0, 0.0);

		glm::vec3 ray_wor = glm::vec3((glm::inverse(viewMatrix) * ray_eye));
		ray_wor = glm::normalize(ray_wor);

		vector<pair<float, pair<int, FaceIter>>> list;
		int i = 0;
		for (FaceIter f = heMesh.facesBegin(); f != heMesh.facesEnd(); f++) {

			glm::mat4 modelMatrixIntersect = glm::mat4(1.0f);
			modelMatrixIntersect = glm::translate(modelMatrixIntersect, glm::vec3(-10.0f, 0.0f, 0.0f));
			modelMatrixIntersect = glm::scale(modelMatrixIntersect, glm::vec3(3.0f));

			HalfedgeIter h = f->halfedge();
			glm::vec3 v0 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

			h = h->next();
			glm::vec3 v1 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

			h = h->next();
			glm::vec3 v2 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

			float u, v;
			float t = rayIntersectsTriangle(cam.position, ray_wor, v0, v1, v2, u, v);


			if (t != -1.0f) {
				std::cout << "t: " << t << " u: " << u << " v: " << v << std::endl;
				list.push_back(pair<float, pair<int, FaceIter>>(t, pair<int, FaceIter>(i, f)));
			}

			i++;
		}



		//If we dont even hit a triangle stop here
		if (!list.empty()) {

			//Sort all triangles hit by ray
			sort(list.begin(), list.end());

			//Find the closest edge to the selected point
			//We need the barycentric coordinates of the intersection ray/triangle ( we could get this smarter then recalculating)
			FaceIter f_sel = list[0].second.second;
			glm::mat4 modelMatrixIntersect = glm::mat4(1.0f);
			modelMatrixIntersect = glm::translate(modelMatrixIntersect, glm::vec3(-10.0f, 0.0f, 0.0f));
			modelMatrixIntersect = glm::scale(modelMatrixIntersect, glm::vec3(3.0f));
			HalfedgeIter h = f_sel->halfedge();
			glm::vec3 v0 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));
			h = h->next();
			glm::vec3 v1 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));
			h = h->next();
			glm::vec3 v2 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

			float u = -1.0f, v = -1.0f;
			float t = rayIntersectsTriangle(cam.position, ray_wor, v0, v1, v2, u, v);
			float w = 1.0f - u - v;

			halfEdgeMeshColors.erase(halfEdgeMeshColors.begin(), halfEdgeMeshColors.end());
			halfEdgeMeshVertices.erase(halfEdgeMeshVertices.begin(), halfEdgeMeshVertices.end());

			for (FaceIter f = heMesh.facesBegin(); f != heMesh.facesEnd(); f++) {
				//1. Half Edge of current face
				HalfedgeIter h = f->halfedge();
				halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
				
				//2. Half Edge of current face
				h = h->next();
				halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));

				//3. Half Edge of current face
				h = h->next();
				halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));

				if (f_sel == f) {
					halfEdgeMeshColors.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
					halfEdgeMeshColors.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
					halfEdgeMeshColors.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
				}
				else {
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
				}
			}

			glGenBuffers(2, vboHalfEdgeMesh);
			glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[0]);
			glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshVertices.size() * sizeof(float) * 3, halfEdgeMeshVertices.data(), GL_STATIC_DRAW);

			glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[1]);
			glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshColors.size() * sizeof(float) * 3, halfEdgeMeshColors.data(), GL_STATIC_DRAW);
		}

	}
}

/*********************************************
Mouse selection on heModel types (iterator linked list)
*********************************************/
void mouseTriangleEdgeSplit() {

	if (leftMouseClick) {

			leftMouseClick = false;

			float x = 2.0f * (float(leftMouseClickX) / float(WIDTH)) - 1.0f;
			float y = 1.0f - 2.0f * (float(leftMouseClickY) / float(HEIGHT));

			glm::vec4 ray_clip = glm::vec4(x, y, -1.0, 1.0);

			glm::vec4 ray_eye = glm::inverse(projMatrix) * ray_clip;
			ray_eye = glm::vec4(ray_eye.x, ray_eye.y, -1.0, 0.0);

			glm::vec3 ray_wor = glm::vec3((glm::inverse(viewMatrix) * ray_eye));
			ray_wor = glm::normalize(ray_wor);

			vector<pair<float, pair<int, FaceIter>>> list;
			int i = 0;
			for (FaceIter f = heMesh.facesBegin(); f != heMesh.facesEnd(); f++) {

				glm::mat4 modelMatrixIntersect = glm::mat4(1.0f);
				modelMatrixIntersect = glm::translate(modelMatrixIntersect, glm::vec3(-10.0f, 0.0f, 0.0f));
				modelMatrixIntersect = glm::scale(modelMatrixIntersect, glm::vec3(3.0f));

				HalfedgeIter h = f->halfedge();
				glm::vec3 v0 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

				h = h->next();
				glm::vec3 v1 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

				h = h->next();
				glm::vec3 v2 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

				float u, v;
				float t = rayIntersectsTriangle(cam.position, ray_wor, v0, v1, v2, u, v);


				if (t != -1.0f) {
					std::cout << "t: " << t << " u: " << u << " v: " << v << std::endl;
					list.push_back(pair<float, pair<int, FaceIter>>(t, pair<int, FaceIter>(i, f)));
				}

				i++;
			}

			

			//If we dont even hit a triangle stop here
			if (!list.empty()) {

				//Sort all triangles hit by ray
				sort(list.begin(), list.end());

				//Find the closest edge to the selected point
				//We need the barycentric coordinates of the intersection ray/triangle ( we could get this smarter then recalculating)
				FaceIter f = list[0].second.second;
				glm::mat4 modelMatrixIntersect = glm::mat4(1.0f);
				modelMatrixIntersect = glm::translate(modelMatrixIntersect, glm::vec3(-10.0f, 0.0f, 0.0f));
				modelMatrixIntersect = glm::scale(modelMatrixIntersect, glm::vec3(3.0f));

				HalfedgeIter h = f->halfedge();
				glm::vec3 v0 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

				h = h->next();
				glm::vec3 v1 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));

				h = h->next();
				glm::vec3 v2 = glm::vec3(modelMatrixIntersect * glm::vec4(h->vertex()->position, 1.0f));


				float u = -1.0f, v = -1.0f;
				float t = rayIntersectsTriangle(cam.position, ray_wor, v0, v1, v2, u, v);
				float w = 1.0f - u - v;
				std::cout << "Chosen result: t: " << t << " u: " << u << " v: " << v << std::endl;
				EdgeIter eToSplit;

				//the POBLEM here...
				//NOT: iterator memorizes no state, it resets itself!
				//WAS: "1-u-v" is a INTEGER 1 --> 1.0f
				if (u <= v && u <= w) {
					std::cout << "u" << std::endl;
					eToSplit = f->halfedge()->next()->next()->edge(); // v0 -> v1 -> v2
					//f->halfedge()->next(); //back to v0 is NOT needed! Iterator memorizes no state.
				}
				else if (v <= u && v <= w) {
					std::cout << "v" << std::endl;
					eToSplit = f->halfedge()->edge(); //v0
				}
				else if (w <= u && w <= v) {
					std::cout << "w" << std::endl;
					eToSplit = f->halfedge()->next()->edge(); // v0 -> v1 
					//f->halfedge()->next()->next(); //back to v2-> v0  is NOT needed! Iterator memorizes no state.

				}

				heMesh.splitEdge(eToSplit);

				halfEdgeMeshColors.erase(halfEdgeMeshColors.begin(), halfEdgeMeshColors.end());
				halfEdgeMeshVertices.erase(halfEdgeMeshVertices.begin(), halfEdgeMeshVertices.end());

				for (FaceIter f = heMesh.facesBegin(); f != heMesh.facesEnd(); f++) {
					//1. Half Edge of current face
					HalfedgeIter h = f->halfedge();
					halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

					//2. Half Edge of current face
					h = h->next();
					halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

					//3. Half Edge of current face
					h = h->next();
					halfEdgeMeshVertices.push_back(glm::vec3(h->vertex()->position.x, h->vertex()->position.y, h->vertex()->position.z));
					halfEdgeMeshColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
				}

				glGenBuffers(2, vboHalfEdgeMesh);
				glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[0]);
				glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshVertices.size() * sizeof(float) * 3, halfEdgeMeshVertices.data(), GL_STATIC_DRAW);

				glBindBuffer(GL_ARRAY_BUFFER, vboHalfEdgeMesh[1]);
				glBufferData(GL_ARRAY_BUFFER, halfEdgeMeshColors.size() * sizeof(float) * 3, halfEdgeMeshColors.data(), GL_STATIC_DRAW);
		}

	}

}

/* *********************************************************************************************************
Display + Main
********************************************************************************************************* */
void display() {
	//Timer
	timer.update();
	//FPS-Counter
	frame++;
	timeCounter = glutGet(GLUT_ELAPSED_TIME);
	if (timeCounter - timebase > 1000) {
		sprintf_s(timeString, "FPS:%4.2f", frame*1000.0 / (timeCounter - timebase));
		timebase = timeCounter;
		frame = 0;
		glutSetWindowTitle(timeString);
	}


	//mouseTriangleSelecction();
	mouseTriangleEdgeSplit();

	//OpenGL Clears
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glClearColor(0.2f, 0.2f, 0.2f, 1);
	
	sponzaStandardScene();
	TwDraw(); //Draw Tweak-Bar

	glutSwapBuffers();
	glutPostRedisplay();

}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_STENCIL);

	glutCreateWindow("Basic Framework");

	setupTweakBar();

	GLenum err = glewInit();
	if (GLEW_OK != err) {
		std::cerr << "Error : " << glewGetErrorString(err) << std::endl;
	}

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMotionFunc(onMouseMove);
	glutMouseFunc(onMouseDown);
	glutReshapeFunc(reshape);
	glutIdleFunc(onIdle);

	glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	TwGLUTModifiersFunc(glutGetModifiers);

	initGL();

	init();

	glutMainLoop();

	TwTerminate();

	return 0;
}










