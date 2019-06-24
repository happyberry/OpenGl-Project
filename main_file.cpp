/*
Niniejszy program jest wolnym oprogramowaniem; możesz go
rozprowadzać dalej i / lub modyfikować na warunkach Powszechnej
Licencji Publicznej GNU, wydanej przez Fundację Wolnego
Oprogramowania - według wersji 2 tej Licencji lub(według twojego
wyboru) którejś z późniejszych wersji.
Niniejszy program rozpowszechniany jest z nadzieją, iż będzie on
użyteczny - jednak BEZ JAKIEJKOLWIEK GWARANCJI, nawet domyślnej
gwarancji PRZYDATNOŚCI HANDLOWEJ albo PRZYDATNOŚCI DO OKREŚLONYCH
ZASTOSOWAŃ.W celu uzyskania bliższych informacji sięgnij do
Powszechnej Licencji Publicznej GNU.
Z pewnością wraz z niniejszym programem otrzymałeś też egzemplarz
Powszechnej Licencji Publicznej GNU(GNU General Public License);
jeśli nie - napisz do Free Software Foundation, Inc., 59 Temple
Place, Fifth Floor, Boston, MA  02110 - 1301  USA
*/

#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "allmodels.h"
#include<iostream>
#include "lodepng.h"
#include "shaderprogram.h"
#include "myCube.h"
#include "myStairs.h"
//#include "nCube.h"
using namespace std;

float speed_x=0;
float speed_y=0;
float aspectRatio=1;

float speed = 0.1;
float mouseLastX;
float mouseLastY;
float sensitivity = 0.2;
float yaw = 90, pitch = 0;
int collisionRecursionDepth;
bool firstMouse = true;
GLuint tex;

glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraPosition = glm::vec3(0.0f, 5.0f, 0.0f);
glm::vec3 cameraRight = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), cameraFront));
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 gravity = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 ellipsoid = glm::vec3(0.5f, 1.0f, 0.5f);

glm::mat3 cbmMatrix = glm::mat3(
    glm::vec3(float(1) / ellipsoid.x, 0, 0),
    glm::vec3(0, float(1) / ellipsoid.y, 0),
    glm::vec3(0, 0, float(1) / ellipsoid.z));




int MazeBase[19][19] = {
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0},
{0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0},
{0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0},
{0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0},
{0,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0},
{0,1,1,1,1,0,1,1,0,0,0,1,0,3,0,0,0,1,0},
{0,1,1,1,1,0,1,1,0,1,1,1,0,2,2,2,0,1,0},
{0,1,1,1,1,0,0,0,0,0,0,0,2,2,0,2,0,0,0},
{0,1,1,0,1,1,1,1,0,1,1,0,2,0,0,2,2,2,0},
{0,0,1,0,0,1,0,0,0,1,1,0,3,0,0,0,0,3,0},
{0,0,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0},
{0,0,0,1,0,1,0,0,0,3,0,0,0,1,0,1,0,0,0},
{0,1,1,1,0,1,0,2,2,2,0,1,1,1,1,1,1,1,0},
{0,1,1,1,0,0,0,2,2,2,0,0,0,1,1,1,1,1,0},
{0,1,1,1,1,1,0,2,2,2,2,2,0,1,1,1,1,1,0},
{0,1,1,1,0,1,0,0,0,0,0,3,0,1,1,1,1,1,0},
{0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}};

float MazeElements[361][3];


//Procedura obsługi błędów
void error_callback(int error, const char* description) {
	fputs(description, stderr);
}


void keyCallback(GLFWwindow* window,int key,int scancode,int action,int mods) {
    /* */
}

void mouseCallback(GLFWwindow* window, double xpos, double ypos) {

    float xOffset = xpos - mouseLastX;
    float yOffset = mouseLastY - ypos;

    //std::cout << cbmMatrix[1].y << std::endl;

    mouseLastX = xpos;
    mouseLastY = ypos;

    if (firstMouse) {
        firstMouse = false;
        return;
    }

    xOffset *= sensitivity;
    yOffset *= sensitivity;

    yaw += xOffset;
    pitch += yOffset;

    pitch = std::min(pitch, float(89));
    pitch = std::max(pitch, float(-89));

    glm::vec3 frontt;
    frontt.x = cos(glm::radians(pitch)) * cos(glm::radians(yaw));
    frontt.y = sin(glm::radians(pitch));
    frontt.z = cos(glm::radians(pitch)) * sin(glm::radians(yaw));
    // printf("%f %f %f\n", cameraFront.x, cameraFront.y, cameraFront.z);
    cameraFront = glm::normalize(frontt);
}

GLuint readTexture(char* filename) {
  GLuint tex;
  glActiveTexture(GL_TEXTURE0);

  //Wczytanie do pamięci komputera
  std::vector<unsigned char> image;   //Alokuj wektor do wczytania obrazka
  unsigned width, height;   //Zmienne do których wczytamy wymiary obrazka
  //Wczytaj obrazek
  unsigned error = lodepng::decode(image, width, height, filename);

  //Import do pamięci karty graficznej
  glGenTextures(1,&tex); //Zainicjuj jeden uchwyt
  glBindTexture(GL_TEXTURE_2D, tex); //Uaktywnij uchwyt
  //Wczytaj obrazek do pamięci KG skojarzonej z uchwytem
  glTexImage2D(GL_TEXTURE_2D, 0, 4, width, height, 0,
    GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char*) image.data());

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  return tex;
}

void windowResizeCallback(GLFWwindow* window,int width,int height) {
    if (height==0) return;
    aspectRatio=(float)width/(float)height;
    glViewport(0,0,width,height);
}

//Procedura inicjująca
void initOpenGLProgram(GLFWwindow* window) {
    initShaders();
    tex=readTexture("bricks.png");
	//************Tutaj umieszczaj kod, który należy wykonać raz, na początku programu************
	glClearColor(0,0,0,1);
	glEnable(GL_DEPTH_TEST);
	glfwSetWindowSizeCallback(window,windowResizeCallback);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetKeyCallback(window,keyCallback);
}


//Zwolnienie zasobów zajętych przez program
void freeOpenGLProgram(GLFWwindow* window) {
    freeShaders();
    glDeleteTextures(1,&tex);
    //************Tutaj umieszczaj kod, który należy wykonać po zakończeniu pętli głównej************
}

void setVertexCoefficients(glm::vec3 point, glm::vec3 basePoint, glm::vec3 vel, float *Aptr, float *Bptr, float *Cptr) {
    *Aptr = glm::dot(vel, vel);
    *Bptr = 2 * (glm::dot(vel, basePoint - point));
    *Cptr = glm::length(point - basePoint) * glm::length(point - basePoint) - 1;
}

void setEdgeCoefficients(glm::vec3 point1, glm::vec3 point2, glm::vec3 basePoint, glm::vec3 vel, float *Aptr, float *Bptr, float *Cptr,
                         float &squaredEdgeLength, float &dotEdgeBTV, float &dotEdgeVel) {
    glm::vec3 edge = point2 - point1;
    glm::vec3 baseToVertex = point1 - basePoint;

    squaredEdgeLength = glm::length(edge) * glm::length(edge);
    dotEdgeBTV = glm::dot(edge, baseToVertex);
    dotEdgeVel = glm::dot(edge, vel);

    *Aptr = squaredEdgeLength * -1 * (glm::length(vel) * glm::length(vel)) + dotEdgeVel * dotEdgeVel;
    *Bptr = squaredEdgeLength * (2 * glm::dot(vel, baseToVertex)) - 2 *(dotEdgeVel * dotEdgeBTV);
    *Cptr = squaredEdgeLength * (1 - glm::length(baseToVertex) * glm::length(baseToVertex)) + dotEdgeBTV * dotEdgeBTV;
}

bool solveQuadraticEquation(float a, float b, float c, float *rootptr, float border) {
    //if (a == 0) std::cout <<"DDD";
    float delta = b * b - 4 * a * c;
    float root1, root2;
    if (delta < 0) return false;
    delta = std::sqrt(delta);
    root1 = (-b - delta)/(2 * a);
    root2 = (-b + delta)/(2 * a);
    //std::cout << root1 << " " << root2 << " " << border << std::endl;
    if (root1 > root2) std::swap(root1, root2);
    *rootptr = root1;
    if (root1 > 0 && root1 < border) return true;
    *rootptr = root2;
    if (root2 > 0 && root2 < border) return true;
    return false;
}



bool isPointInTriangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 point) {

    glm::vec3 ab = b - a;
    glm::vec3 ac = c - a;
    glm::vec3 apoint = point - a;

    float dotacac = glm::dot(ac, ac);
    float dotabab = glm::dot(ab, ab);
    float dotabac = glm::dot(ab, ac);
    float dotapointac = glm::dot(apoint, ac);
    float dotapointab = glm::dot(apoint, ab);

    float inversedDen = 1 / (dotabab * dotacac - dotabac * dotabac);
    float u = (dotabab * dotapointac - dotabac * dotapointab) * inversedDen;
    float v = (dotacac * dotapointab - dotabac * dotapointac) * inversedDen;

    return (u >= 0) && (v >= 0) && (u + v < 1);
}

bool isCollision(glm::vec3 playerPosition, glm::vec3 velocity, float finalT, float triangles[][12], int triangleIndex) {

    // DETECTION

    //std::cout << "DEPTH: " << collisionRecursionDepth << std::endl;
    float t0, t1, signedDistance, den, planeConstant, root, A, B, C, f, squaredEdgeLength, dotEdgeBTV, dotEdgeVel;
    float *rootptr = &root, *Aptr = &A, *Bptr = &B, *Cptr = &C, maxT, finalDistance;
    bool collisionFound = false, embeddedInPlane = false;

    glm::vec3 ellipsoidVelocity = cbmMatrix * velocity;
    glm::vec3 basePoint = cbmMatrix * playerPosition;
    glm::vec3 normalToPlane;
    glm::vec3 intersectionPoint;
    glm::vec3 resultingPoint;

    for (int i = 0; i < triangleIndex; i++) {

        glm::vec3 pointA = glm::vec3(triangles[i][0], triangles[i][1], triangles[i][2]);
        glm::vec3 pointB = glm::vec3(triangles[i][4], triangles[i][5], triangles[i][6]);
        glm::vec3 pointC = glm::vec3(triangles[i][8], triangles[i][9], triangles[i][10]);

        pointA = cbmMatrix * pointA;
        pointB = cbmMatrix * pointB;
        pointC = cbmMatrix * pointC;

        normalToPlane = -glm::cross((pointC - pointA), (pointB - pointA));
        normalToPlane = glm::normalize(normalToPlane);

        planeConstant = -glm::dot(normalToPlane, pointA);
        signedDistance = glm::dot(normalToPlane, basePoint) + planeConstant;
        den = glm::dot(normalToPlane, ellipsoidVelocity);

        if (den == 0.0f) {
            if (std::abs(signedDistance) < 1) {
                t0 = 0;
                t1 = 1;
                embeddedInPlane = true;
            }
            else continue;
        }
        else {
            t0 = (-1 - signedDistance) / den;
            t1 = (1 - signedDistance) / den;
        }

        if (t0 > t1) std::swap(t0, t1);
        if (t0 > 1 || t1 < 0) continue;

        t0 = glm::clamp(t0, 0.0f, 1.0f);
        t1 = glm::clamp(t1, 0.0f, 1.0f);
        intersectionPoint = basePoint - normalToPlane + t0 * ellipsoidVelocity;
        maxT = 1.0;
        if (isPointInTriangle(pointA, pointB, pointC, intersectionPoint) && t0 < finalT) {
            finalT = t0;
            resultingPoint = intersectionPoint;
            finalDistance = t0 * glm::length(ellipsoidVelocity);
            collisionFound = true;
            cout << "INTERSECTION POINT: " << resultingPoint.x << " " << resultingPoint.y << " " << resultingPoint.z << endl;
        }
        else {
            setVertexCoefficients(pointA, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                finalT = *rootptr;
                resultingPoint = pointA;
                finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                collisionFound = true;
                cout << "V" << endl;
            }

            setVertexCoefficients(pointB, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                finalT = *rootptr;
                resultingPoint = pointB;
                finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                collisionFound = true;
                cout << "V" << endl;
            }

            setVertexCoefficients(pointC, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                finalT = *rootptr;
                resultingPoint = pointC;
                finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                collisionFound = true;
                cout << "V" << endl;
            }

            setEdgeCoefficients(pointA, pointB, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.0) {
                    intersectionPoint = pointA + f * (pointB - pointA);
                    finalT = *rootptr;
                    resultingPoint = intersectionPoint;
                    finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                    collisionFound = true;
                    cout << "e" << endl;
                    cout << "INTERSECTION POINT: " << resultingPoint.x << " " << resultingPoint.y << " " << resultingPoint.z << endl;
                }
            }

            setEdgeCoefficients(pointB, pointC, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.0) {
                    intersectionPoint = pointB + f * (pointC - pointB);
                    finalT = *rootptr;
                    resultingPoint = intersectionPoint;
                    finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                    collisionFound = true;
                    cout << "e" << endl;
                    cout << "INTERSECTION POINT: " << resultingPoint.x << " " << resultingPoint.y << " " << resultingPoint.z << endl;
                }
            }

            setEdgeCoefficients(pointC, pointA, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT) && *rootptr < finalT) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.0) {
                    intersectionPoint = pointC + f * (pointA - pointC);
                    finalT = *rootptr;
                    resultingPoint = intersectionPoint;
                    finalDistance = *rootptr * glm::length(ellipsoidVelocity);
                    collisionFound = true;
                    cout << "e" << endl;
                    cout << "INTERSECTION POINT: " << resultingPoint.x << " " << resultingPoint.y << " " << resultingPoint.z << endl;
                }
            }
        }
    }

    // RESPONSE

    if (collisionFound) {
        std::cout << "Camera point: " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
        /*std::cout << "Intersection point: " << resultingPoint.x << " " << resultingPoint.y << " " << resultingPoint.z << std::endl;
        std::cout << "Intersection distance: " << finalDistance << std::endl;
        std::cout << "Time before collision:" << finalT << std::endl;*/
        /*glm::vec3 newPosition = basePoint + finalT * ellipsoidVelocity;
        glm::vec3 planeNormal = newPosition - resultingPoint;
        glm::vec3 destination = basePoint + ellipsoidVelocity;
        glm::normalize(planeNormal);
        float distance = glm::dot(planeNormal, destination) - (planeNormal.x * resultingPoint.x + planeNormal.y * resultingPoint.y + planeNormal.z * resultingPoint.z);
        //std::cout << "Distance: " << distance << std::endl;
        //std::cout << "Sliding Plane Normal:" << planeNormal.x << " " <<  planeNormal.y << " " << planeNormal.z << std::endl;
        glm::vec3 newResultingPoint = destination - distance * planeNormal;
        glm::vec3 finalVelocity = newResultingPoint - resultingPoint;
        finalVelocity = inversedCbmMatrix * finalVelocity;
        //std::cout << "Velocity: " << finalVelocity.x << " " << finalVelocity.y << " " << finalVelocity.z << std::endl << std::endl;
        cameraPosition += finalVelocity;
        if (++collisionRecursionDepth < 5) {
            return isCollision(inversedCbmMatrix * cameraPosition, inversedCbmMatrix * finalVelocity, maxT - finalT, triangles, triangleIndex);
        }*/
    }
    else cameraPosition += velocity;
    return collisionFound;
}

void CreateMazeVec(int Check, int Check2){
    int MazeIndex = 0;
    for(int i = 0;i < 19; i++){
        for(int j = 0;j < 19; j++){
            if(MazeBase[i][j] == Check || MazeBase[i][j] == Check2){
                MazeElements[MazeIndex][0] = float(18 - 2*j);
                MazeElements[MazeIndex][1] = 0.0f;
                MazeElements[MazeIndex][2] = float(37 - 2*i);
            }
            else{
                MazeElements[MazeIndex][0] = 0.0f;
                MazeElements[MazeIndex][1] = 0.0f;
                MazeElements[MazeIndex][2] = 0.0f;
            }
            MazeIndex++;
        }
    }
}

//Procedura rysująca zawartość sceny
void drawScene(GLFWwindow* window,float angle_x,float angle_y) {
	//************Tutaj umieszczaj kod rysujący obraz******************l
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glm::mat4 V=glm::lookAt(
         cameraPosition,
         cameraPosition + cameraFront,
         cameraUp); //Wylicz macierz widoku


    glm::mat4 P=glm::perspective(50.0f*PI/180.0f, aspectRatio, 0.01f, 50.0f); //Wylicz macierz rzutowania

    glm::mat4 M=glm::mat4(1.0f);

    //glm::mat4 M1 = glm::translate(M,glm::vec3(-2.0f,0.0f,-1.0f));
    //glUniformMatrix4fv(spLambert->u("M"),1,false,glm::value_ptr(M1));

    spTextured->use();
    glUniformMatrix4fv(spTextured->u("P"),1,false,glm::value_ptr(P));
    glUniformMatrix4fv(spTextured->u("V"),1,false,glm::value_ptr(V));
    glUniform3fv(spTextured->u("cameraPosition"), 1, glm::value_ptr(cameraPosition));

    glEnableVertexAttribArray(spTextured->a("vertex"));
    glVertexAttribPointer(spTextured->a("vertex"),4,GL_FLOAT,false,0,myCubeVertices);

    glEnableVertexAttribArray(spTextured->a("texCoord"));
    glVertexAttribPointer(spTextured->a("texCoord"),2,GL_FLOAT,false,0,texCoords);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D,tex);
    glUniform1i(spLambertTextured->u("tex"),0);

    CreateMazeVec(0,0);
    for(int i = 0; i <361; i++){//mur w ziemi nizszy
        if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){

            for(float j = -2.0f; j < 5; j+= 2.0f){//4 warstwy muru

                glm::mat4 M1 = glm::translate(M,glm::vec3(MazeElements[i][0],j,MazeElements[i][2]));
                glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
                glDrawArrays( GL_TRIANGLES, 0, myCubeVertexCount );
            }

        }
    }

    glDisableVertexAttribArray(spTextured->a("vertex"));
    glDisableVertexAttribArray(spTextured->a("texCoord"));


    spColored->use();
    glUniformMatrix4fv(spColored->u("P"),1,false,glm::value_ptr(P));
    glUniformMatrix4fv(spColored->u("V"),1,false,glm::value_ptr(V));

    glEnableVertexAttribArray(spColored->a("vertex"));
    glVertexAttribPointer(spColored->a("vertex"),4,GL_FLOAT,false,0,myCubeVertices);

    CreateMazeVec(1,3);
    glEnableVertexAttribArray(spColored->a("color"));
    glVertexAttribPointer(spColored->a("color"),4,GL_FLOAT,false,0,myCubeColors1);

    for(int i = 0; i <361; i++){//podloga 1 poziom
        if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){
            glm::mat4 M1 = glm::translate(M,glm::vec3(MazeElements[i][0],-2.0f,MazeElements[i][2]));
            glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
            glDrawArrays( GL_TRIANGLES, 0, myCubeVertexCount );
        }
    }

    CreateMazeVec(2,2);
    for(int i = 0; i <361; i++){//podloga 1 poziom
        if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){
            glm::mat4 M1 = glm::translate(M,glm::vec3(MazeElements[i][0],-2.0f,MazeElements[i][2]));
            glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
            glDrawArrays( GL_TRIANGLES, 0, myCubeVertexCount );
            M1 = glm::translate(M,glm::vec3(MazeElements[i][0],0.0f,MazeElements[i][2]));
            glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
            glDrawArrays( GL_TRIANGLES, 0, myCubeVertexCount );
        }
    }

    glEnableVertexAttribArray(spColored->a("vertex"));
    glVertexAttribPointer(spColored->a("vertex"),4,GL_FLOAT,false,0,myStairsVertices);

    glEnableVertexAttribArray(spColored->a("color"));
    glVertexAttribPointer(spColored->a("color"),4,GL_FLOAT,false,0,myStairsColors);

    //PIERWSZE SCHODY
    glm::mat4 M1 = glm::translate(M,glm::vec3(-4.0f, 0.0f, 5.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(-6.0f, 0.0f, 17.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(-16.0f, 0.0f, 17.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    //glEnableVertexAttribArray(spColored->a("vertex"));
    //glVertexAttribPointer(spColored->a("vertex"),4,GL_FLOAT,false,0,myStairsVerticesRotated);

    M1 = glm::translate(M,glm::vec3(-8.0f, 0.0f, 25.0f));
    M1 = glm::rotate(M1, 3.14f, glm::vec3(0.0f, 1.0f, 0.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(0.0f, 0.0f, 13.0f));
    M1 = glm::rotate(M1, 3.14f, glm::vec3(0.0f, 1.0f, 0.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);


    glfwSwapBuffers(window); //Przerzuć tylny bufor na przedni
}


int main(void)
{
        CreateMazeVec(0,0);
        static float triangles[100000][12];
        int triangleindex = 0;
        for(int i = 0; i <361; i++){//mur w ziemi nizszy
            if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){
                for(float k=-2.0f; k <5.0f; k+= 2.0f){
                    for(int j=0;j<12;j++){
                        triangles[triangleindex][0] = myCubeVertices[12*j]+MazeElements[i][0];
                        triangles[triangleindex][4] = myCubeVertices[12*j+4]+MazeElements[i][0];
                        triangles[triangleindex][8] = myCubeVertices[12*j+8]+MazeElements[i][0];
                        triangles[triangleindex][1] = myCubeVertices[12*j+1]+k;
                        triangles[triangleindex][5] = myCubeVertices[12*j+5]+k;
                        triangles[triangleindex][9] = myCubeVertices[12*j+9]+k;
                        triangles[triangleindex][2] = myCubeVertices[12*j+2]+MazeElements[i][2];
                        triangles[triangleindex][6] = myCubeVertices[12*j+6]+MazeElements[i][2];
                        triangles[triangleindex][10] = myCubeVertices[12*j+10]+MazeElements[i][2];
                        triangles[triangleindex][3] = 1.0f;
                        triangles[triangleindex][7] = 1.0f;
                        triangles[triangleindex][11] = 1.0f;
                        triangleindex++;
                    }
                }
            }

        }

    CreateMazeVec(1,3);
    for(int i = 0; i <361; i++){//podloga 1 poziom
        if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){
             for(int j=0;j<12;j++){
                triangles[triangleindex][0] = myCubeVertices[12*j]+MazeElements[i][0];
                triangles[triangleindex][4] = myCubeVertices[12*j+4]+MazeElements[i][0];
                triangles[triangleindex][8] = myCubeVertices[12*j+8]+MazeElements[i][0];
                triangles[triangleindex][1] = myCubeVertices[12*j+1]-2.0f;
                triangles[triangleindex][5] = myCubeVertices[12*j+5]-2.0f;
                triangles[triangleindex][9] = myCubeVertices[12*j+9]-2.0f;
                triangles[triangleindex][2] = myCubeVertices[12*j+2]+MazeElements[i][2];
                triangles[triangleindex][6] = myCubeVertices[12*j+6]+MazeElements[i][2];
                triangles[triangleindex][10] = myCubeVertices[12*j+10]+MazeElements[i][2];
                triangles[triangleindex][3] = 1.0f;
                triangles[triangleindex][7] = 1.0f;
                triangles[triangleindex][11] = 1.0f;
                triangleindex++;
            }
        }
    }

    CreateMazeVec(2,2);
    for(int i = 0; i <361; i++){//podloga 1 poziom
        if(!(MazeElements[i][0] == 0.0f && MazeElements[i][1] == 0.0f && MazeElements[i][2] == 0.0f)){
            for(float j =-2.0; j < 1.0f; j += 2.0f){
                for(int j=0;j<12;j++){
                    triangles[triangleindex][0] = myCubeVertices[12*j]+MazeElements[i][0];
                    triangles[triangleindex][4] = myCubeVertices[12*j+4]+MazeElements[i][0];
                    triangles[triangleindex][8] = myCubeVertices[12*j+8]+MazeElements[i][0];
                    triangles[triangleindex][1] = myCubeVertices[12*j+1]+j;
                    triangles[triangleindex][5] = myCubeVertices[12*j+5]+j;
                    triangles[triangleindex][9] = myCubeVertices[12*j+9]+j;
                    triangles[triangleindex][2] = myCubeVertices[12*j+2]+MazeElements[i][2];
                    triangles[triangleindex][6] = myCubeVertices[12*j+6]+MazeElements[i][2];
                    triangles[triangleindex][10] = myCubeVertices[12*j+10]+MazeElements[i][2];
                    triangles[triangleindex][3] = 1.0f;
                    triangles[triangleindex][7] = 1.0f;
                    triangles[triangleindex][11] = 1.0f;
                    triangleindex++;
                }
            }
        }
    }
    float StairsTranslations[5][3] = {
    {-4,0,5},
    {-6,0,17},
    {-16,0,17},
    {-8,0,25},
    {0,0,13}};

    for(int k = 0; k < 3; k++){
        for(int j = 0; j < 24; j++)
        {
            int sign;
            if(k<3)sign=1; else sign = -1;
            triangles[triangleindex][0] = sign*myStairsVertices[24*j]+StairsTranslations[k][0];
            triangles[triangleindex][4] = sign*myStairsVertices[24*j+4]+StairsTranslations[k][0];
            triangles[triangleindex][8] = sign*myStairsVertices[24*j+8]+StairsTranslations[k][0];
            triangles[triangleindex][1] = myStairsVertices[24*j+1]+StairsTranslations[k][1];
            triangles[triangleindex][5] = myStairsVertices[24*j+5]+StairsTranslations[k][1];
            triangles[triangleindex][9] = myStairsVertices[24*j+9]+StairsTranslations[k][1];
            triangles[triangleindex][2] = sign*myStairsVertices[24*j+2]+StairsTranslations[k][2];
            triangles[triangleindex][6] = sign*myStairsVertices[24*j+6]+StairsTranslations[k][2];
            triangles[triangleindex][10] = sign*myStairsVertices[24*j+10]+StairsTranslations[k][2];
            triangles[triangleindex][3] = 1.0f;
            triangles[triangleindex][7] = 1.0f;
            triangles[triangleindex][11] = 1.0f;
            triangleindex++;
        }
    }

    /*StairsTranslations2[2][3] = {
    {-8,0,25},
    {0,0,13}};
    for(int k = 0; k < 2; k++){
        for(int j = 0; j < 24; j++)
        {
            triangles[triangleindex][0] = myStairsVertices[24*j]+StairsTranslations2[k][0];
            triangles[triangleindex][4] = myStairsVertices[24*j+4]+StairsTranslations2[k][0];
            triangles[triangleindex][8] = myStairsVertices[24*j+8]+StairsTranslations2[k][0];
            triangles[triangleindex][1] = myStairsVertices[24*j+1]+StairsTranslations2[k][1];
            triangles[triangleindex][5] = myStairsVertices[24*j+5]+StairsTranslations2[k][1];
            triangles[triangleindex][9] = myStairsVertices[24*j+9]+StairsTranslations2[k][1];
            triangles[triangleindex][2] = myStairsVertices[24*j+2]+StairsTranslations2[k][2];
            triangles[triangleindex][6] = myStairsVertices[24*j+6]+StairsTranslations2[k][2];
            triangles[triangleindex][10] = myStairsVertices[24*j+10]+StairsTranslations2[k][2];
            triangles[triangleindex][3] = 1.0f;
            triangles[triangleindex][7] = 1.0f;
            triangles[triangleindex][11] = 1.0f;
            triangleindex++;
        }
    }*/




    printf("%d\n",triangleindex);

	GLFWwindow* window; //Wskaźnik na obiekt reprezentujący okno

	glfwSetErrorCallback(error_callback);//Zarejestruj procedurę obsługi błędów

	if (!glfwInit()) { //Zainicjuj bibliotekę GLFW
		fprintf(stderr, "Nie można zainicjować GLFW.\n");
		exit(EXIT_FAILURE);
	}

	window = glfwCreateWindow(500, 500, "OpenGL", NULL, NULL);  //Utwórz okno 500x500 o tytule "OpenGL" i kontekst OpenGL.

	if (!window) //Jeżeli okna nie udało się utworzyć, to zamknij program
	{
		fprintf(stderr, "Nie można utworzyć okna.\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window); //Od tego momentu kontekst okna staje się aktywny i polecenia OpenGL będą dotyczyć właśnie jego.
	glfwSwapInterval(1); //Czekaj na 1 powrót plamki przed pokazaniem ukrytego bufora

	if (glewInit() != GLEW_OK) { //Zainicjuj bibliotekę GLEW
		fprintf(stderr, "Nie można zainicjować GLEW.\n");
		exit(EXIT_FAILURE);
	}

	initOpenGLProgram(window); //Operacje inicjujące

	//Główna pętla
	float angle_x=0; //Aktualny kąt obrotu obiektu
	float angle_y=0; //Aktualny kąt obrotu obiektu
	float accel;
	bool spacePressed;
	int framesSinceLastJump;
	glm::vec3 velocity;
	glfwSetTime(0); //Zeruj timer
	while (!glfwWindowShouldClose(window)) //Tak długo jak okno nie powinno zostać zamknięte
	{
        glfwSetTime(0); //Zeruj timer

        velocity = glm::vec3(0.0f, 0.0f, 0.0f);
        collisionRecursionDepth = 0;
        spacePressed = false;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
            //velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * -speed;
            //isCollision(cameraPosition, velocity, 1);

            velocity -= glm::normalize(glm::cross(cameraFront, cameraUp));
        }
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
            //velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * speed;
            //isCollision(cameraPosition, velocity, 1);

            velocity += glm::normalize(glm::cross(cameraFront, cameraUp));
        }
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
            //velocity = cameraFront * speed;
            //isCollision(cameraPosition, velocity, 1);

            velocity += cameraFront;
        }
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
            //velocity = cameraFront * -speed;
            //isCollision(cameraPosition, velocity, 1);

            velocity -= cameraFront;
        }
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && framesSinceLastJump >= 10) {
            //velocity = cameraFront * -speed;
            //isCollision(cameraPosition, velocity, 1);
            framesSinceLastJump = 0;
            velocity *= speed;
            spacePressed = true;
            velocity.y += 0.4;
        }
        else framesSinceLastJump++;
        if (velocity.y < 0) velocity.y = 0;
        //if (velocity.y > 0) velocity.y *= 2;
        if (velocity.x != 0 || velocity.y != 0 || velocity.z != 0) {
            if (!spacePressed) velocity *= speed;
            isCollision(cameraPosition, velocity, 1, triangles, triangleindex);
        }
        if(!isCollision(cameraPosition, gravity, 1, triangles, triangleindex)) {
            gravity.y -= 0.01;
        }
        else gravity.y = 0;

		drawScene(window,angle_x,angle_y); //Wykonaj procedurę rysującą
		glfwPollEvents(); //Wykonaj procedury callback w zalezności od zdarzeń jakie zaszły.
	}

	freeOpenGLProgram(window);

	glfwDestroyWindow(window); //Usuń kontekst OpenGL i okno
	glfwTerminate(); //Zwolnij zasoby zajęte przez GLFW
	exit(EXIT_SUCCESS);
}
