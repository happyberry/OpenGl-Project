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
#define GLM_FORCE_SWIZZLE

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "constants.h"
#include "lodepng.h"
#include "shaderprogram.h"
#include "myCube.h"
#include "myTeapot.h"

float speed_x=0;
float speed_y=0;
float aspectRatio=1;
float speed = 0.1;
float mouseLastX;
float mouseLastY;
float sensitivity = 0.2;
float yaw = 90, pitch = 0;
bool firstMouse = true;

glm::vec3 zeroVec3 = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraPosition = glm::vec3(0.0f, 0.0f, -5.0f);
glm::vec3 cameraRight = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), cameraFront));
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 ellipsoid = glm::vec3(1.0f, 1.0f, 1.0f);

glm::mat3 cbmMatrix = glm::mat3(
    glm::vec3(float(1) / ellipsoid.x, 0, 0),
    glm::vec3(0, float(1) / ellipsoid.y, 0),
    glm::vec3(0, 0, float(1) / ellipsoid.z));

ShaderProgram *sp;

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
    if (a == 0) std::cout <<"DDD";
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

bool isCollision(glm::vec3 playerPosition) {

    float t0, t1, signedDistance, den, planeConstant, root, A, B, C, f, squaredEdgeLength, dotEdgeBTV, dotEdgeVel;
    float *rootptr = &root, *Aptr = &A, *Bptr = &B, *Cptr = &C, maxT;
    bool collisionFound = false, embeddedInPlane = false;

    glm::vec3 ellipsoidVelocity = cbmMatrix * velocity;
    glm::vec3 basePoint = cbmMatrix * playerPosition;
    glm::vec3 normalToPlane;
    glm::vec3 intersectionPoint;

    for (int i = 0; i < myCubeVertexCount * 4; i += 12) {

        glm::vec3 pointA = glm::vec3(myCubeVertices[i], myCubeVertices[i + 1], myCubeVertices[i + 2]);
        glm::vec3 pointB = glm::vec3(myCubeVertices[i + 4], myCubeVertices[i + 5], myCubeVertices[i + 6]);
        glm::vec3 pointC = glm::vec3(myCubeVertices[i + 8], myCubeVertices[i + 9], myCubeVertices[i + 10]);

        pointA = cbmMatrix * pointA;
        pointB = cbmMatrix * pointB;
        pointC = cbmMatrix * pointC;

        normalToPlane = glm::cross((pointC - pointA), (pointB - pointA));
        normalToPlane = glm::normalize(normalToPlane);

        planeConstant = -(normalToPlane.x * pointA.x + normalToPlane.y * pointA.y + normalToPlane.z * pointA.z);
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
        if (isPointInTriangle(pointA, pointB, pointC, intersectionPoint)) {
            std::cout << "INSIDE " << t0 << std::endl;
            std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
            std::cout << "Vector normal to triangle plane: " << normalToPlane.x << " " <<   normalToPlane.y << " " << normalToPlane.z << std::endl;
            std::cout << "Velocity vector: " << ellipsoidVelocity.x << " " <<   ellipsoidVelocity.y << " " << ellipsoidVelocity.z << std::endl;
            std::cout << "Intersection point: " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << std::endl;
            std::cout << "Intersection distance: " << t0 * glm::length(ellipsoidVelocity) << std::endl << std::endl;
            collisionFound = true;
        }
        else {
            setVertexCoefficients(pointA, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                std::cout << "VERTEX " << t0 << " " << *rootptr << std::endl;
                std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                std::cout << "Intersection point: " << pointA.x << " " << pointA.y << " " << pointA.z << std::endl;
                std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                collisionFound = true;
            }

            setVertexCoefficients(pointB, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                std::cout << "VERTEX " << t0 << " " << *rootptr << std::endl;
                std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                std::cout << "Intersection point: " << pointB.x << " " << pointB.y << " " << pointB.z << std::endl;
                std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                collisionFound = true;
            }

            setVertexCoefficients(pointC, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                std::cout << "VERTEX " << t0 << " " << *rootptr << std::endl;
                std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                std::cout << "Intersection point: " << pointC.x << " " << pointC.y << " " << pointC.z << std::endl;
                std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                collisionFound = true;
            }

            setEdgeCoefficients(pointA, pointB, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.0) {
                    intersectionPoint = pointA + f * (pointB - pointA);
                    std::cout << "EDGE: " << std::endl;
                    std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                    std::cout << "Intersection point: " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << std::endl << std::endl;
                    std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                    collisionFound = true;
                }
            }

            setEdgeCoefficients(pointB, pointC, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.1) {
                    intersectionPoint = pointB + f * (pointC - pointB);
                    std::cout << "EDGE: " << std::endl;
                    std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                    std::cout << "Intersection point: " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << std::endl << std::endl;
                    std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                    collisionFound = true;
                }
            }

            setEdgeCoefficients(pointC, pointA, basePoint, ellipsoidVelocity, Aptr, Bptr, Cptr, squaredEdgeLength, dotEdgeBTV, dotEdgeVel);
            if (solveQuadraticEquation(A, B, C, rootptr, maxT)) {
                f = (dotEdgeVel * (*rootptr) - dotEdgeBTV) / squaredEdgeLength;
                if (f >= 0.0 && f <= 1.1) {
                    intersectionPoint = pointC + f * (pointA - pointC);
                    std::cout << "EDGE: " << std::endl;
                    std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
                    std::cout << "Intersection point: " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << std::endl << std::endl;
                    std::cout << "Intersection distance: " << *rootptr * glm::length(ellipsoidVelocity) << std::endl << std::endl;
                    collisionFound = true;
                }
            }
        }
    }
    //std::cout << "CAMERA POSITION: " <<  playerPosition.x << " " <<   playerPosition.y << " " << playerPosition.z << std::endl;
    return collisionFound;
}

//Procedura obsługi błędów
void error_callback(int error, const char* description) {
	fputs(description, stderr);
}


void keyCallback(GLFWwindow* window,int key,int scancode,int action,int mods) {
    /*if (action==GLFW_RELEASE) {
        if (key==GLFW_KEY_LEFT) velocity = zeroVec3;
        if (key==GLFW_KEY_RIGHT) velocity = zeroVec3;
        if (key==GLFW_KEY_UP) velocity = zeroVec3;
        if (key==GLFW_KEY_DOWN) velocity = zeroVec3;
    }
    if (action==GLFW_PRESS) {
        if (key==GLFW_KEY_LEFT) velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * -speed;
        if (key==GLFW_KEY_RIGHT) velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * speed;
        if (key==GLFW_KEY_UP) velocity = cameraFront * speed;
        if (key==GLFW_KEY_DOWN) velocity = cameraFront * -speed;
    }*/


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

void windowResizeCallback(GLFWwindow* window,int width,int height) {
    if (height==0) return;
    aspectRatio=(float)width/(float)height;
    glViewport(0,0,width,height);
}

//Procedura inicjująca
void initOpenGLProgram(GLFWwindow* window) {
	//************Tutaj umieszczaj kod, który należy wykonać raz, na początku programu************
	glClearColor(0,0,0,1);
	glEnable(GL_DEPTH_TEST);
	glfwSetWindowSizeCallback(window,windowResizeCallback);
	glfwSetKeyCallback(window,keyCallback);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	sp=new ShaderProgram("v_simplest.glsl",NULL,"f_simplest.glsl");
}


//Zwolnienie zasobów zajętych przez program
void freeOpenGLProgram(GLFWwindow* window) {
    //************Tutaj umieszczaj kod, który należy wykonać po zakończeniu pętli głównej************

    delete sp;
}




//Procedura rysująca zawartość sceny
void drawScene(GLFWwindow* window,float angle_x,float angle_y) {
	//************Tutaj umieszczaj kod rysujący obraz******************l
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * -speed;
        if (!isCollision(cameraPosition)) cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * speed;
        if (!isCollision(cameraPosition)) cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        velocity = cameraFront * speed;
        if (!isCollision(cameraPosition)) cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        velocity = cameraFront * -speed;
        if (!isCollision(cameraPosition)) cameraPosition += velocity;
    }
	glm::mat4 V=glm::lookAt(
         cameraPosition,
         cameraPosition + cameraFront,
         cameraUp); //Wylicz macierz widoku

    glm::mat4 P=glm::perspective(50.0f*PI/180.0f, aspectRatio, 0.01f, 50.0f); //Wylicz macierz rzutowania

    glm::mat4 M=glm::mat4(1.0f);
	M=glm::rotate(M,angle_y,glm::vec3(1.0f,0.0f,0.0f)); //Wylicz macierz modelu
	M=glm::rotate(M,angle_x,glm::vec3(0.0f,1.0f,0.0f)); //Wylicz macierz modelu

    sp->use();//Aktywacja programu cieniującego
    //Przeslij parametry programu cieniującego do karty graficznej
    glUniformMatrix4fv(sp->u("P"),1,false,glm::value_ptr(P));
    glUniformMatrix4fv(sp->u("V"),1,false,glm::value_ptr(V));
    glUniformMatrix4fv(sp->u("M"),1,false,glm::value_ptr(M));

    glEnableVertexAttribArray(sp->a("vertex"));  //Włącz przesyłanie danych do atrybutu vertex
    glEnableVertexAttribArray(sp->a("color"));
    glVertexAttribPointer(sp->a("vertex"),4,GL_FLOAT,false,0,myCubeVertices); //Wskaż tablicę z danymi dla atrybutu vertex
    glVertexAttribPointer(sp->a("color"),4,GL_FLOAT,false,0,myCubeColors);

    glDrawArrays(GL_TRIANGLES,0,myCubeVertexCount); //Narysuj obiekt

    //glVertexAttribPointer(sp->a("vertex"),4,GL_FLOAT,false,0,myTeapotVertices); //Wskaż tablicę z danymi dla atrybutu vertex
    //glVertexAttribPointer(sp->a("color"),4,GL_FLOAT,false,0,myTeapotColors);

    //glDrawArrays(GL_TRIANGLES,0,myTeapotVertexCount); //Narysuj obiekt

    glDisableVertexAttribArray(sp->a("vertex"));  //Wyłącz przesyłanie danych do atrybutu vertex

    glfwSwapBuffers(window); //Przerzuć tylny bufor na przedni
}


int main(void)
{
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
	glfwSetTime(0); //Zeruj timer
	while (!glfwWindowShouldClose(window)) //Tak długo jak okno nie powinno zostać zamknięte
	{
        angle_x+=speed_x*glfwGetTime(); //Zwiększ/zmniejsz kąt obrotu na podstawie prędkości i czasu jaki upłynał od poprzedniej klatki
        angle_y+=speed_y*glfwGetTime(); //Zwiększ/zmniejsz kąt obrotu na podstawie prędkości i czasu jaki upłynał od poprzedniej klatki
        glfwSetTime(0); //Zeruj timer
		drawScene(window,angle_x,angle_y); //Wykonaj procedurę rysującą
		glfwPollEvents(); //Wykonaj procedury callback w zalezności od zdarzeń jakie zaszły.
	}

	freeOpenGLProgram(window);

	glfwDestroyWindow(window); //Usuń kontekst OpenGL i okno
	glfwTerminate(); //Zwolnij zasoby zajęte przez GLFW
	exit(EXIT_SUCCESS);
}
