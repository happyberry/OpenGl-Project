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

void windowResizeCallback(GLFWwindow* window,int width,int height) {
    if (height==0) return;
    aspectRatio=(float)width/(float)height;
    glViewport(0,0,width,height);
}

//Procedura inicjująca
void initOpenGLProgram(GLFWwindow* window) {
    initShaders();
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
    //************Tutaj umieszczaj kod, który należy wykonać po zakończeniu pętli głównej************
}

void CreateMazeVec(int Check, int Check2)
{
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

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * -speed;
        //if (!isCollision(cameraPosition))
        cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        velocity = glm::normalize(glm::cross(cameraFront, cameraUp)) * speed;
        //if (!isCollision(cameraPosition))
        cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        velocity = cameraFront * speed;
        //if (!isCollision(cameraPosition))
        cameraPosition += velocity;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        velocity = cameraFront * -speed;
        //if (!isCollision(cameraPosition))
        cameraPosition += velocity;
    }
	glm::mat4 V=glm::lookAt(
         cameraPosition,
         cameraPosition + cameraFront,
         cameraUp); //Wylicz macierz widoku


    glm::mat4 P=glm::perspective(50.0f*PI/180.0f, aspectRatio, 0.01f, 50.0f); //Wylicz macierz rzutowania

    spLambert->use();//Aktywacja programu cieniującego
    //Przeslij parametry programu cieniującego do karty graficznej
    glUniform4f(spLambert->u("color"),0,1,0,1);
    glUniformMatrix4fv(spLambert->u("P"),1,false,glm::value_ptr(P));
    glUniformMatrix4fv(spLambert->u("V"),1,false,glm::value_ptr(V));


    glm::mat4 M=glm::mat4(1.0f);

    //glm::mat4 M1 = glm::translate(M,glm::vec3(-2.0f,0.0f,-1.0f));
    //glUniformMatrix4fv(spLambert->u("M"),1,false,glm::value_ptr(M1));

    spColored->use();
    glUniformMatrix4fv(spColored->u("P"),1,false,glm::value_ptr(P));
    glUniformMatrix4fv(spColored->u("V"),1,false,glm::value_ptr(V));

    glEnableVertexAttribArray(spColored->a("vertex"));
    glVertexAttribPointer(spColored->a("vertex"),4,GL_FLOAT,false,0,myCubeVertices);

    glEnableVertexAttribArray(spColored->a("color"));
    glVertexAttribPointer(spColored->a("color"),4,GL_FLOAT,false,0,myCubeColors);

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

    CreateMazeVec(1,3);
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

    M1 = glm::translate(M,glm::vec3(-8.0f, 0.0f, 25.0f));
    M1 = glm::rotate(M1, 3.14f, glm::vec3(0.0f, 1.0f, 0.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(0.0f, 0.0f, 13.0f));
    M1 = glm::rotate(M1, 3.14f, glm::vec3(0.0f, 1.0f, 0.0f));
    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);




    /*M1 = glm::translate(M,glm::vec3(-1.0f, -1.0f, 0.0f));//DRUGIE SCHODY
    M1 = glm::translate(M1,glm::vec3(2.0f, 0.0f, 14.0f));
    M1 = glm::scale(M1, glm::vec3(2.0f, 2.0f, 2.0f));
    M1 = glm::rotate(M1,3.14f, glm::vec3(0.0f, 1.0f, 0.0f));

    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(-1.0f, -1.0f, 0.0f));//TRZECIE
    M1 = glm::translate(M1,glm::vec3(-6.0f, 0.0f, 16.0f));
    M1 = glm::scale(M1, glm::vec3(2.0f, 2.0f, 2.0f));

    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(-1.0f, -1.0f, 0.0f));//CZWARTE
    M1 = glm::translate(M1,glm::vec3(-6.0f, 0.0f, 26.0f));
    M1 = glm::scale(M1, glm::vec3(2.0f, 2.0f, 2.0f));
    M1 = glm::rotate(M1,3.14f, glm::vec3(0.0f, 1.0f, 0.0f));

    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    M1 = glm::translate(M,glm::vec3(-1.0f, -1.0f, 0.0f));//CZWARTE
    M1 = glm::translate(M1,glm::vec3(-14.0f, 0.0f, 24.0f));
    M1 = glm::scale(M1, glm::vec3(2.0f, 2.0f, 2.0f));
    M1 = glm::rotate(M1,3.14f/2, glm::vec3(0.0f, 1.0f, 0.0f));

    glUniformMatrix4fv(spColored->u("M"),1,false,glm::value_ptr(M1));
    glDrawArrays( GL_TRIANGLES, 0, myStairsVertexCount);

    glDisableVertexAttribArray(spColored->a("vertex"));
    glDisableVertexAttribArray(spColored->a("color"));
    */


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
    }
    float StairsTranslations[3][3] = {
    {-4,0,5},
    {-6,0,17},
    {-16,0,17}};
    for(int k = 0; k < 3; k++){
        for(int j = 0; j < 24; j++)
        {
            triangles[triangleindex][0] = myStairsVertices[24*j]+StairsTranslations[k][0];
            triangles[triangleindex][4] = myStairsVertices[24*j+4]+StairsTranslations[k][0];
            triangles[triangleindex][8] = myStairsVertices[24*j+8]+StairsTranslations[k][0];
            triangles[triangleindex][1] = myStairsVertices[24*j+1]+StairsTranslations[k][1];
            triangles[triangleindex][5] = myStairsVertices[24*j+5]+StairsTranslations[k][1];
            triangles[triangleindex][9] = myStairsVertices[24*j+9]+StairsTranslations[k][1];
            triangles[triangleindex][2] = myStairsVertices[24*j+2]+StairsTranslations[k][2];
            triangles[triangleindex][6] = myStairsVertices[24*j+6]+StairsTranslations[k][2];
            triangles[triangleindex][10] = myStairsVertices[24*j+10]+StairsTranslations[k][2];
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
