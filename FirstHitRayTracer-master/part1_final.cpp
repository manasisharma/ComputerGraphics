#include <GL/glut.h>
#include<iostream>

struct rgb_values
{
    float r; float g; float b;
};

void display()
{
 
    int i,j;
    float m,n,o;
    rgb_values pixels[512][512];
    for(i=0;i<512;i++)
    {
        if(i>=1&&i<100)
        { m=1.0; n=0.0; o=0.0; }
        else if(i>=100&&i<200)
        {   n=1.0; m=0.0; o=0.0; }
        else if(i>=200&&i<300)
        {   o=1.0; m=0.0; n=0.0; }
        else if(i>300&&i<400)
        { m=1.0; n=1.0; o=0.0; }
        else
        { n=0.0; m=1.0; o=1.0; }
        for(j=0;j<512;j++)
        {
            pixels[i][j].r=m;
            pixels[i][j].g=n;
            pixels[i][j].b=o;
        }
    }
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(512,512,GL_RGB,GL_FLOAT,pixels);
    glutSwapBuffers();
}

int main(int argc,char**argv)
{
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
    glutInitWindowSize(512,512);
    glutCreateWindow("glut1");
    glutDisplayFunc( display );
    glutMainLoop();
    return 0;
}
