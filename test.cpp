/*
Copyright (c) 2011-2013 Gerhard Reitmayr, TU Graz

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "kfusion.h"
#include "helpers.h"

#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef __APPLE__
#include <GLUT/glut.h>
#elif defined(WIN32)
#define GLUT_NO_LIB_PRAGMA
#include <glut.h>
#else
#include <GL/glut.h>
#endif

#include "perfstats.h"

using namespace std;
using namespace TooN;

SE3<float> preTrans(makeVector(0.0, 0, -0.9, 0, 0, 0));
SE3<float> rot(makeVector(0.0, 0, 0, 0, 0, 0));
SE3<float> trans(makeVector(0.5, 0.5, 0.5, 0, 0, 0));

KFusion kfusion;
Volume reference;
Image<float3, HostDevice> vertex, normal;
Image<float, HostDevice> depth;
Image<uchar4, HostDevice> rgb;


int2 outputSize = make_int2(1280, 800);
int counter = 0;
bool benchmark = false;

void display(void) {

    static bool integrate = true;

    const uint2 imageSize = kfusion.configuration.inputSize;

    const double start = Stats.start();
    renderInput(vertex.getDeviceImage(), normal.getDeviceImage(), depth.getDeviceImage(), reference, toMatrix4( trans * rot * preTrans ) * getInverseCameraMatrix(kfusion.configuration.camera), kfusion.configuration.nearPlane, kfusion.configuration.farPlane, kfusion.configuration.stepSize(), 0.01, outputSize);
    cudaDeviceSynchronize();
    Stats.sample("ground raycast");
    Stats.sample("ground copy");

    glRasterPos2i(0,0);
    glDrawPixels(vertex);
    glRasterPos2i(imageSize.x, 0);
    glDrawPixels(normal);
    glRasterPos2i(imageSize.x * 2, 0);
    glDrawPixels(depth);
    Stats.sample("ground draw");

    kfusion.setDepth( depth.getDeviceImage() );
    cudaDeviceSynchronize();
    const double track_start = Stats.sample("process depth");

    if(counter > 1){
        integrate = kfusion.Track();
        cudaDeviceSynchronize();
        Stats.sample("track");
    }

    renderTrackResult(rgb.getDeviceImage(), kfusion.reduction);
    cudaDeviceSynchronize();
    Stats.sample("track render");
    Stats.sample("track copy");

    if(integrate){
        kfusion.Integrate();
        cudaDeviceSynchronize();
        Stats.sample("integration");
        kfusion.Raycast();
        cudaDeviceSynchronize();
        Stats.sample("raycast");
        vertex = kfusion.vertex;
        normal = kfusion.normal;
        Stats.sample("raycast get");
    }

    glRasterPos2i(0,imageSize.y * 1);
    glDrawPixels(vertex);
    glRasterPos2i(imageSize.x, imageSize.y * 1);
    glDrawPixels(normal);
    glRasterPos2i(2 * imageSize.x, imageSize.y * 1);
    glDrawPixels(rgb);
    Stats.sample("track draw");

    Stats.sample("total track", Stats.get_time() - track_start, PerfStats::TIME);

    renderInput(vertex.getDeviceImage(), normal.getDeviceImage(), depth.getDeviceImage(), kfusion.integration,  kfusion.pose * getInverseCameraMatrix(kfusion.configuration.camera), kfusion.configuration.nearPlane, kfusion.configuration.farPlane, kfusion.configuration.stepSize(), 0.7 * kfusion.configuration.mu, outputSize);
    cudaDeviceSynchronize();
    Stats.sample("view raycast");
    Stats.sample("view copy");

    glRasterPos2i(0,imageSize.y * 2);
    glDrawPixels(vertex);
    glRasterPos2i(imageSize.x, imageSize.y * 2);
    glDrawPixels(normal);
    glRasterPos2i(imageSize.x * 2, imageSize.y * 2);
    glDrawPixels(depth);
    Stats.sample("view draw");

    Stats.sample("events");
    Stats.sample("total all", Stats.get_time() - start, PerfStats::TIME);

    if(counter % 30 == 0){
        Stats.print();
        Stats.reset();
        cout << endl;
    }

    ++counter;

    printCUDAError();

    glutSwapBuffers();
}

void keys(unsigned char key, int x, int y) {
    switch(key){
    case 'r':
        kfusion.setPose( toMatrix4( trans * rot * preTrans ));
        break;
    case 'c':
        kfusion.Reset();
        kfusion.setPose( toMatrix4( trans * rot * preTrans ));
        break;
    case 'd':
        cout << kfusion.pose << endl;
        break;
    case 'q':
        exit(0);
        break;
    }
    glutPostRedisplay();
}

void specials(int key, int x, int y){
    switch(key){
    case GLUT_KEY_LEFT:
        rot *= SE3<float>(makeVector(0.0, 0, 0, 0, 0.1, 0));
        break;
    case GLUT_KEY_RIGHT:
        rot *= SE3<float>(makeVector(0.0, 0, 0, 0, -0.1, 0));
        break;
    case GLUT_KEY_UP:
        rot *= SE3<float>(makeVector(0.0, 0, 0, -0.1, 0, 0));
        break;
    case GLUT_KEY_DOWN:
        rot *= SE3<float>(makeVector(0.0, 0, 0, 0.1, 0, 0));
        break;
    }
    glutPostRedisplay();
}

void idle(void) {
    if(counter > 100 && benchmark)
        exit(0);

    if(benchmark)
        rot *= SE3<float>(makeVector(0.0, 0, 0, 0, 0.02, 0));

    glutPostRedisplay();
}

void reshape(int width, int height){
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glColor3f(1.0f,1.0f,1.0f);
    glRasterPos2f(-1, 1);
    glOrtho(-0.375, width-0.375, height-0.375, -0.375, -1 , 1); //offsets to make (0,0) the top left pixel (rather than off the display)
    glPixelZoom(1,-1);
}

int main(int argc, char ** argv) {

    benchmark = argc > 1 && string(argv[1]) == "-b";

    KFusionConfig config;
    config.volumeSize = make_uint3(128);

    config.combinedTrackAndReduce = false;

    config.iterations[0] = 10;
    config.iterations[1] = 5;
    config.iterations[2] = 5;

    config.inputSize = make_uint2(320, 240);
    config.camera = make_float4(100, 100, 160, 120);
    config.nearPlane = 0.001;

    config.maxweight = 100;
    config.mu = 0.1;

    config.dist_threshold = 0.2f;
    config.normal_threshold = 0.8f;

    kfusion.Init(config);
    if(printCUDAError()){
        cudaDeviceReset();
        exit(1);
    }

    reference.init(config.volumeSize, config.volumeDimensions);

    initVolumeWrap(reference, 1.0f);
    setBoxWrap(reference, make_float3(0.1f,0.1f,0.8f), make_float3(0.9f, 0.9f, 0.9f), -1.0f);
    setBoxWrap(reference, make_float3(0.1f,0.8f,0.1f), make_float3(0.9f, 0.9f, 0.9f), -1.0f);
    setBoxWrap(reference, make_float3(0.8f,0.1f,0.1f), make_float3(0.9f, 0.9f, 0.9f), -1.0f);
    setSphereWrap(reference, make_float3(0.5f), 0.2f, -1.0f);

    kfusion.setPose( toMatrix4( trans * rot * preTrans ));

    vertex.alloc(config.inputSize);
    normal.alloc(config.inputSize);
    depth.alloc(config.inputSize);
    rgb.alloc(config.inputSize);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE );
    glutInitWindowSize(config.inputSize.x * 3, config.inputSize.y * 3);
    glutCreateWindow("kfusion test");

    glutDisplayFunc(display);
    glutKeyboardFunc(keys);
    glutSpecialFunc(specials);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);

    glutMainLoop();

    cudaDeviceReset();

    return 0;
}
