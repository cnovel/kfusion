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
#include "interface.h"
#include "perfstats.h"

// OVR include
#include "LibOVR/Src/OVR_CAPI.h"
#include "LibOVR/Src/OVR_SensorFusion.h"
#include "LibOVR/Src/Kernel/OVR_Math.h"
#include "LibOVR/Include/OVR.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <memory>

#ifdef __APPLE__
#include <GLUT/glut.h>
#elif defined(WIN32)
#define GLUT_NO_LIB_PRAGMA
#include <glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;
using namespace TooN;

KFusion kfusion;
Image<uchar4, HostDevice> lightScene, trackModel, lightModel, viewLeft, viewRight, barrel;
Image<uint16_t, HostDevice> depthImage[2];
Image<uchar3, HostDevice> rgbImage;

const float3 light = make_float3(1, 1, -1.0);
const float3 lightbis = make_float3(-1, -1, 1.0);
const float3 ambient = make_float3(0.1, 0.1, 0.1);
const float3 ambientbis = make_float3(1, 1, 1);

SE3<float> initPose;

int counter = 0;
int integration_rate = 2;

int2 outputSize = make_int2(1280, 800);
bool reset = true;
bool should_integrate = true;
bool render_texture = false;
bool ovr_sensor_tracking = false;
bool ovr_mode = false;

Image<float3, Device> posRight, normalsRight, posLeft, normalsLeft;
Image<float, Device> depRight, depLeft;

SE3<float> preTrans, trans, rot(makeVector(0.0, 0, 0, 0, 0, 0));

// OVR Object
ovrHmd hmd;
ovrHmdDesc hmdDesc;
float offset = 0.0f; //0.0325f

void display(void){

    // OVR Sensor set up
    ovrSensorState ss = ovrHmd_GetSensorState(hmd, 0.0);
    float yYaw = 0;
    float xEyePitch = 0;
    float zEyeRoll = 0;

    if (ss.StatusFlags & (ovrStatus_OrientationTracked)) {
        // Obtain the orientation of the HMD
        OVR::Transformf pose = ss.Recorded.Pose;
        pose.Rotation.GetEulerAngles<OVR::Axis_Y, OVR::Axis_X, OVR::Axis_Z>(&yYaw, &xEyePitch, &zEyeRoll);
        yYaw = -yYaw;
        zEyeRoll = -zEyeRoll;
    }

    static bool integrate = true;

    glClear( GL_COLOR_BUFFER_BIT );
    const double startFrame = Stats.start();
    const double startProcessing = Stats.sample("kinect");

    kfusion.setKinectDeviceDepth(depthImage[GetKinectFrame()].getDeviceImage());

    integrate = kfusion.Track();

    if((should_integrate && integrate && ((counter % integration_rate) == 0)) || reset){
        kfusion.Integrate();
        kfusion.Raycast();
        Stats.sample("integrate");
        if(counter > 2) // use the first two frames to initialize
            reset = false;
    }

    Matrix4 cameraView = getInverseCameraMatrix(kfusion.configuration.camera * 2);
    Matrix4 ovrPose = kfusion.pose;
    
    if (ovr_sensor_tracking)
        viewMatrixUpdate(ovrPose, yYaw, zEyeRoll, xEyePitch);

    Matrix4 leftEye = cameraView;
    Matrix4 rightEye = cameraView;

    leftEye.data[0].w -= offset;
    rightEye.data[0].w += offset;

    float3 lux = make_float3(ovrPose.data[0].w, ovrPose.data[1].w, ovrPose.data[2].w);

    if (!ovr_mode) {
        renderInput( posLeft, normalsLeft, depLeft, kfusion.integration, ovrPose * leftEye, kfusion.configuration.nearPlane, kfusion.configuration.farPlane, kfusion.configuration.stepSize(), 0.75 * kfusion.configuration.mu, outputSize);
        renderInput( posRight, normalsRight, depRight, kfusion.integration, ovrPose * rightEye, kfusion.configuration.nearPlane, kfusion.configuration.farPlane, kfusion.configuration.stepSize(), 0.75 * kfusion.configuration.mu, outputSize);

        if(render_texture) {
            renderTexture(viewLeft.getDeviceImage(), posLeft, normalsLeft, rgbImage.getDeviceImage(), getCameraMatrix(2*kfusion.configuration.camera) * inverse(kfusion.pose), lux);
            renderTexture(viewRight.getDeviceImage(), posRight, normalsRight, rgbImage.getDeviceImage(), getCameraMatrix(2*kfusion.configuration.camera) * inverse(kfusion.pose), light);
        } else {
            renderLight(viewLeft.getDeviceImage(), posLeft, normalsLeft, lux, ambient);
            renderLight(viewRight.getDeviceImage(), posRight, normalsRight, lux, ambient);
        }
    } else {
        renderOculusCam(barrel.getDeviceImage(), kfusion.integration, rgbImage.getDeviceImage(), ovrPose, getCameraMatrix(2*kfusion.configuration.camera) * inverse(kfusion.pose), kfusion.configuration.nearPlane, kfusion.configuration.farPlane, kfusion.configuration.stepSize(), 0.75 * kfusion.configuration.mu, lux, ambient);
    }

    cudaDeviceSynchronize();

    glClear(GL_COLOR_BUFFER_BIT);
    glRasterPos2i(0, 0);
    if (ovr_mode){
        glDrawPixels(barrel);
    } else {
        glDrawPixels(viewLeft);
        //glRasterPos2i(640, 0);
        //glDrawPixels(viewRight);
    }

    if(printCUDAError()) {
        exit(1);
    }

    ++counter;

    glutSwapBuffers();    
}

void idle(void){
    if(KinectFrameAvailable())
        glutPostRedisplay();
}

void keys(unsigned char key, int x, int y){
    switch(key){
    case 'c':
        kfusion.Reset();
        kfusion.setPose(toMatrix4(initPose));
        reset = true;
        break;
    case 'q':
        exit(0);
        break;
    case 'i':
        should_integrate = !should_integrate;
        break;
    case 't':
        render_texture = !render_texture;
        break;
    case 'o':
        ovr_sensor_tracking = !ovr_sensor_tracking;
        break;
    case 'r':
        ovr_mode = !ovr_mode;
        break;//*
    case '+':
        offset += .01f;
        break;
    case '-':
        if (offset > .01f)
            offset -= .02f;
        break;//*/
    }

}

void specials(int key, int x, int y){
    switch(key){
    case GLUT_KEY_LEFT:
        rot = SE3<float>(makeVector(0.0, 0, 0, 0, 0.1, 0)) * rot;
        break;
    case GLUT_KEY_RIGHT:
        rot = SE3<float>(makeVector(0.0, 0, 0, 0, -0.1, 0)) * rot;
        break;
    case GLUT_KEY_UP:
        rot *= SE3<float>(makeVector(0.0, 0, 0, -0.1, 0, 0));
        break;
    case GLUT_KEY_DOWN:
        rot *= SE3<float>(makeVector(0.0, 0, 0, 0.1, 0, 0));
        break;
    }
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

void exitFunc(void){
    CloseKinect();
    kfusion.Clear();
    cudaDeviceReset();
}

int main(int argc, char ** argv) {
    const float size = (argc > 1) ? atof(argv[1]) : 5.f;
    int sizeX = outputSize.x;
    int sizeY = outputSize.y;

    KFusionConfig config;

    // it is enough now to set the volume resolution once.
    // everything else is derived from that.
    // config.volumeSize = make_uint3(64);
    // config.volumeSize = make_uint3(128);
    config.volumeSize = make_uint3(640);

    // these are physical dimensions in meters
    config.volumeDimensions = make_float3(size);
    config.nearPlane = 0.2f;
    config.farPlane = 10.0f;
    config.mu = 0.1;
    config.combinedTrackAndReduce = false;

    // change the following parameters for using 640 x 480 input images
    config.inputSize = make_uint2(320,240);
    config.camera =  make_float4(531.15/2, 531.15/2, 160, 120);

    // config.iterations is a vector<int>, the length determines
    // the number of levels to be used in tracking
    // push back more then 3 iteraton numbers to get more levels.
    config.iterations[0] = 10;
    config.iterations[1] = 5;
    config.iterations[2] = 4;

    config.dist_threshold = (argc > 2 ) ? atof(argv[2]) : config.dist_threshold;
    config.normal_threshold = (argc > 3 ) ? atof(argv[3]) : config.normal_threshold;

    initPose = SE3<float>(makeVector(size/2, size/2, 0, 0, 0, 0));

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE );
    glutInitWindowSize(1280, 800);
    glutCreateWindow("OFusion");

    kfusion.Init(config);

    // input buffers
    depthImage[0].alloc(make_uint2(640, 480));
    depthImage[1].alloc(make_uint2(640, 480));
    rgbImage.alloc(make_uint2(640, 480));

    // render buffers
    lightScene.alloc(config.inputSize), trackModel.alloc(config.inputSize), lightModel.alloc(config.inputSize);
    viewLeft.alloc(make_uint2(sizeX, sizeY)), viewRight.alloc(make_uint2(sizeX, sizeY));
    barrel.alloc(make_uint2(1280, 800));

    posLeft.alloc(make_uint2(1280, 800)), normalsLeft.alloc(make_uint2(1280, 800)), depLeft.alloc(make_uint2(640, 800));
    posRight.alloc(make_uint2(1280, 800)), normalsRight.alloc(make_uint2(1280, 800)), depRight.alloc(make_uint2(640, 800));


    if(printCUDAError()) {
        cudaDeviceReset();
        return 1;
    }

    memset(depthImage[0].data(), 0, depthImage[0].size.x*depthImage[0].size.y * sizeof(uint16_t));
    memset(depthImage[1].data(), 0, depthImage[1].size.x*depthImage[1].size.y * sizeof(uint16_t));
    memset(rgbImage.data(), 0, rgbImage.size.x*rgbImage.size.y * sizeof(uchar3));

    uint16_t * buffers[2] = {depthImage[0].data(), depthImage[1].data()};
    if(InitKinect(buffers, (unsigned char *)rgbImage.data())){
        cudaDeviceReset();
        return 1;
    }

    kfusion.setPose(toMatrix4(initPose));

    // model rendering parameters
    preTrans = SE3<float>::exp(makeVector(0.0, 0, -size, 0, 0, 0));
    trans = SE3<float>::exp(makeVector(0.5, 0.5, 0.5, 0, 0, 0) * size);


    // OVR Part
    // Initializes LibOVR

    ovr_Initialize();
    
    hmd = ovrHmd_Create(0);

    if (hmd) {
        ovrHmd_GetDesc(hmd, &hmdDesc);
    }

    // Start the sensors (don't need position)
    ovrHmd_StartSensor(hmd, ovrSensorCap_Orientation | ovrSensorCap_YawCorrection, ovrSensorCap_Orientation);
    //*/

    atexit(exitFunc);
    glutDisplayFunc(display);
    glutKeyboardFunc(keys);
    glutSpecialFunc(specials);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);

    glutMainLoop();

    // Close Kinect & Rift
    CloseKinect();
    ovrHmd_Destroy(hmd);
    ovr_Shutdown();

    return 0;
}
