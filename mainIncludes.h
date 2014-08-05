#ifndef MAININCLUDES_H
#define MAININCLUDES_H

#include "kfusion.h"
#include "helpers.h"
#include "interface.h"
#include "perfstats.h"
#include "track.h"
#include "gridHelpers.h"

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
#include <queue>

#include <thrust/host_vector.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#elif defined(WIN32)
#define GLUT_NO_LIB_PRAGMA
#include <glut.h>
#else
#include <GL/glut.h>
#endif

#endif // MAININCLUDES_H