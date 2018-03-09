/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

MGLpoly_mode poly_mode;
MGLmatrix_mode matrix_mode;

vec3 color;

mat4 Identity = {{1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1}};

vector<mat4> projectionStack = {Identity};
vector<mat4> modelviewStack = {Identity};

vector<vector<MGLfloat>> min_zBulk;

// Stack Helper Functions
mat4& topOfActiveMatrixStack()
{
    if(matrix_mode == MGL_PROJECTION)
        return projectionStack.back();
    else
        return modelviewStack.back();    
}

struct Vertex
{
    vec4 position;
    vec3 col;

    Vertex() : position(0, 0, 0, 0), col(0, 0, 0) { }
    Vertex(MGLfloat x, MGLfloat y, MGLfloat z, MGLfloat w, vec3 c) 
          : position(x, y, z, w), col(c) { }
};

vector<Vertex> vertices;

struct Triangle
{
    Vertex A, B, C;
    
    Triangle(Vertex a, Vertex b, Vertex c) : A(a), B(b), C(c) { }
};

vector<Triangle> triangles;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}



MGLfloat calcArea(vec2 a, vec2 b, vec2 c)
{
    return a[0] * (b[1] - c[1]) + a[1] * (c[0] - b[0]) + (b[0] * c[1] - b[1] * c[0]);
}

Triangle findNPC_Coordinates(const Triangle& tri)
{
    Triangle t = tri;

    t.A.position[0] = tri.A.position[0] / tri.A.position[3];
    t.A.position[1] = tri.A.position[1] / tri.A.position[3];
    t.A.position[2] = tri.A.position[2] / tri.A.position[3];

    t.B.position[0] = tri.B.position[0] / tri.B.position[3];
    t.B.position[1] = tri.B.position[1] / tri.B.position[3];
    t.B.position[2] = tri.B.position[2] / tri.B.position[3];

    t.C.position[0] = tri.C.position[0] / tri.C.position[3];
    t.C.position[1] = tri.C.position[1] / tri.C.position[3];
    t.C.position[2] = tri.C.position[2] / tri.C.position[3];

    return t;
}

bool isInClippingBox(const Triangle& t)
{
    return 
           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3])) &&
           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3])) &&

           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3])) &&
           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3])) &&
 
           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3])) &&
           ((t.A.position[3] >= 0 && t.A.position[0] <= t.A.position[3]) ||
           (t.A.position[3] <= 0 && t.A.position[0] >= t.A.position[3]));
}
// Calculates perspective correct color interpolation, based on input, 
// index, which is alpha, beta or gamma
vector<MGLfloat> getColors(MGLfloat alpha, MGLfloat beta,  MGLfloat gamma, Triangle t)
{
    vector<MGLfloat> colors;

    MGLfloat div = (alpha / t.A.position[3]) 
                    + (beta / t.B.position[3]) 
                    + (gamma / t.C.position[3]);
    
    colors.push_back((alpha / t.A.position[3]) / div);
    colors.push_back((beta / t.B.position[3]) / div);
    colors.push_back((gamma / t.C.position[3]) / div);
    
    return colors;
}

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data)
{
    Triangle t = findNPC_Coordinates(tri);
    
    vec2 a, b, c;
    a[0] = ((t.A.position[0] + 1) * width / 2) - .5; 
    a[1] = ((t.A.position[1] + 1) * height / 2) - .5; 

    b[0] = ((t.B.position[0] + 1) * width / 2) - .5; 
    b[1] = ((t.B.position[1] + 1) * height / 2) - .5; 

    c[0] = ((t.C.position[0] + 1) * width / 2) - .5; 
    c[1] = ((t.C.position[1] + 1) * height / 2) - .5; 

    MGLfloat div, alphaInter, betaInter, gammaInter, red, green, blue;

    // Bounding box, to make test 25 faster
    int xMin = min(a[0], b[0]), 
        xMax = max(a[0], b[0]), 
        yMin = min(a[1], b[1]), 
        yMax = max(a[1], b[1]);
    xMin = min(xMin, static_cast<int>(c[0]));
    xMax = max(xMax, static_cast<int>(c[0]));
    yMin = min(yMin, static_cast<int>(c[1]));
    yMax = max(yMax, static_cast<int>(c[1]));
      
    //std::cout << xMin << " " << xMax << " " << yMin << " " << yMax << std::endl;
  //  std::cout << a[0] << " " << a[1] << " " << b[0] << " " << b[1] << " " << c[0] << " " << c[1] << std::endl;
    
    int iMin = max(xMin, 0), iMax = min(width - 1, xMax), 
        jMin = max(yMin, 0), jMax = min(height - 1, yMax);

    //std::cout << iMin << " " << iMax << " " << jMin << " " << jMax << std::endl;
    for(int i = iMin; i <= iMax; ++i)
    {
        for(int j = jMin; j <= jMax; ++j)
        {
            vec2 p = {i, j};
            MGLfloat area = calcArea(a, b, c);
            MGLfloat alpha = calcArea(p, b, c) / area;
            MGLfloat beta = calcArea(a, p, c) / area;
            MGLfloat gamma = calcArea(a, b, p) / area;

           // if(isInClippingBox(tri))
            {
            if(alpha >= 0 && beta >= 0 && gamma >= 0)
            {
                MGLfloat z_depth = alpha * t.A.position[2] + beta * t.B.position[2]
                                   + gamma * t.C.position[2];
                if(z_depth <= 1 && z_depth >= -1)
                {
                    if(z_depth < min_zBulk[i][j])
                    {
                        div = (alpha / t.A.position[3]) 
                                        + (beta / t.B.position[3]) 
                                        + (gamma / t.C.position[3]);
    
                        alphaInter = ((alpha / t.A.position[3]) / div);
                        betaInter = ((beta / t.B.position[3]) / div);
                        gammaInter = ((gamma / t.C.position[3]) / div);
                       
                        red = t.A.col[0] * alphaInter 
                                     + t.B.col[0] * betaInter 
                                     + t.C.col[0] * gammaInter;
                        green = t.A.col[1] * alphaInter 
                                     + t.B.col[1] * betaInter 
                                     + t.C.col[1] * gammaInter;
                        blue = t.A.col[2] * alphaInter 
                                     + t.B.col[2] * betaInter 
                                     + t.C.col[2] * gammaInter;

                        // Fix for perspective correct color interpolation
            //            if(red != 0 && green != 0 && blue != 0)
                            data[i+j*width] = Make_Pixel(255 * red, 255 * green, 255 * blue);
                        min_zBulk[i][j] = z_depth;
                    }
                }
            }
            }
        }
    }
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    Make_Pixel(0, 0, 0);
    
    // Allocate size
    min_zBulk.resize(width);
    for(size_t i = 0; i < width; i++)
        min_zBulk[i].resize(height);

    // Initialize
    for(size_t i = 0; i < width; i++)
        for(size_t j = 0; j < height; j++)
            min_zBulk[i][j] = 2;

    for(size_t i = 0; i < triangles.size(); i++)
       Rasterize_Triangle(triangles[i], width, height, data);

    triangles.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    poly_mode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    if(poly_mode == MGL_TRIANGLES)
    {
        for(size_t i = 0; i < vertices.size(); i += 3)
        {
            Triangle t(vertices[i], vertices[i+1], vertices[i+2]);
            triangles.push_back(t);
        }
    }
    else
    {
        for(size_t i = 0; i < vertices.size(); i += 4)
        {
            Triangle t1(vertices[i], vertices[i+1], vertices[i+2]);
            triangles.push_back(t1);

            Triangle t2(vertices[i], vertices[i+2], vertices[i+3]);
            triangles.push_back(t2);
        }
    }
    vertices.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x, y, 0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 pos = { x, y, z, 1};
    pos = projectionStack.back() * modelviewStack.back() * pos; 
    
    Vertex v(0, 0, 0, 0, color);
    v.position = pos;
    vertices.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    matrix_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    if(matrix_mode == MGL_MODELVIEW)
        modelviewStack.push_back(modelviewStack.back());
    else
        projectionStack.push_back(projectionStack.back());
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    if(matrix_mode == MGL_MODELVIEW)
        modelviewStack.pop_back();
    else
        projectionStack.pop_back();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    mat4& currentM = topOfActiveMatrixStack();
    currentM = Identity;
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    mat4 tempM;

    for(int i = 0; i < 16; i++)
       tempM.values[i] = matrix[i];

    mat4& currentM = topOfActiveMatrixStack();
    currentM = tempM;    
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mat4 tempM;

    for(int i = 0; i < 16; i++)
       tempM.values[i] = matrix[i];

    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * tempM;    
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4 temp_matrix = {{1, 0, 0, 0,
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                         x, y, z, 1}};

    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * temp_matrix; 
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    vec3 v = {x, y, z};
    v = v.normalized();

    double radAngle = (angle * M_PI) / 180;
    MGLfloat c = cos(radAngle), s = sin(radAngle);

    mat4 temp_matrix = {{v[0] * v[0] * (1 - c) + c, 
                         v[0] * v[1] * (1 - c) + v[2] * s,
                         v[0] * v[2] * (1 - c) - v[1] * s, 
                         0,
                         v[0] * v[1] * (1 - c) - v[2] * s,
                         v[1] * v[1] * (1 - c) + c,
                         v[1] * v[2] * (1 - c) + v[0] * s,
                         0,
                         v[0] * v[2] * (1 - c) + v[1] * s,
                         v[1] * v[2] * (1 - c) - v[0] * s,
                         v[2] * v[2] * (1 - c) + c,
                         0,
                         0, 0, 0, 1}};	 

    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * temp_matrix; 
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 temp_matrix = {{x, 0, 0, 0,
                         0, y, 0, 0,
                         0, 0, z, 0,
                         0, 0, 0, 1}};
 
    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * temp_matrix; 
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    MGLfloat A, B, C, D;
    A = (right + left) / (right - left);
    B = (top + bottom) / (top - bottom);
    C = (far + near) / (far - near);
    D = (2 * far * near) / (far - near);

    mat4 temp_matrix = {{ (2 * near) / (right - left), 0, 0, 0, 
                          0, (2 * near) / (top - bottom), 0, 0, 
                          A, B, C, -1, 
                          0, 0, D, 0}};
    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * temp_matrix; 
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    MGLfloat t_x, t_y, t_z;
    t_x = -1 * (right + left) / (right - left);
    t_y = -1 * (top + bottom) / (top - bottom);
    t_z = -1 * (far + near) / (far - near);

    mat4 temp_matrix = {{ 2 / (right - left), 0, 0, 0, 
                          0, 2 / (top - bottom), 0, 0, 
                          0, 0, -2 / (far - near), 0, 
                          t_x, t_y, t_z, 1}};
    mat4& currentM = topOfActiveMatrixStack();
    currentM = currentM * temp_matrix; 
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
   color = {red, green, blue};
}
