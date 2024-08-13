/*
REMEBER: In this setup the observer is at the origin (0,0,0).


// 2D ASCII projection
Dimmest to brightest .,-~:;=!*#$@
Rather than render the whole surface, we only render points on the surface
and since the ASCII chars are so large if looks as if the surface is smooth.
Distence from the observer to view port z'
Distence from the observer to object z
points on the 3D object (x,y,z)
points on the 2D plane of projection (x',y')
points on the 3D plane centered at the 2D plane of projection (x',y',z')
Projected x coordinate: x' = x * z'/z
Projected y coordinate: y' = y * z'/z
Coordinates on the projection (x',y') = (x * z'/z, y * z'/z)
The article defines K1 := z'
z' is an arbitrary constant based on the size of the view port
If z' is too large or small the object will be too close or too far to see.

Its useful to calculate 1/z since 1/z = 0 corresponds to infinite depth.
1. If you init our z-buffer to 0 everything is infinitly far way and therefore 
    invisible.
2. we can reuse 1/z when calculating x' and y'.

// Torus
A torus is a solid of revolution so draw a circle and sweep it in a circle (theta from 0 to 2pi).
Let r1 be the raduis of the tube and r2 be the radius of the sweeping circle.

Radius and Angle corresponding to the tube: r1, theta
Radius and Angle corresponding to the donut: r2, phi


Forming all points on a circle:
We use the equation for a circle in terms or an angle theta.
  (x,y,z)
    = (r2,0,0) + (r1*cos(theta), r1*sin(theta), 0) 
    = (r2 + r1*cos(theta), r1*sin(theta), 0)

Rotating the circle to form the Torus T. 
  In 3D there are 3 rotation matrices, Rx, Ry, and Rz, one for each axis.
    Rx rotates about the x-axis with angle A
    Ry rotates about the y-axis with angle B
    Rz rotates about the z-axis with angle C
  These rotations can be composed
  T = (r2 + r1*cos(theta), r1*sin(theta), 0) * [Rotation Matrix(phi)]

Rotating the Torus (for fun!)
  T * [one or more rotation matrices] 
  = ((r2 + r1*cos(theta), r1*sin(theta), 0) * [Rotation Matrix(phi)]) * [one or more rotation matrices]

// Projecting
This step actually has two parts
1. move the torus in from of the observer.
2. project from 3D to 2D.

Move the torus backward with some constant K2, giving
  (x', y') = ( x * z'/(K2+z) , y * z'/(K2+z) )

// Illumination
To find the illumination of a given point take the dot product of surface normal and the light direction.
Let N be the surface normal, i.e. the direction perpendicular to the surface. We use the standard unit vector
in terms of our angle theta and multiple it by the same rotations we applied to the Torus. WE can use theta 
because we can find the normal for the circle then rotate those normals around with the circle to find the rest.
  N = (Nx, Ny, Nz) = (cos(mu), sin(mu), 0) * [one or more rotation matrices]

Since the observer is at the origin let's have the light source be behind and above they
  L = (0,1,-1)

Then take the dot prodoct of the two 
  N * L

*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#define PI 3.14159265358979323846

#define R_ONE 1
#define R_TWO 2
#define K_TWO 5
#define SCREEN_WIDTH 160
#define SCREEN_HEIGHT 25
#define SCREEN_DIM SCREEN_WIDTH*SCREEN_HEIGHT



const float theta_resolution = 0.07f;
const float phi_resolution = 0.02f;
const float R1 = R_ONE;
const float R2 = R_TWO;
const float K2 = K_TWO;


/* Calculate K1 base on screen size: the max x-distence occurse roughly at the edge of the torus 
* which is x=R1+R2, z=0. We want to be displaced 3/8th of the width of the screen, which is 3/4ths 
* of the way from the center to the side of the screen */
// SCREEN_WIDTH * 3/8 = K1(R1+R2)/(K2+0)
const float K1 = SCREEN_WIDTH * K_TWO * 3 / (8 * (R_ONE + R_TWO));

void render_frame(float A, float B) {
  float cosA = cos(A);
  float cosB = cos(B);
  float sinA = sin(A);
  float sinB = sin(B);

  char output[SCREEN_DIM];
  float zbuffer[SCREEN_DIM];
  memset(output, ' ', SCREEN_DIM);
  memset(zbuffer, 0, SCREEN_DIM);

  // theta goes around the cross-sectional circle of a torus
  for (float theta = 0.0f; theta < 2.0f*PI; theta += theta_resolution) {
    float costheta = cos(theta);
    float sintheta = sin(theta);

    for (float phi = 0.0f; phi < 2.0f*PI; phi += phi_resolution) {
      float cosphi = cos(phi);
      float sinphi = sin(phi);

      // the x,y coordinate of the circle, before revolving (factored out of 
      // the above equations)
      float circlex = R2 + R1*costheta;
      float circley = R1*sintheta;

      // final 3D (x,y,z) coordinate after rotations, directly from our math above
      float x = circlex*(cosB*cosphi + sinA*sinB*sinphi) - circley*cosA*cosB;
      float y = circlex*(sinB*cosphi - sinA*cosB*sinphi) + circley*cosA*cosB;
      float z = K2 + cosA*circlex*sinphi + circley*sinA;
      float z_inverse = 1.0f / z;

      // x and y projection. note that y is negate here bcause y goes up in 3D space but down in 2D displays.
      int xp = (int) (SCREEN_WIDTH/2 + K1*z_inverse*x);
      int yp = (int) (SCREEN_HEIGHT/2 - K1*z_inverse*y);

      // calculate luminance. ugly, but correct.
      float luminance = cosphi*costheta*sinB - cosA*costheta*sinphi - sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);
      // luminance ranges from -sqrt(2) to +sqrt(2). If it's < 0, the surface is pointing away from us, so we won't bother trying to plot it.
      if (luminance > 0) {
        // test against z-buffer. A larger z_inverse means the pixel is closer to the viewer then what's already plotted.
        if(z_inverse > zbuffer[xp*yp]) {
          zbuffer[xp*yp] = z_inverse; 
          int luminance_index = 8.0f*luminance;
          // luminance_index is now in range 0 through 11 (8*sqrt(2) = 11.3) now we lookup the char 
          // corresponding to the luminance and plot it in our output:
          output[xp * yp] = ".,-~:;=!*#$@"[luminance_index];
        }
      }
    }
  }

  // now, dump output[] to the screen.
  // bring the cursor to "home" location, in just about any currently-used terminal emulation mode.
  printf("\x1b[H");
  for (int k = 0; k < SCREEN_HEIGHT * SCREEN_WIDTH; k++) {
    putchar(k % SCREEN_WIDTH ? output[k] : 10);
  }
}


int main() {
  render_frame(0.5f, 0.5f);
  return 0;
}
