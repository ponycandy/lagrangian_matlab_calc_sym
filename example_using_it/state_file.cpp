//
// File: state_file.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 18-Apr-2023 16:23:32
//

// Include Files
#include "state_file.h"
#include <math.h>

// Function Definitions
//
// STATE_FILE
//     D_STATE = STATE_FILE(DTHETA,DX,F,G,L,M1,M2,THETA)
//
// Arguments    : double Dtheta
//                double Dx
//                double F
//                double g
//                double l
//                double m1
//                double m2
//                double theta
//                double d_state[4]
// Return Type  : void
//
void state_file(double Dtheta, double Dx, double F, double g, double l,
                double m1, double m2, double theta, double d_state[4])
{
  double b_d_state_tmp;
  double d_state_tmp;
  double t2;
  double t3;
  double t4;
  double t9;
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     18-Apr-2023 16:06:41
  t2 = cos(theta);
  t3 = sin(theta);
  t4 = Dtheta * Dtheta;
  t9 = 1.0 / ((m1 + m2) + -(m2 * (t2 * t2)));
  d_state[0] = Dtheta;
  d_state_tmp = g * m2;
  b_d_state_tmp = l * m2;
  d_state[1] = -(t9 * (((F * t2 + g * m1 * t3) + d_state_tmp * t3) +
                       b_d_state_tmp * t2 * t3 * t4)) /
               l;
  d_state[2] = Dx;
  d_state[3] = t9 * ((F + d_state_tmp * t2 * t3) + b_d_state_tmp * t3 * t4);
}

//
// File trailer for state_file.cpp
//
// [EOF]
//
