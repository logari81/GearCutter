// --------------------------------------------------------------------------
//
// GearCutter 0.2
// Generation of involute gear tooth profiles
// Copyright (C) 2009-2024 Konstantinos Poulios
//
// --------------------------------------------------------------------------
//
// This file is part of GearCutter.
//
// GearCutter is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// --------------------------------------------------------------------------

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <algorithm>
#include "gear_cutter.h"

using std::min;
using std::max;
const double PI = M_PI;

///////////////////////////////////////////////////////////////////////
// Structures

struct CSpos {
  double x,y,fi;
};

struct ContourPoint {
  int id_min, id_max;
  double r, theta_min, theta_max, x_min, y_min, x_max, y_max;
};


///////////////////////////////////////////////////////////////////////
// Procedures

void tool_to_part_coords(const CSpos &tool_pos, const CSpos &part_pos,
                         const double x_tool, const double y_tool,
                         double &x_part, double &y_part) {
  double cc = cos(tool_pos.fi);
  double ss = sin(tool_pos.fi);

  double x_orig = x_tool*cc - y_tool*ss + tool_pos.x - part_pos.x;
  double y_orig = x_tool*ss + y_tool*cc + tool_pos.y - part_pos.y;

  cc = cos(-part_pos.fi);
  ss = sin(-part_pos.fi);

  x_part = x_orig*cc - y_orig*ss;
  y_part = x_orig*ss + y_orig*cc;
}

double normalized_angle(double angle) {
  while (angle < 0.) angle += 2*PI;
  while (angle > 2*PI) angle -= 2*PI;
  return angle;
}

void normalize_angles(double &angle1, double &angle2) {
  while (angle1 < 0.) { angle1 += 2*PI;  angle2 += 2*PI; }
  while (angle1 > 2*PI) { angle1 -= 2*PI; angle2 -= 2*PI; }
}

double norm(const double a, const double b) {
  return sqrt(a*a+b*b);
}

void tool_params_internal_gear(const int z, const int z_0,
                               const double m_t, const double alfa_t,
                               const double xm_n, const double r_a0,
                               const double rho_a0, const double h_f1,
                               const double s_P0,
//                               double &r_w1, double &r_f1,
                               double &r_Ff1, // double &h_FaP0,
                               double &a, double &x_m, double &y_m) {
  const double cos_alfa_t = cos(alfa_t),
               sin_alfa_t = sin(alfa_t),
               tan_alfa_t = tan(alfa_t);
  const double AB = abs(z)*m_t*cos_alfa_t/2,
               DE = z_0*m_t*cos_alfa_t/2,
               AC = abs(z)*m_t/2;
//  r_w1 = AC;

  const double DF = r_a0 - rho_a0,
               EF = sqrt(DF*DF - DE*DE),
               EG = EF + rho_a0,
               AK = AC - xm_n + h_f1;
//  r_f1 = AK;

  const double AD = AK - r_a0;
  a = AD;

  const double BE = sqrt(AD*AD - (AB-DE)*(AB-DE)),
               BC = AC*sin_alfa_t,
               EC = BC - BE,
               CG = EG - EC;
  const double h_FaP0 = CG*sin_alfa_t + xm_n;
  r_Ff1 = norm(AB, BE + EG);

  const double HI = CG*cos_alfa_t + h_FaP0*tan_alfa_t + s_P0/2;
  const double theta = acos((AB-DE)/AD) - alfa_t,
               phi = acos(DE/DF) - alfa_t - theta,
               phi0 = phi - (HI/AC - theta)*abs(z) / double(z_0);
  x_m = DF*sin(phi0);
  y_m = DF*cos(phi0);
}

void cutting_process(const std::vector<Segment> &tool_contour,
                     const std::vector<CSpos> &tool_pos,
                     const std::vector<CSpos> &part_pos,
                     std::vector<ContourPoint> &pts,
                     const int mode) {

  const int n_sim = tool_pos.size(),
            n_profile = pts.size();
  for (int it=0; it < n_sim; ++it) {
    const double Dfi = tool_pos[it].fi - part_pos[it].fi;
    for (const Segment &cur_segm : tool_contour) {
      if (cur_segm.arc) {
        Arc &cur_arc = *(cur_segm.arc);
        tool_to_part_coords(tool_pos[it], part_pos[it],
                            cur_arc.x_tool_m, cur_arc.y_tool_m,
                            cur_arc.x_part_m, cur_arc.y_part_m);
        cur_arc.r_m = norm(cur_arc.x_part_m, cur_arc.y_part_m);
        cur_arc.fi_m = atan2(cur_arc.y_part_m, cur_arc.x_part_m);
        cur_arc.fi_part_1 =  cur_arc.fi_tool_1 + Dfi;
        cur_arc.fi_part_2 =  cur_arc.fi_tool_2 + Dfi;
        normalize_angles(cur_arc.fi_part_1, cur_arc.fi_part_2);
      } else {
        LineSeg &cur_line = *(cur_segm.line);
        tool_to_part_coords(tool_pos[it], part_pos[it],
                            cur_line.x_tool_1, cur_line.y_tool_1,
                            cur_line.x_part_1, cur_line.y_part_1);
        tool_to_part_coords(tool_pos[it], part_pos[it],
                            cur_line.x_tool_2, cur_line.y_tool_2,
                            cur_line.x_part_2, cur_line.y_part_2);
        cur_line.r1 = norm(cur_line.x_part_1, cur_line.y_part_1);
        cur_line.r2 = norm(cur_line.x_part_2, cur_line.y_part_2);
        cur_line.fi_part_12  = normalized_angle(cur_line.fi_tool_12 + Dfi);
      }
    }

    for (int it1=0; it1 < n_profile; ++it1) {
      double r = pts[it1].r;

      int id_min=0, id_max=0;
      double x_min=0, y_min=0, x_max=0, y_max=0,
             theta_min=1e10,
             theta_max=-1e10;
      for (const Segment &cur_segm : tool_contour) {
        if (cur_segm.arc) { // arc
          Arc &cur_arc = *(cur_segm.arc);
          const double RCS = cur_arc.rho,
                       ROC = cur_arc.r_m,
                       fiCO = cur_arc.fi_m + PI,
                       cc = (RCS*RCS + ROC*ROC - r*r) / (2*ROC*RCS);
          if (fabs(cc) <= 1.) {
            for (const double root : {fiCO - acos(cc), fiCO + acos(cc)}) {
              double Dfi = normalized_angle(root - cur_arc.fi_part_1);
              if (Dfi > 0. && Dfi < cur_arc.fi_part_2 - cur_arc.fi_part_1) {
                double my_x = cur_arc.x_part_m + RCS*cos(root),
                       my_y = cur_arc.y_part_m + RCS*sin(root),
                       my_theta = normalized_angle(atan2(my_y,my_x));
                if (my_theta < theta_min) {
                  id_min = cur_segm.id;
                  theta_min = my_theta;
                  x_min = my_x;
                  y_min = my_y;
                }
                if (my_theta > theta_max) {
                  id_max = cur_segm.id;
                  theta_max = my_theta;
                  x_max = my_x;
                  y_max = my_y;
                }
              }
            }
          }
        } else { // straight line
          LineSeg &cur_line = *(cur_segm.line);
          const double x1 = cur_line.x_part_1,
                       y1 = cur_line.y_part_1,
                       x2 = cur_line.x_part_2,
                       y2 = cur_line.y_part_2,
                       dx = x2 - x1,
                       dy = y2 - y1,
                       L = norm(dx,dy);
          const double tmp = (x1*dy - y1*dx) / (L*L),
                       xP = tmp*dy,
                       yP = -tmp*dx,
                       d = norm(xP, yP);
//                       fiP = atan2(yP,xP);std::
          double lambda1, lambda2;
          if (fabs(dx) > fabs(dy)) {
            lambda1 = (x1 - xP) / dx;
            lambda2 = (x2 - xP) / dx;
          } else {
            lambda1 = (y1 - yP) / dy;
            lambda2 = (y2 - yP) / dy;
          }
          const double minlim = min(lambda1,lambda2)*L,
                       maxlim = max(lambda1,lambda2)*L;
          if (r >= d-1e-12) {
            const double lambda = sqrt(std::max(r*r - d*d, 0.));
            for (double root : {-lambda,lambda}) {
              if (root < maxlim && root > minlim) {
                const double my_x = xP + root*dx/L,
                             my_y = yP + root*dy/L,
                             my_theta = normalized_angle(atan2(my_y,my_x));
                if (my_theta < theta_min) {
                   id_min = cur_segm.id;
                   theta_min = my_theta;
                   x_min = my_x;
                   y_min = my_y;
                }
                if (my_theta > theta_max) {
                   id_max = cur_segm.id;
                   theta_max = my_theta;
                   x_max = my_x;
                   y_max = my_y;
                }
              }
            }
          }
        }
      } // for (const Segment &cur_segm : tool_contour)

      if (mode == 0 && theta_min <= theta_max) {
        if (theta_min < pts[it1].theta_min) {
          pts[it1].id_min = id_min;
          pts[it1].theta_min = theta_min;
          pts[it1].x_min = x_min;
          pts[it1].y_min = y_min;
        }
        if (theta_max > pts[it1].theta_max) {
          pts[it1].id_max = id_max;
          pts[it1].theta_max = theta_max;
          pts[it1].x_max = x_max;
          pts[it1].y_max = y_max;
        }
      } else if (mode==-1 && theta_min <= theta_max) {
        if (theta_max > pts[it1].theta_min) {
          pts[it1].id_min = id_max;
          pts[it1].theta_min = theta_max;
          pts[it1].x_min = x_max;
          pts[it1].y_min = y_max;
        }
      } else if (mode==1 && theta_min <= theta_max) {
        if (theta_min < pts[it1].theta_max) {
          pts[it1].id_max = id_min;
          pts[it1].theta_max = theta_min;
          pts[it1].x_max = x_min;
          pts[it1].y_max = y_min;
        }
      }
    } //it1=1,n_profile
  } //it=1,n_sim

} // cutting_process


void build_tool_contour(double m_n,     // in mm
                        double beta,    // in radians
                        double h_aP0,   // in mm
                        double alfa_P0, // in radians
                        double rho_aP0, // in mm
                        double h_fP0,   // in mm
                        std::vector<Segment> &contour) {
  // account for helix angle
  double m_t = m_n/cos(beta);
  alfa_P0 = atan(tan(alfa_P0)/cos(beta));
  // the impact of the helix angle on rho_aP0 is neglected

  // auxiliary variables
  double s_P0 = PI*m_t/2,
         h_FaP0 = h_aP0 - (1. - sin(alfa_P0))*rho_aP0;

  // build tool contour
  //1st segment (line)
  contour.emplace_back(1, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = 2*h_fP0;
    cur_line.x_tool_2   = cur_line.x_tool_1; // - s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_2   = h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //2nd segment (line)
  contour.emplace_back(2, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = h_fP0;
    cur_line.x_tool_2   = -s_P0/2. + h_FaP0*tan(alfa_P0);
    cur_line.y_tool_2   = -h_FaP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //3rd segment (arc)
  contour.emplace_back(3, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contour.back().arc);
    cur_arc.x_tool_m   = -s_P0/2. + h_FaP0*tan(alfa_P0) + rho_aP0*cos(alfa_P0);
    cur_arc.y_tool_m   = -h_aP0 + rho_aP0;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = PI + alfa_P0;
    cur_arc.fi_tool_2  = 3.*PI/2.;
  }
  //4th segment (line)
  contour.emplace_back(4, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -s_P0/2. + h_FaP0*tan(alfa_P0) + rho_aP0*cos(alfa_P0);
    cur_line.y_tool_1   = -h_aP0;
    cur_line.x_tool_2   = -cur_line.x_tool_1;
    cur_line.y_tool_2   = cur_line.y_tool_1;
    cur_line.fi_tool_12 = 0.;
  }
  //5th segment (arc)
  contour.emplace_back(5, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contour.back().arc);
    cur_arc.x_tool_m   = s_P0/2. - h_FaP0*tan(alfa_P0) - rho_aP0*cos(alfa_P0);
    cur_arc.y_tool_m   = -h_aP0 + rho_aP0;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = 3.*PI/2.;
    cur_arc.fi_tool_2  = 2.*PI - alfa_P0;
  }
  //6th segment (line)
  contour.emplace_back(6, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = s_P0/2. - h_FaP0*tan(alfa_P0);
    cur_line.y_tool_1   = -h_FaP0;
    cur_line.x_tool_2   = s_P0/2. + h_fP0*tan(alfa_P0);
    cur_line.y_tool_2   = h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //7th segment (line)
  contour.emplace_back(7, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = s_P0/2. + h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = h_fP0;
    cur_line.x_tool_2   = s_P0/2. + h_fP0*tan(alfa_P0);
    cur_line.y_tool_2   = 2*h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
}


void build_tool_with_protuberance_contour(double m_n,       // in mm
                                          double beta,      // in radians
                                          double h_prP0,    // in mm
                                          double alfa_P0,   // in radians
                                          double rho_aP0,   // in mm
                                          double h_fP0,     // in mm
                                          double h_FaP0,    // in mm
                                          double alfa_KP0,  // in radians
                                          double alfa_prP0, // in radians
                                          std::vector<Segment> &contour) {
  // account for helix angle
  double m_t = m_n/cos(beta);
  alfa_P0 = atan(tan(alfa_P0)/cos(beta));
  alfa_KP0 = atan(tan(alfa_KP0)/cos(beta));
  alfa_prP0 = atan(tan(alfa_prP0)/cos(beta));
  // the impact of the helix angle on rho_aP0 is neglected

  // auxiliary variables
  double s_P0 = PI*m_t/2.,
         h_FprP0 = h_prP0 - (1. - sin(alfa_prP0))*rho_aP0,
         x_orig = -s_P0/2. + h_FaP0*tan(alfa_P0);

  // build tool contour
  //1st segment (line)
  contour.emplace_back(1, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = 2*h_fP0;
    cur_line.x_tool_2   = cur_line.x_tool_1; // - s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_2   = h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //2nd segment (line)
  contour.emplace_back(2, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -s_P0/2. - h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = h_fP0;
    cur_line.x_tool_2   = -s_P0/2. + h_FaP0*tan(alfa_P0);
    cur_line.y_tool_2   = -h_FaP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //3rd segment (line)
  contour.emplace_back(3, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = x_orig;
    cur_line.y_tool_1   = -h_FaP0;
    cur_line.x_tool_2   = x_orig + h_FprP0*tan(alfa_prP0);
    cur_line.y_tool_2   = -h_FaP0 - h_FprP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //4th segment (arc)
  contour.emplace_back(4, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contour.back().arc);
    cur_arc.x_tool_m   = x_orig + h_FprP0*tan(alfa_prP0) + rho_aP0*cos(alfa_prP0);
    cur_arc.y_tool_m   = -h_FaP0 - h_prP0 + rho_aP0;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = PI + alfa_P0;
    cur_arc.fi_tool_2  = 3.*PI/2.;
  }
  //5th segment (line)
  contour.emplace_back(5, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = x_orig + h_FprP0*tan(alfa_prP0) + rho_aP0*cos(alfa_prP0);
    cur_line.y_tool_1   = -h_FaP0 - h_prP0;
    cur_line.x_tool_2   = -cur_line.x_tool_1;
    cur_line.y_tool_2   = cur_line.y_tool_1;
    cur_line.fi_tool_12 = 0.;
  }
  //6th segment (arc)
  contour.emplace_back(6, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contour.back().arc);
    cur_arc.x_tool_m   = -x_orig - h_FprP0*tan(alfa_prP0) - rho_aP0*cos(alfa_prP0);
    cur_arc.y_tool_m   = -h_FaP0 - h_prP0 + rho_aP0;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = 3.*PI/2.;
    cur_arc.fi_tool_2  = 2.*PI - alfa_prP0;
  }
  //7th segment (line)
  contour.emplace_back(7, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = -x_orig;
    cur_line.y_tool_1   = -h_FaP0;
    cur_line.x_tool_2   = -x_orig - h_FprP0*tan(alfa_prP0);
    cur_line.y_tool_2   = -h_FaP0 - h_FprP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //8th segment (line)
  contour.emplace_back(8, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = s_P0/2. + h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = h_fP0;
    cur_line.x_tool_2   = s_P0/2. - h_FaP0*tan(alfa_P0);
    cur_line.y_tool_2   = -h_FaP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //9th segment (line)
  contour.emplace_back(9, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contour.back().line);
    cur_line.x_tool_1   = s_P0/2. + h_fP0*tan(alfa_P0);
    cur_line.y_tool_1   = 2*h_fP0;
    cur_line.x_tool_2   = cur_line.x_tool_1;
    cur_line.y_tool_2   = h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
}


void build_internal_gear_tool_contours(double m_n,       // in mm
                                       double beta,      // in radians
                                       double alfa_P0,   // in degrees
                                       double h_aP0,     // in mm
                                       double rho_aP0,   // in mm
                                       double x,         // dimensionless
                                       int z,
                                       int z_0,
                                       double r_a0,      // in mm
                                       std::vector<std::vector<Segment> > &contours,
                                       double &r_Ff,
                                       double &a_tool) {
  // account for helix angle
  double m_t = m_n/cos(beta);
  alfa_P0 = atan(tan(alfa_P0)/cos(beta));
  // the impact of the helix angle on rho_aP0 is neglected

  // auxiliary variables
  double s_P0 = PI*m_t/2.,
         h_fP0 = s_P0/(2*tan(alfa_P0));
  double x_m, y_m;
  tool_params_internal_gear(z, z_0, m_t, alfa_P0, x*m_n, r_a0,
                            rho_aP0, h_aP0, s_P0,
                            r_Ff, a_tool, x_m, y_m);
  // define a contour for the auxiliary tool (rack) that will generate
  // the actual cutting tool (wheel)

  contours.resize(3); // 3 tools needed for simulation of internal gear cutting

  //1st segment (line)
  contours[0].clear();
  contours[0].emplace_back(1, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contours[0].back().line);
    cur_line.x_tool_1   = -s_P0/2.-2*h_aP0*tan(alfa_P0);
    cur_line.y_tool_1   = 2*h_aP0;
    cur_line.x_tool_2   = 0.;
    cur_line.y_tool_2   = -h_fP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }
  //2nd segment (line)
  contours[0].emplace_back(2, std::make_unique<LineSeg>());
  {
    LineSeg &cur_line = *(contours[0].back().line);
    cur_line.x_tool_1   = 0.;
    cur_line.y_tool_1   = -h_fP0;
    cur_line.x_tool_2   = s_P0/2.+2*h_aP0*tan(alfa_P0);
    cur_line.y_tool_2   = 2*h_aP0;
    cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                cur_line.x_tool_2 - cur_line.x_tool_1);
  }

  //3rd segment (arc)
  contours[1].clear();
  contours[1].emplace_back(3, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contours[1].back().arc);
    cur_arc.x_tool_m   = x_m;
    cur_arc.y_tool_m   = y_m;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = 0.;
    cur_arc.fi_tool_2  = PI/2.; //atan2(y_m,x_m)
  }

  //4th segment (arc)
  contours[2].clear();
  contours[2].emplace_back(4, std::make_unique<Arc>());
  {
    Arc &cur_arc = *(contours[2].back().arc);
    cur_arc.x_tool_m   = -x_m;
    cur_arc.y_tool_m   = y_m;
    cur_arc.rho        = rho_aP0;
    cur_arc.fi_tool_1  = PI/2.; //atan2(y_m,x_m)
    cur_arc.fi_tool_2  = PI;
  }
} // build_internal_gear_tool_contours



// Generating the tooth contour from a gear specification given in fs_in
void generate_profile(ToolType tool_type,
                      const double m_n,
                      const double beta,
                      const double h_aP0,
                      const int z,
                      const double x,
                      const double r_a,
                      const int z_0,       // only for internal gears
                      const double r_Ff,   // only for internal gears
                      const double a_tool, // only for internal gears
                      const int n_sim,
                      const int n_profile,
                      const int n_tip,
                      const std::vector<std::vector<Segment> > &tool_contours,
                      std::vector<std::vector<double> > &contour_xy, // output
                      std::vector<int> &contour_xy_id)               // output
{
  double r_f;
  std::vector<ContourPoint> pts(n_profile, {0,0,0.,0.,0.,0.,0.,0.,0.});
  {
    double r_w = abs(z)*m_n/cos(beta)/2.;
    r_f = r_w + (z > 0 ? 1 : -1) * (x*m_n - h_aP0);
    std::vector<CSpos> tool_pos(n_sim, {0., r_w + (z > 0 ? x : -x)*m_n, 0.}),
                       part_pos(n_sim, {0., 0., 0.});
    for (int it=0; it < n_sim; ++it) {
      part_pos[it].fi = PI * (0.5*double(it) / double(n_sim-1) - 0.25);
      tool_pos[it].x = -r_w*part_pos[it].fi;
    }

    // initialize and calculate tooth profile
    const double incr = (r_a-r_f) / double(n_profile-1);
    for (int it=0; it < n_profile; ++it) {
      pts[it].r = r_f + it*incr;
      pts[it].theta_min = 1e8;
      pts[it].theta_max = -1e8;
    }

    cutting_process(tool_contours[0], tool_pos, part_pos, pts, 0);
    if (tool_type == ToolType::INTERNAL) {
      // redefine the tool trajectory
      for (int it=0; it < n_sim; ++it) {
        tool_pos[it].x = 0.;
        tool_pos[it].y = a_tool;
        tool_pos[it].fi = (abs(z)*part_pos[it].fi) / double(z_0);
      }
      for (ContourPoint &pt : pts) {
        if (pt.r <= r_Ff)
          break;
        pt.theta_min = -1e8;
        pt.theta_max = 1e8;
        pt.id_min = -1;
        pt.id_max = -1;
      }
      cutting_process(tool_contours[1], tool_pos, part_pos, pts, 1);
      cutting_process(tool_contours[2], tool_pos, part_pos, pts, -1);
    }
  }

  // locate the first point of the profile (corresponds to the tooth root)
  int n_start = 0, n_root = 0;
  if (z > 0) {
    while (pts[n_start].id_max == -1)
      ++n_start;
    n_root = int(ceil((pts[n_start].theta_max - PI/2.) /
                      (pts[n_start+1].theta_max - pts[n_start].theta_max)));
  } else {
    while (pts[n_start].theta_min < PI/2 - PI/abs(z) ||
           pts[n_start].id_min == -1)
      ++n_start;
    n_root = int(ceil((pts[n_start].theta_min - PI/2 + PI/abs(z)) /
                      (pts[n_start+1].theta_min - pts[n_start].theta_min)));
  }
  // completion of the tooth root profile with an arc
  contour_xy.assign(n_root + n_tip + n_profile - n_start, std::vector<double>(2,0));
  contour_xy_id.assign(n_root + n_tip + n_profile - n_start, 0);
  int N_contour=0;
  if (n_root > 0) {
    double incr = (z > 0)
                ? (pts[n_start].theta_max - PI/2) / double(n_root)
                : (pts[n_start].theta_min - PI/2 + PI/abs(z)) / double(n_root);
    for (int it=0; it < n_root; ++it) {
      contour_xy_id[it] = -1;
      const double my_theta = PI/2 - PI/abs(z) + incr*it;
      contour_xy[it][0] = r_f*cos(my_theta);
      contour_xy[it][1] = r_f*sin(my_theta);
    }
    N_contour = n_root;
  }

  // transfer data from pts to contour_xy with rotation (for
  // external gears) by half pitch
  if (z > 0) {
    const double dx = sin(PI/z),
                 dy = cos(PI/z);
    for (int it=n_start; it < n_profile; ++it) {
      if (pts[it].theta_min < pts[it].theta_max) {
        contour_xy_id[N_contour] = pts[it].id_max;
        contour_xy[N_contour][0] = pts[it].x_max*dy + pts[it].y_max*dx;
        contour_xy[N_contour][1] = -pts[it].x_max*dx + pts[it].y_max*dy;
        ++N_contour;
      } else
        std::cout << "Error in the calculated profile at r="
                  << pts[it].r << std::endl;
    }
  } else {
    for (int it=n_start; it < n_profile; ++it) {
      if (pts[it].theta_min < pts[it].theta_max) {
        contour_xy_id[N_contour] = pts[it].id_min;
        contour_xy[N_contour][0] = pts[it].x_min;
        contour_xy[N_contour][1] = pts[it].y_min;
        ++N_contour;
      } else
        std::cout << "Error in the calculated profile at r="
                  << pts[it].r << std::endl;
    }
  }

  // completion of the tooth tip profile with an arc
  const double r = pts[n_profile-1].r,
               my_theta = (z > 0) ? pts[n_profile-1].theta_max - PI/double(z)
                                   : pts[n_profile-1].theta_min,
               incr = (PI/2 - my_theta) / double(n_tip);
  for (int it=0; it < n_tip; ++it) {
    contour_xy_id[N_contour] = -2;
    contour_xy[N_contour][0] = r*cos(my_theta + incr*(it+1));
    contour_xy[N_contour][1] = r*sin(my_theta + incr*(it+1));
    ++N_contour;
  }
}


void generate_external_gear_tooth_profile
  (double m_n, double beta, double h_aP0, double alfa_P0, double rho_aP0, // tooth shape parameters
   int z, double x, double r_a,                                           // gear parameters
   double h_fP0,                                                          // tool specific parameters
   int n_sim, int n_profile, int n_tip,                                   // simulation parameters
   std::vector<std::vector<Segment> > &tool_contours,  // output
   std::vector<std::vector<double> > &contour_xy,      // output
   std::vector<int> &contour_xy_id)                    // output
{
  build_tool_contour(m_n, beta, h_aP0, alfa_P0, rho_aP0,
                     h_fP0,
                     tool_contours[0]);
  generate_profile(ToolType::EXTERNAL,
                   m_n, beta, h_aP0, z, x, r_a, 0, 0., 0.,
                   n_sim, n_profile, n_tip, tool_contours,
                   contour_xy, contour_xy_id);
}

void generate_external_gear_tooth_profile // with protuberance
  (double m_n, double beta, double h_aP0, double alfa_P0, double rho_aP0, // tooth shape parameters
   int z, double x, double r_a,                                           // gear parameters
   double h_fP0, double h_FaP0, double alfa_KP0, double alfa_prP0,        // tool specific parameters
   int n_sim, int n_profile, int n_tip,                                   // simulation parameters
   std::vector<std::vector<Segment> > &tool_contours,  // output
   std::vector<std::vector<double> > &contour_xy,      // output
   std::vector<int> &contour_xy_id)                    // output
{
  const double h_prP0 = h_aP0;
  build_tool_with_protuberance_contour(m_n, beta, h_prP0, alfa_P0, rho_aP0,
                                       h_fP0, h_FaP0, alfa_KP0, alfa_prP0,
                                       tool_contours[0]);
  generate_profile(ToolType::EXTERNAL_WITH_PROTUBERANCE,
                   m_n, beta, h_aP0, z, x, r_a, 0, 0., 0.,
                   n_sim, n_profile, n_tip, tool_contours,
                   contour_xy, contour_xy_id);
}


void generate_internal_gear_tooth_profile // internal gear cutting
  (double m_n, double beta, double h_aP0, double alfa_P0, double rho_aP0, // tooth shape parameters
   int z, double x, double r_a,                                           // gear parameters
   int z_0, double r_a0,                                                  // tool specific parameters
   int n_sim, int n_profile, int n_tip,                                   // simulation parameters
   std::vector<std::vector<Segment> > &tool_contours,     // output
   std::vector<std::vector<double> > &contour_xy,         // output
   std::vector<int> &contour_xy_id)                       // output
{
  double r_Ff, a_tool;
  build_internal_gear_tool_contours(m_n, beta, h_aP0, alfa_P0, rho_aP0,
                                    z, x, z_0, r_a0,
                                    tool_contours, r_Ff, a_tool);
  generate_profile(ToolType::INTERNAL,
                   m_n, beta, h_aP0, z, x, r_a,
                   z_0, r_Ff, a_tool, // args only needed for internal gears
                   n_sim, n_profile, n_tip, tool_contours,
                   contour_xy, contour_xy_id);
}
