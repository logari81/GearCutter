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

void write_to_file(const std::string &fname, const std::vector<Segment> &contour) {

  std::fstream fs(fname, std::fstream::out);
  for (const Segment &cur_segm : contour) {
    if (cur_segm.arc) {
      const Arc &cur_arc = *(cur_segm.arc);
      fs << "Arc" << std::endl << "x_M=" << cur_arc.x_tool_m
                               << " y_M=" << cur_arc.y_tool_m
                               << " rho=" << cur_arc.rho;
    } else {
      const LineSeg &cur_line = *(cur_segm.line);
      fs << "LineSeg" << std::endl << "x_A=" << cur_line.x_tool_1
                                   << " y_A=" << cur_line.y_tool_1
                                   << " x_E=" << cur_line.x_tool_2
                                   << " y_E=" << cur_line.y_tool_2;
    }
    fs << std::endl;
  }
} // write_to_file

void read_polygon_tool_contour(double &m_n,
                               const double beta,
                               std::ifstream &fs_in,
                               std::vector<Segment> &contour,
                               double &h_aP0) {
  std::string tmp;
  fs_in >> tmp >> m_n;

  h_aP0 = 0.;
  contour.clear();
  int it = 0;
  double x0, y0, x, y;
  fs_in >> x0 >> y0;
  while (fs_in >> x >> y) {
    contour.emplace_back(++it, std::make_unique<LineSeg>());
    {
      LineSeg &cur_line = *(contour.back().line);
      cur_line.x_tool_1   = x0;
      cur_line.y_tool_1   = y0/cos(beta);
      cur_line.x_tool_2   = x;
      cur_line.y_tool_2   = y/cos(beta);
      cur_line.fi_tool_12 = atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                 cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    x0 = x;
    y0 = y;
    if (h_aP0 < -y) h_aP0 = -y;
  }
}

///////////////////////////////////////////////////////////////////////
// Main Programm
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

  int n_sim = 1000,
      n_profile = 1200,
      n_tip = 7;
  std::string input_file("Gear_data.txt"),
              output_file("tooth_profile.dat"),
              output_file2("tool_contour.dat");

  // commandline parameters
  for (int i = 1; i < argc; ++i) {
    bool missing_argument(false);
    std::string option(argv[i]);
    if (option == "-i") {
      if (i + 1 < argc) input_file = argv[++i];
      else              missing_argument = true; 
    } else if (option == "-o") {
      if (i + 1 < argc) {
        output_file = output_file2 = argv[++i];
        output_file2.append(".1");
      } else
        missing_argument = true; 
    } else if (option == "-n_sim") {
      if (i + 1 < argc)  n_sim = std::stoi(argv[++i]);
      else               missing_argument = true; 
    } else if (option == "-n_profile") {
      if (i + 1 < argc)  n_profile = std::stoi(argv[++i]);
      else               missing_argument = true; 
    } else if (option == "-n_tip") {
      if (i + 1 < argc)  n_tip = std::stoi(argv[++i]);
      else               missing_argument = true; 
    } else {
      std::cerr << "Unknown option " << option << std::endl;
      return 1;
    }  

    if (missing_argument) {
      std::cerr << option << " option requires one argument." << std::endl;
      return 1;
    }  
  }

  std::ifstream fs_in(input_file, std::ifstream::in);
  if (!fs_in.is_open()) {
    std::cout << "Input file \"" << input_file << "\" not found." << std::endl;
    return 1;
  }

  std::vector<std::vector<double> > contour_xy;
  std::vector<int> contour_xy_id;
  std::vector<std::vector<Segment> > tool_contours(1); // just 1 tool contour by default

  {
    ToolType tool_type;
    int z, z_0, itmp;
    double r_a, x, beta, m_n, h_aP0;
    std::string tmp;
    fs_in >> itmp;              tool_type = static_cast<ToolType>(itmp);
    fs_in >> tmp >> z;
    fs_in >> tmp >> r_a;        r_a = fabs(r_a)/2;
    fs_in >> tmp >> x;
    fs_in >> tmp >> beta;       beta *= PI/180.;

    std::vector<Segment> contour1, contour2;
    if (tool_type == ToolType::EXTERNAL) { // without protuberance
      double alfa_P0, h_fP0, rho_aP0;
      fs_in >> tmp >> alfa_P0;   alfa_P0 *= PI/180.;
      fs_in >> tmp >> m_n;
      fs_in >> tmp >> h_fP0;     h_fP0    *= m_n;
      fs_in >> tmp >> h_aP0;     h_aP0    *= m_n;
      fs_in >> tmp >> rho_aP0;   rho_aP0  *= m_n;
      generate_external_gear_tooth_profile
        (m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and gear parameters
         h_fP0,                                         // tool specific parameters
         n_sim, n_profile, n_tip,                       // simulation parameters
         tool_contours, contour_xy, contour_xy_id);
    } else if (tool_type == ToolType::EXTERNAL_WITH_PROTUBERANCE) {
      double alfa_P0, alfa_KP0, alfa_prP0, h_fP0, h_FaP0, h_prP0, rho_aP0;
      fs_in >> tmp >> alfa_P0;   alfa_P0 *= PI/180.;
      fs_in >> tmp >> alfa_KP0;  alfa_KP0 *= PI/180.;
      fs_in >> tmp >> alfa_prP0; alfa_prP0 *= PI/180.;
      fs_in >> tmp >> m_n;
      fs_in >> tmp >> h_fP0;     h_fP0   *= m_n;
      fs_in >> tmp >> h_FaP0;    h_FaP0  *= m_n;
      fs_in >> tmp >> h_prP0;    h_prP0  *= m_n;
      fs_in >> tmp >> rho_aP0;   rho_aP0 *= m_n;
      generate_external_gear_tooth_profile
        (m_n, beta, h_prP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and parameters
         h_fP0, h_FaP0, alfa_KP0, alfa_prP0,             // tool specific parameters
         n_sim, n_profile, n_tip,                        // simulation parameters
         tool_contours, contour_xy, contour_xy_id);
    } else if (tool_type == ToolType::INTERNAL) { // internal gear
      double alfa_P0, rho_aP0, r_a0;
      fs_in >> tmp >> alfa_P0;   alfa_P0 *= PI/180.;
      fs_in >> tmp >> m_n;
      fs_in >> tmp >> h_aP0;     h_aP0   *= m_n;
      fs_in >> tmp >> rho_aP0;   rho_aP0 *= m_n;
      fs_in >> tmp >> r_a0;      r_a0    /= 2.;  // cutting wheel outer diameter/radius
      fs_in >> tmp >> z_0;   // cutting wheel number of teeth
      generate_internal_gear_tooth_profile
        (m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and gear parameters
         z_0, r_a0,                                     // tool specific parameters
         n_sim, n_profile, n_tip,                       // simulation parameters
         tool_contours, contour_xy, contour_xy_id);
    } else if (tool_type == ToolType::POLYLINE) { // polygon from file
      read_polygon_tool_contour(m_n, beta, fs_in, tool_contours[0], h_aP0);
      generate_profile(tool_type, m_n, beta, h_aP0, z, x, r_a, 0, 0., 0.,
                       n_sim, n_profile, n_tip, tool_contours,
                       contour_xy, contour_xy_id);
    }
  }

  
  std::fstream fs_out(output_file, std::fstream::out);
  fs_out << "#             x               y  id" << std::endl;
  fs_out << std::fixed << std::setprecision(10);
  for (int it=0; it < int(contour_xy.size()); ++it)
    fs_out << std::setw(15) << contour_xy[it][0] << ","
           << std::setw(15) << contour_xy[it][1] << ","
           << std::setw(3) << contour_xy_id[it] << std::endl;
  contour_xy.clear();
  contour_xy_id.clear();

  write_to_file(output_file2, tool_contours[0]);
  return 0;
} // main

