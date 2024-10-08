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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <gear_cutter.h>

const double PI = M_PI;

typedef struct {
  PyObject_HEAD
  double m_n;
  double h_aP0;
  double alfa_P0;
  double rho_aP0;
  double h_fP0;
} cutting_tool_rack;

typedef struct {
  PyObject_HEAD
  double m_n;
  double h_prP0;
  double alfa_P0;
  double rho_aP0;
  double h_fP0;
  double h_FaP0;
  double alfa_KP0;
  double alfa_prP0;
} cutting_tool_rack_with_protuberance;

typedef struct {
  PyObject_HEAD
  double m_n;
  double h_aP0;
  double alfa_P0;
  double rho_aP0;
  double h_fP0;
  double z_0;
  double d_a0;
} cutting_tool_pinion; // cutting tool for internal gear
      
// getters for cutting_tool_rack
static PyObject* cutting_tool_rack_get_m_n(PyObject* self, void* closure) 
{ return PyFloat_FromDouble(((cutting_tool_rack*)self)->m_n); }
static PyObject* cutting_tool_rack_get_h_aP0(PyObject* self, void* closure) 
{ return PyFloat_FromDouble(((cutting_tool_rack*)self)->h_aP0); }
static PyObject* cutting_tool_rack_get_alfa_P0(PyObject* self, void* closure) 
{ return PyFloat_FromDouble(((cutting_tool_rack*)self)->alfa_P0); }
static PyObject* cutting_tool_rack_get_rho_aP0(PyObject* self, void* closure) 
{ return PyFloat_FromDouble(((cutting_tool_rack*)self)->rho_aP0); }
static PyObject* cutting_tool_rack_get_h_fP0(PyObject* self, void* closure) 
{ return PyFloat_FromDouble(((cutting_tool_rack*)self)->h_fP0); }
// getters for cutting_tool_rack_with_protuberance
static PyObject* cutting_tool_rack_with_protuberance_get_m_n(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->m_n); }
static PyObject* cutting_tool_rack_with_protuberance_get_h_prP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->h_prP0); }
static PyObject* cutting_tool_rack_with_protuberance_get_alfa_P0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->alfa_P0); }
static PyObject* cutting_tool_rack_with_protuberance_get_rho_aP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->rho_aP0); }
static PyObject* cutting_tool_rack_with_protuberance_get_h_fP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->h_fP0); }
static PyObject* cutting_tool_rack_with_protuberance_get_h_FaP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->h_FaP0); }
static PyObject* cutting_tool_rack_with_protuberance_get_alfa_KP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->alfa_KP0); }
static PyObject* cutting_tool_rack_with_protuberance_get_alfa_prP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_rack_with_protuberance*)self)->alfa_prP0); }
// getters for cutting_tool_pinion
static PyObject* cutting_tool_pinion_get_m_n(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->m_n); }
static PyObject* cutting_tool_pinion_get_h_aP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->h_aP0); }
static PyObject* cutting_tool_pinion_get_alfa_P0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->alfa_P0); }
static PyObject* cutting_tool_pinion_get_rho_aP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->rho_aP0); }
static PyObject* cutting_tool_pinion_get_h_fP0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->h_fP0); }
static PyObject* cutting_tool_pinion_get_z_0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->z_0); }
static PyObject* cutting_tool_pinion_get_d_a0(PyObject* self, void* closure)
{ return PyFloat_FromDouble(((cutting_tool_pinion*)self)->d_a0); }

// setters for cutting_tool_rack
static int cutting_tool_rack_set_m_n(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack*)self)->m_n = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_set_h_aP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack*)self)->h_aP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_set_alfa_P0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack*)self)->alfa_P0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_set_rho_aP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack*)self)->rho_aP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_set_h_fP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack*)self)->h_fP0 = PyFloat_AsDouble(value);
  return 0;
}
// setters for cutting_tool_rack_with_protuberance
static int cutting_tool_rack_with_protuberance_set_m_n(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->m_n = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_h_prP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->h_prP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_alfa_P0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->alfa_P0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_rho_aP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->rho_aP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_h_fP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->h_fP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_h_FaP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->h_FaP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_alfa_KP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->alfa_KP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_rack_with_protuberance_set_alfa_prP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_rack_with_protuberance*)self)->alfa_prP0 = PyFloat_AsDouble(value);
  return 0;
}
// setters for cutting_tool_pinion
static int cutting_tool_pinion_set_m_n(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->m_n = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_h_aP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->h_aP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_alfa_P0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->alfa_P0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_rho_aP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->rho_aP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_h_fP0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->h_fP0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_z_0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->z_0 = PyFloat_AsDouble(value);
  return 0;
}
static int cutting_tool_pinion_set_d_a0(PyObject* self, PyObject* value, void* closure) {
  ((cutting_tool_pinion*)self)->d_a0 = PyFloat_AsDouble(value);
  return 0;
}

static PyGetSetDef cutting_tool_rack_getset[] = {
  {"m_n", cutting_tool_rack_get_m_n, cutting_tool_rack_set_m_n,
   "The normal module", NULL},
  {"h_aP0", cutting_tool_rack_get_h_aP0, cutting_tool_rack_set_h_aP0,
   "The height factor of the cutting tool teeth tip above the reference line", NULL},
  {"alfa_P0", cutting_tool_rack_get_alfa_P0, cutting_tool_rack_set_alfa_P0,
   "The pressure angle of the cutting tool teeth", NULL},
  {"rho_aP0", cutting_tool_rack_get_rho_aP0, cutting_tool_rack_set_rho_aP0,
   "The fillet radius factor at the tip of the cutting tool teeth", NULL},
  {"h_fP0", cutting_tool_rack_get_h_fP0, cutting_tool_rack_set_h_fP0,
   "The depth factor of the cutting tool teeth root below the reference line", NULL},
  {NULL}  // Sentinel
};
static PyGetSetDef cutting_tool_rack_with_protuberance_getset[] = {
  {"m_n", cutting_tool_rack_with_protuberance_get_m_n, cutting_tool_rack_with_protuberance_set_m_n,
   "The normal module", NULL},
  {"h_prP0", cutting_tool_rack_with_protuberance_get_h_prP0, cutting_tool_rack_with_protuberance_set_h_prP0,
   "The height factor of the cutting tool teeth tip above the reference line", NULL},
  {"alfa_P0", cutting_tool_rack_with_protuberance_get_alfa_P0, cutting_tool_rack_with_protuberance_set_alfa_P0,
   "The pressure angle of the cutting tool teeth", NULL},
  {"rho_aP0", cutting_tool_rack_with_protuberance_get_rho_aP0, cutting_tool_rack_with_protuberance_set_rho_aP0,
   "The fillet radius factor at the tip of the cutting tool teeth", NULL},
  {"h_fP0", cutting_tool_rack_with_protuberance_get_h_fP0, cutting_tool_rack_with_protuberance_set_h_fP0,
   "The depth factor of the cutting tool teeth root below the reference line", NULL},
  {"h_FaP0", cutting_tool_rack_with_protuberance_get_h_FaP0, cutting_tool_rack_with_protuberance_set_h_FaP0,
   "The depth factor for the end of the straight flank of the protuberance tool rack teeth", NULL},
  {"alfa_KP0", cutting_tool_rack_with_protuberance_get_alfa_KP0, cutting_tool_rack_with_protuberance_set_alfa_KP0,
   "The pressure angle at the chamfer root of the protuberance tool rack", NULL},
  {"alfa_prP0", cutting_tool_rack_with_protuberance_get_alfa_prP0, cutting_tool_rack_with_protuberance_set_alfa_prP0,
   "The pressure angle of the protuberance section of the tool teeth", NULL},
  {NULL}  // Sentinel
};
static PyGetSetDef cutting_tool_pinion_getset[] = {
  {"m_n", cutting_tool_pinion_get_m_n, cutting_tool_pinion_set_m_n,
   "The normal module", NULL},
  {"h_aP0", cutting_tool_pinion_get_h_aP0, cutting_tool_pinion_set_h_aP0,
   "The height factor of the cutting tool teeth tip above the reference line", NULL},
  {"alfa_P0", cutting_tool_pinion_get_alfa_P0, cutting_tool_pinion_set_alfa_P0,
   "The pressure angle of the cutting tool teeth", NULL},
  {"rho_aP0", cutting_tool_pinion_get_rho_aP0, cutting_tool_pinion_set_rho_aP0,
   "The fillet radius factor at the tip of the cutting tool teeth", NULL},
  {"h_fP0", cutting_tool_pinion_get_h_fP0, cutting_tool_pinion_set_h_fP0,
   "The depth factor of the cutting tool teeth root below the reference line", NULL},
  {"z_0", cutting_tool_pinion_get_z_0, cutting_tool_pinion_set_z_0,
   "The depth of the pinion", NULL},
  {"d_a0", cutting_tool_pinion_get_d_a0, cutting_tool_pinion_set_d_a0,
   "The tip diameter of the cutting tool pinion", NULL},
  {NULL}  // Sentinel
};


// Initialize method for cutting_tool_rack object
static int
cutting_tool_rack_init(cutting_tool_rack *self, PyObject *args, PyObject *kwds)
{
  char* keywords[6];
  int i = 0;
  for (const std::string kw : {"m_n", "h_aP0", "alfa_P0", "rho_aP0", "h_fP0"}) {
    keywords[i] = new char[kw.length() + 1];
    strlcpy(keywords[i++], kw.c_str(), kw.length() + 1);
  }
  keywords[5] = NULL;
  self->h_aP0 = 1.25;
  self->alfa_P0 = 20.;
  self->rho_aP0 = 0.38;
  self->h_fP0 = 1.;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|dddd", keywords,
                                   &self->m_n,
                                   &self->h_aP0,
                                   &self->alfa_P0,
                                   &self->rho_aP0,
                                   &self->h_fP0))
    return -1;
  return 0;
}static int
// Initialize method for cutting_tool_rack_with_protuberance object
cutting_tool_rack_with_protuberance_init(cutting_tool_rack_with_protuberance *self, PyObject *args, PyObject *kwds)
{
  char* keywords[9];
  int i = 0;
  for (const std::string kw : {"m_n", "h_prP0", "alfa_P0", "rho_aP0", "h_fP0", "h_FaP0", "alfa_KP0", "alfa_prP0"}) {
    keywords[i] = new char[kw.length() + 1];
    strlcpy(keywords[i++], kw.c_str(), kw.length() + 1);
  }
  keywords[8] = NULL;
  self->h_prP0 = 1.25;
  self->alfa_P0 = 20.;
  self->rho_aP0 = 0.38;
  self->h_fP0 = 1.;
  self->h_FaP0 = 0.8;
  self->alfa_KP0 = 37.5;
  self->alfa_prP0 = 12.5;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|dddddd", keywords,
                                   &self->m_n,
                                   &self->h_prP0,
                                   &self->alfa_P0,
                                   &self->rho_aP0,
                                   &self->h_fP0,
                                   &self->h_FaP0,
                                   &self->alfa_KP0,
                                   &self->alfa_prP0))
    return -1;
  return 0;
}
// Initialize method for cutting_tool_pinion object
static int
cutting_tool_pinion_init(cutting_tool_pinion *self, PyObject *args, PyObject *kwds)
{
  char* keywords[8];
  int i = 0;
  for (const std::string kw : {"m_n", "h_aP0", "alfa_P0", "rho_aP0", "h_fP0", "z_0", "d_a0"}) {
    keywords[i] = new char[kw.length() + 1];
    strlcpy(keywords[i++], kw.c_str(), kw.length() + 1);
  }
  keywords[7] = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd", keywords,
                                   &self->m_n,
                                   &self->h_aP0,
                                   &self->alfa_P0,
                                   &self->rho_aP0,
                                   &self->h_fP0,
                                   &self->z_0,
                                   &self->d_a0))
    return -1;
  return 0;
}

// Type definition for cutting_tool_rack object
static PyTypeObject cutting_tool_rack_type = {
  .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "gear_cutter.cutting_tool_rack",
  .tp_basicsize = sizeof(cutting_tool_rack),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "Type containing the parameters of gear cutting tool rack.\n"
            "It has five float members: m_n, h_aP0, alfa_P0, rho_aP0 and h_fP0.\n"
            "It can be initialized with these five float members as the constructor arguments.",
  .tp_getset = cutting_tool_rack_getset,
  .tp_init = (initproc)cutting_tool_rack_init,
  .tp_new = PyType_GenericNew,
};
// Type definition for cutting_tool_rack_with_protuberance object
static PyTypeObject cutting_tool_rack_with_protuberance_type = {
  .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "gear_cutter.cutting_tool_rack_with_protuberance",
  .tp_basicsize = sizeof(cutting_tool_rack_with_protuberance),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "Type containing the parameters of gear cutting tool rack with protuberance.\n"
            "It has eight float members: m_n, h_prP0, alfa_P0, rho_aP0, h_fP0, h_FaP0, alfa_KP0, alfa_prP0.\n"
            "It can be initialized with these eight float members as the constructor arguments.",
  .tp_getset = cutting_tool_rack_with_protuberance_getset,
  .tp_init = (initproc)cutting_tool_rack_with_protuberance_init,
};
// Type definition for cutting_tool_pinion object
static PyTypeObject cutting_tool_pinion_type = {
  .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "gear_cutter.cutting_tool_pinion",
  .tp_basicsize = sizeof(cutting_tool_pinion),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "Type containing the parameters of gear cutting tool pinion for cutting internal gears.\n"
            "It has seven float members: m_n, h_aP0, alfa_P0, rho_aP0, h_fP0, z_0, d_a0.\n"
            "It can be initialized with these seven float members as the constructor arguments.",
  .tp_getset = cutting_tool_pinion_getset,
  .tp_init = (initproc)cutting_tool_pinion_init,
};


// Generate tooth profile function
static PyObject *
gear_cutter_generate_tooth_profile(PyObject *self, PyObject *args,
                                   PyObject* kwargs) {
  /* Parse arguments */
//  static char* keywords[] = {"tool", "z", "beta", "x", "d_a",
//                             "N_sim", "N_profile", "N_tip", NULL};
  char* keywords[9];
  int i = 0;
  for (const std::string kw : { "tool", "z", "beta", "x", "d_a",
                                "N_sim", "N_profile", "N_tip"}) {
    keywords[i] = new char[kw.length() + 1];
    strlcpy(keywords[i++], kw.c_str(), kw.length() + 1);
  }
  keywords[8] = NULL;

  PyObject *tool = NULL;
  int z;
  double beta, x, d_a;
  int N_sim=1000, N_profile=1200, N_tip=7;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                   "Oiddd|iii", keywords,
                                   &tool,
                                   &z,
                                   &beta,
                                   &x,
                                   &d_a,
                                   &N_sim,
                                   &N_profile,
                                   &N_tip))
    return NULL;

  for (int i = 0; i < 9; ++i)
    delete[] keywords[i];

  beta *= PI/180.;
  double r_a = abs(d_a)/2;

  std::vector<std::vector<double> > contourXY;
  std::vector<int> contourXY_id;
  std::vector<std::vector<Segment> > tool_contours(1);
  if (PyObject_TypeCheck(tool, &cutting_tool_rack_type)) { // without protuberance
    double m_n     = ((cutting_tool_rack*)tool)->m_n;
    double h_aP0   = m_n * ((cutting_tool_rack*)tool)->h_aP0;
    double alfa_P0 = PI/180 * ((cutting_tool_rack*)tool)->alfa_P0;
    double rho_aP0 = m_n * ((cutting_tool_rack*)tool)->rho_aP0;
    double h_fP0   = m_n * ((cutting_tool_rack*)tool)->h_fP0;
    generate_external_gear_tooth_profile
      (m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and gear parameters
       h_fP0,                                         // tool specific parameters
       N_sim, N_profile, N_tip,                       // simulation parameters
       tool_contours, contourXY, contourXY_id);
  } else if (PyObject_TypeCheck(tool, &cutting_tool_rack_with_protuberance_type)) { // with protuberance
    double m_n       = ((cutting_tool_rack_with_protuberance*)tool)->m_n;
    double h_prP0    = m_n * ((cutting_tool_rack_with_protuberance*)tool)->h_prP0;
    double alfa_P0   = PI/180 * ((cutting_tool_rack_with_protuberance*)tool)->alfa_P0;
    double rho_aP0   = m_n * ((cutting_tool_rack_with_protuberance*)tool)->rho_aP0;
    double h_fP0     = m_n * ((cutting_tool_rack_with_protuberance*)tool)->h_fP0;
    double h_FaP0    = m_n * ((cutting_tool_rack_with_protuberance*)tool)->h_FaP0;
    double alfa_KP0  = PI/180 * ((cutting_tool_rack_with_protuberance*)tool)->alfa_KP0;
    double alfa_prP0 = PI/180 * ((cutting_tool_rack_with_protuberance*)tool)->alfa_prP0;
    generate_external_gear_tooth_profile
      (m_n, beta, h_prP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and gear parameters
       h_fP0, h_FaP0, alfa_KP0, alfa_prP0,             // tool specific parameters
       N_sim, N_profile, N_tip,                        // simulation parameters
       tool_contours, contourXY, contourXY_id);
  } else if (PyObject_TypeCheck(tool, &cutting_tool_pinion_type)) { // internal gear
    double m_n     = ((cutting_tool_pinion*)tool)->m_n;
    double h_aP0   = m_n * ((cutting_tool_pinion*)tool)->h_aP0;
    double alfa_P0 = PI/180 * ((cutting_tool_pinion*)tool)->alfa_P0;
    double rho_aP0 = m_n * ((cutting_tool_pinion*)tool)->rho_aP0;
    double z_0     = ((cutting_tool_pinion*)tool)->z_0;
    double r_a0    = 0.5 * ((cutting_tool_pinion*)tool)->d_a0;
    generate_internal_gear_tooth_profile
      (m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a, // tooth and gear parameters
       z_0, r_a0,                                     // tool specific parameters
       N_sim, N_profile, N_tip,                       // simulation parameters
       tool_contours, contourXY, contourXY_id);
//  } else if (PyObject_TypeCheck(tool, &cutting_tool_polyline_type)) { // polyline
//    read_polygon_tool_contour(m_n, beta, fs_in, tool_contours[0], h_aP0);
//    generate_profile(ToolType::POLYLINE, m_n, beta, h_aP0, z, x, r_a, 0, 0., 0.,
//                     N_sim, N_profile, N_tip, tool_contours,
//                     contourXY, contourXY_id);
  } else {
    std::string errmsg("Unknown tool type. Expected: cutting_tool_rack, "
                       "cutting_tool_rack_with_protuberance or cutting_tool_pinion");
    PyErr_SetString(PyExc_ValueError, errmsg.c_str());
    return NULL;
  }
  // need to allocate the memory for the returned objects
  double *X = new double[contourXY.size()];
  double *Y = new double[contourXY.size()];
  int *ID = new int[contourXY.size()];

  for (int it=0; it < int(contourXY.size()); ++it) {
    X[it] = contourXY[it][0];
    Y[it] = contourXY[it][1];
    ID[it] = contourXY_id[it];
  }

  npy_intp dims[] = { int(contourXY.size()) };
  PyObject * Xarray = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, X);
  PyObject * Yarray = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, Y);
  PyObject * IDarray =PyArray_SimpleNewFromData(1, dims, NPY_INT, ID);
 return PyTuple_Pack(3, Xarray, Yarray, IDarray);
}


// Methods of the module gear_cutter
static PyMethodDef gear_cutter_methods[] = {
    // tooth profile generation function
    {"generate_tooth_profile",
     (PyCFunction)gear_cutter_generate_tooth_profile,
     METH_VARARGS | METH_KEYWORDS,
     "Simulates the gear cutting process and returns an array"
     " of points describing the contour of half tooth\n"
     "The input arguments are expected in this order:\n"
     "  tool: the object describing the cutting tool,\n"
     "        it should be one of\n"
     "        - cutting_tool_rack\n"
     "        - cutting_tool_rack_with_protuberance\n"
     "        - cutting_tool_pinion\n"
     "  z: the number of teeth of the gear\n"
     "  m_n: the normal module of the gear teeth\n"
     "  beta: the helix angle of the gear teeth\n"
     "  x: the profile shift of the gear\n"
     "  d_a: the tip diameter of the gear\n"
     "  N_sim: the number of points to be generated in the simulation\n"
     "  N_profile: the number of points to be generated for the tooth profile\n"
     "  N_tip: the number of points to be generated for the tip arc\n"
     "The function returns a tuple of three arrays\n"
     "containing the x, y coordinates and the id of the points of the tooth contour"},
    {NULL, NULL, 0, NULL}};



static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "gear_cutter",       /* m_name */
  "Module for simulation of the cutting of involute gears",  /* m_doc */
  -1,                     /* m_size */
  gear_cutter_methods,    /* m_methods */
  NULL,                   /* m_reload */
  NULL,                   /* m_traverse */
  NULL,                   /* m_clear */
  NULL,                   /* m_free */
};


// "PyInit_" is the prefix and "gear_cutter" is the module name
PyMODINIT_FUNC
PyInit_gear_cutter()
{
  import_array(); // needs to be called in order to work with numpy arrays
  PyObject *mod = PyModule_Create(&moduledef);
  if (mod == NULL)
    return NULL;

  if (PyType_Ready(&cutting_tool_rack_type) < 0) {
    Py_DECREF(mod);
    return NULL;
  }
  Py_INCREF(&cutting_tool_rack_type);
  PyModule_AddObject(mod, "cutting_tool_rack",
                     (PyObject *)&cutting_tool_rack_type);

  if (PyType_Ready(&cutting_tool_rack_with_protuberance_type) < 0) {
    Py_DECREF(mod);
    return NULL;
  }
  Py_INCREF(&cutting_tool_rack_with_protuberance_type);
  PyModule_AddObject(mod, "cutting_tool_rack_with_protuberance",
                     (PyObject *)&cutting_tool_rack_with_protuberance_type);

  if (PyType_Ready(&cutting_tool_pinion_type) < 0) {
    Py_DECREF(mod);
    return NULL;
  }
  Py_INCREF(&cutting_tool_pinion_type);
  PyModule_AddObject(mod, "cutting_tool_pinion",
                     (PyObject *)&cutting_tool_pinion_type);

  return mod;
}

