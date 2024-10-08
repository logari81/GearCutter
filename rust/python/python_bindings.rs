// --------------------------------------------------------------------------
//
// GearCutter 0.2
// Generation of involute gear tooth profiles
// Copyright (C) 2024-2024 Konstantinos Poulios
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

use pyo3::prelude::*;
use pyo3::types::PyAny;
use pyo3::wrap_pyfunction;
use gear_cutter::*;

#[pyclass]
#[allow(non_snake_case)]
struct CuttingToolRack {
    #[pyo3(get)]
    m_n: f64,
    #[pyo3(get)]
    h_aP0: f64,
    #[pyo3(get)]
    alfa_P0: f64,
    #[pyo3(get)]
    rho_aP0: f64,
    #[pyo3(get)]
    h_fP0: f64,
}

#[pymethods]
#[allow(non_snake_case)]
impl CuttingToolRack {
    #[new]
    #[pyo3(signature = (m_n, h_aP0=None, alfa_P0=None, rho_aP0=None, h_fP0=None))]
    fn new(m_n: f64,
           h_aP0: Option<f64>,
           alfa_P0: Option<f64>,
           rho_aP0: Option<f64>,
           h_fP0: Option<f64>
    ) -> Self {
        Self {
            m_n,
            h_aP0: h_aP0.unwrap_or(1.25),
            alfa_P0: alfa_P0.unwrap_or(20.),
            rho_aP0: rho_aP0.unwrap_or(0.38),
            h_fP0: h_fP0.unwrap_or(1.),
        }
    }
}

#[pyclass]
#[allow(non_snake_case)]
struct CuttingToolRackWithProtuberance {
    #[pyo3(get)]
    m_n: f64,
    #[pyo3(get)]
    h_prP0: f64,
    #[pyo3(get)]
    alfa_P0: f64,
    #[pyo3(get)]
    rho_aP0: f64,
    #[pyo3(get)]
    h_fP0: f64,
    #[pyo3(get)]
    h_FaP0: f64,
    #[pyo3(get)]
    alfa_KP0: f64,
    #[pyo3(get)]
    alfa_prP0: f64,
}

#[pymethods]
#[allow(non_snake_case)]
impl CuttingToolRackWithProtuberance {
    #[new]
    #[pyo3(signature = (m_n, h_prP0=None, alfa_P0=None, rho_aP0=None, h_fP0=None, h_FaP0=None, alfa_KP0=None, alfa_prP0=None))]
    fn new(
        m_n: f64,
        h_prP0: Option<f64>,
        alfa_P0: Option<f64>,
        rho_aP0: Option<f64>,
        h_fP0: Option<f64>,
        h_FaP0: Option<f64>,
        alfa_KP0: Option<f64>,
        alfa_prP0: Option<f64>,
    ) -> Self {
        Self {
            m_n,
            h_prP0: h_prP0.unwrap_or(1.25),
            alfa_P0: alfa_P0.unwrap_or(20.0),
            rho_aP0: rho_aP0.unwrap_or(0.38),
            h_fP0: h_fP0.unwrap_or(1.0),
            h_FaP0: h_FaP0.unwrap_or(0.8),
            alfa_KP0: alfa_KP0.unwrap_or(37.5),
            alfa_prP0: alfa_prP0.unwrap_or(12.5),
        }
    }
}

#[pyclass]
#[allow(non_snake_case)]
struct CuttingToolPinion {
    #[pyo3(get)]
    m_n: f64,
    #[pyo3(get)]
    h_aP0: f64,
    #[pyo3(get)]
    alfa_P0: f64,
    #[pyo3(get)]
    rho_aP0: f64,
    #[pyo3(get)]
    h_fP0: f64,
    #[pyo3(get)]
    z_0: i32,
    #[pyo3(get)]
    d_a0: f64,
}

#[pymethods]
#[allow(non_snake_case)]
impl CuttingToolPinion {
    #[new]
    fn new(
        m_n: f64,
        h_aP0: f64,
        alfa_P0: f64,
        rho_aP0: f64,
        h_fP0: f64,
        z_0: i32,
        d_a0: f64,
    ) -> Self {
        Self {
            m_n,
            h_aP0,
            alfa_P0,
            rho_aP0,
            h_fP0,
            z_0,
            d_a0,
        }
    }
}

#[pyfunction]
fn generate_tooth_profile(
    tool: &Bound<'_, PyAny>,
    z: i32,
    beta: f64,
    x: f64,
    d_a: f64,
    n_sim: usize,
    n_profile: usize,
    n_tip: usize,
) -> (Vec<f64>, Vec<f64>, Vec<i32>) {
    let mut tool_contours: Vec<Vec<Segment>> = Vec::new();
    tool_contours.push(vec![]);
    let mut contour_xy: Vec<Vec<f64>> = Vec::new();
    let mut contour_xy_id: Vec<i32> = Vec::new();
    if let Ok(tool1) = tool.extract::<PyRef<CuttingToolRack>>() {
        gear_cutter::generate_external_gear_tooth_profile(
            tool1.m_n,
            beta,
            tool1.h_aP0 * tool1.m_n,
            tool1.alfa_P0.to_radians(),
            tool1.rho_aP0 * tool1.m_n,
            z,
            x,
            d_a/2.0,
            tool1.h_fP0 * tool1.m_n,
            n_sim,
            n_profile,
            n_tip,
            &mut tool_contours,
            &mut contour_xy,
            &mut contour_xy_id,
        );
    } else if let Ok(tool2) = tool.extract::<PyRef<CuttingToolRackWithProtuberance>>() {
        gear_cutter::generate_external_gear_tooth_profile_with_protuberance(
            tool2.m_n,
            beta,
            tool2.h_prP0 * tool2.m_n,
            tool2.alfa_P0.to_radians(),
            tool2.rho_aP0 * tool2.m_n,
            z,
            x,
            d_a/2.0,
            tool2.h_fP0 * tool2.m_n,
            tool2.h_FaP0 * tool2.m_n,
            tool2.alfa_KP0.to_radians(),
            tool2.alfa_prP0.to_radians(),
            n_sim,
            n_profile,
            n_tip,
            &mut tool_contours,
            &mut contour_xy,
            &mut contour_xy_id,
        )
    } else if let Ok(tool3) = tool.extract::<PyRef<CuttingToolPinion>>() {
        gear_cutter::generate_internal_gear_tooth_profile(
            tool3.m_n,
            beta,
            tool3.h_aP0 * tool3.m_n,
            tool3.alfa_P0.to_radians(),
            tool3.rho_aP0 * tool3.m_n,
            z,
            x,
            d_a/2.0,
            tool3.z_0,
            tool3.d_a0 / 2.0,
            n_sim,
            n_profile,
            n_tip,
            &mut tool_contours,
            &mut contour_xy,
            &mut contour_xy_id,
        );
    } else {
        panic!("Unknown tool type. Expected: External, ExternalWithProtuberance or Internal");
    }
    let xvec: Vec<f64> = contour_xy.iter().map(|p| p[0]).collect();
    let yvec: Vec<f64> = contour_xy.iter().map(|p| p[1]).collect();
    (xvec, yvec, contour_xy_id)
}

#[pymodule]
fn gear_cutter_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(generate_tooth_profile, m)?)?;
    m.add_class::<CuttingToolRack>()?;
    m.add_class::<CuttingToolRackWithProtuberance>()?;
    m.add_class::<CuttingToolPinion>()?;
    Ok(())
}

