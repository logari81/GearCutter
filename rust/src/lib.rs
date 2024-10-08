use std::f64::consts::PI;
use std::vec::Vec;
use std::fs::File;
use std::io::{BufReader};
use std::io::Write;
use std::io::Lines;

// Structures

#[derive(Default,Clone)]
struct Arc {
    x_tool_m: f64,
    y_tool_m: f64,
    fi_tool_1: f64,
    fi_tool_2: f64,
    rho: f64,
    x_part_m: f64,
    y_part_m: f64,
    fi_part_1: f64,
    fi_part_2: f64,
    r_m: f64,
    fi_m: f64,
}

#[derive(Default,Clone)]
struct LineSeg {
    x_tool_1: f64,
    y_tool_1: f64,
    x_tool_2: f64,
    y_tool_2: f64,
    fi_tool_12: f64,
    x_part_1: f64,
    y_part_1: f64,
    x_part_2: f64,
    y_part_2: f64,
    fi_part_12: f64,
    r1: f64,
    r2: f64,
}

#[derive(Clone)]
pub struct Segment {
    id: i32,
    arc: Option<Arc>,
    line: Option<LineSeg>,
}

#[derive(PartialEq)]
pub enum ToolType {
    External = 0,
    ExternalWithProtuberance = 1,
    Internal = 2,
    Polyline = 3,
}



#[derive(Clone)]
struct CSpos {
    x: f64,
    y: f64,
    fi: f64,
}

#[derive(Default,Clone)]
struct ContourPoint {
    id_min: i32,
    id_max: i32,
    r: f64,
    theta_min: f64,
    theta_max: f64,
    x_min: f64,
    y_min: f64,
    x_max: f64,
    y_max: f64,
}

// Procedures
fn tool_to_part_coords(tool_pos: &CSpos, part_pos: &CSpos, x_tool: f64, y_tool: f64) -> (f64, f64) {
    let cc = tool_pos.fi.cos();
    let ss = tool_pos.fi.sin();

    let x_orig = x_tool * cc - y_tool * ss + tool_pos.x - part_pos.x;
    let y_orig = x_tool * ss + y_tool * cc + tool_pos.y - part_pos.y;

    let cc = (-part_pos.fi).cos();
    let ss = (-part_pos.fi).sin();

    let x_part = x_orig * cc - y_orig * ss;
    let y_part = x_orig * ss + y_orig * cc;

    (x_part, y_part)
}

fn normalized_angle(angle: f64) -> f64 {
    let mut normalized = angle;
    while normalized < 0.0 {
        normalized += 2.0 * PI;
    }
    while normalized > 2.0 * PI {
        normalized -= 2.0 * PI;
    }
    normalized
}

fn normalize_angles(angle1: &mut f64, angle2: &mut f64) {
    while *angle1 < 0.0 {
        *angle1 += 2.0 * PI;
        *angle2 += 2.0 * PI;
    }
    while *angle1 > 2.0 * PI {
        *angle1 -= 2.0 * PI;
        *angle2 -= 2.0 * PI;
    }
}

fn norm(a: f64, b: f64) -> f64 {
    (a * a + b * b).sqrt()
}

#[allow(non_snake_case)]
fn tool_params_internal_gear(
    z: i32,
    z_0: i32,
    m_t: f64,
    alfa_t: f64,
    xm_n: f64,
    r_a0: f64,
    rho_a0: f64,
    h_f1: f64,
    s_P0: f64,
) -> (f64, f64, f64, f64) {
    let cos_alfa_t = alfa_t.cos();
    let sin_alfa_t = alfa_t.sin();
    let tan_alfa_t = alfa_t.tan();
    let AB = (z.abs() as f64) * m_t * cos_alfa_t / 2.0;
    let DE = z_0 as f64 * m_t * cos_alfa_t / 2.0;
    let AC = (z.abs() as f64) * m_t / 2.0;
    let DF = r_a0 - rho_a0;
    let EF = (DF * DF - DE * DE).sqrt();
    let EG = EF + rho_a0;
    let AK = AC - xm_n + h_f1;
    let AD = AK - r_a0;
    let a = AD;
    let BE = (AD * AD - (AB - DE) * (AB - DE)).sqrt();
    let BC = AC * sin_alfa_t;
    let EC = BC - BE;
    let CG = EG - EC;
    let h_FaP0 = CG * sin_alfa_t + xm_n;
    let r_Ff1 = norm(AB, BE + EG);
    let HI = CG * cos_alfa_t + h_FaP0 * tan_alfa_t + s_P0 / 2.0;
    let theta = (AB - DE).acos() - alfa_t;
    let phi = DE.acos() - alfa_t - theta;
    let phi0 = phi - (HI / AC - theta) * (z.abs() as f64) / (z_0 as f64);
    let x_m = DF * phi0.sin();
    let y_m = DF * phi0.cos();
    (r_Ff1, a, x_m, y_m)
}

fn cutting_process(
    tool_contour: &Vec<Segment>,
    tool_pos: &Vec<CSpos>,
    part_pos: &Vec<CSpos>,
    pts: &mut Vec<ContourPoint>,
    mode: i32,
) {
    let n_sim = tool_pos.len();
    let n_profile = pts.len();

     // copy so that it can be modified for each position during the simulation
    let mut tool_contour_copy = tool_contour.clone();
    for it in 0..n_sim {
        let dfi = tool_pos[it].fi - part_pos[it].fi;

        for cur_segm in tool_contour_copy.iter_mut() {
            if let Some(cur_arc) = &mut cur_segm.arc {
                (cur_arc.x_part_m, cur_arc.y_part_m)
                    = tool_to_part_coords(&tool_pos[it],
                                          &part_pos[it],
                                          cur_arc.x_tool_m,
                                          cur_arc.y_tool_m);
                cur_arc.r_m = norm(cur_arc.x_part_m, cur_arc.y_part_m);
                cur_arc.fi_m = f64::atan2(cur_arc.y_part_m, cur_arc.x_part_m);
                cur_arc.fi_part_1 = cur_arc.fi_tool_1 + dfi;
                cur_arc.fi_part_2 = cur_arc.fi_tool_2 + dfi;
                normalize_angles(&mut cur_arc.fi_part_1, &mut cur_arc.fi_part_2);
            } else if let Some(cur_line) = &mut cur_segm.line {
                (cur_line.x_part_1, cur_line.y_part_1)
                    = tool_to_part_coords(&tool_pos[it],
                                          &part_pos[it],
                                          cur_line.x_tool_1,
                                          cur_line.y_tool_1);
                (cur_line.x_part_2, cur_line.y_part_2)
                    = tool_to_part_coords(&tool_pos[it],
                                          &part_pos[it],
                                          cur_line.x_tool_2,
                                          cur_line.y_tool_2);
                cur_line.r1 = norm(cur_line.x_part_1, cur_line.y_part_1);
                cur_line.r2 = norm(cur_line.x_part_2, cur_line.y_part_2);
                cur_line.fi_part_12 = normalized_angle(cur_line.fi_tool_12 + dfi);
            }
        }

        for it1 in 0..n_profile {
            let r = pts[it1].r;
            let mut id_min = 0;
            let mut id_max = 0;
            let mut x_min = 0.0;
            let mut y_min = 0.0;
            let mut x_max = 0.0;
            let mut y_max = 0.0;
            let mut theta_min = 1e10;
            let mut theta_max = -1e10;

            for cur_segm in tool_contour_copy.iter() {
                #[allow(non_snake_case)]
                if let Some(cur_arc) = &cur_segm.arc {
                    let RCS = cur_arc.rho;
                    let ROC = cur_arc.r_m;
                    let fiCO = cur_arc.fi_m + PI;
                    let cc = (RCS * RCS + ROC * ROC - r * r) / (2.0 * ROC * RCS);

                    if cc.abs() <= 1.0 {
                        for root in [fiCO - cc.acos(), fiCO + cc.acos()].iter() {
                            let dfi = normalized_angle(root - cur_arc.fi_part_1);

                            if dfi > 0.0 && dfi < cur_arc.fi_part_2 - cur_arc.fi_part_1 {
                                let my_x = cur_arc.x_part_m + RCS * root.cos();
                                let my_y = cur_arc.y_part_m + RCS * root.sin();
                                let my_theta = normalized_angle(f64::atan2(my_y,my_x));

                                if my_theta < theta_min {
                                    id_min = cur_segm.id;
                                    theta_min = my_theta;
                                    x_min = my_x;
                                    y_min = my_y;
                                }
                                if my_theta > theta_max {
                                    id_max = cur_segm.id;
                                    theta_max = my_theta;
                                    x_max = my_x;
                                    y_max = my_y;
                                }
                            }
                        }
                    }
                } else if let Some(cur_line) = &cur_segm.line {
                    let x1 = cur_line.x_part_1;
                    let y1 = cur_line.y_part_1;
                    let x2 = cur_line.x_part_2;
                    let y2 = cur_line.y_part_2;
                    let dx = x2 - x1;
                    let dy = y2 - y1;
                    let l = norm(dx, dy);
                    let tmp = (x1 * dy - y1 * dx) / (l * l);
                    let x0 = tmp * dy;
                    let y0 = -tmp * dx;
                    let d = norm(x0, y0);

                    let lambda1;
                    let lambda2;

                    if dx.abs() > dy.abs() {
                        lambda1 = (x1 - x0) / dx;
                        lambda2 = (x2 - x0) / dx;
                    } else {
                        lambda1 = (y1 - y0) / dy;
                        lambda2 = (y2 - y0) / dy;
                    }

                    let minlim = lambda1.min(lambda2) * l;
                    let maxlim = lambda1.max(lambda2) * l;

                    if r >= d - 1e-12 {
                        let lambda = (r * r - d * d).max(0.0).sqrt();

                        for &root in [-lambda, lambda].iter() {
                            if root < maxlim && root > minlim {
                                let my_x = x0 + root * dx / l;
                                let my_y = y0 + root * dy / l;
                                let my_theta = normalized_angle(f64::atan2(my_y,my_x));

                                if my_theta < theta_min {
                                    id_min = cur_segm.id;
                                    theta_min = my_theta;
                                    x_min = my_x;
                                    y_min = my_y;
                                }
                                if my_theta > theta_max {
                                    id_max = cur_segm.id;
                                    theta_max = my_theta;
                                    x_max = my_x;
                                    y_max = my_y;
                                }
                            }
                        }
                    }
                }
            }

            if mode == 0 && theta_min <= theta_max {
                if theta_min < pts[it1].theta_min {
                    pts[it1].id_min = id_min;
                    pts[it1].theta_min = theta_min;
                    pts[it1].x_min = x_min;
                    pts[it1].y_min = y_min;
                }
                if theta_max > pts[it1].theta_max {
                    pts[it1].id_max = id_max;
                    pts[it1].theta_max = theta_max;
                    pts[it1].x_max = x_max;
                    pts[it1].y_max = y_max;
                }
            } else if mode == -1 && theta_min <= theta_max {
                if theta_max > pts[it1].theta_min {
                    pts[it1].id_min = id_max;
                    pts[it1].theta_min = theta_max;
                    pts[it1].x_min = x_max;
                    pts[it1].y_min = y_max;
                }
            } else if mode == 1 && theta_min <= theta_max {
                if theta_min < pts[it1].theta_max {
                    pts[it1].id_max = id_min;
                    pts[it1].theta_max = theta_min;
                    pts[it1].x_max = x_min;
                    pts[it1].y_max = y_min;
                }
            }
        }
    }
}



#[allow(non_snake_case)]
fn build_tool_contour(
    m_n: f64,
    beta: f64,
    h_aP0: f64,
    alfa_P0: f64,
    rho_aP0: f64,
    h_fP0: f64,
    contour: &mut Vec<Segment>,
) {
    // Account for helix angle
    let m_t = m_n / beta.cos();
    let alfa_P0 = (alfa_P0.tan() / beta.cos()).atan();
    // The impact of the helix angle on rho_aP0 is neglected

    // Auxiliary variables
    let s_P0 = PI * m_t / 2.0;
    let h_FaP0 = h_aP0 - (1.0 - alfa_P0.sin()) * rho_aP0;

    // Build tool contour
    // 1st segment (line)
    contour.push(Segment { id: 1, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 - h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = 2.0 * h_fP0;
        cur_line.x_tool_2 = cur_line.x_tool_1;
        cur_line.y_tool_2 = h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 2nd segment (line)
    contour.push(Segment { id: 2, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 - h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = h_fP0;
        cur_line.x_tool_2 = -s_P0 / 2.0 + h_FaP0 * alfa_P0.tan();
        cur_line.y_tool_2 = -h_FaP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 3rd segment (arc)
    contour.push(Segment { id: 3, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contour.last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = -s_P0 / 2.0 + h_FaP0 * alfa_P0.tan() + rho_aP0 * alfa_P0.cos();
        cur_arc.y_tool_m = -h_aP0 + rho_aP0;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = PI + alfa_P0;
        cur_arc.fi_tool_2 = 3.0 * PI / 2.0;
    }
    // 4th segment (line)
    contour.push(Segment { id: 4, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 + h_FaP0 * alfa_P0.tan() + rho_aP0 * alfa_P0.cos();
        cur_line.y_tool_1 = -h_aP0;
        cur_line.x_tool_2 = -cur_line.x_tool_1;
        cur_line.y_tool_2 = cur_line.y_tool_1;
        cur_line.fi_tool_12 = 0.0;
    }
    // 5th segment (arc)
    contour.push(Segment { id: 5, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contour.last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = s_P0 / 2.0 - h_FaP0 * alfa_P0.tan() - rho_aP0 * alfa_P0.cos();
        cur_arc.y_tool_m = -h_aP0 + rho_aP0;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = 3.0 * PI / 2.0;
        cur_arc.fi_tool_2 = 2.0 * PI - alfa_P0;
    }
    // 6th segment (line)
    contour.push(Segment { id: 6, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = s_P0 / 2.0 - h_FaP0 * alfa_P0.tan();
        cur_line.y_tool_1 = -h_FaP0;
        cur_line.x_tool_2 = s_P0 / 2.0 + h_fP0 * alfa_P0.tan();
        cur_line.y_tool_2 = h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 7th segment (line)
    contour.push(Segment { id: 7, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = s_P0 / 2.0 + h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = h_fP0;
        cur_line.x_tool_2 = s_P0 / 2.0 + h_fP0 * alfa_P0.tan();
        cur_line.y_tool_2 = 2.0 * h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
}


#[allow(non_snake_case)]
fn build_tool_with_protuberance_contour(
    m_n: f64,
    beta: f64,
    h_prP0: f64,
    alfa_P0: f64,
    rho_aP0: f64,
    h_fP0: f64,
    h_FaP0: f64,
    _alfa_KP0: f64,
    alfa_prP0: f64,
    contour: &mut Vec<Segment>,
) {
    // Account for helix angle
    let m_t = m_n / beta.cos();
    let alfa_P0 = (alfa_P0.tan() / beta.cos()).atan();
    //let alfa_KP0 = (alfa_KP0.tan() / beta.cos()).atan();
    let alfa_prP0 = (alfa_prP0.tan() / beta.cos()).atan();
    // The impact of the helix angle on rho_aP0 is neglected

    // Auxiliary variables
    let s_P0 = PI * m_t / 2.0;
    let h_FprP0 = h_prP0 - (1.0 - alfa_prP0.sin()) * rho_aP0;
    let x_orig = -s_P0 / 2.0 + h_FaP0 * alfa_P0.tan();

    // Build tool contour
    // 1st segment (line)
    contour.push(Segment { id: 1, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 - h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = 2.0 * h_fP0;
        cur_line.x_tool_2 = cur_line.x_tool_1;
        cur_line.y_tool_2 = h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 2nd segment (line)
    contour.push(Segment { id: 2, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 - h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = h_fP0;
        cur_line.x_tool_2 = -s_P0 / 2.0 + h_FaP0 * alfa_P0.tan();
        cur_line.y_tool_2 = -h_FaP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 3rd segment (line)
    contour.push(Segment { id: 3, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = x_orig;
        cur_line.y_tool_1 = -h_FaP0;
        cur_line.x_tool_2 = x_orig + h_FprP0 * alfa_prP0.tan();
        cur_line.y_tool_2 = -h_FaP0 - h_FprP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 4th segment (arc)
    contour.push(Segment { id: 4, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contour.last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = x_orig + h_FprP0 * alfa_prP0.tan() + rho_aP0 * alfa_prP0.cos();
        cur_arc.y_tool_m = -h_FaP0 - h_prP0 + rho_aP0;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = PI + alfa_P0;
        cur_arc.fi_tool_2 = 3.0 * PI / 2.0;
    }
    // 5th segment (line)
    contour.push(Segment { id: 5, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = x_orig + h_FprP0 * alfa_prP0.tan() + rho_aP0 * alfa_prP0.cos();
        cur_line.y_tool_1 = -h_FaP0 - h_prP0;
        cur_line.x_tool_2 = -cur_line.x_tool_1;
        cur_line.y_tool_2 = cur_line.y_tool_1;
        cur_line.fi_tool_12 = 0.0;
    }
    // 6th segment (arc)
    contour.push(Segment { id: 6, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contour.last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = -x_orig - h_FprP0 * alfa_prP0.tan() - rho_aP0 * alfa_prP0.cos();
        cur_arc.y_tool_m = -h_FaP0 - h_prP0 + rho_aP0;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = 3.0 * PI / 2.0;
        cur_arc.fi_tool_2 = 2.0 * PI - alfa_prP0;
    }
    // 7th segment (line)
    contour.push(Segment { id: 7, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -x_orig;
        cur_line.y_tool_1 = -h_FaP0;
        cur_line.x_tool_2 = -x_orig - h_FprP0 * alfa_prP0.tan();
        cur_line.y_tool_2 = -h_FaP0 - h_FprP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 8th segment (line)
    contour.push(Segment { id: 8, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = s_P0 / 2.0 + h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = h_fP0;
        cur_line.x_tool_2 = s_P0 / 2.0 - h_FaP0 * alfa_P0.tan();
        cur_line.y_tool_2 = -h_FaP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
    // 9th segment (line)
    contour.push(Segment { id: 9, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contour.last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = s_P0 / 2.0 + h_fP0 * alfa_P0.tan();
        cur_line.y_tool_1 = 2.0 * h_fP0;
        cur_line.x_tool_2 = cur_line.x_tool_1;
        cur_line.y_tool_2 = h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }
}


#[allow(non_snake_case)]
fn build_internal_gear_tool_contours(
    m_n: f64,
    beta: f64,
    alfa_P0: f64,
    h_aP0: f64,
    rho_aP0: f64,
    x: f64,
    z: i32,
    z_0: i32,
    r_a0: f64,
    contours: &mut Vec<Vec<Segment>>,
) -> (f64, f64) {
    // Account for helix angle
    let m_t = m_n / beta.cos();
    let alfa_P0 = (alfa_P0.to_radians().tan() / beta.cos()).atan();
    // The impact of the helix angle on rho_aP0 is neglected

    // Auxiliary variables
    let s_P0 = PI * m_t / 2.0;
    let h_fP0 = s_P0 / (2.0 * alfa_P0.tan());
    let (r_Ff, a_tool, x_m, y_m) = tool_params_internal_gear(z,
                                                             z_0,
                                                             m_t,
                                                             alfa_P0,
                                                             x * m_n,
                                                             r_a0,
                                                             rho_aP0,
                                                             h_aP0,
                                                             s_P0);

    // Define a contour for the auxiliary tool (rack) that will generate
    // the actual cutting tool (wheel)
    contours.clear();
    contours.resize(3, Vec::<Segment>::new());

    // 1st segment (line)
    contours[0].push(Segment { id: 1, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contours[0].last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = -s_P0 / 2.0 - 2.0 * h_aP0 * alfa_P0.tan();
        cur_line.y_tool_1 = 2.0 * h_aP0;
        cur_line.x_tool_2 = 0.0;
        cur_line.y_tool_2 = -h_fP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }

    // 2nd segment (line)
    contours[0].push(Segment { id: 2, arc: None, line: Some(LineSeg::default()) });
    {
        let cur_line = &mut contours[0].last_mut().unwrap().line.as_mut().unwrap();
        cur_line.x_tool_1 = 0.0;
        cur_line.y_tool_1 = -h_fP0;
        cur_line.x_tool_2 = s_P0 / 2.0 + 2.0 * h_aP0 * alfa_P0.tan();
        cur_line.y_tool_2 = 2.0 * h_aP0;
        cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                         cur_line.x_tool_2 - cur_line.x_tool_1);
    }

    // 3rd segment (arc)
    contours[1].push(Segment { id: 3, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contours[1].last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = x_m;
        cur_arc.y_tool_m = y_m;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = 0.0;
        cur_arc.fi_tool_2 = f64::atan2(y_m,x_m);
    }

    // 4th segment (arc)
    contours[2].push(Segment { id: 4, arc: Some(Arc::default()), line: None });
    {
        let cur_arc = &mut contours[2].last_mut().unwrap().arc.as_mut().unwrap();
        cur_arc.x_tool_m = -x_m;
        cur_arc.y_tool_m = y_m;
        cur_arc.rho = rho_aP0;
        cur_arc.fi_tool_1 = f64::atan2(y_m,x_m);
        cur_arc.fi_tool_2 = PI;
    }
    (r_Ff, a_tool)
}



#[allow(non_snake_case)]
pub fn read_polygon_tool_contour(
    m_n: &mut f64,
    beta: f64,
    lines: &mut Lines<BufReader<File>>,
    contour: &mut Vec<Segment>,
    h_aP0: &mut f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let next_line = lines.next().unwrap_or_else(|| {panic!("Unexpected end of input file");})?;
    if next_line.starts_with("m_n") {
        if let Some(valstr) = next_line.split(":").last().map(|s| s.trim()) {
            *m_n = valstr.parse().unwrap_or_else(|e| {
                   panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                          next_line, e);
                   });
        } else {
            panic!("Syntax error in input file, m_n entry does not contain a colon: {}",
                   next_line);
        }
    } else {
        panic!("Syntax error in input file, expected variable m_n not found");
    }

    *h_aP0 = 0.0;
    contour.clear();
    let mut it = 0;
    let mut x0 = 0.0;
    let mut y0 = 0.0;
        
    for line in lines {
        if let Ok(line_str) = line {
            let mut values = line_str.split_whitespace();
            if let (Some(x_str), Some(y_str)) = (values.next(), values.next()) {
                if let (Ok(x), Ok(y)) = (x_str.parse::<f64>(), y_str.parse::<f64>()) {
                    contour.push(Segment { id: it + 1, arc: None, line: Some(LineSeg::default()) });
                    let cur_line = contour.last_mut().unwrap().line.as_mut().unwrap();
                    cur_line.x_tool_1 = x0;
                    cur_line.y_tool_1 = y0 / beta.cos();
                    cur_line.x_tool_2 = x;
                    cur_line.y_tool_2 = y / beta.cos();
                    cur_line.fi_tool_12 = f64::atan2(cur_line.y_tool_2 - cur_line.y_tool_1,
                                                     cur_line.x_tool_2 - cur_line.x_tool_1);
                    x0 = x;
                    y0 = y;
                    if *h_aP0 < -y {
                        *h_aP0 = -y;
                    }
                    it += 1;
                }
            }
        }
    }
    return Ok(());
}


pub fn write_to_file(fname: &str, contour: &Vec<Segment>) {
    if let Ok(mut fs) = File::create(fname) {
        for cur_segm in contour {
            if cur_segm.arc.is_some() {
                let cur_arc = cur_segm.arc.as_ref().unwrap();
                write!(
                    fs,
                    "Arc\nx_M={:.5} y_M={:.5} rho={:.5}",
                    cur_arc.x_tool_m, cur_arc.y_tool_m, cur_arc.rho
                )
                .unwrap();
            } else if cur_segm.line.is_some() {
                let cur_line = cur_segm.line.as_ref().unwrap();
                write!(
                    fs,
                    "LineSeg\nx_A={:.5} y_A={:.5} x_E={:.5} y_E={:.5}",
                    cur_line.x_tool_1, cur_line.y_tool_1, cur_line.x_tool_2, cur_line.y_tool_2
                )
                .unwrap();
            }
            write!(fs, "\n").unwrap();
        }
    }
}


#[allow(non_snake_case)]
pub fn generate_profile(
    tool_type: ToolType,
    m_n: f64,
    beta: f64,
    h_aP0: f64,
    z: i32,
    x: f64,
    r_a: f64,
    z_0: i32,
    r_Ff: f64,
    a_tool: f64,
    n_sim: usize,
    n_profile: usize,
    n_tip: usize,
    tool_contours: &Vec<Vec<Segment>>,
    contour_xy: &mut Vec<Vec<f64>>,
    contour_xy_id: &mut Vec<i32>,
) {
    let r_f: f64;
    let mut pts: Vec<ContourPoint> = vec![ContourPoint::default(); n_profile];

    let r_w = (z.abs() as f64) * m_n / beta.cos() / 2.0;
    r_f = r_w + if z > 0 { x * m_n - h_aP0 } else { -x * m_n - h_aP0 };

    let mut tool_pos: Vec<CSpos> = vec![CSpos {
        x: 0.0,
        y: r_w + if z > 0 { x } else { -x } * m_n,
        fi: 0.0,
    }; n_sim];
    let mut part_pos: Vec<CSpos> = vec![CSpos {
        x: 0.0,
        y: 0.0,
        fi: 0.0,
    }; n_sim];

    for it in 0..n_sim {
        part_pos[it].fi = PI * (0.5 * (it as f64) / (n_sim as f64 - 1.0) - 0.25);
        tool_pos[it].x = -r_w * part_pos[it].fi;
    }

    let incr = (r_a - r_f) / (n_profile - 1) as f64;
    for it in 0..n_profile {
        pts[it].r = r_f + (it as f64) * incr;
        pts[it].theta_min = 1e8;
        pts[it].theta_max = -1e8;
    }

    cutting_process(&tool_contours[0], &tool_pos, &part_pos, &mut pts, 0);

    if tool_type == ToolType::Internal {
        for it in 0..n_sim {
            tool_pos[it].x = 0.0;
            tool_pos[it].y = a_tool;
            tool_pos[it].fi = (z.abs() as f64) * part_pos[it].fi / (z_0 as f64);
        }

        for pt in pts.iter_mut() {
            if pt.r <= r_Ff {
                break;
            }
            pt.theta_min = -1e8;
            pt.theta_max = 1e8;
            pt.id_min = -1;
            pt.id_max = -1;
        }

        cutting_process(&tool_contours[1], &tool_pos, &part_pos, &mut pts, 1);
        cutting_process(&tool_contours[2], &tool_pos, &part_pos, &mut pts, -1);
    }

    let mut n_start = 0 as usize;
    let n_root : usize;

    if z > 0 {
        while pts[n_start].id_max == -1 {
            n_start += 1;
        }

        n_root = ((pts[n_start].theta_max - PI / 2.0) /
            (pts[n_start + 1].theta_max - pts[n_start].theta_max)).ceil() as usize;
    } else {
        while pts[n_start].theta_min < (PI / 2.0 - PI / z.abs() as f64) ||
              pts[n_start].id_min == -1 {
            n_start += 1;
        }

        n_root = ((pts[n_start].theta_min - (PI / 2.0 - PI / z.abs() as f64)) /
            (pts[n_start + 1].theta_min - pts[n_start].theta_min)).ceil() as usize;
    }

    contour_xy.clear();
    contour_xy_id.clear();
    contour_xy.resize(n_root + n_tip + n_profile - n_start, vec![0 as f64; 2]);
    contour_xy_id.resize(n_root + n_tip + n_profile - n_start, 0 as i32);
    let mut n_contour = 0 as usize;

    if n_root > 0 {
        let incr = if z > 0 {
            (pts[n_start].theta_max - PI / 2.0) / (n_root as f64)
        } else {
            (pts[n_start].theta_min - (PI / 2.0 - PI / z.abs() as f64)) / (n_root as f64)
        };

        for it in 0..n_root {
            contour_xy_id[it] = -1;
            let my_theta = PI / 2.0 - PI / z.abs() as f64 + incr * (it as f64);
            contour_xy[it][0] = r_f * my_theta.cos();
            contour_xy[it][1] = r_f * my_theta.sin();
        }

        n_contour = n_root;
    }

    if z > 0 {
        let dx = (PI / z.abs() as f64).sin();
        let dy = (PI / z.abs() as f64).cos();

        for it in n_start..n_profile {
            if pts[it].theta_min < pts[it].theta_max {
                contour_xy_id[n_contour] = pts[it].id_max;
                contour_xy[n_contour][0] = pts[it].x_max * dy + pts[it].y_max * dx;
                contour_xy[n_contour][1] = -pts[it].x_max * dx + pts[it].y_max * dy;
                n_contour += 1;
            } else {
                println!("Error in the calculated profile at r = {}", pts[it].r);
            }
        }
    } else {
        for it in n_start..n_profile {
            if pts[it].theta_min < pts[it].theta_max {
                contour_xy_id[n_contour] = pts[it].id_min;
                contour_xy[n_contour][0] = pts[it].x_min;
                contour_xy[n_contour][1] = pts[it].y_min;
                n_contour += 1;
            } else {
                println!("Error in the calculated profile at r = {}", pts[it].r);
            }
        }
    }

    let r = pts[n_profile - 1].r;
    let my_theta = if z > 0 {
        pts[n_profile - 1].theta_max - PI / (z as f64)
    } else {
        pts[n_profile - 1].theta_min
    };
    let incr = (PI / 2.0 - my_theta) / (n_tip as f64);

    for it in 0..n_tip {
        contour_xy_id[n_contour + it] = -2;
        contour_xy[n_contour + it][0] = r * (my_theta + incr * ((it+1) as f64)).cos();
        contour_xy[n_contour + it][1] = r * (my_theta + incr * ((it+1) as f64)).sin();
    }
}




#[allow(non_snake_case)]
pub fn generate_external_gear_tooth_profile(
    m_n: f64,
    beta: f64,
    h_aP0: f64,
    alfa_P0: f64,
    rho_aP0: f64,
    z: i32,
    x: f64,
    r_a: f64,
    h_fP0: f64,
    n_sim: usize,
    n_profile: usize,
    n_tip: usize,
    tool_contours: &mut Vec<Vec<Segment>>,
    contour_xy: &mut Vec<Vec<f64>>,
    contour_xy_id: &mut Vec<i32>,
) {
    build_tool_contour(
        m_n,
        beta,
        h_aP0,
        alfa_P0,
        rho_aP0,
        h_fP0,
        &mut tool_contours[0],
    );

    generate_profile(
        ToolType::External,
        m_n,
        beta,
        h_aP0,
        z,
        x,
        r_a,
        0,
        0.0,
        0.0,
        n_sim,
        n_profile,
        n_tip,
        &tool_contours,
        contour_xy,
        contour_xy_id,
    );
}

#[allow(non_snake_case)]
pub fn generate_external_gear_tooth_profile_with_protuberance(
    m_n: f64,
    beta: f64,
    h_aP0: f64,
    alfa_P0: f64,
    rho_aP0: f64,
    z: i32,
    x: f64,
    r_a: f64,
    h_fP0: f64,
    h_FaP0: f64,
    alfa_KP0: f64,
    alfa_prP0: f64,
    n_sim: usize,
    n_profile: usize,
    n_tip: usize,
    tool_contours: &mut Vec<Vec<Segment>>,
    contour_xy: &mut Vec<Vec<f64>>,
    contour_xy_id: &mut Vec<i32>,
) {
    let h_prP0 = h_aP0;
    build_tool_with_protuberance_contour(
        m_n,
        beta,
        h_prP0,
        alfa_P0,
        rho_aP0,
        h_fP0,
        h_FaP0,
        alfa_KP0,
        alfa_prP0,
        &mut tool_contours[0],
    );

    generate_profile(
        ToolType::ExternalWithProtuberance,
        m_n,
        beta,
        h_aP0,
        z,
        x,
        r_a,
        0,
        0.0,
        0.0,
        n_sim,
        n_profile,
        n_tip,
        &tool_contours,
        contour_xy,
        contour_xy_id,
    );
}

#[allow(non_snake_case)]
pub fn generate_internal_gear_tooth_profile(
    m_n: f64,
    beta: f64,
    h_aP0: f64,
    alfa_P0: f64,
    rho_aP0: f64,
    z: i32,
    x: f64,
    r_a: f64,
    z_0: i32,
    r_a0: f64,
    n_sim: usize,
    n_profile: usize,
    n_tip: usize,
    tool_contours: &mut Vec<Vec<Segment>>,
    contour_xy: &mut Vec<Vec<f64>>,
    contour_xy_id: &mut Vec<i32>,
) {
    let (r_Ff, a_tool) = build_internal_gear_tool_contours(
        m_n,
        beta,
        alfa_P0,
        h_aP0,
        rho_aP0,
        x,
        z,
        z_0,
        r_a0,
        tool_contours,
    );
    generate_profile(
        ToolType::Internal,
        m_n,
        beta,
        h_aP0,
        z,
        x,
        r_a,
        z_0,
        r_Ff,
        a_tool,
        n_sim,
        n_profile,
        n_tip,
        &tool_contours,
        contour_xy,
        contour_xy_id,
    );
}

