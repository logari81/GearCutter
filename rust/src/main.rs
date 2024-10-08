use std::f64::consts::PI;
use std::env;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::io::Result;
use std::string::String;
use std::fs::File;

//mod gear_cutter;
use gear_cutter::Segment;
use gear_cutter::ToolType;
use gear_cutter::read_polygon_tool_contour;
use gear_cutter::generate_profile;
use gear_cutter::generate_external_gear_tooth_profile;
use gear_cutter::generate_external_gear_tooth_profile_with_protuberance;
use gear_cutter::generate_internal_gear_tooth_profile;
use gear_cutter::write_to_file;

fn main() -> Result<()> {
    let mut n_sim = 1000;
    let mut n_profile = 1200;
    let mut n_tip = 7;
    let mut input_file = String::from("Gear_data.txt");
    let mut output_file = String::from("tooth_profile.dat");
    let mut output_file2 = String::from("tool_contour.dat.1");

    let args: Vec<String> = env::args().collect();
    let mut i = 1;

    while i < args.len() {
        let option = &args[i];
        let mut missing_argument = false;

        match option.as_str() {
            "-i" => {
                if i + 1 < args.len() {
                    input_file = args[i + 1].clone();
                    i += 1;
                } else {
                    missing_argument = true;
                }
            }
            "-o" => {
                if i + 1 < args.len() {
                    output_file = args[i + 1].clone();
                    output_file2 = format!("{}.1", &output_file);
                    i += 1;
                } else {
                    missing_argument = true;
                }
            }
            "-N_sim" => {
                if i + 1 < args.len() {
                    n_sim = args[i + 1].parse().unwrap();
                    i += 1;
                } else {
                    missing_argument = true;
                }
            }
            "-N_profile" => {
                if i + 1 < args.len() {
                    n_profile = args[i + 1].parse().unwrap();
                    i += 1;
                } else {
                    missing_argument = true;
                }
            }
            "-N_tip" => {
                if i + 1 < args.len() {
                    n_tip = args[i + 1].parse().unwrap();
                    i += 1;
                } else {
                    missing_argument = true;
                }
            }
            _ => {
                eprintln!("Unknown option {}", option);
                return Ok(());
            }
        }

        if missing_argument {
            eprintln!("{} option requires one argument.", option);
            return Ok(());
        }

        i += 1;
    }

    let mut tool_contours: Vec<Vec<Segment>> = Vec::new();
    tool_contours.push(vec![]);
    let mut contour_xy: Vec<Vec<f64>> = Vec::new();
    let mut contour_xy_id: Vec<i32> = Vec::new();
    {
        let fs_in = File::open(&input_file)?;
        let reader = BufReader::new(fs_in);
        let mut lines = reader.lines();
        let next_line = lines.next().unwrap()?;
        let itmp: i32 = next_line.trim().parse().unwrap();
        let tool_type = match itmp {
            0 => ToolType::External,
            1 => ToolType::ExternalWithProtuberance,
            2 => ToolType::Internal,
            3 => ToolType::Polyline,
            _ => {
                eprintln!("Invalid tool type.");
                return Ok(());
            }
        };

        let mut z: i32 = 0;
        let mut r_a: f64 = 0.;
        let mut x: f64 = 0.;
        let mut beta: f64 = 0.;
        for &name in ["z", "d_a", "x", "beta"].iter() {
            let next_line = lines.next().unwrap_or_else(|| {panic!("Unexpected end of input file");})?;
            if next_line.starts_with(name) {
                if let Some(valstr) = next_line.split(":").last().map(|s| s.trim()) {
                    if name == "z" { // i32
                        z = valstr.parse().unwrap_or_else(|e| {
                            panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                   next_line, e);
                            });
                    } else { // f64
                        let val: f64 =
                            valstr.parse().unwrap_or_else(|e| {
                            panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                   next_line, e);
                            });
                        match name {
                            "d_a"  => { r_a = val/2.; }
                            "x"    => { x = val; }
                            "beta" => { beta = val; }
                            &_ => { panic!("Should never happen, internal error"); }
                        }
                    }
                } else {
                    panic!("Syntax error in input file, {} entry does not contain a colon:\n{}",
                           name, next_line);
                }
            } else {
                panic!("Syntax error in input file, expected variable {} not found", name);
            }
        }
    
        #[allow(non_snake_case)]
        match tool_type {
            ToolType::External => {
                let mut alfa_P0: f64 = 0.;
                let mut m_n: f64 = 0.;
                let mut h_fP0: f64 = 0.;
                let mut h_aP0: f64 = 0.;
                let mut rho_aP0: f64 = 0.;
                for &name in ["alfa_P0", "m_n", "h_fP0", "h_aP0", "rho_aP0"].iter() {
                    let next_line = lines.next().unwrap_or_else(|| {panic!("Unexpected end of input file");})?;
                    if next_line.starts_with(name) {
                        if let Some(valstr) = next_line.split(":").last().map(|s| s.trim()) {
                            let val: f64 =
                                valstr.parse().unwrap_or_else(|e| {
                                panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                       next_line, e);
                                });
                            match name {
                                "alfa_P0" => { alfa_P0 = PI/180. * val; }
                                "m_n"     => { m_n = val; }
                                "h_fP0"   => { h_fP0 = m_n*val; }
                                "h_aP0"   => { h_aP0 = m_n*val; }
                                "rho_aP0" => { rho_aP0 = m_n*val; }
                                &_ => { panic!("Should never happen, internal error"); }
                            }
                        } else {
                            panic!("Syntax error in input file, {} entry does not contain a colon:\n{}",
                                   name, next_line);
                        }
                    } else {
                        panic!("Syntax error in input file, expected variable {} not found", name);
                    }
                }
                generate_external_gear_tooth_profile(
                    m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a,
                    h_fP0, n_sim, n_profile, n_tip,
                    &mut tool_contours, &mut contour_xy, &mut contour_xy_id,
                );
            }
            ToolType::ExternalWithProtuberance => {
                let mut alfa_P0: f64 = 0.;
                let mut alfa_KP0: f64 = 0.;
                let mut alfa_prP0: f64 = 0.;
                let mut m_n: f64 = 0.;
                let mut h_fP0: f64 = 0.;
                let mut h_FaP0: f64 = 0.;
                let mut h_prP0: f64 = 0.;
                let mut rho_aP0: f64 = 0.;
                for &name in ["alfa_P0", "m_n", "h_fP0", "h_aP0", "rho_aP0"].iter() {
                    let next_line = lines.next().unwrap_or_else(|| {panic!("Unexpected end of input file");})?;
                    if next_line.starts_with(name) {
                        if let Some(valstr) = next_line.split(":").last().map(|s| s.trim()) {
                            let val: f64 =
                                valstr.parse().unwrap_or_else(|e| {
                                panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                       next_line, e);
                                });
                            match name {
                                "alfa_P0"   => { alfa_P0 = PI/180. * val; }
                                "alfa_KP0"  => { alfa_KP0 = PI/180. * val; }
                                "alfa_prP0" => { alfa_prP0 = PI/180. * val; }
                                "m_n"       => { m_n = val; }
                                "h_fP0"     => { h_fP0 = m_n*val; }
                                "h_FaP0"    => { h_FaP0 = m_n*val; }
                                "h_prP0"    => { h_prP0 = m_n*val; }
                                "rho_aP0"   => { rho_aP0 = m_n*val; }
                                &_ => { panic!("Should never happen, internal error"); }
                            }
                        } else {
                            panic!("Syntax error in input file, {} entry does not contain a colon:\n{}",
                                   name, next_line);
                        }
                    } else {
                        panic!("Syntax error in input file, expected variable {} not found", name);
                    }
                }
                generate_external_gear_tooth_profile_with_protuberance(
                    m_n, beta, h_prP0, alfa_P0, rho_aP0, z, x, r_a,
                    h_fP0, h_FaP0, alfa_KP0, alfa_prP0,
                    n_sim, n_profile, n_tip,
                    &mut tool_contours, &mut contour_xy, &mut contour_xy_id,
                );
            }
            ToolType::Internal => {
                let mut alfa_P0: f64 = 0.;
                let mut m_n: f64 = 0.;
                let mut h_aP0: f64 = 0.;
                let mut rho_aP0: f64 = 0.;
                let mut r_a0: f64 = 0.;
                let mut z_0: i32 = 0;
                for &name in ["alfa_P0", "m_n", "h_aP0", "rho_aP0", "d_a0", "z_0"].iter() {
                    let next_line = lines.next().unwrap_or_else(|| {panic!("Unexpected end of input file");})?;
                    if next_line.starts_with(name) {
                        if let Some(valstr) = next_line.split(":").last().map(|s| s.trim()) {
                            if name == "z_0" { // i32
                                z_0 = valstr.parse().unwrap_or_else(|e| {
                                      panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                             next_line, e);
                                    });
                            } else { // f64
                                let val: f64 =
                                    valstr.parse().unwrap_or_else(|e| {
                                    panic!("Syntax error in input file, could not parse line\n{}\nError: {}",
                                           next_line, e);
                                        });
                                 match name {
                                    "alfa_P0" => { alfa_P0 = PI/180. * val; }
                                    "m_n"     => { m_n = val; }
                                    "h_aP0"   => { h_aP0 = m_n*val; }
                                    "rho_aP0" => { rho_aP0 = m_n*val; }
                                    "d_a0"    => { r_a0 = val/2.; }
                                    &_ => { panic!("Should never happen, internal error"); }
                                }
                            }
                        } else {
                            panic!("Syntax error in input file, {} entry does not contain a colon:\n{}",
                                   name, next_line);
                        }
                    } else {
                        panic!("Syntax error in input file, expected variable {} not found", name);
                    }
                }
                generate_internal_gear_tooth_profile(
                    m_n, beta, h_aP0, alfa_P0, rho_aP0, z, x, r_a,
                    z_0, r_a0,
                    n_sim, n_profile, n_tip,
                    &mut tool_contours, &mut contour_xy, &mut contour_xy_id,
                );
            }
            ToolType::Polyline => {
                let mut m_n: f64 = 0.;
                let mut h_aP0: f64 = 0.;
                match read_polygon_tool_contour(&mut m_n, beta, &mut lines, &mut tool_contours[0], &mut h_aP0) {
                    Ok(()) => {
                        generate_profile(
                            tool_type,
                            m_n, beta, h_aP0, z, x, r_a, 0, 0., 0.,
                            n_sim, n_profile, n_tip,
                            &mut tool_contours, &mut contour_xy, &mut contour_xy_id,
                        );
                    },
                    Err(e) => {
                        eprintln!("Error: {}", e);
                    }
                }
            }
        }
    }

    let mut fs_out = File::create(&output_file)?;
    writeln!(fs_out, "#             x               y  id")?;
    for it in 0..contour_xy.len() {
        writeln!(
            fs_out,
            "{:>15.10},{:>15.10},{:>3}",
            contour_xy[it][0],
            contour_xy[it][1],
            contour_xy_id[it]
        )?;
    }

    write_to_file(&output_file2, &tool_contours[0]);

    Ok(())
}

