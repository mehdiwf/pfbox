const PFBOX_VERSION: &str = env!("CARGO_PKG_VERSION");

use chrono;
use std::io::Write; // to use "write_all" method
use std::path::Path;
use std::fs; // to read/write a file contents
use rgsl::logarithm::log;
use rgsl::exponential::exp;
use crate::maths::fcts::*;
use crate::maths::mystructs::*;
use crate::maths::fcts::*;
use crate::configfile::cfg_io;
use crate::configfile::cfg_struct;

pub fn do_sim(configinput: cfg_struct::ConfigInput,
              overwrite: bool) {

    println!("Starting vdw simulation");

    let mut logproblem_counts = 0;
    let print_logproblems = false;
    let mut output_path = "./src/testoutput";
    let mut config = cfg_io::read_cfg_file(
        "./src/procedures/vdw_default_cfg.toml"); // the default config file
    
    match configinput {
        cfg_struct::ConfigInput::Empty => {
            config = cfg_io::read_cfg_file(
                "./src/procedures/vdw_default_cfg.toml");
        },
        cfg_struct::ConfigInput::Path(input_path) =>
        {
            config = cfg_io::read_cfg_file(&input_path);
        }
    }
    
    let output_dir = format!("./src/testoutput/{}",
                             &config.save_config.directory_name);

    let dir_exist = Path::new(&output_dir).is_dir();
    if (dir_exist && overwrite)
    {
        fs::remove_dir_all(&output_dir);
        println!("deleted already existing directory");
        fs::create_dir(&output_dir);
        println!("created fresh directory");
    }
    else
    {fs::create_dir(&output_dir);}

    let simdata_dir = format!("{}/sim_data",output_dir);
    fs::create_dir(&simdata_dir);
    let mut readmefile = fs::File::create(
        &format!("{}/README.txt", output_dir))
        .expect("couldn't create readme file");

    let mut logfile = fs::File::create(
        &format!("{}/logfile", output_dir))
        .expect("couldn't create logfile file");
    
    let simulation_utc_start_time = chrono::offset::Utc::now();
    
    let readme_message = format!(
        " ----------------------------------------\n\
         | Van Der Waals simulation readme file   |\n\
         | Pfbox version: {:>21}   |\n \
         ----------------------------------------\n\
         \n\
         user comment: {}\n\n\
         [TIME]\n\
         Simulation started at: {}\n",
        &PFBOX_VERSION,
        &config.save_config.user_comment,
        chrono::offset::Utc::now());
    readmefile.write_all(readme_message.as_bytes()).unwrap();

    cfg_io::write_cfg_file(&config,
                           &format!("{}/config_file.toml", output_dir));

    let cfg_struct::SimCfg{
        physics_config: 
        cfg_struct::PhyParam{
            dx,
            dy,
            dt,
            max_sim_time: max_time_step,
            temper0,
            temper1,
            rho_liq,
            rho_vap,
            Tc,
            aa,
            w,
            b,
            zeta0,
            eta0,
            m,
            lambda0,
            nb_col: x_max,
            nb_row: y_max
                
        },
        initial_time_config: 
        cfg_struct::InitCfg{
            n_liq,
        },
        save_config: 
        cfg_struct::SaveCfg{
            directory_name,
            user_comment,
            histo_freq,
            histo_save,
            num_bin,
        },
        physics_constants: 
        cfg_struct::PhyConstants{
            kB,
            DeBroglie0,
            dim,
            lambda_grad: lambda}
    } = config;

    let output_dir = format!("./src/testoutput/{}",
                             directory_name);
    // System initialisation
    let mut step = 0;
    let step_count_before_save = histo_save;

    let print_frequency = 10.;
    let mut print_percertage_threshold = 100./print_frequency;
    
    // let dt = 1e-2;
    let mut time = 0.;

    // let ncol_size = 100usize;
    // let nrow_size = 2usize;

    let rho_liq0 = rho_liq;
    let ln_rho_liq0 = log(rho_liq0);
    let rho_vap0 = rho_vap;
    let ln_rho_vap0 = log(rho_vap0);
    let temp0 = temper0;

    let box_info = BoxInfo{x_max: x_max as i32,
                           y_max: y_max as i32,
                           dx,
                           dy};

    //auie physics quantities definition
    ////////////////////////////////////////////////////////////////////////////
    // Physics quantities definition
    ////////////////////////////////////////////////////////////////////////////
    // GD for grid
    let mut GD_rho = ScalarField2D::new(x_max, y_max);
    let mut GD_temp = ScalarField2D::new(x_max, y_max);
    let mut GD_pressure = TensorField2D::new(x_max, y_max);
    // momentum, also known as J = rho*velocity
    let mut GD_J = VectorField2D::new(x_max, y_max);
    // velocity
    let mut GD_v = VectorField2D::new(x_max, y_max);

    //auie quantities used for the computations definition
    ////////////////////////////////////////////////////////////////////////////
    // Quantities used for the computations definition
    ////////////////////////////////////////////////////////////////////////////
    
    let mut GD_ln_rho = ScalarField2D::new(x_max, y_max);
    let mut GD_grad_rho = VectorField2D::new(x_max, y_max);
    let mut GD_lap_rho = ScalarField2D::new(x_max, y_max);

    let mut GD_vJ = TensorField2D::new(x_max, y_max);

    let mut GD_grad_v = TensorField2D::new(x_max, y_max);

    let mut GD_div_v = ScalarField2D::new(x_max, y_max);

    let mut GD_traceless_grad_v = TensorField2D::new(x_max, y_max);

    let mut GD_lap_v = VectorField2D::new(x_max, y_max);

    let mut GD_div_vJ = VectorField2D::new(x_max, y_max);
    let mut GD_grad_div_v = VectorField2D::new(x_max, y_max);

    let mut GD_div_press = VectorField2D::new(x_max, y_max);

    let mut GD_ln_rho_traceless_grad_v = VectorField2D::new(x_max, y_max);

    let mut GD_traceless_grad_v_dyadic_grad_v = ScalarField2D::new(x_max, y_max);

    let mut GD_grad_ln_rho_scalar_grad_T = ScalarField2D::new(x_max, y_max);

    let mut GD_grad_ln_rho = VectorField2D::new(x_max, y_max);

    let mut GD_v_scalar_grad_ln_rho = ScalarField2D::new(x_max, y_max);

    let mut GD_grad_ln_rho_traceless_grad_v = VectorField2D::new(x_max, y_max);
    
    let mut GD_grad_T = VectorField2D::new(x_max, y_max);

    let mut GD_lap_T = ScalarField2D::new(x_max, y_max);

    let mut GD_v_scalar_grad_T = ScalarField2D::new(x_max, y_max);

    let inv_cv = 1.0/(1.5*kB);

    //auie fluid initial state
    ////////////////////////////////////////////////////////////////////////////
    // Fluid initial state
    ////////////////////////////////////////////////////////////////////////////
    
    for yi in 0..y_max {
        for xi in 0..x_max {
            // putting liquid in the first half
            if (xi < x_max/2){
                GD_rho.set_pos(xi, yi,
                               rho_liq0);
                GD_ln_rho.set_pos(xi, yi,
                                  ln_rho_liq0);}
            else {GD_rho.set_pos(xi, yi,
                                 rho_vap0);
                  GD_ln_rho.set_pos(xi, yi,
                                    ln_rho_vap0);}

            // setting initial temperature
            GD_temp.set_pos(xi, yi, temp0);
        }}
    

    //auie computation variables update
    ////////////////////////////////////////////////////////////////////////////
    // Computations variables update
    ////////////////////////////////////////////////////////////////////////////
    println!("\nBeginning time loop\n");
    //auie time loop
    for i_time_step in 0..max_time_step {

        step = i_time_step;

        if (i_time_step <= 2)
            {
        println!("before update");
        println!("GD_div_press_x_std = {:.8}",
                 GD_div_press.x.s.std(1.));}
        
    // update of computations variables
    for yi in 0..y_max {
        for xi in 0..x_max {

            let yi_i32 = yi as i32;
            let xi_i32 = xi as i32;

            // -------------------------------------------------------
            // GD_lap_rho begin update
            GD_lap_rho
                .set_pos(xi, yi,
                         laplacian(&GD_rho,
                                   xi_i32, yi_i32,
                                   lambda,
                                   &box_info));
            // GD_lap_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_T begin update
            GD_lap_T
                .set_pos(xi, yi,
                         laplacian(&GD_temp,
                                   xi_i32, yi_i32,
                                   lambda,
                                   &box_info));
            // GD_lap_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_T begin update
            GD_grad_T
                .set_pos(xi, yi,
                         &gradient(&GD_temp,
                                   xi_i32, yi_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_rho begin update
            GD_grad_rho
                .set_pos(xi, yi,
                         &gradient(&GD_rho,
                                   xi_i32, yi_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_ln_rho begin update
            let rho = GD_rho.get_pos(xi, yi);
            if (rho < 0.) {
                let str_to_append = format!("IMPOSSIBLE COMPUTATION: NEGATIVE LOG\n\
                                             step {}, col={}, row={}\n\
                                             negative log of rho avoided!\n\
                                             rho value: {}\n\
                                             ------------\n",
                                            &step, &yi, &xi, &rho);
                logproblem_counts += 1;
                if print_logproblems {
                    println!("negative log detected ! check log for more info")};
                logfile.write_all(&str_to_append.as_bytes())
                    .expect("write failed");
                println!("error step {}:\n\
                          negative rho: rho = {}", step, rho);
            GD_ln_rho
                .set_pos(xi, yi,
                         0.);}
            else {
                let ln_rho = log(rho);
                GD_ln_rho
                    .set_pos(xi, yi,
                             ln_rho);}
            // GD_ln_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_v begin update
            GD_lap_v
                .set_pos(xi, yi,
                         &laplacian_vector(&GD_v,
                                           xi_i32, yi_i32,
                                           lambda,
                                           &box_info));
            // GD_lap_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_div_v begin update
            GD_div_v
                .set_pos(xi, yi,
                         div_vector(&GD_v,
                                    xi_i32, yi_i32,
                                    lambda,
                                    &box_info));
            // GD_div_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_v begin update
            GD_grad_v
                .set_pos(xi, yi,
                         &gradient_vector(&GD_v,
                                          xi_i32, yi_i32,
                                          lambda,
                                          &box_info));
            // GD_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho begin update
            GD_grad_ln_rho
                .set_pos(xi, yi,
                         &gradient(&GD_ln_rho,
                                   xi_i32, yi_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_ln_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_div_v begin update
            GD_grad_div_v
                .set_pos(xi, yi,
                         &grad_div_vel(&GD_v,
                                       xi_i32, yi_i32,
                                       lambda,
                                       &box_info));
            // GD_grad_div_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_traceless_grad_v begin update
            {
                let grad_v = GD_grad_v.get_pos(xi, yi);
                let div_v = GD_div_v.get_pos(xi, yi);
                
                let traceless_grad_v = tens2D {
                    xx: 2.*grad_v.xx - (2./(1.*dim as f64)) * div_v,
                    xy: grad_v.xy + grad_v.yx,
                    yx: grad_v.xy + grad_v.yx,
                    yy: 2.*grad_v.yy - (2./(1.*dim as f64)) * div_v};
                
                GD_traceless_grad_v.set_pos(xi, yi,
                                            &traceless_grad_v);
            }
            // GD_traceless_grad_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_vJ begin update
            {
                let v = GD_v.get_pos(xi, yi);
                let J = GD_J.get_pos(xi, yi);
                let tens_vJ = tens2D{
                    xx: v.x * J.x,
                    xy: v.x * J.y,
                    yx: v.y * J.x,
                    yy: v.y * J.y
                };
                GD_vJ
                    .set_pos(xi, yi,
                             &tens_vJ)
            }
            // GD_vJ end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_div_vJ begin update            
            GD_div_vJ
                .set_pos(xi, yi,
                         &div_tensor(&GD_vJ,
                                     xi_i32, yi_i32,
                                     lambda,
                                     &box_info));
            // GD_div_vJ end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_v_scal_grad_T begin update            
            GD_v_scalar_grad_T
                .set_pos(xi, yi,
                         scal_product(&GD_v.get_pos(xi, yi),
                                      &GD_grad_T.get_pos(xi, yi)));
            // GD_v_scal_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_traceless_grad_v_dyadic_grad_v begin update            
            GD_traceless_grad_v_dyadic_grad_v
                .set_pos(xi, yi,
                         dyadic_product(&GD_traceless_grad_v.get_pos(xi, yi),
                                        &GD_grad_v.get_pos(xi, yi)));
            // GD_traceless_grad_v_dyadic_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_v_scalar_grad_ln_rho begin update            
            GD_v_scalar_grad_ln_rho
                .set_pos(xi, yi,
                         scal_product(&GD_v.get_pos(xi, yi),
                                      &GD_grad_ln_rho.get_pos(xi, yi)));
            // GD_v_scalar_grad_ln_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_pressure begin update            
            GD_pressure
                .set_pos(xi, yi,
                         &pressure(GD_rho.get_pos(xi, yi),
                                   &GD_grad_rho.get_pos(xi, yi),
                                   GD_lap_rho.get_pos(xi, yi),
                                   GD_temp.get_pos(xi, yi),
                                   kB, aa, b, w, step + 10*xi as i32 + 10*yi as i32));
            // GD_pressure end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho_traceless_grad_v begin update
            
            GD_grad_ln_rho_traceless_grad_v
                .set_pos(xi, yi,
                         &tens_product_vec(
                             &GD_traceless_grad_v.get_pos(xi, yi),
                             &GD_grad_ln_rho.get_pos(xi, yi)));
            // GD_grad_ln_rho_traceless_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho_scalar_grad_T begin update
            GD_grad_ln_rho_scalar_grad_T
                .set_pos(xi, yi,
                         scal_product(&GD_grad_ln_rho.get_pos(xi, yi),
                                      &GD_grad_T.get_pos(xi, yi)));
            // GD_grad_ln_rho_scalar_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // div_press begin update
            GD_div_press
                .set_pos(xi, yi,
                         &div_tensor(&GD_pressure, xi_i32, yi_i32,
                                     lambda,
                                     &box_info));
            // div_press end update
            // -------------------------------------------------------

        }} // updating computations values end parenthesis


        if (i_time_step <= 2)
            {        
        println!("after update");
        
        println!("GD_div_press_x_std = {:.8}",
                 GD_div_press.x.s.std(1.));}
        
                let percentage_done = 100.*(step as f64/max_time_step as f64);
        if (percentage_done > print_percertage_threshold)
        {
            print_percertage_threshold += 100./print_frequency;
            println!("completed {percentage_done:.1}%");
            // println!("div_vJ.y = {:.8}\n\
            //           lap_v.y = {:.8}\n\
            //           v_y= {:.8}\n\
            //           grad_ln_rho_traceless_grad_v.y = {:.8}\n\
            //           div_press.y = {:.8}",
            //          GD_div_vJ.get_pos(0,25).y,
            //          GD_lap_v.get_pos(0,25).y,
            //          GD_v.get_pos(0,25).y,
            //          GD_ln_rho_traceless_grad_v.get_pos(0,25).y,
            //          GD_div_press.get_pos(0,25).y);
            // println!("-------");
            // println!("ln_rho = {:.8}\n\
            //           div_v = {:.8}\n\
            //           v_x= {:.8}\n\
            //           v_y= {:.8}",
            //          GD_ln_rho.get_pos(0,25),
            //          GD_div_v.get_pos(0,25),
            //          GD_v.get_pos(0,25).x,
            //          GD_v.get_pos(0,25).y);
            // println!("-------");            
        }
        let a = 0..=2;
        if (a.contains(&step))
            // (- div_vJ.x
            //  + eta0 * rho * lap_v.x
            //  + eta0 * (1.-2./(1.*dim as f64) + zeta0)
            //  * rho * grad_div_v.x
            //  + eta0 * rho * grad_ln_rho_traceless_grad_v.x
            //  + zeta0 * rho * grad_ln_rho.x * div_v
            //  - div_press.x)
            // printing init doing
        {   let print_pos = 0;
            println!("-----------------[I]\n\
                      STEP {step}");
            // println!("GD_div_vJ.x = {:.8}\n\
            //          GD_div_vJ.y = {:.8}",
            //          GD_div_vJ.x.get_pos(print_pos, 0),
            //          GD_div_vJ.y.get_pos(print_pos, 0));
            // println!("GD_lap_v.x = {:.8}\n\
            //          GD_lap_v.y = {:.8}",
            //          GD_lap_v.x.get_pos(print_pos, 0),
            //          GD_lap_v.y.get_pos(print_pos, 0));
            // println!("GD_div_v = {:.8}",
            //          GD_div_v.get_pos(print_pos, 0));
            // println!("GD_v.x = {:.8}\n\
            //           GD_v.y = {:.8}",
            //          GD_v.x.get_pos(print_pos, 0),
            //          GD_v.y.get_pos(print_pos, 0));
            // println!("GD_J.x = {:.8}\n\
            //           GD_J.y = {:.8}",
            //          GD_J.x.get_pos(print_pos, 0),
            //          GD_J.y.get_pos(print_pos, 0));
            println!("GD_rho_size: {:?}", GD_rho.s.raw_dim());
            println!("x_max: {} y_max: {}", x_max, y_max);
            println!("pos: {} {}", print_pos, 0);
            // println!("GD_rho = {:.8}",
            //          GD_rho.get_pos(print_pos, 0));
            // println!("GD_ln_rho = {:.8}",
            //          GD_ln_rho.get_pos(print_pos, 0));
            // println!("GD_grad_ln_rho.x = {:.8}\n\
            //           GD_grad_ln_rho.y = {:.8}",
            //          GD_grad_ln_rho.x.get_pos(print_pos, 0),
            //          GD_grad_ln_rho.y.get_pos(print_pos, 0));
            // println!("GD_div_press.x = {:.8}\n\
            //           GD_div_press.y = {:.8}",
            //          GD_div_press.x.get_pos(print_pos, 0),
            //          GD_div_press.y.get_pos(print_pos, 0));
            // println!("GD_div_press.x = {:?}",
            //          GD_div_press.x.s);

            //////////// DOING

            // println!("GD_div_press.xx = {:?}",
            //          GD_div_press.x.s);
            // println!("GD_div_press.xy = {:?}",
            //          GD_div_press.x.s);

            // println!("GD_div_vJ.xx = {:?}",
            //          GD_div_vJ.x.s);
            // println!("GD_div_vJ.xy = {:?}",
            //          GD_div_vJ.x.s);

            println!("GD_div_v = {:?}",
                     GD_div_v.s);

            println!("GD_v.x = {:?}",
                     GD_v.x.s);
            println!("GD_v.y = {:?}",
                     GD_v.y.s);
           
            // println!("GD_pressure.xx = {:?}",
            //          GD_pressure.xx.s);
            // println!("GD_pressure.xy = {:?}",
            //          GD_pressure.xy.s);
            // println!("GD_pressure.yx = {:?}",
            //          GD_pressure.yx.s);
            // println!("GD_pressure.yy = {:?}",
            //          GD_pressure.yy.s);

            // println!("GD_pressure.xx = {:.8}\n\
            //           GD_pressure.xy = {:.8}\n\
            //           GD_pressure.yx = {:.8}\n\
            //           GD_pressure.yy = {:.8}",
            //          GD_pressure.xx.get_pos(print_pos, 0),
            //          GD_pressure.xy.get_pos(print_pos, 0),
            //          GD_pressure.yx.get_pos(print_pos, 0),
            //          GD_pressure.yy.get_pos(print_pos, 0));
        }

    //auie main loop
    ////////////////////////////////////////////////////////////////////////////
    // Main loop
    ////////////////////////////////////////////////////////////////////////////
    
    for xi in 0..x_max {
        for yi in 0..y_max {

            let div_vJ = GD_div_vJ.get_pos(xi, yi);
            let rho = GD_rho.get_pos(xi, yi);
            let lap_v = GD_lap_v.get_pos(xi, yi);
            let grad_div_v = GD_grad_div_v.get_pos(xi, yi);
            let grad_ln_rho_traceless_grad_v =
                GD_grad_ln_rho_traceless_grad_v.get_pos(xi, yi);
            let grad_ln_rho = GD_grad_ln_rho.get_pos(xi, yi);
            let div_v = GD_div_v.get_pos(xi, yi);
            let div_press = GD_div_press.get_pos(xi, yi);
            let ln_rho = GD_ln_rho.get_pos(xi, yi);
            let v_grad_ln_rho = GD_v_scalar_grad_ln_rho.get_pos(xi, yi);
            let temp = GD_temp.get_pos(xi, yi);            
            let traceless_grad_v_dyadic_grad_v = GD_traceless_grad_v_dyadic_grad_v.get_pos(xi, yi);
            let grad_ln_rho_scalar_grad_T = GD_grad_ln_rho_scalar_grad_T.get_pos(xi, yi);
            let lap_T = GD_lap_T.get_pos(xi, yi);
            let v_scalar_grad_T = GD_v_scalar_grad_T.get_pos(xi, yi);
            let J = GD_J.get_pos(xi, yi);
            
            //bépo MOMENTUM conservation

            let mut new_J = vec2D
            {
                x: J.x +
                    (- div_vJ.x
	             + eta0 * rho * lap_v.x
                     + eta0 * (1.-2./(1.*dim as f64) + zeta0)
                     * rho * grad_div_v.x
	             + eta0 * rho * grad_ln_rho_traceless_grad_v.x
                     + zeta0 * rho * grad_ln_rho.x * div_v
                     - div_press.x)
                    * dt,
                y: J.y + 
                    (- div_vJ.y
	             + eta0 * rho * lap_v.y
                     + eta0 * (1.-2./(1.*dim as f64) + zeta0)
                     * rho * grad_div_v.y
	             + eta0 * rho * grad_ln_rho_traceless_grad_v.y
                     + zeta0 * rho * grad_ln_rho.y * div_v
                     - div_press.y)
                    * dt
            };

            // if you want gravity
            // J.y += -rho * gravity * dt;

            GD_J.set_pos(xi, yi, &new_J);

            //bépo MASS conservation

            // without ln_rho :
            //rho[i][j] -= div_J[i][j]*dt;
            
            let mut new_ln_rho = ln_rho -
                (div_v + v_grad_ln_rho) * dt;
            let mut new_rho = exp(new_ln_rho);
            
            GD_ln_rho.set_pos(xi, yi, new_ln_rho);
            GD_rho.set_pos(xi, yi, new_rho);

            //bépo VELOCITY from momentum
            GD_v.set_pos(xi, yi,
                         &vec2D{x: new_J.x/new_rho,
                                y: new_J.y/new_rho});
            
            //bépo TEMPERATURE ENERGY conservation
            
            // term l div_v

            let mut new_T = temp +
                inv_cv *
                (
                    // term l div_v
                    -kB * temp * (1. + rho * b/(1.-rho * b)) * div_v 
                    // term dissipative_stress_grad_v
                    + eta0 * traceless_grad_v_dyadic_grad_v
                    + zeta0 * div_v * div_v
                    // term laplacian T
                        + lambda0 * (grad_ln_rho_scalar_grad_T + lap_T)
                ) * dt
                - v_scalar_grad_T * dt;
            GD_temp.set_pos(xi, yi, new_T);
            
        }} // i, j loop closing parenthesis
            //bépo WRITING part

        if (step % step_count_before_save == 0) {
            
            let filename = format!("{}/step_{}",
                                   simdata_dir, i_time_step);
            let mut file = fs::File::create(&filename)
                .expect("couldn't create log file");
        
            file.write_all(
                "# column density temperature v_x v_y \
                 P_xx P_xy P_yx P_yy\n".as_bytes())
                .expect("write failed");

            let rho_profile = GD_rho.x_profile();
            let temp_profile = GD_temp.x_profile();
            let v_profile = GD_v.x_profile();
            let P_profile = GD_pressure.x_profile();


        for x_index in 0..x_max
        {
            let str_to_append = format!("{} {} {} {} {} {} {} {} {}\n",
                                        &x_index,
                                        &rho_profile[x_index],
                                        &temp_profile[x_index],
                                        &v_profile.x[x_index],
                                        &v_profile.y[x_index],
                                        &P_profile.xx[x_index],
                                        &P_profile.xy[x_index],
                                        &P_profile.yx[x_index],
                                        &P_profile.yy[x_index],
            );

            file.write_all(&str_to_append.as_bytes())
                .expect("write failed");
        }}
        
    } // time step closing parenthesis

    let simulation_utc_end_time = chrono::offset::Utc::now();
    let simulation_utc_duration = simulation_utc_end_time - simulation_utc_start_time;
    readmefile.write_all(
        &format!("Simulation ended at:   {}\n\
                  Duration: {}\n\n",
                 chrono::offset::Utc::now(),
                 simulation_utc_duration).as_bytes()).unwrap();

    readmefile.write_all(
        b"* MORE INFO *\n\
          You can access the simulation configuration\n\
          in config_file.toml, written in Tom's Obvious Minimal\n\
          Language. The simulation data is stored in sim_data directory, and\n\
          the logfile file contains possible warnings\n\n\
          The data stored in sim_data is the mean value on each column\n\
          for different parameters, which written in the first\n\
          commented line of each data file. The first numbers are the\n\
          column number");
    println!("\nSimulation finished, problems logged in log: {}", logproblem_counts);    
    
} // main definition closing parenthesis
