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

fn create_scalar_grid(ncol_size: i32,
                      nrow_size: i32) -> ScalarField2D
{
    return ScalarField2D {
        s: vec![vec![0.; ncol_size as usize];
                nrow_size as usize],
    };
}

fn create_vector_grid(ncol_size: i32,
                      nrow_size: i32) -> VectorField2D
{
    let result = VectorField2D
    {
        x: vec![vec![0.; ncol_size as usize];
                nrow_size as usize],
        y: vec![vec![0.; ncol_size as usize];
                nrow_size as usize]
    };
    return result;
}

fn create_tensor_grid(ncol_size: i32,
                      nrow_size: i32) -> TensorField2D
{
    let result = TensorField2D
        {
            xx: vec![vec![0.; ncol_size as usize];
                     nrow_size as usize],
            xy: vec![vec![0.; ncol_size as usize];
                     nrow_size as usize],
            yx: vec![vec![0.; ncol_size as usize];
                     nrow_size as usize],
            yy: vec![vec![0.; ncol_size as usize];
                     nrow_size as usize]
        };
    return result;
}

pub fn do_sim(configinput: cfg_struct::ConfigInput,
              overwrite: bool) {
    let mut logproblem_counts = 0;
    let print_logproblems = false;
    let mut output_path = "./src/testoutput";
    let mut config = cfg_io::read_cfg_file(
        "src/procedures/vdw_default_cfg.toml"); // the default config file
    
    match configinput {
        cfg_struct::ConfigInput::Empty => {
            config = cfg_io::read_cfg_file(
                "src/procedures/vdw_default_cfg.toml");
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
    {
        fs::create_dir(&output_dir);
    }

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
            lambda0
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

    let print_frequency = 20.;
    let mut print_percertage_threshold = 100./print_frequency;
    
    // let dt = 1e-2;
    let mut time = 0.;

    let ncol_size = 100;
    let nrow_size = 2;

    let rho_liq0 = rho_liq;
    let ln_rho_liq0 = log(rho_liq0);
    let rho_vap0 = rho_vap;
    let ln_rho_vap0 = log(rho_vap0);
    let temp0 = temper0;

    let box_info = BoxInfo{col_max: ncol_size,
                           row_max: nrow_size,
                           col_dx: dx,
                           row_dx: dy};

    //auie physics quantities definition
    ////////////////////////////////////////////////////////////////////////////
    // Physics quantities definition
    ////////////////////////////////////////////////////////////////////////////
    // GD for grid
    let mut GD_rho = create_scalar_grid(ncol_size, nrow_size);
    let mut GD_temp = create_scalar_grid(ncol_size, nrow_size);
    let mut GD_pressure = create_tensor_grid(ncol_size, nrow_size);
    // momentum, also known as J = rho*velocity
    let mut GD_J = create_vector_grid(ncol_size, nrow_size);
    // velocity
    let mut GD_v = create_vector_grid(ncol_size, nrow_size);

    //auie quantities used for the computations definition
    ////////////////////////////////////////////////////////////////////////////
    // Quantities used for the computations definition
    ////////////////////////////////////////////////////////////////////////////
    
    let mut GD_ln_rho = create_scalar_grid(ncol_size, nrow_size);
    let mut GD_grad_rho = create_vector_grid(ncol_size, nrow_size);
    let mut GD_lap_rho = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_vJ = create_tensor_grid(ncol_size, nrow_size);

    let mut GD_grad_v = create_tensor_grid(ncol_size, nrow_size);

    let mut GD_div_v = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_traceless_grad_v = create_tensor_grid(ncol_size, nrow_size);

    let mut GD_lap_v = create_vector_grid(ncol_size, nrow_size);

    let mut GD_div_vJ = create_vector_grid(ncol_size, nrow_size);
    let mut GD_grad_div_v = create_vector_grid(ncol_size, nrow_size);

    let mut GD_div_press = create_vector_grid(ncol_size, nrow_size);

    let mut GD_ln_rho_traceless_grad_v = create_vector_grid(ncol_size, nrow_size);
    
    let inv_cv = 1.0/(1.5*kB);

    let mut GD_traceless_grad_v_dyadic_grad_v = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_grad_ln_rho_scalar_grad_T = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_grad_ln_rho = create_vector_grid(ncol_size, nrow_size);

    let mut GD_v_scalar_grad_ln_rho = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_grad_ln_rho_traceless_grad_v = create_vector_grid(ncol_size, nrow_size);
    
    let mut GD_grad_T = create_vector_grid(ncol_size, nrow_size);

    let mut GD_lap_T = create_scalar_grid(ncol_size, nrow_size);

    let mut GD_v_scalar_grad_T = create_scalar_grid(ncol_size, nrow_size);

    //auie fluid initial state
    ////////////////////////////////////////////////////////////////////////////
    // Fluid initial state
    ////////////////////////////////////////////////////////////////////////////
    
    for col in 0usize..ncol_size as usize {
        for row in 0usize..nrow_size as usize {
            // putting liquid in the first half
            if ((col as i32) < ncol_size/2){
                GD_rho.set_pos(row, col,
                               &rho_liq0);
                GD_ln_rho.set_pos(row, col,
                                  &ln_rho_liq0);}
            else {GD_rho.set_pos(row, col,
                                 &rho_vap0);
                  GD_ln_rho.set_pos(row, col,
                                    &ln_rho_vap0);}

            // setting initial temperature
            GD_temp.set_pos(row, col, &temp0);
        }}
    

    //auie computation variables update
    ////////////////////////////////////////////////////////////////////////////
    // Computations variables update
    ////////////////////////////////////////////////////////////////////////////
    println!("\nBeginning time loop\n");
    //auie time loop
    for i_time_step in 0..max_time_step {

        step = i_time_step;
        let percentage_done = 100.*(step as f64/max_time_step as f64);
        if (percentage_done > print_percertage_threshold)
        {
            print_percertage_threshold += 100./print_frequency;
            println!("completed {percentage_done:.1}%");
        }

    // update of computations variables
    for col in 0usize..ncol_size as usize {
        for row in 0usize..nrow_size as usize {

            let col_i32 = col as i32;
            let row_i32 = row as i32;

            // -------------------------------------------------------
            // GD_lap_rho begin update
            GD_lap_rho
                .set_pos(row, col,
                         &laplacian(&GD_rho,
                                    row_i32, col_i32,
                                    lambda,
                                    &box_info));
            // GD_lap_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_T begin update
            GD_lap_T
                .set_pos(row, col,
                         &laplacian(&GD_temp,
                                    row_i32, col_i32,
                                    lambda,
                                    &box_info));
            // GD_lap_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_T begin update
            GD_grad_T
                .set_pos(row, col,
                         &gradient(&GD_temp,
                                   row_i32, col_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_rho begin update
            GD_grad_rho
                .set_pos(row, col,
                         &gradient(&GD_rho,
                                   row_i32, col_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_ln_rho begin update
            // :todo:log:
            let rho = GD_rho.get_pos(row, col);
            if (rho < 0.) {
                let str_to_append = format!("step {}, col={}, row={}\n\
                                             neg log: {}\n\
                                             ------------\n",
                                            &step, &col, &row, &rho);
                // appending the string to 
                file.write_all(&str_to_append.as_bytes())
                    .expect("write failed");
                println!("error step {}:\n\
                          negative rho: rho = {}", step, rho);
            GD_ln_rho
                .set_pos(row, col,
                         &0.);}
            else {
                let ln_rho = log(rho);
                GD_ln_rho
                    .set_pos(row, col,
                             &ln_rho);}
            // GD_ln_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_lap_v begin update
            GD_lap_v
                .set_pos(row, col,
                         &laplacian_vector(&GD_v,
                                           row_i32, col_i32,
                                           lambda,
                                           &box_info));
            // GD_lap_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_div_v begin update
            GD_div_v
                .set_pos(row, col,
                         &div_vector(&GD_v,
                                     row_i32, col_i32,
                                     lambda,
                                     &box_info));
            // GD_div_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_v begin update
            GD_grad_v
                .set_pos(row, col,
                         &gradient_vector(&GD_v,
                                          row_i32, col_i32,
                                          lambda,
                                          &box_info));
            // GD_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho begin update
            GD_grad_ln_rho
                .set_pos(row, col,
                         &gradient(&GD_ln_rho,
                                   row_i32, col_i32,
                                   lambda,
                                   &box_info));
            // GD_grad_ln_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_div_v begin update
            GD_grad_div_v
                .set_pos(row, col,
                         &grad_div_vel(&GD_v,
                                       row_i32, col_i32,
                                       lambda,
                                       &box_info));
            // GD_grad_div_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_traceless_grad_v begin update
            {
                let grad_v = GD_grad_v.get_pos(row, col);
                let div_v = GD_div_v.get_pos(row, col);
                
                let traceless_grad_v = tens2D {
                    xx: 2.*grad_v.xx - (2./(1.*dim as f64)) * div_v,
                    xy: grad_v.xy + grad_v.yx,
                    yx: grad_v.xy + grad_v.yx,
                    yy: 2.*grad_v.yy - (2./(1.*dim as f64)) * div_v};
                
                GD_traceless_grad_v.set_pos(row, col,
                                            &traceless_grad_v);
            }
            // GD_traceless_grad_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_vJ begin update
            {
                let v = GD_v.get_pos(row, col);
                let J = GD_J.get_pos(row, col);
                let tens_vJ = tens2D{
                    xx: v.x * J.x,
                    xy: v.x * J.y,
                    yx: v.y * J.x,
                    yy: v.y * J.y
                };
                GD_vJ
                    .set_pos(row, col,
                             &tens_vJ)
            }
            // GD_vJ end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_div_vJ begin update            
            GD_div_vJ
                .set_pos(row, col,
                         &div_tensor(&GD_vJ,
                                     row_i32, col_i32,
                                     lambda,
                                     &box_info));
            // GD_div_vJ end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_v_scal_grad_T begin update            
            GD_v_scalar_grad_T
                .set_pos(row, col,
                         &scal_product(&GD_v.get_pos(row, col),
                                       &GD_grad_T.get_pos(row, col)));
            // GD_v_scal_grad_T end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_traceless_grad_v_dyadic_grad_v begin update            
            GD_traceless_grad_v_dyadic_grad_v
                .set_pos(row, col,
                         &dyadic_product(&GD_traceless_grad_v.get_pos(row, col),
                                         &GD_grad_v.get_pos(row, col)));
            // GD_traceless_grad_v_dyadic_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_v_scalar_grad_ln_rho begin update            
            GD_v_scalar_grad_ln_rho
                .set_pos(row, col,
                         &scal_product(&GD_v.get_pos(row, col),
                                       &GD_grad_ln_rho.get_pos(row, col)));
            // GD_v_scalar_grad_ln_rho end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_pressure begin update            
            GD_pressure
                .set_pos(row, col,
                         &pressure(GD_rho.get_pos(row, col),
                                   &GD_grad_rho.get_pos(row, col),
                                   GD_lap_rho.get_pos(row, col),
                                   GD_temp.get_pos(row, col),
                                   kB, aa, b, w));
            // GD_pressure end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho_traceless_grad_v begin update
            
            GD_grad_ln_rho_traceless_grad_v
                .set_pos(row, col,
                         &tens_product_vec(
                             &GD_traceless_grad_v.get_pos(row, col),
                             &GD_grad_ln_rho.get_pos(row, col)));
            // GD_grad_ln_rho_traceless_grad_v end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho_scalar_grad_T begin update
            GD_grad_ln_rho_scalar_grad_T
                .set_pos(row, col,
                         &scal_product(&GD_grad_ln_rho.get_pos(row, col),
                                       &GD_grad_T.get_pos(row, col)));
            // GD_grad_ln_rho_scalar_grad_T end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // div_press begin update
            GD_div_press
                .set_pos(row, col,
                         &div_tensor(&GD_pressure, row_i32, col_i32,
                                     lambda,
                                     &box_info));
            // div_press end update
            // -------------------------------------------------------


        }} // updating computations values end parenthesis

    //bépo WRITING part

        if (step % step_count_before_save == 0) {
        
        let filename = format!("{}/step_{}",
                               simdata_dir, i_time_step);
        let mut file = fs::File::create(&filename)
            .expect("couldn't create log file");
        
        file.write_all(
            "# column density temperature\n".as_bytes())
            .expect("write failed");

        let rho_profile = GD_rho.x_profile();
        let temp_profile = GD_temp.x_profile();
        
        for col_index in 0usize..ncol_size as usize
        {
            let str_to_append = format!("{} {} {}\n",
                                        &col_index,
                                        &rho_profile[col_index],
                                        &temp_profile[col_index]);

            file.write_all(&str_to_append.as_bytes())
                .expect("write failed");
        }}
    // let str_to_append = format!("step {}, i={}, j={}\n\
    //                              neg log {}\n\
    //                              ------------\n",
    //                             &step, &i, &j, &rho);
    //     // appending the string to 
    //     file.write_all(&str_to_append.as_bytes())
    //         .expect("write failed");
        

    //auie main loop
    ////////////////////////////////////////////////////////////////////////////
    // Main loop
    ////////////////////////////////////////////////////////////////////////////
    
    for row in 0usize..nrow_size as usize {
        for col in 0usize..ncol_size as usize {

            let row_i32 = row as i32;
            let col_i32 = col as i32;

            let div_vJ = GD_div_vJ.get_pos(row, col);
            let rho = GD_rho.get_pos(row, col);
            let lap_v = GD_lap_v.get_pos(row, col);
            let grad_div_v = GD_grad_div_v.get_pos(row, col);
            let grad_ln_rho_traceless_grad_v =
                GD_grad_ln_rho_traceless_grad_v.get_pos(row, col);
            let grad_ln_rho = GD_grad_ln_rho.get_pos(row, col);
            let div_v = GD_div_v.get_pos(row, col);
            let div_press = GD_div_press.get_pos(row, col);
            let ln_rho = GD_ln_rho.get_pos(row, col);
            let v_grad_ln_rho = GD_v_scalar_grad_ln_rho.get_pos(row, col);
            let temp = GD_temp.get_pos(row, col);            
            let traceless_grad_v_dyadic_grad_v = GD_traceless_grad_v_dyadic_grad_v.get_pos(row, col);
            let grad_ln_rho_scalar_grad_T = GD_grad_ln_rho_scalar_grad_T.get_pos(row, col);
            let lap_T = GD_lap_T.get_pos(row, col);
            let v_scalar_grad_T = GD_v_scalar_grad_T.get_pos(row, col);
            let J = GD_J.get_pos(row, col);
            
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

            GD_J.set_pos(row, col, &new_J);

            //bépo MASS conservation

            // without ln_rho :
            //rho[i][j] -= div_J[i][j]*dt;
            
            let mut new_ln_rho = ln_rho -
                (div_v + v_grad_ln_rho) * dt;
            let mut new_rho = exp(new_ln_rho);
            
            GD_ln_rho.set_pos(row, col, &new_ln_rho);
            GD_rho.set_pos(row, col, &new_rho);

            //bépo VELOCITY from momentum
            GD_v.set_pos(row, col,
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
            GD_temp.set_pos(row, col, &new_T);
            
        }} // i, j loop closing parenthesis
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
