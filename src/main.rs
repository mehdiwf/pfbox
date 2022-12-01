#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused)]
#![allow(dead_code)]

mod maths;

use rgsl::logarithm::log;
use maths::fcts::*;
use maths::mystructs::*;

fn create_scalar_grid(x_size: i32,
                      y_size: i32) -> ScalarField2D
    // Vec<Vec<f64>>
{
    return ScalarField2D {
        s: vec![vec![0.; y_size as usize];
                x_size as usize],
    };
}

fn create_vector_grid(x_size: i32,
                      y_size: i32) -> VectorField2D
{
    let result = VectorField2D
    {
        x: vec![vec![0.; y_size as usize];
                x_size as usize],
        y: vec![vec![0.; y_size as usize];
                x_size as usize]
    };
    return result;
}

fn create_tensor_grid(x_size: i32,
                      y_size: i32) -> TensorField2D
{
    let result = TensorField2D
        {
            xx: vec![vec![0.; y_size as usize];
                     x_size as usize],
            xy: vec![vec![0.; y_size as usize];
                     x_size as usize],
            yx: vec![vec![0.; y_size as usize];
                     x_size as usize],
            yy: vec![vec![0.; y_size as usize];
                     x_size as usize]
        };
    return result;
}

fn main() {

    // System initialisation
    let mut step = 0;
    let max_step = 3;

    let dt = 0.1;
    let mut time = 0.;

    let x_size = 8;
    let y_size = 2;

    let rho_liq0 = 0.8;
    let ln_rho_liq0 = log(rho_liq0);
    let rho_vap0 = 0.02;
    let ln_rho_vap0 = log(rho_vap0);
    let temp0 = 0.7;

    let box_info = BoxInfo{imax: x_size,
                           jmax: y_size};

    ////////////////////////////////////////////////////////////////////////////
    // Physics quantities
    ////////////////////////////////////////////////////////////////////////////
    // GD for grid
    let mut GD_rho = create_scalar_grid(x_size, y_size);
    let mut GD_temp = create_scalar_grid(x_size, y_size);
    let mut GD_pressure = create_tensor_grid(x_size, y_size);
    // momentum, also known as J = rho*velocity
    let mut GD_J = create_vector_grid(x_size, y_size);
    // velocity
    let mut GD_v = create_vector_grid(x_size, y_size);

    ////////////////////////////////////////////////////////////////////////////
    // Quantities used for the computations
    ////////////////////////////////////////////////////////////////////////////
    
    let mut GD_ln_rho = create_scalar_grid(x_size, y_size);
    let mut GD_grad_rho = create_vector_grid(x_size, y_size);
    let mut GD_lap_rho = create_scalar_grid(x_size, y_size);

    let mut GD_vJ = create_tensor_grid(x_size, y_size);

    let mut GD_grad_v = create_tensor_grid(x_size, y_size);

    let mut GD_div_v = create_scalar_grid(x_size, y_size);

    let mut GD_traceless_grad_v = create_tensor_grid(x_size, y_size);

    let mut GD_lap_v = create_vector_grid(x_size, y_size);

    let mut GD_div_vJ = create_vector_grid(x_size, y_size);
    let mut GD_grad_div_v = create_vector_grid(x_size, y_size);

    let mut GD_div_press = create_vector_grid(x_size, y_size);

    let mut GD_ln_rho_traceless_grad_v = create_vector_grid(x_size, y_size);
    
    let inv_cv = 1.0/(1.5*kB);

    let mut GD_traceless_grad_v_dyadic_grad_v = create_scalar_grid(x_size, y_size);

    let mut GD_grad_ln_rho_scalar_grad_T = create_scalar_grid(x_size, y_size);

    let mut GD_grad_ln_rho = create_vector_grid(x_size, y_size);

    let mut GD_v_scal_grad_ln_rho = create_scalar_grid(x_size, y_size);

    let mut GD_grad_ln_rho_traceless_grad_v = create_vector_grid(x_size, y_size);
    
    let mut GD_grad_T = create_vector_grid(x_size, y_size);

    let mut GD_lap_T = create_scalar_grid(x_size, y_size);

    let mut GD_v_scal_grad_T = create_scalar_grid(x_size, y_size);

    // the grid used to compute stuff
    let mut GD_buffer = create_scalar_grid(x_size, y_size);

    let mut GD_buffer_vec = create_vector_grid(x_size, y_size);

    let mut GD_buffer_tens = create_tensor_grid(x_size, y_size);

    ////////////////////////////////////////////////////////////////////////////
    // Fluid initial state
    ////////////////////////////////////////////////////////////////////////////
    
    // initializing fluid with first half: vapor
    for i in 0..x_size {
        for j in 0..y_size {
            // putting liquid in the first half
            if (i < x_size/2){
                GD_rho.set_pos(i as usize, j as usize, &rho_liq0);
                GD_ln_rho.set_pos(i as usize, j as usize, &ln_rho_liq0);}}}

    // for i in 0..x_size {
    //     for j in 0..y_size {}}

    ////////////////////////////////////////////////////////////////////////////
    // Computations variables update
    ////////////////////////////////////////////////////////////////////////////

    // update of computations variables
    for i in 0usize..x_size as usize {
        for j in 0usize..y_size as usize {

            let i_i32 = i as i32;
            let j_i32 = j as i32;


            // -------------------------------------------------------
            // div_press begin update
            GD_div_press
                .set_pos(i, j,
                         &div_tensor(&GD_pressure, &i_i32, &j_i32,
                                     &box_info));
            // div_press end update
            // -------------------------------------------------------


            // -------------------------------------------------------
            // GD_grad_ln_rho_scalar_grad_T begin update
            GD_grad_ln_rho_scalar_grad_T
                .set_pos(i, j,
                         &scal_product(&GD_grad_ln_rho.get_pos(i, j),
                                       &GD_grad_T.get_pos(i, j)));
            // GD_grad_ln_rho_scalar_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_ln_rho_traceless_grad_v begin update
            
            GD_grad_ln_rho_traceless_grad_v
                .set_pos(i, j,
                         &tens_product_vec(
                             &GD_traceless_grad_v.get_pos(i, j),
                             &GD_grad_ln_rho.get_pos(i, j)));
            // GD_grad_ln_rho_scalar_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_ln_rho_traceless_grad_v begin update
            
            GD_grad_ln_rho_traceless_grad_v
                .set_pos(i, j,
                         &tens_product_vec(
                             &GD_traceless_grad_v.get_pos(i, j),
                             &GD_grad_ln_rho.get_pos(i, j)));
            // GD_grad_ln_rho_scalar_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_pressure begin update            
            GD_pressure
                .set_pos(i, j,
                         &pressure(&GD_rho.get_pos(i,j),
                                   &GD_grad_rho.get_pos(i,j),
                                   &GD_lap_rho.get_pos(i,j),
                                   &GD_lap_rho.get_pos(i,j)));
            // GD_pressure end update
            // -------------------------------------------------------
            
            // -------------------------------------------------------
            // GD_v_grad_ln_rho begin update            
            GD_v_scal_grad_ln_rho
                .set_pos(i, j,
                         &scal_product(&GD_v.get_pos(i, j),
                                       &GD_grad_ln_rho.get_pos(i, j)));
            // GD_v_grad_ln_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_traceless_grad_v_dyadic_grad_v begin update            
            GD_traceless_grad_v_dyadic_grad_v
                .set_pos(i, j,
                         &dyadic_product(&GD_traceless_grad_v.get_pos(i, j),
                                         &GD_grad_v.get_pos(i, j)));
            // GD_traceless_grad_v_dyadic_grad_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_v_scal_grad_T begin update            
            GD_v_scal_grad_T
                .set_pos(i, j,
                         &scal_product(&GD_v.get_pos(i, j),
                                       &GD_grad_T.get_pos(i, j)));
            // GD_v_scal_grad_T end update
            // -------------------------------------------------------
            
            // -------------------------------------------------------
            // GD_div_vJ begin update            
            GD_div_vJ
                .set_pos(i, j,
                         &div_tensor(&GD_vJ,
                         &i_i32, &j_i32, &box_info));
            // GD_div_vJ end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_vJ begin update
            {
                let v = GD_v.get_pos(i, j);
                let J = GD_J.get_pos(i, j);
                let tens_vJ = tens2D{
                    xx: v.x * J.x,
                    xy: v.x * J.y,
                    yx: v.y * J.x,
                    yy: v.y * J.y
                };
                GD_vJ
                    .set_pos(i, j,
                             &tens_vJ)
            }
            // GD_vJ end update
            // -------------------------------------------------------
            
            
            // -------------------------------------------------------
            // GD_grad_div_v begin update
            GD_grad_div_v
                .set_pos(i, j,
                         &grad_div_vel(&GD_v,
                                       &i_i32, &j_i32,
                                       &box_info));
            // GD_grad_div_v end update
            // -------------------------------------------------------
            
            // -------------------------------------------------------
            // GD_grad_div_v begin update
            GD_grad_ln_rho
                .set_pos(i, j,
                         &gradient(&GD_ln_rho,
                                   &i_i32, &j_i32,
                                   &box_info));
            // GD_grad_div_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_v begin update
            GD_grad_v
                .set_pos(i, j,
                         &gradient_vector(&GD_v,
                                          &i_i32, &j_i32,
                                          &box_info));
            // GD_grad_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_div_v begin update
            GD_div_v
                .set_pos(i, j,
                         &div_vector(&GD_v,
                                     &i_i32, &j_i32,
                                     &box_info));
            // GD_div_v end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_v begin update
            GD_lap_v
                .set_pos(i, j,
                         &laplacian_vector(&GD_v,
                                           &i_i32, &j_i32,
                                           &box_info));
            // GD_lap_v end update
            // -------------------------------------------------------
            
            // -------------------------------------------------------
            // GD_ln_rho begin update
            // :todo:log:
            let rho = GD_rho.get_pos(i, j);
            if (rho <= 1.) {
            println!("neg log: {}", rho);
            GD_ln_rho
                .set_pos(i, j,
                         &0.);}
            else {
                let ln_rho = log(GD_rho.get_pos(i, j));
                GD_ln_rho
                    .set_pos(i, j,
                             &ln_rho);}
            // GD_ln_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_rho begin update
            GD_grad_rho
                .set_pos(i, j,
                         &gradient(&GD_rho,
                                   &i_i32, &j_i32,
                                   &box_info));
            // GD_grad_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_rho begin update
            GD_lap_rho
                .set_pos(i, j,
                         &laplacian(&GD_rho,
                                    &i_i32, &j_i32,
                                    &box_info));
            // GD_lap_rho end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_lap_T begin update
            GD_lap_T
                .set_pos(i, j,
                         &laplacian(&GD_temp,
                                    &i_i32, &j_i32,
                                    &box_info));
            // GD_lap_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_grad_T begin update
            GD_grad_T
                .set_pos(i, j,
                         &gradient(&GD_temp,
                                   &i_i32, &j_i32,
                                   &box_info));
            // GD_grad_T end update
            // -------------------------------------------------------

            // -------------------------------------------------------
            // GD_traceless_grad_v begin update
            {
                let grad_v = GD_grad_v.get_pos(i, j);
                let div_v = GD_div_v.get_pos(i, j);
                
                let traceless_grad_v = tens2D {
                    xx: 2.*grad_v.xx - (2./(1.*dim as f64)) * div_v,
                    xy: grad_v.xy + grad_v.yx,
                    yx: grad_v.xy + grad_v.yx,
                    yy: 2.*grad_v.yy - (2./(1.*dim as f64)) * div_v};
                
                GD_traceless_grad_v.set_pos(i, j,
                                            &traceless_grad_v);
            }
            //---------------------------------------------    
            // GD_traceless_grad_v.xx[i][j] =
            //     2.*grad_v.xx - (2./(1.*dim as f64)) * div_v;
	    
            // GD_traceless_grad_v.xy[i][j] = 
            //     grad_v.xy + grad_v.yx;
            // GD_traceless_grad_v.yx[i][j] = 
            //     grad_v.xy + grad_v.yx;

	    // GD_traceless_grad_v.yy[i][j] =
            //     2.*grad_v.yy - (2./(1.*dim as f64)) * div_v;
            //---------------------------------------------
            // GD_traceless_grad_v end update
            // -------------------------------------------------------
            

            // let grad_v = grad_vector(&GD_vel,
            //                          &i_i32, &j_i32, 
            //                          &box_info);
            
            // // :doing:
            // GD_grad_v.x[i][j] = grad_v.x;
            // GD_div_v = div_vector(&GD_vel,
            //                       &i_i32, &j_i32, 
            //                       &box_info);

        }}
    

    // main loop
    for i in 0usize..x_size as usize {
        for j in 0usize..y_size as usize {

            let i_i32 = i as i32;
            let j_i32 = j as i32;

            // momentum equation: variables needed
            let div_vJ = div_tensor(&GD_vJ, 
                                    &i_i32, &j_i32, 
                                    &box_info);
            let lap_v = laplacian_vector(&GD_v, 
                                         &i_i32, &j_i32, 
                                         &box_info);

            
            
            
        }}

    // v = tens_product_vec(&t, &v);
    // f = dyadic_product(&t, &t2);
    // t = grad_vector(&vecfield, &2, &2, 
    //                 &box_info);
    // f = div_vector(&vecfield, &2, &2, 
    //                &box_info);
    // v = div_tensor(&tensfield, &2, &2, 
    //                &box_info);
    // f = laplacian(&grid_y2, &2, &2, 
    //               &box_info);
    // v = laplacian_vector(&vecfield, &2, &2, &box_info);

    // v = grad_div_vel(&vecfield, &2, &2, &box_info);

    // println!("{:?}", v);
}
