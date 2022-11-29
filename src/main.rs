#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused)]
#![allow(dead_code)]

mod maths;

use maths::fcts::*;
use maths::mystructs::*;

fn main() {

    // System initialisation
    let mut step = 0;
    let max_step = 3;

    let dt = 0.1;
    let mut time = 0.;

    let x_size = 8;
    let y_size = 2;

    let rho_liq0 = 0.8;
    let rho_vap0 = 0.02;
    let temp0 = 0.7;

    let box_info = BoxInfo{imax: x_size,
                           jmax: y_size};

    let mut grid_rho = vec![vec![rho_vap0; y_size as usize]; 
                            x_size as usize];
    let mut grid_temp = vec![vec![temp0; y_size as usize]; 
                             x_size as usize];
    let mut grid_pressure = vec![vec![0.; y_size as usize];
                                 x_size as usize];
    // momentum, also known as J = rho*velocity
    let mut grid_momentum = vec![vec![0.; y_size as usize];
                                 x_size as usize];
    let mut grid_vel = vec![vec![0.; y_size as usize];
                            x_size as usize];

    // initializing fluid with step

    for i in 0..x_size {
        for j in 0..y_size {
            if (i < x_size/2){
            grid_rho[i as usize][j as usize] = 
                rho_liq0;}}}

    // for i in 0..x_size {
    //     for j in 0..y_size {}}
    



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
