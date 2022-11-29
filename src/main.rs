#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused)]
#![allow(dead_code)]

mod maths;
use maths::fcts::*;
use maths::mystructs::*;

fn main() {
    // an uniform grid 
    let mut uni_grid = vec![vec![1.; 4]; 4];
    // a grid growing in the x axis
    let mut grid_x = vec![vec![1., 2., 3., 4.]; 4];
    // a grid growing in the y axis

    let grid_y = vec![vec![1., 1., 1., 1.],
                      vec![2., 2., 2., 2.],
                      vec![3., 3., 3., 3.],
                      vec![4., 4., 4., 4.]];

    let grid_y2 = vec![vec![1., 1., 1., 1.],
                      vec![2., 2., 2., 2.],
                      vec![4., 4., 4., 4.],
                      vec![8., 8., 8., 8.]];
    let mut x = 1.5;
    let mut i = 1;
    let mut f = 0.;
    let mut v = vec2D{x: 1., y: 2.};
    let mut t = tens2D{xx: 1., xy: 1.,
                   yx: 1., yy: 1.};
    let mut t2 = tens2D{xx:2., xy:2.,
                        yx:2., yy:2.};
    let mut y = 1.;
    let mut z = 0.;
    let box_info = BoxInfo {
        imax: 4,
        jmax: 4};

    let vecfield = VectorField2D{
        x: vec![vec![1., 2., 3., 4.]; 4],
        // y: uni_grid
        y: vec![vec![1., 2., 3., 4.]; 4]
   };

    let tensfield = TensorField2D{
        xx: vec![vec![1., 2., 3., 4.]; 4],
        xy: vec![vec![1., 2., 3., 4.]; 4],
        yx: vec![vec![1., 2., 3., 4.]; 4],
        yy: vec![vec![1., 2., 3., 4.]; 4]
   };
    
    // z = shear_viscosity(&x, &x, &x);
    // t2 = dissipative_stress(&x, &x, &t, &x, &t2);
    // x = scal_product(&v, &v);
    // v = v_nabla_v(&v, &v, &t);
    v = tens_product_vec(&t, &v);
    f = dyadic_product(&t, &t2);
    t = grad_vector(&vecfield, &2, &2, 
                    &box_info);
    f = div_vector(&vecfield, &2, &2, 
                   &box_info);
    v = div_tensor(&tensfield, &2, &2, 
                   &box_info);
    f = laplacian(&grid_y2, &2, &2, 
                  &box_info);
    println!("{:?}", f);
}
