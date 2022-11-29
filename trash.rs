/*

// const Tc : f32 = 1.0;			   // critical temperature
// const rhom_c : f32 = 1.0/3.0; // critical mass density
// const Pc : f32 = 1.0; // critical pressure (not used)

// const aa : f32 = 1.0;
// const b : f32 = 1.0/(3.0*rhom_c);
// // square gradient coefficient 
// const w : f32 = 1.0;
// // Cahn Hilliard coefficient
// const G : f32 = 0.0;


// // chosen so that kBTc=1!
// const kB : f32 = 8./27.;

// const DeBroglie0 : f32 = 0.1;

// // molecule mass
// const m : f32 = 1.0;
// const inv_m : f32 = 1./m ;
// // shear viscosity
// const eta0 : f32 = 1.0;
// // bulk viscosity 
// const zeta0 : f32 = 1.0;
// // thermal conductivity
// const lambda0 : f32 = 0.2;
// // coefficient in front of the evaporative flux
// const Jev : f32 = 0.0;
// // liquid/vapor latent heat
// const hlv : f32 = 1.0;

// // coefficient to evaluate the gradients
// const lambda : f32 = 0.66;

// // space dimensionality 
// const dim : i32 = 2;

// // elementary unit cell length in the x direction
// const dx : f32 = 0.5;
// const dy : f32 = 0.5;

// // force density in the x_direction 
// const forcex : f32 = 0.00;
// // constant heat flux at the bottom of the crenel
// const flux : f32 = 0.00;

// // dimensions of the simulation box 
// // :mifzmodif:
// const NX : i32 = 100;
// const NY : i32 = 2 ;

// const j_wall_bot : i32 =  0;
// const j_wall_top : i32 =  NY-1;

// const rho_wall : f32 = 0.2/b;

// const Tw : f32 = 0.9;

// // number of time step of equilibration without energy transfer 
// const nsteps_eq_heat : i32 = 100000;
 
// const rho_min : f32 = 0.01/b;

// // a small number 
// const EPS : f32 = 0.004;

// const PI : f32 = 3.14159265358;

// #[derive(Debug)]
// struct vec2D {x: f32, y: f32}
// #[derive(Debug)]
// struct int2D {i: i32, j: i32}
// #[derive(Debug)]
// struct tens2D {
//     xx: f32,
//     xy: f32,
//     yx: f32,
//     yy: f32}

// fn add(x: &f32, y: &f32) -> f32
// {return x+y;}

// fn shear_viscosity(rho: &f32,
//                    temperature: &f32,
//                    shear_visc: &f32) -> f32
// {
//     eta0 * m * rho
// }

// fn bulk_viscosity(rho: &f32) -> f32
// {
//     zeta0 * m * rho
// }

// fn dissipative_stress(shear_visc: &f32,
//                       bulk_visc: &f32,
//                       grad_v: &tens2D, 
// 		      div_v: &f32,
//                       diss_stress: &tens2D) -> tens2D
// {
//     let new_diss_stress = tens2D
//     {
//         xx: (bulk_visc - shear_visc / 4.0) * div_v
//             + 2. * shear_visc * grad_v.xx
//             + shear_visc * ( grad_v.xy + grad_v.yx ),
//         yy: ( bulk_visc - shear_visc / 4.0 ) * div_v
//             + 2. * shear_visc * grad_v.yy,
//         xy: diss_stress.xy,
//         yx: diss_stress.yx,
//     };
//     return new_diss_stress
// }

// fn v_nabla_v(v: &vec2D, v_imp: &vec2D,
//              grad_v: &tens2D) -> vec2D
// {
//     let mut v_grad_v = vec2D{x: 0., y:0.};
//     v_grad_v.x = (v.x + v_imp.x) * grad_v.xx
//                  + (v.y + v_imp.y) * (grad_v.yx);
    
//     v_grad_v.y = (v.x + v_imp.x) * (grad_v.xy )
//                  + (v.y + v_imp.y) * grad_v.yy;
//     return v_grad_v;
// }

// fn scal_product(va: &vec2D,
//                 vb: &vec2D) -> f32
// {
//     return va.x * vb.x
//            + va.y * vb.y;
// }

// fn tens_product_vec(tens_a: &tens2D, 
//                     vec_b: &vec2D) -> vec2D
// {
//     let mut vec_ab = vec2D{x: 0., y: 0.};
//     vec_ab.x = tens_a.xx * vec_b.x + tens_a.yx * vec_b.y;
//     vec_ab.y = tens_a.xy * vec_b.x + tens_a.yy * vec_b.y;
//     return vec_ab;
// }

// fn dyadic_product(tens_a: &tens2D, tens_b: &tens2D) -> f32
// {
//     let dy_product: f32 = 
//         tens_a.xx * tens_b.xx + tens_a.yx * tens_b.yx
//         + tens_a.xy * tens_b.xy + tens_a.yy * tens_b.yy;
//     return dy_product;
// }


// if needed?
pub fn vector_x_access(vector: &vec2D) -> &f32
{
    return &vector.x;
}

pub fn vector_y_access(vector: &vec2D) -> &f32
{
    return &vector.y;
}

pub fn tensor_xx_access(tensor: &tens2D) -> &f32
{
    return &tensor.xx;
}

pub fn tensor_xy_access(tensor: &tens2D) -> &f32
{
    return &tensor.xy;
}

pub fn tensor_yx_access(tensor: &tens2D) -> &f32
{
    return &tensor.yx;
}

pub fn tensor_yy_access(tensor: &tens2D) -> &f32
{
    return &tensor.yy;
}


// pub fn test_access(tensor: &Vec<Vec<f32>>) -> &f32
// {
//     return &tensor[0][0];
// }
*/
