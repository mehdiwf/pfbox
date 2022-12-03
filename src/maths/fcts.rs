use super::mystructs::*;

pub const Tc : f64 = 1.0;			   // critical temperature
pub const rhom_c : f64 = 1.0/3.0; // critical mass density
pub const Pc : f64 = 1.0; // critical pressure (not used)

pub const aa : f64 = 1.0;
pub const b : f64 = 1.0/(3.0*rhom_c);
// square gradient coefficient 
pub const w : f64 = 1.0;
// Cahn Hilliard coefficient
pub const G : f64 = 0.0;


// chosen so that kBTc=1!
pub const kB : f64 = 8./27.;

pub const DeBroglie0 : f64 = 0.1;

// molecule mass
pub const m : f64 = 1.0;
pub const inv_m : f64 = 1./m ;
// shear viscosity
pub const eta0 : f64 = 1.0;
// bulk viscosity 
pub const zeta0 : f64 = 1.0;
// thermal conductivity
pub const lambda0 : f64 = 0.2;
// coefficient in front of the evaporative flux
pub const Jev : f64 = 0.0;
// liquid/vapor latent heat
pub const hlv : f64 = 1.0;

// coefficient to evaluate the gradients
pub const lambda : f64 = 0.66;

// space dimensionality 
pub const dim : i32 = 2;

// elementary unit cell length in the x direction
pub const dx : f64 = 0.02;
pub const dy : f64 = 0.02;

// force density in the x_direction 
pub const forcex : f64 = 0.00;
// constant heat flux at the bottom of the crenel
pub const flux : f64 = 0.00;

// dimensions of the simulation box 
// :mifzmodif:
pub const NX : i32 = 100;
pub const NY : i32 = 2 ;

pub const j_wall_bot : i32 =  0;
pub const j_wall_top : i32 =  NY-1;

pub const rho_wall : f64 = 0.2/b;

pub const Tw : f64 = 0.9;

// number of time step of equilibration without energy transfer 
pub const nsteps_eq_heat : i32 = 100000;
 
pub const rho_min : f64 = 0.01/b;

// a small number 
pub const EPS : f64 = 0.004;

pub const PI : f64 = 3.14159265358;

pub fn shear_viscosity(rho: f64) -> f64
{
    eta0 * m * rho
}

pub fn bulk_viscosity(rho: f64) -> f64
{
    zeta0 * m * rho
}

pub fn dissipative_stress(shear_visc: f64,
                          bulk_visc: f64,
                          grad_v: &tens2D, 
		          div_v: f64,
                          diss_stress: &tens2D) -> tens2D
{
    let new_diss_stress = tens2D
    {
        xx: (bulk_visc - shear_visc / 4.0) * div_v
            + 2. * shear_visc * grad_v.xx
            + shear_visc * ( grad_v.xy + grad_v.yx ),
        yy: ( bulk_visc - shear_visc / 4.0 ) * div_v
            + 2. * shear_visc * grad_v.yy,
        xy: diss_stress.xy,
        yx: diss_stress.yx,
    };
    return new_diss_stress
}

pub fn v_nabla_v(v: &vec2D, v_imp: &vec2D,
                 grad_v: &tens2D) -> vec2D
{
    let v_grad_v = vec2D{
        x: (v.x + v_imp.x) * grad_v.xx
            + (v.y + v_imp.y) * grad_v.yx,
        y: (v.x + v_imp.x) * grad_v.xy
            + (v.y + v_imp.y) * grad_v.yy};

    return v_grad_v;
}

pub fn scal_product(va: &vec2D,
                    vb: &vec2D) -> f64
{
    return va.x * vb.x
           + va.y * vb.y;
}

pub fn tens_product_vec(tens_a: &tens2D, 
                    vec_b: &vec2D) -> vec2D
{
    let vec_ab = vec2D{
        x: tens_a.xx * vec_b.x + tens_a.yx * vec_b.y,
        y: tens_a.xy * vec_b.x + tens_a.yy * vec_b.y};
    return vec_ab;
}

pub fn dyadic_product(tens_a: &tens2D, tens_b: &tens2D) -> f64
{
    let dy_product: f64 = 
        tens_a.xx * tens_b.xx + tens_a.yx * tens_b.yx
        + tens_a.xy * tens_b.xy + tens_a.yy * tens_b.yy;
    return dy_product;
}

/// returns the partial derivative of a scalar field
/// in i,j and in the direction 0 (j axis) or 1 (i axis)
/// WARNING!!!: partial derivative in the x direction is 1, not 0
/// WARNING!!!: partial derivative in the y direction is 0, not 1
pub fn partial_deriv(scal_field: &Vec<Vec<f64>>,
                     i: i32, j: i32,
                     direction: i32,
                     box_info: &BoxInfo) -> f64
{
    let BoxInfo { col_max: box_col_max,
                  row_max: box_row_max } = box_info;
    // i+1 with Periodic Boundaries
    let ip = ((i+1) % box_row_max) as usize;
    // i-1 with Periodic Boundaries
    let im = ((i - 1 + box_row_max) % box_row_max) as usize;
    // j+1 with Periodic Boundaries
    let jp = ((j+1) % box_col_max) as usize;
    // j-1 with Periodic Boundaries
    let jm = ((j - 1 + box_col_max) % box_col_max) as usize;
    let (i, j) = (i as usize, j as usize);

    match direction
        {
            // on the x axis
            0 => {
   let derivative = 
          lambda*(scal_field[ip][j] - scal_field[im][j])/(2.*dx)
	  + 0.25*lambda*(scal_field[ip][jp] - scal_field[im][jp])/(2.*dx)
          + 0.25*lambda*(scal_field[ip][jm] - scal_field[im][jm])/(2.*dx);
                return derivative;
            },
            // on the y axis
            1 => {
   let derivative = 
          lambda*(scal_field[i][jp] - scal_field[i][jm])/(2.*dx)
	  + 0.25*lambda*(scal_field[ip][jp] - scal_field[ip][jm])/(2.*dx)	
          + 0.25*lambda*(scal_field[im][jp] - scal_field[im][jm])/(2.*dx);
                return derivative;
            },
            // if different than 0/1, panic
            _ => panic!("direction should be 0 or 1")
        }
    
}

pub fn grad_scalar(scalar_field: &Vec<Vec<f64>>,
                   i: i32, j: i32,
                   box_info: &BoxInfo) -> vec2D
{
    let grad = vec2D {
        // x is index 1 because it's the columns in the simulation
        x: partial_deriv(&scalar_field, i, j, 1, &box_info),
        y: partial_deriv(&scalar_field, i, j, 0, &box_info),};
        
    return grad;
}

pub fn gradient(scalar_field: &ScalarField2D,
                i: i32, j: i32,
                box_info: &BoxInfo) -> vec2D
{
    let field = &scalar_field.s;
    return grad_scalar(field, i, j,
                       &box_info);
}

pub fn gradient_vector(vector_field: &VectorField2D,
                    i: i32, j: i32,
                    box_info: &BoxInfo) -> tens2D
{
    let tens = tens2D {
        xx: partial_deriv(&vector_field.x, i, j, 1, &box_info),
        xy: partial_deriv(&vector_field.x, i, j, 0, &box_info),
        yx: partial_deriv(&vector_field.y, i, j, 1, &box_info),
        yy: partial_deriv(&vector_field.y, i, j, 0, &box_info)};
    return tens;
}

pub fn div_vector(vector_field: &VectorField2D,
                  i: i32, j: i32,
                  box_info: &BoxInfo) -> f64
{
    let dVx_dx = partial_deriv(&vector_field.x, i, j, 1, &box_info);
    let dVy_dy = partial_deriv(&vector_field.y, i, j, 0, &box_info);

    return dVx_dx + dVy_dy;
}

pub fn div_tensor(tensor_field: &TensorField2D,
                  i: i32, j: i32,
                  box_info: &BoxInfo) -> vec2D
{
    let vector = vec2D{
        x: partial_deriv(&tensor_field.xx, i, j, 1, &box_info)
            + partial_deriv(&tensor_field.yx, i, j, 0, &box_info),
        y: partial_deriv(&tensor_field.xy, i, j, 1, &box_info)
            + partial_deriv(&tensor_field.yy, i, j, 0, &box_info)};
    return vector;
}

pub fn lap(scalar_field: &Vec<Vec<f64>>,
                 i: i32, j: i32,
                 box_info: &BoxInfo) -> f64
{

    let BoxInfo { col_max: box_col_max,
                  row_max: box_row_max } = box_info;
    let ip = (i+1) % box_row_max;
    let im = (i - 1 + box_row_max) % box_row_max;
    let jp = (j+1) % box_col_max;
    let jm = (j - 1 + box_col_max) % box_col_max;
            // on the x axis
    let laplacian_value = 
        (
            2.0*(scalar_field[ip as usize][j as usize]  
                 + scalar_field[im as usize][j as usize]  
                 + scalar_field[i as usize][jp as usize]  
                 + scalar_field[i as usize][jm as usize])
            + scalar_field[ip as usize][jp as usize] 
            + scalar_field[im as usize][jm as usize] 
            + scalar_field[im as usize][jp as usize] 
            + scalar_field[ip as usize][jm as usize]
            - 12.0*scalar_field[i as usize][j as usize]
        )
        /(3.0*dx*dy);
                return laplacian_value;
}

pub fn laplacian(scalar_field: &ScalarField2D,
                 i: i32, j: i32,
                 box_info: &BoxInfo) -> f64
{
    let field = &scalar_field.s;
    let laplacian_value = lap(field,
                              i, j, &box_info);
    return laplacian_value;
}

pub fn laplacian_vector(vector_field: &VectorField2D,
                        i: i32, j: i32,
                        box_info: &BoxInfo) -> vec2D
{
    let vec = vec2D {
        x: lap(&vector_field.x,
                     i, j,
                     &box_info),
        y: lap(&vector_field.y,
                     i, j,
                     &box_info)};

                return vec;
}

pub fn grad_div_vel(vector_field: &VectorField2D,
                    i: i32, j: i32,
                    box_info: &BoxInfo) -> vec2D
{
    let BoxInfo { col_max: box_col_max,
                  row_max: box_row_max } = box_info;
    let ip = (i+1) % box_row_max;
    let im = (i - 1 + box_row_max) % box_row_max;
    let jp = (j+1) % box_col_max;
    let jm = (j - 1 + box_col_max) % box_col_max;

    let vec = vec2D{
        x: 
        (10.* vector_field.x[ip as usize][j as usize]
         + 10.* vector_field.x[im as usize][j as usize]
         + vector_field.x[ip as usize][jp as usize] 
         + vector_field.x[im as usize][jp as usize] 
         + vector_field.x[im as usize][jm as usize] 
         + vector_field.x[ip as usize][jm as usize]
	 -2.* vector_field.x[i as usize][jp as usize] 
         - 2.* vector_field.x[i as usize][jm as usize]      
         - 20.* vector_field.x[i as usize][j as usize])
            /(12.*dx*dy)	
        + ((vector_field.y[ip as usize][jp as usize]  
            - vector_field.y[im as usize][jp as usize])  
           + (vector_field.y[im as usize][jm as usize] 
              - vector_field.y[ip as usize][jm as usize]))
            /(4.0*dx*dy),


        y: 
        ( 10. * vector_field.y[i as usize][jp as usize]
          + 10. * vector_field.y[i as usize][jm as usize]
          + vector_field.y[ip as usize][jp as usize] 
          + vector_field.y[ip as usize][jm as usize]
          + vector_field.y[im as usize][jp as usize]
          + vector_field.y[im as usize][jm as usize]
          -2.* vector_field.y[ip as usize][j as usize] 
          -2.* vector_field.y[im as usize][j as usize]
          - 20.0*vector_field.y[i as usize][j as usize])
            /(12.0*dx*dy)
            + ((vector_field.x[ip as usize][jp as usize]
                - vector_field.x[im as usize][jp as usize])  
            + (vector_field.x[im as usize][jm as usize]
               - vector_field.x[ip as usize][jm as usize]))/(4.0*dx*dy)};

    return vec;
}

pub fn pressure(rho: f64, grad_rho: &vec2D, 
                lap_rho: f64, temp: f64) -> tens2D
{
    let p_thermo = rho * kB * temp/(1. - b * rho)
                   - aa * rho * rho;
    let p_iso = p_thermo - 0.5 * w * 
        (grad_rho.x * grad_rho.x + grad_rho.y*grad_rho.y)
        -  w*rho*lap_rho;
    
    let pressure = tens2D {
        xx: p_iso + w*grad_rho.x*grad_rho.x,
        xy: w * grad_rho.x * grad_rho.y,
	yx: w * grad_rho.x * grad_rho.y,
        yy: p_iso  + w * grad_rho.y * grad_rho.y
};
    return pressure;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scal_product_test() {
        let va = vec2D {x: 2., y: 1.};
        let vb = vec2D {x: 3., y: -2.};        
        assert_eq!(4., 
                   scal_product(&vb, &va));
    }
    #[test]
    fn tens_product_vec_test() {
        let t = tens2D {xx: 3., xy: 1.,
                        yx: -1., yy: 2.};
        let v = vec2D {x: 3., y: -2.};  

        let result = vec2D{x: 11., y: -1.};
        println!("{:?}", result);
        println!("{:?}", tens_product_vec(&t, &v));
        assert_eq!(result,
                   tens_product_vec(&t, &v));
    }

    #[test]
    fn grad_scalar_test() {
        // :todo:
        let mut v = vec2D{x: 0., y: 0.};

        let box_info = BoxInfo {
            col_max: 4,
            row_max: 4};
        
        // three tests:
        // first one with an uniform field

        let uni_grid = vec![vec![1.; 4]; 4];
        
        v = grad_scalar(&uni_grid, 2, 2, &box_info);
        assert_eq!(v.x, 0., "testing the gradient is zero \
                             with a uniform field");
        assert_eq!(v.y, 0., "testing the gradient is zero \
                             with a uniform field");

        // second one with a x growing field

        let grid_x = vec![vec![1., 2., 3., 4.]; 4];

        v = grad_scalar(&grid_x, 2, 2, &box_info);
        assert!((v.x > 0.), "testing the gradient.x is positive \
                            for a field growing with respect of x");
        assert_eq!(v.y, 0., "testing the gradient.y is zero \
                            for a field growing with respect of x");

        // second one with a y growing field

        let grid_y = vec![vec![1., 1., 1., 1.],
                          vec![2., 2., 2., 2.],
                          vec![3., 3., 3., 3.],
                          vec![4., 4., 4., 4.]];

        v = grad_scalar(&grid_y, 2, 2, &box_info);
        println!("y {:?}", v);
        assert!((v.y > 0.), "testing the gradient.y is positive \
                            for a field growing with respect of y");
        assert_eq!(v.x, 0., "testing the gradient.x is zero \
                            for a field growing with respect of y");
}
    // todo: test gradient fct
}
