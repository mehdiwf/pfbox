use super::mystructs::*;
use ndarray::prelude::*;

// an enum used to indicate the direction of a derivative calculation
// in 2D
pub enum DerivDirection {
    Rows,
    Columns}

pub fn shear_viscosity(rho: f64, mass: f64, dyn_visc: f64) -> f64
{
    dyn_visc * mass * rho
}

pub fn bulk_viscosity(rho: f64, mass: f64, cin_visc: f64) -> f64
{
    cin_visc * mass * rho
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

pub fn tens_product_vec(tens: &tens2D, 
                    vec: &vec2D) -> vec2D
{
    let product = vec2D{
        x: tens.xx * vec.x + tens.yx * vec.y,
        y: tens.xy * vec.x + tens.yy * vec.y};
    return product;
}

pub fn dyadic_product(tens_a: &tens2D, tens_b: &tens2D) -> f64
{
    let dy_product: f64 = 
        tens_a.xx * tens_b.xx + tens_a.yx * tens_b.yx
        + tens_a.xy * tens_b.xy + tens_a.yy * tens_b.yy;
    return dy_product;
}

/// returns the partial derivative of a scalar field
/// in i,j and in one direction
/// partial derivative in the x direction is DerivDirection::Rows
/// partial derivative in the y direction is DerivDirection::Columns
pub fn partial_deriv(a: &ScalarField2D,
                     i: i32, j: i32,
                     direction: DerivDirection,
                     lambda: f64,
                     box_info: &BoxInfo) -> f64
{
    let &BoxInfo { col_max: box_col_max,
                  row_max: box_row_max,
                  col_dx: dx,
                  row_dx: dy} = box_info;
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
            // on the x axis on samy's code
            DerivDirection::Rows => {
   let derivative = 
          lambda*(a.get_pos(ip, j) - a.get_pos(im, j))/(2.*dx)
	  + 0.25*lambda*(a.get_pos(ip, jp) - a.get_pos(im, jp))/(2.*dx)
          + 0.25*lambda*(a.get_pos(ip, jm) - a.get_pos(im, jm))/(2.*dx);
                return derivative;
            },
            // on the y axis on samy's code
            DerivDirection::Columns => {
   let derivative = 
          lambda*(a.get_pos(i, jp) - a.get_pos(i, jm))/(2.*dx)
	  + 0.25*lambda*(a.get_pos(ip, jp) - a.get_pos(ip, jm))/(2.*dx)	
          + 0.25*lambda*(a.get_pos(im, jp) - a.get_pos(im, jm))/(2.*dx);
                return derivative;
            }}}

pub fn grad_scalar(scalar_field: &ScalarField2D,
                   i: i32, j: i32,
                   lambda: f64,
                   box_info: &BoxInfo) -> vec2D
{
    let grad = vec2D {
        x: partial_deriv(&scalar_field, i, j,
                         DerivDirection::Columns, lambda, &box_info),
        y: partial_deriv(&scalar_field, i, j,
                         DerivDirection::Rows, lambda, &box_info),};
        
    return grad;
}

pub fn gradient(scalar_field: &ScalarField2D,
                i: i32, j: i32,
                lambda: f64,
                box_info: &BoxInfo) -> vec2D
{
    return grad_scalar(&scalar_field, i, j,
                       lambda, &box_info);
}

pub fn gradient_vector(vector_field: &VectorField2D,
                       i: i32, j: i32,
                       lambda: f64,
                       box_info: &BoxInfo) -> tens2D
{
    let tens = tens2D {
        xx: partial_deriv(&vector_field.x, i, j,
                          DerivDirection::Columns, lambda, &box_info),
        xy: partial_deriv(&vector_field.y, i, j,
                          DerivDirection::Columns, lambda, &box_info),
        yx: partial_deriv(&vector_field.x, i, j,
                          DerivDirection::Rows, lambda, &box_info),
        yy: partial_deriv(&vector_field.y, i, j,
                          DerivDirection::Rows, lambda, &box_info)};
    return tens;
}

pub fn div_vector(vector_field: &VectorField2D,
                  i: i32, j: i32,
                  lambda: f64,
                  box_info: &BoxInfo) -> f64
{
    let dVx_dx = partial_deriv(&vector_field.x, i, j,
                               DerivDirection::Columns, lambda, &box_info);
    let dVy_dy = partial_deriv(&vector_field.y, i, j,
                               DerivDirection::Rows, lambda, &box_info);

    return dVx_dx + dVy_dy;
}

pub fn div_tensor(tensor_field: &TensorField2D,
                  i: i32, j: i32,
                  lambda: f64,
                  box_info: &BoxInfo) -> vec2D
{
    let vector = vec2D{
        x: partial_deriv(&tensor_field.xx, i, j,
                         DerivDirection::Columns, lambda, &box_info)
            + partial_deriv(&tensor_field.yx, i, j,
                            DerivDirection::Rows, lambda, &box_info),
        
        y: partial_deriv(&tensor_field.xy, i, j,
                         DerivDirection::Columns, lambda,
                         &box_info)
            + partial_deriv(&tensor_field.yy, i, j,
                            DerivDirection::Rows, lambda,
                            &box_info)};
    return vector;
}

pub fn lap_scalar(a: &ScalarField2D,
                  i: i32, j: i32,
                  lambda: f64,
                  box_info: &BoxInfo) -> f64
{

    let BoxInfo { col_max: box_col_max,
                  row_max: box_row_max,
                  col_dx: dx,
                  row_dx: dy} = *box_info;
    let ip = ((i+1) % box_row_max) as usize;
    let im = ((i - 1 + box_row_max) % box_row_max) as usize;
    let jp = ((j+1) % box_col_max) as usize;
    let jm = ((j - 1 + box_col_max) % box_col_max) as usize;
    let (i, j) = (i as usize, j as usize);
            // on the x axis
    let laplacian_value = 
        (
            2.0*(a.get_pos(ip, j) + a.get_pos(im, j)  
                 + a.get_pos(i, jp) + a.get_pos(i, jm))
            + a.get_pos(ip, jp) + a.get_pos(im, jm) 
            + a.get_pos(im, jp) + a.get_pos(ip, jm)
            - 12.0*a.get_pos(i, j)
        )
        /(3.0*dx*dy);
                return laplacian_value;
}

pub fn laplacian(scalar_field: &ScalarField2D,
                 i: i32, j: i32,
                 lambda: f64,
                 box_info: &BoxInfo) -> f64
{
    let field = &scalar_field;
    let laplacian_value = lap_scalar(field,
                                     i, j, lambda, &box_info);
    return laplacian_value;
}

pub fn laplacian_vector(vector_field: &VectorField2D,
                        i: i32, j: i32,
                        lambda: f64,
                        box_info: &BoxInfo) -> vec2D
{
    let vec = vec2D {
        x: lap_scalar(&vector_field.x,
                      i, j,
                      lambda, &box_info),
        y: lap_scalar(&vector_field.y,
                      i, j,
                      lambda, &box_info)};

                return vec;
}

pub fn grad_div_vel(v: &VectorField2D,
                    i: i32, j: i32,
                    lambda: f64,
                    box_info: &BoxInfo) -> vec2D
{
    let BoxInfo { col_max: box_col_max,
                  row_max: box_row_max,
                  col_dx: dx,
                  row_dx: dy} = *box_info;
    
    let ip = ((i+1) % box_row_max) as usize;
    let im = ((i - 1 + box_row_max) % box_row_max) as usize;
    let jp = ((j+1) % box_col_max) as usize;
    let jm = ((j - 1 + box_col_max) % box_col_max) as usize;
    let (i, j) = (i as usize, j as usize);
    
    let vec = vec2D{
        x: 
        (10.* v.x.get_pos(ip, j) + 10.* v.x.get_pos(im, j)
         + v.x.get_pos(ip, jp) + v.x.get_pos(im, jp) 
         + v.x.get_pos(im, jm) + v.x.get_pos(ip, jm)
	 - 2.* v.x.get_pos(i, jp) - 2.* v.x.get_pos(i, jm)
         - 20.* v.x.get_pos(i, j))
            /(12.*dx*dy)
        + ((v.y.get_pos(ip, jp) - v.y.get_pos(im, jp))  
           + (v.y.get_pos(im, jm) 
           - v.y.get_pos(ip, jm)))
            /(4.0*dx*dy),


        y: 
        (10.* v.y.get_pos(i, jp) + 10.* v.y.get_pos(i, jm)
         + v.y.get_pos(ip, jp) + v.y.get_pos(ip, jm)
         + v.y.get_pos(im, jp) + v.y.get_pos(im, jm)
         - 2.* v.y.get_pos(ip, j) -2.* v.y.get_pos(im, j)
         - 20.0*v.y.get_pos(i, j))
            /(12.0*dx*dy)
            
          + ((v.x.get_pos(ip, jp) - v.x.get_pos(im, jp))
             + (v.x.get_pos(im, jm) - v.x.get_pos(ip, jm)))
            /(4.0*dx*dy)};

    return vec;
}

pub fn pressure(rho: f64, grad_rho: &vec2D, 
                lap_rho: f64, temp: f64,
                kB: f64, aa: f64, b: f64,
                w: f64) -> tens2D
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
            row_max: 4,
            col_dx: 0.1,
            row_dx: 0.1};
        
        // three tests:
        // first one with an uniform field

        let uni_grid = ScalarField2D { s: Array::<f64, Ix2>::ones((4, 4).f()) };
        
        v = grad_scalar(&uni_grid, 2, 2, 0.66, &box_info);
        assert_eq!(v.x, 0., "testing the gradient is zero \
                             with a uniform field");
        assert_eq!(v.y, 0., "testing the gradient is zero \
                             with a uniform field");

        // second one with a x growing field

        let grid_x = ScalarField2D {
            s:array![[1., 2., 3., 4.],
                     [1., 2., 3., 4.],
                     [1., 2., 3., 4.],
                     [1., 2., 3., 4.]]};

        v = grad_scalar(&grid_x, 2, 2, 0.66, &box_info);
        assert!((v.y > 0.), "testing the gradient.y is positive \
                            for a field growing with respect of y");
        assert_eq!(v.x, 0., "testing the gradient.x is zero \
                            for a field growing with respect of y");

        // second one with a y growing field

        let grid_y = ScalarField2D{
            s:array![[1., 1., 1., 1.],
                     [2., 2., 2., 2.],
                     [3., 3., 3., 3.],
                     [4., 4., 4., 4.]]};

        v = grad_scalar(&grid_y, 2, 2, 0.66, &box_info);
        println!("y {:?}", v);
        assert!((v.x > 0.), "testing the gradient.x is positive \
                            for a field growing with respect of x");
        assert_eq!(v.y, 0., "testing the gradient.y is zero \
                            for a field growing with respect of x");
}
    // todo: test gradient fct
}
