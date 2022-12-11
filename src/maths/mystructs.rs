use ndarray::prelude::*;

// an enum used to indicate the direction of a derivative calculation
// in 2D
pub enum DerivDirection {
    X_axis,
    Y_axis}

// a function able to produce a profile of a scalar field
pub fn x_profile(scalar_field: &ScalarField2D) -> Array1<f64>
{
    let profile = scalar_field.s.mean_axis(Axis(0))
        .expect("could not compute a scalar field profile");
    return profile;
}

#[derive(Debug, PartialEq)]
pub struct vec2D {pub x: f64, pub y: f64}

#[derive(Debug, PartialEq)]
pub struct int2D {pub i: i32, pub j: i32}

#[derive(Debug, PartialEq)]
pub struct tens2D {
    pub xx: f64,
    pub xy: f64,
    pub yx: f64,
    pub yy: f64}

#[derive(Debug, PartialEq)]
pub struct ScalarField2D {
    pub s: Array2<f64>,
}    

#[derive(Debug, PartialEq)]
pub struct VectorField2D {
    pub x: ScalarField2D,
    pub y: ScalarField2D
}

#[derive(Debug, PartialEq)]
pub struct TensorField2D {
    pub xx: ScalarField2D,
    pub xy: ScalarField2D,
    pub yx: ScalarField2D,
    pub yy: ScalarField2D
}

#[derive(Debug, PartialEq)]
pub struct BoxInfo {
    pub x_max: i32,
    pub y_max: i32,
    pub dx: f64,
    pub dy: f64}

#[derive(Debug, PartialEq)]
pub struct VectorProfile2D{
    pub x: Array1<f64>,
    pub y: Array1<f64>
}

#[derive(Debug, PartialEq)]
pub struct TensorProfile2D{
    pub xx: Array1<f64>,
    pub xy: Array1<f64>,
    pub yx: Array1<f64>,
    pub yy: Array1<f64>
}

impl ScalarField2D {
    pub fn new(x_size: usize, y_size: usize) -> ScalarField2D {
        return ScalarField2D {
        s: Array::<f64, Ix2>::zeros((y_size, x_size).f())};}

    pub fn new_growx(x_size: usize, y_size: usize, 
                     minvalue: f64, maxvalue: f64) -> ScalarField2D {
        let ones = Array::<f64, Ix2>::ones((y_size, x_size));
        let mult_array = Array::linspace(minvalue, maxvalue, x_size);
        let growx_array = ones*mult_array;
        return ScalarField2D {s: growx_array};}

    pub fn new_growy(x_size: usize, y_size: usize, 
                     minvalue: f64, maxvalue: f64) -> ScalarField2D {
        
        let growx = ScalarField2D::new_growx(x_size, y_size,
                                             minvalue, maxvalue);
        let growy = growx.s.t().to_owned();
        return ScalarField2D {s: growy};}

    pub fn get_pos(&self, x: usize, y: usize) -> f64
    {
        // println!("get x = {}, y = {}", x, y);
        return self.s[[y,x]];
    }

    pub fn set_pos(&mut self, x: usize, y: usize, f: f64)
    {
        // println!("set x = {}, y = {}", x, y);
        self.s[[y,x]] = f;
    }

    pub fn x_profile(&self) -> Array1<f64>
    {return x_profile(&self);}
}

impl VectorField2D {

    pub fn new(nrow: usize, ncol: usize) -> VectorField2D {
        return VectorField2D {
            x: ScalarField2D::new(nrow, ncol),
            y: ScalarField2D::new(nrow, ncol)}}

    pub fn get_pos(&self, x: usize, y: usize) -> vec2D
    {
        let vec_at_pos = vec2D {
            x: self.x.get_pos(x, y),
            y: self.y.get_pos(x, y)};
        return vec_at_pos;
    }

    pub fn set_pos(&mut self, x: usize, y: usize, vec: &vec2D)
    {
        self.x.set_pos(x, y, vec.x);
        self.y.set_pos(x, y, vec.y);
    }

    pub fn x_profile(&self) -> VectorProfile2D
    {
        let profile = VectorProfile2D
        {
            x: x_profile(&self.x),
            y: x_profile(&self.y)
        };
        return profile;
    }
}

impl TensorField2D {
    pub fn new(nrow: usize, ncol: usize) -> TensorField2D {
        return TensorField2D {
            xx: ScalarField2D::new(nrow, ncol),
            xy: ScalarField2D::new(nrow, ncol),
            yx: ScalarField2D::new(nrow, ncol),
            yy: ScalarField2D::new(nrow, ncol)}}

    pub fn get_pos(&self, x: usize, y: usize) -> tens2D
    {
        let tens_at_pos = tens2D {
            xx: self.xx.get_pos(x, y),
            xy: self.xy.get_pos(x, y),
            yx: self.yx.get_pos(x, y),
            yy: self.yy.get_pos(x, y)};
        return tens_at_pos;
    }

    pub fn set_pos(&mut self, x: usize, y: usize, tens: &tens2D)
    {
        self.xx.set_pos(x, y, tens.xx);
        self.xy.set_pos(x, y, tens.xy);
        self.yx.set_pos(x, y, tens.yx);
        self.yy.set_pos(x, y, tens.yy);
    }
    
    pub fn x_profile(&self) -> TensorProfile2D
    {
        let profile = TensorProfile2D
        {
            xx: x_profile(&self.xx),
            xy: x_profile(&self.xy),
            yx: x_profile(&self.yx),
            yy: x_profile(&self.yy)
        };
        return profile;
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x_profile_test() {
        let grid_y = ScalarField2D {
            s: array![[1., 1., 0., -1.],
                      [0., 0., 0., -2.],
                      [-1., 2., 0., 0.]]};
        let profile = array![0., 1., 0., -1.];
        assert_eq!(profile, 
                   x_profile(&grid_y));
        }

}
