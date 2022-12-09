use ndarray::prelude::*;

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
    pub col_max: i32,
    pub row_max: i32,
    pub col_dx: f64,
    pub row_dx: f64}

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
    pub fn new(nrow: usize, ncol: usize) -> ScalarField2D {
        return ScalarField2D {
        s: Array::<f64, Ix2>::zeros((nrow, ncol).f())};}

    pub fn get_pos(&self, i: usize, j: usize) -> f64
    {
        return self.s[[i,j]];
    }

    pub fn set_pos(&mut self, i: usize, j: usize, f: f64)
    {
        self.s[[i,j]] = f;
    }

    pub fn x_profile(&self) -> Array1<f64>
    {return x_profile(&self);}
}

impl VectorField2D {

    pub fn new(nrow: usize, ncol: usize) -> VectorField2D {
        return VectorField2D {
            x: ScalarField2D::new(nrow, ncol),
            y: ScalarField2D::new(nrow, ncol)}}

    pub fn get_pos(&self, i: usize, j: usize) -> vec2D
    {
        let vec_at_pos = vec2D {
            x: self.x.get_pos(i, j),
            y: self.y.get_pos(i, j)};
        return vec_at_pos;
    }

    pub fn set_pos(&mut self, i: usize, j: usize, vec: &vec2D)
    {
        self.x.set_pos(i, j, vec.x);
        self.y.set_pos(i, j, vec.y);
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

    pub fn get_pos(&self, i: usize, j: usize) -> tens2D
    {
        let tens_at_pos = tens2D {
            xx: self.xx.get_pos(i, j),
            xy: self.xy.get_pos(i, j),
            yx: self.yx.get_pos(i, j),
            yy: self.yy.get_pos(i, j)};
        return tens_at_pos;
    }

    pub fn set_pos(&mut self, i: usize, j: usize, tens: &tens2D)
    {
        self.xx.set_pos(i, j, tens.xx);
        self.xy.set_pos(i, j, tens.xy);
        self.yx.set_pos(i, j, tens.yx);
        self.yy.set_pos(i, j, tens.yy);
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
                      [-1., 2., 0., 0.]] };
        let profile = array![0., 1., 0., -1.];
        assert_eq!(profile, 
                   x_profile(&grid_y));
        }

}
