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
    pub s: Vec<Vec<f64>>,
}    

#[derive(Debug, PartialEq)]
pub struct VectorField2D {
    pub x: Vec<Vec<f64>>,
    pub y: Vec<Vec<f64>>
}

#[derive(Debug, PartialEq)]
pub struct TensorField2D {
    pub xx: Vec<Vec<f64>>,
    pub xy: Vec<Vec<f64>>,
    pub yx: Vec<Vec<f64>>,
    pub yy: Vec<Vec<f64>>
}

#[derive(Debug, PartialEq)]
pub struct BoxInfo {
    pub imax: i32,
    pub jmax: i32}

impl ScalarField2D {
    pub fn get_pos(&self, i: usize, j: usize) -> f64
    {
        return self.s[i][j];
    }

    pub fn set_pos(&mut self, i: usize, j: usize, f: &f64)
    {
        self.s[i][j] = *f;
    }
}

impl VectorField2D {
    pub fn get_pos(&self, i: usize, j: usize) -> vec2D
    {
        let vec_at_pos = vec2D {
            x: self.x[i][j],
            y: self.y[i][j]};
        return vec_at_pos;
    }

    pub fn set_pos(&mut self, i: usize, j: usize, vec: &vec2D)
    {
        self.x[i][j] = vec.x;
        self.y[i][j] = vec.y;
    }
}

impl TensorField2D {
    pub fn get_pos(&self, i: usize, j: usize) -> tens2D
    {
        let tens_at_pos = tens2D {
            xx: self.xx[i][j],
            xy: self.xy[i][j],
            yx: self.yx[i][j],
            yy: self.yy[i][j]};
        return tens_at_pos;
    }

    pub fn set_pos(&mut self, i: usize, j: usize, tens: &tens2D)
    {
        self.xx[i][j] = tens.xx;
        self.xy[i][j] = tens.xy;
        self.yx[i][j] = tens.yx;
        self.yy[i][j] = tens.yy;
    }
    
}
