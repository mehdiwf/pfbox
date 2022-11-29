#[derive(Debug, PartialEq)]
pub struct vec2D {pub x: f32, pub y: f32}

#[derive(Debug, PartialEq)]
pub struct int2D {pub i: i32, pub j: i32}

#[derive(Debug, PartialEq)]
pub struct tens2D {
    pub xx: f32,
    pub xy: f32,
    pub yx: f32,
    pub yy: f32}

// a scalar field is just "Vec<Vec<f32>>"... therefore not defined

#[derive(Debug, PartialEq)]
pub struct VectorField2D {
    pub x: Vec<Vec<f32>>,
    pub y: Vec<Vec<f32>>
}

#[derive(Debug, PartialEq)]
pub struct TensorField2D {
    pub xx: Vec<Vec<f32>>,
    pub xy: Vec<Vec<f32>>,
    pub yx: Vec<Vec<f32>>,
    pub yy: Vec<Vec<f32>>
}

#[derive(Debug, PartialEq)]
pub struct BoxInfo {
    pub imax: i32,
    pub jmax: i32}
