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
    pub col_max: i32,
    pub row_max: i32}

pub struct VectorProfile2D{
    x: Vec<f64>,
    y: Vec<f64>
}

pub fn x_profile(scalar_field: &Vec<Vec<f64>>) -> Vec<f64>
{
    let row_len = scalar_field.len();
    let col_len = scalar_field[0].len();
    let mut profile = vec![0.; col_len];
    for col in 0usize..col_len {
        let mut col_value = 0.;
        for row in 0usize..row_len
        {col_value += scalar_field[row][col];}
        profile[col] = col_value/(row_len as f64);}
    
    return profile;
}

impl ScalarField2D {
    pub fn get_pos(&self, i: usize, j: usize) -> f64
    {
        return self.s[i][j];
    }

    pub fn set_pos(&mut self, i: usize, j: usize, f: &f64)
    {
        self.s[i][j] = *f;
    }

    pub fn x_profile(&self) -> Vec<f64>
    {return x_profile(&self.s);}
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x_profile_test() {
        let grid_y = vec![vec![1., 1., 0., -1.],
                          vec![0., 0., 0., -2.],
                          vec![-1., 2., 0., 0.]];
        let profile = vec![0., 1., 0., -1.];
        assert_eq!(profile, 
                   x_profile(&grid_y));
        }

}
