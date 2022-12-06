use serde_derive::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct PhyParam {
    pub dx: f64,
    pub dy: f64,
    pub dt: f64,
    pub max_sim_time: i32,
    pub temper0: f64,
    pub temper1: f64,
    pub rho_liq: f64,
    pub rho_vap: f64,
    pub Tc: f64,
    pub aa: f64,
    pub w: f64,
    pub b: f64,
    pub zeta0: f64,
    pub eta0: f64,
    pub m: f64,
    pub lambda0: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct InitCfg {
    pub n_liq: i32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SaveCfg {
    pub directory_name: String,
    pub user_comment: String,
    pub histo_freq: i32,
    pub histo_save: i32,
    pub num_bin: i32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PhyConstants {
    pub kB: f64,
    pub DeBroglie0: f64,
    pub dim: i32,
    pub lambda_grad: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SimCfg {
    pub physics_config: PhyParam,
    pub initial_time_config: InitCfg,
    pub save_config: SaveCfg,
    pub physics_constants: PhyConstants
}

#[derive(Debug)]
pub enum ConfigInput {
    Empty,
    Path(String)}

// pub fn get_default_config() -> SimCfg
// {
//     let config_template = SimCfg {
//         physics_config: PhyCfg
//         { x_size: 20,
//           dx: 0.1,
//           temper0: 0.5,
//           temper1: 0.7,
//           rho_liq: 0.8,
//           rho_vap: 0.02,
//           dt: 0.01 },
//         initial_time_config: InitCfg
//         { n_liq: 40.0 },
//         save_config: SaveCfg
//         { directory_name: "test_simu".to_string(),
//           user_comment: "".to_string(),
//           histo_freq: 10000,
//           histo_save: 10000,
//           num_bin: 100 }};

//     return config_template;
// }

