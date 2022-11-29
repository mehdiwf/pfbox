use serde_derive::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct PhyCfg {
    pub x_size: i32,
    pub dx: f32,
    pub temper0: f32,
    pub temper1: f32,
    pub rho_liq: f32,
    pub rho_vap: f32,
    pub dt: f32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct InitCfg {
    pub n_liq: f32
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
pub struct SimCfg {
    pub physics_config: PhyCfg,
    pub initial_time_config: InitCfg,
    pub save_config: SaveCfg
}

pub fn get_default_config() -> SimCfg
{
    let config_template = SimCfg {
        physics_config: PhyCfg
        { x_size: 20,
          dx: 0.1,
          temper0: 0.5,
          temper1: 0.7,
          rho_liq: 0.8,
          rho_vap: 0.02,
          dt: 0.01 },
        initial_time_config: InitCfg
        { n_liq: 40.0 },
        save_config: SaveCfg
        { directory_name: "test_simu".to_string(),
          user_comment: "".to_string(),
          histo_freq: 10000,
          histo_save: 10000,
          num_bin: 100 }};

    return config_template;
}

