use toml;

use std::fs; // to read/write a file contents
use std::fs::File;
use std::io::Write;
// use std::path::Path; // for path

use crate::configfile::cfg_struct::SimCfg;

pub fn read_cfg_file(path : &str) -> SimCfg
{
    let contents = fs::read_to_string(&path)
        .expect("Problem reading the file");

    let config: SimCfg = toml::from_str(&contents)
        .expect("Config file is not valid");
    return config;
}

pub fn write_cfg_file(config: &SimCfg, path: &str)
{
    let toml_export = toml::to_string(&config).unwrap();
    fs::write(&path, toml_export)
        .expect("Unable to write file");
}
