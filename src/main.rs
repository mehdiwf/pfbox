#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused)]
#![allow(dead_code)]

const PFBOX_VERSION: &str = env!("CARGO_PKG_VERSION");

mod maths;
mod procedures;
mod configfile;

use std::env; // to have a backtrace + input arguments
use std::io::Write; // to use "write_all" method
use std::fs; // to read/write a file contents

use procedures::vdw_simulation;
use configfile::cfg_io;
use configfile::cfg_struct;

// ConfigInput

fn main() {
    env::set_var("RUST_BACKTRACE", "1");

    println!("pfbox version {}\n", PFBOX_VERSION);

    let args: Vec<String> = env::args().collect();
    let mut config_input = cfg_struct::ConfigInput::Empty;
    let mut output_path = "";
    let mut overwrite = true;
    
    if args.len() < 2 {

        println!{"No input file given as argument, default config file\
        used. Type -h (:todo:) for more information\n"};

    }
    else
    {
        let inputfile = &args[1];
        config_input = cfg_struct::ConfigInput::Path(inputfile.to_string());
    }
    vdw_simulation::do_sim(config_input, overwrite);
} // main definition closing parenthesis
