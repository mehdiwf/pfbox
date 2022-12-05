#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused)]
#![allow(dead_code)]

mod maths;
mod procedures;

use std::env; // to have a backtrace (search backtrace in this file)
use std::io::Write; // to use "write_all" method
use std::fs; // to read/write a file contents

use procedures::vdw_simulation;

fn main() {
    
    vdw_simulation::do_sim();
} // main definition closing parenthesis
