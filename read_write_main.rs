// to test it:
////
// cargo run -- ./data/tomlfile.toml
////

use std::collections::HashMap;
use toml;

use std::env; // to get cli arguments
use std::fs; // to read/write a file contents
use std::fs::File;
use std::io::Write;
use std::fmt; // format!
use std::path::Path; // for path

use crate::configfile::cfg_struct;
use crate::configfile::cfg_io;
pub mod configfile;

#[allow(unused_variables)]
#[allow(unused_assignments)]
#[allow(unused_imports)]
#[allow(dead_code)]

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {panic!(
        "Please enter an input file, \
         or create a default one using (:todo:)")}
    // dbg!(args);

    let inputfile = &args[1];
    let testpath = "/home/mehdi/workdir/dossiers/ilm/these/code_simulations/rust_implementation/pfbox/src/data/tomlfile.toml";

    let config = cfg_io::read_cfg_file(inputfile);

    println!("test print\n {}", config.physics_config.dx);

    cfg_io::write_cfg_file(&config, 
                           "./testoutput/testtomlexp.toml");

    println!("{:?}", config);

    let path = Path::new("./testoutput/append_test.txt");
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why),
        Ok(file) => file,
    };

    // open a file to append stuff to it
    let mut file = fs::OpenOptions::new()
        .write(true)
        .append(true) // This is needed to append to file
        .open(&path)
        .unwrap();

    // default config template
    let config_template = cfg_struct::get_default_config();
    

    println!("starting test");
    for i in 0..10
    {
        // formatting the string to append to the file
        let str_to_append = format!("to app {}\n", &i);
        // appending the string to 
        file.write_all(&str_to_append.as_bytes())
            .expect("write failed");
}
}
