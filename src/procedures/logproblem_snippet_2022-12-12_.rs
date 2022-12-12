                let str_to_append = format!("IMPOSSIBLE COMPUTATION: NEGATIVE LOG\n\
                                             step {}, col={}, row={}\n\
                                             negative log of rho avoided!\n\
                                             rho value: {}\n\
                                             ------------\n",
                                            &step, &yi, &xi, &rho);
                logproblem_counts += 1;
                if print_logproblems {
                    println!("negative log detected ! check log for more info")};
                logfile.write_all(&str_to_append.as_bytes())
                    .expect("write failed");
                println!("error step {}:\n\
                          negative rho: rho = {}", step, rho);
