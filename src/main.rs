use std::env;
use std::time::Instant;
use rdm_twin_density::api::*;

fn print_usage(prog: &str) {
    println!("Usage: {} <command> [options]", prog);
    println!();
    println!("Commands:");
    println!("  constant              Show the computed 2C2 constant (Euler product)");
    println!("  predict <N>           Predict pi_2(N) using the density formula");
    println!("  verify <N>            Predict AND count actual twin primes using LEI scan");
    println!("  verify-range <N1> <N2> Count actual primes and twin primes in a specific range");
    println!("  export <N> <file>     Export the exact API(k) bit-collapse trace to a CSV file");
    println!();
    println!("Examples:");
    println!("  {} constant", prog);
    println!("  {} predict 1000000", prog);
    println!("  {} verify 100000", prog);
    println!("  {} verify-range 1000 10000", prog);
    println!("  {} export 100000 trace.csv", prog);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let prog = &args[0];

    if args.len() < 2 {
        print_usage(prog);
        return;
    }

    match args[1].as_str() {
        "constant" => {
            println!("============================================================");
            println!("  RDM(TM) TWIN PRIME DENSITY ENGINE -- 2C2 Euler Product");
            println!("============================================================");
            let start = Instant::now();
            let (two_c2, pcount) = compute_twin_prime_constant(1_000_000);
            println!("Computed 2C2 Constant       : {:.16}", two_c2);
            println!("Euler Product Primes Used   : {}", pcount);
            println!();
            println!("Step-by-step (first 15 primes):");
            let small_primes: [u64; 15] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
            let mut running = 1.0_f64;
            for &p in &small_primes {
                let pf = p as f64;
                let factor = pf * (pf - 2.0) / ((pf - 1.0) * (pf - 1.0));
                running *= factor;
                println!("  p={:<4} factor={:.10}  running C2={:.12}  2C2={:.12}", p, factor, running, 2.0 * running);
            }
            println!();
            println!("Final 2C2 (converged)       : {:.16}", two_c2);
            println!("Proof Basis                 : Trace-to-Lane Splitting -- Uchechukwu Second Law (LEI)");
            println!("Execution Time              : {:?}", start.elapsed());
            println!("============================================================");
        }
        "predict" => {
            if args.len() < 3 {
                println!("Error: predict requires N. Usage: {} predict <N>", prog);
                return;
            }
            let n: f64 = match args[2].parse() {
                Ok(v) if v >= 5.0 => v,
                _ => { println!("Error: N must be >= 5"); return; }
            };
            let start = Instant::now();
            let (two_c2, pcount) = compute_twin_prime_constant(1_000_000);
            let predicted = predict_twin_density(n, two_c2);
            let ln_n = n.ln();

            println!("============================================================");
            println!("  RDM(TM) TWIN PRIME DENSITY ENGINE -- Density Prediction");
            println!("============================================================");
            println!("Computed 2C2 Constant       : {:.16}", two_c2);
            println!("(via Euler product over {} primes)", pcount);
            println!("Input N                     : {}", n as u64);
            println!("ln(N)                       : {:.6}", ln_n);
            println!("Predicted pi_2(N)           : {:.2}", predicted);
            println!("Formula                     : 2C2 * N / (ln N)^2");
            println!("                            : {:.10} * {} / {:.6}", two_c2, n as u64, ln_n * ln_n);
            println!("Lane Configuration          : L- x L+ (RDI spectral projection)");
            println!("Proof Basis                 : Trace-to-Lane Splitting -- Uchechukwu Second Law (LEI)");
            println!("Execution Time              : {:?}", start.elapsed());
            println!("============================================================");
        }
        "verify" => {
            if args.len() < 3 {
                println!("Error: verify requires N. Usage: {} verify <N>", prog);
                return;
            }
            let n: f64 = match args[2].parse() {
                Ok(v) if v >= 5.0 => v,
                _ => { println!("Error: N must be >= 5"); return; }
            };
            let n_u64 = n as u64;
            // The engine is fully unrestricted. The user assumes the execution time for massive N.

            let start = Instant::now();
            let (two_c2, pcount) = compute_twin_prime_constant(1_000_000);
            let predicted = predict_twin_density(n, two_c2);

            let vstart = Instant::now();
            let actual = count_twin_primes_lei(n_u64);
            let vtime = vstart.elapsed();

            println!("============================================================");
            println!("  RDM(TM) TWIN PRIME DENSITY ENGINE -- Verified Prediction");
            println!("============================================================");
            println!("Computed 2C2 Constant       : {:.16}", two_c2);
            println!("(via Euler product over {} primes)", pcount);
            println!("Input N                     : {}", n_u64);
            println!("Predicted pi_2(N)           : {:.2}", predicted);
            println!("Actual pi_2(N) (LEI count)  : {}", actual);
            println!("Prediction Error            : {:.2}", (actual as f64 - predicted).abs());
            println!("Prediction Ratio            : {:.4}", actual as f64 / predicted);
            println!("Lane Configuration          : L- x L+ (RDI spectral projection)");
            println!("Proof Basis                 : Trace-to-Lane Splitting -- Uchechukwu Second Law (LEI)");
            println!("LEI Verification Time       : {:?}", vtime);
            println!("Total Execution Time        : {:?}", start.elapsed());
            println!("============================================================");
        }
        "verify-range" => {
            if args.len() < 4 {
                println!("Error: verify-range requires <N1> <N2>. Usage: {} verify-range <N1> <N2>", prog);
                return;
            }
            let n1: f64 = match args[2].parse() {
                Ok(v) if v >= 0.0 => v,
                _ => { println!("Error: N1 must be >= 0"); return; }
            };
            let n2: f64 = match args[3].parse() {
                Ok(v) if v >= n1 => v,
                _ => { println!("Error: N2 must be >= N1"); return; }
            };
            
            let start = Instant::now();
            let config = SelectedCounts {
                minus: false,
                plus: false,
                twin: true,
                cumulative: true,
            };
            
            // Run reconstruction for N2
            let res2 = run_reconstruction(n2 as u128, config);
            
            // Run reconstruction for N1 - 1
            let n1_minus_1 = if n1 as u128 > 0 { n1 as u128 - 1 } else { 0 };
            let res1 = run_reconstruction(n1_minus_1, config);
            
            let range_primes = res2.cumulative_primes.unwrap_or(0) - res1.cumulative_primes.unwrap_or(0);
            let range_twins = res2.twin_primes.unwrap_or(0) - res1.twin_primes.unwrap_or(0);
            
            println!("============================================================");
            println!("  RDM(TM) TWIN PRIME DENSITY ENGINE -- Range Verification");
            println!("============================================================");
            println!("Range                       : {} to {}", n1 as u128, n2 as u128);
            println!("Total Primes in Range       : {}", range_primes);
            println!("Twin Primes in Range        : {}", range_twins);
            println!("Execution Time              : {:?}", start.elapsed());
            println!("============================================================");
        }
        "export" => {
            if args.len() < 4 {
                println!("Error: export requires <N> <filename>. Usage: {} export <N> <filename>", prog);
                return;
            }
            let n: f64 = match args[2].parse() {
                Ok(v) if v >= 5.0 => v,
                _ => { println!("Error: N must be >= 5"); return; }
            };
            let filename = &args[3];
            let n_u64 = n as u128;
            
            let start = Instant::now();
            match export_trace_to_csv(n_u64, filename) {
                Ok(_) => {
                    println!("============================================================");
                    println!("  RDM(TM) TWIN PRIME DENSITY ENGINE -- Trace Export");
                    println!("============================================================");
                    println!("Input N                     : {}", n_u64);
                    println!("Export File                 : {}", filename);
                    println!("Status                      : SUCCESS ✅");
                    println!("Execution Time              : {:?}", start.elapsed());
                    println!("============================================================");
                }
                Err(e) => println!("Error exporting to CSV: {}", e),
            }
        }
        _ => {
            println!("Error: Unknown command '{}'", args[1]);
            print_usage(prog);
        }
    }
}
