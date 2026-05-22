use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelectedCounts {
    pub minus: bool,
    pub plus: bool,
    pub twin: bool,
    pub cumulative: bool,
}

pub struct APIResults {
    pub limit_x: u128,
    pub k_limit: u128,
    pub minus_lane_primes: Option<u128>,
    pub plus_lane_primes: Option<u128>,
    pub twin_primes: Option<u128>,
    pub cumulative_primes: Option<u128>,
    pub exec_time: std::time::Duration,
}

pub fn run_reconstruction(x: u128, config: SelectedCounts) -> APIResults {
    let start_time = Instant::now();
    let k_limit = (x / 6) as usize;
    let a_limit = ((( (x as f64).sqrt() + 1.0 ) / 6.0) as usize) + 2;

    let need_minus = config.minus || config.twin || config.cumulative;
    let need_plus = config.plus || config.twin || config.cumulative;

    let mut xi_minus = if need_minus { vec![false; k_limit + 1] } else { Vec::new() };
    let mut xi_plus = if need_plus { vec![false; k_limit + 1] } else { Vec::new() };

    for a in 1..=a_limit {
        if need_plus {
            let denom1 = 6 * a + 1;
            if k_limit >= a {
                let b_max1 = (k_limit - a) / denom1;
                for b in a..=b_max1 {
                    let k = 6 * a * b + a + b;
                    if k <= k_limit { xi_plus[k] = true; }
                }
            }

            let denom2 = 6 * a - 1;
            if denom2 > 0 {
                let b_max2 = (k_limit + a) / denom2;
                for b in a..=b_max2 {
                    let k = 6 * a * b - a - b;
                    if k <= k_limit { xi_plus[k] = true; }
                }
            }
        }

        if need_minus {
            let denom3 = 6 * a + 1;
            let b_max3 = (k_limit + a) / denom3;
            for b in 1..=b_max3 {
                let k = 6 * a * b - a + b;
                if k <= k_limit { xi_minus[k] = true; }
            }

            let denom4 = 6 * a - 1;
            if denom4 > 0 && k_limit >= a {
                let b_max4 = (k_limit - a) / denom4;
                for b in 1..=b_max4 {
                    let k = 6 * a * b + a - b;
                    if k <= k_limit { xi_minus[k] = true; }
                }
            }
        }
    }

    let minus_lane_primes = if config.minus || config.cumulative {
        let mut comp_count = 0;
        for k in 1..=k_limit {
            if xi_minus[k] { comp_count += 1; }
        }
        Some((k_limit as u128) - comp_count)
    } else { None };

    let plus_lane_primes = if config.plus || config.cumulative {
        let mut comp_count = 0;
        let mut valid_k_limit = k_limit;
        if 6 * (k_limit as u128) + 1 > x {
            if valid_k_limit > 0 { valid_k_limit -= 1; }
        }
        for k in 1..=valid_k_limit {
            if xi_plus[k] { comp_count += 1; }
        }
        Some((valid_k_limit as u128) - comp_count)
    } else { None };

    let twin_primes = if config.twin {
        let mut comp_count = 0;
        let mut valid_k_limit = k_limit;
        if 6 * (k_limit as u128) + 1 > x {
            if valid_k_limit > 0 { valid_k_limit -= 1; }
        }
        for k in 1..=valid_k_limit {
            if xi_minus[k] || xi_plus[k] { comp_count += 1; }
        }
        let twin_count = if x >= 5 {
            (valid_k_limit as u128) - comp_count + 1
        } else { 0 };
        Some(twin_count)
    } else { None };

    let cumulative_primes = if config.cumulative {
        let m_primes = minus_lane_primes.unwrap();
        let p_primes = plus_lane_primes.unwrap();
        let count = if x < 2 { 0 } else if x < 3 { 1 } else if x < 5 { 2 } else { m_primes + p_primes + 2 };
        Some(count)
    } else { None };

    APIResults {
        limit_x: x,
        k_limit: k_limit as u128,
        minus_lane_primes,
        plus_lane_primes,
        twin_primes,
        cumulative_primes,
        exec_time: start_time.elapsed(),
    }
}

pub fn count_twin_primes_lei(limit: u64) -> u64 {
    let config = SelectedCounts {
        minus: false,
        plus: false,
        twin: true,
        cumulative: false,
    };
    let res = run_reconstruction(limit as u128, config);
    res.twin_primes.unwrap_or(0) as u64
}

pub fn compute_twin_prime_constant(limit: u64) -> (f64, usize) {
    let mut is_prime = vec![true; limit as usize + 1];
    is_prime[0] = false;
    if limit >= 1 { is_prime[1] = false; }
    let sqrt_limit = (limit as f64).sqrt() as usize;
    for i in 2..=sqrt_limit {
        if is_prime[i] {
            let mut j = i * i;
            while j <= limit as usize {
                is_prime[j] = false;
                j += i;
            }
        }
    }
    let mut primes = Vec::new();
    for i in 3..=limit as usize {
        if is_prime[i] {
            primes.push(i as u64);
        }
    }
    let mut running_c2 = 1.0_f64;
    for &p in &primes {
        let pf = p as f64;
        running_c2 *= (pf * (pf - 2.0)) / ((pf - 1.0) * (pf - 1.0));
    }
    (2.0 * running_c2, primes.len())
}

pub fn predict_twin_density(n: f64, two_c2: f64) -> f64 {
    let ln_n = n.ln();
    two_c2 * n / (ln_n * ln_n)
}

pub fn export_trace_to_csv(x: u128, filename: &str) -> Result<(), std::io::Error> {
    use std::fs::File;
    use std::io::Write;
    
    let k_limit = (x / 6) as usize;
    let a_limit = ((( (x as f64).sqrt() + 1.0 ) / 6.0) as usize) + 2;

    let mut xi_minus = vec![false; k_limit + 1];
    let mut xi_plus = vec![false; k_limit + 1];

    for a in 1..=a_limit {
        let denom1 = 6 * a + 1;
        if k_limit >= a {
            let b_max1 = (k_limit - a) / denom1;
            for b in a..=b_max1 {
                let k = 6 * a * b + a + b;
                if k <= k_limit { xi_plus[k] = true; }
            }
        }

        let denom2 = 6 * a - 1;
        if denom2 > 0 {
            let b_max2 = (k_limit + a) / denom2;
            for b in a..=b_max2 {
                let k = 6 * a * b - a - b;
                if k <= k_limit { xi_plus[k] = true; }
            }
        }

        let denom3 = 6 * a + 1;
        let b_max3 = (k_limit + a) / denom3;
        for b in 1..=b_max3 {
            let k = 6 * a * b - a + b;
            if k <= k_limit { xi_minus[k] = true; }
        }

        let denom4 = 6 * a - 1;
        if denom4 > 0 && k_limit >= a {
            let b_max4 = (k_limit - a) / denom4;
            for b in 1..=b_max4 {
                let k = 6 * a * b + a - b;
                if k <= k_limit { xi_minus[k] = true; }
            }
        }
    }

    let mut file = File::create(filename)?;
    writeln!(file, "k,L_minus_prime,L_plus_prime,Is_Twin")?;
    
    for k in 1..=k_limit {
        let l_minus_prime = !xi_minus[k];
        let l_plus_prime = !xi_plus[k];
        let is_twin = l_minus_prime && l_plus_prime;
        writeln!(file, "{},{},{},{}", k, l_minus_prime, l_plus_prime, is_twin)?;
    }
    
    Ok(())
}

