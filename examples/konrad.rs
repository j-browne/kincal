extern crate kincal;
use std::env;
use std::process::exit;
use kincal::{ReactionKinematics, Value, Outgoing};

fn main() {
    let args: Vec<String> = env::args().collect();

    let e_in = match args[1].parse::<f64>() {
        Ok(e) => e,
        Err(_) => {
            eprintln!("failed to parse e_in");
            exit(1)
        }
    };
    let theta_p = match args[2].parse::<f64>() {
        Ok(t) => t,
        Err(_) => {
            eprintln!("failed to parse theta_p");
            exit(1)
        }
    };
    let e_ex = match args[3].parse::<f64>() {
        Ok(e) => e,
        Err(_) => {
            eprintln!("failed to parse e_ex");
            exit(1)
        }
    };
    let beam = &args[4];

    let m = match beam.as_str() {
        "34S" => [33.967867012 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.965902584 * 931.494095 + e_ex],
        "34Cl" => [33.973762491 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.966776314 * 931.494095 + e_ex],
        "34Ar" => [33.970409992 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.962969140 * 931.494095 + e_ex],
        _ => panic!("Unknown beam species"),
    };
    let r = ReactionKinematics::new(m, e_in);

    use Value::{Zero, One, Two};
    use Outgoing::{Ejectile, Recoil};
    match (r.thi_to_thj(theta_p, Ejectile, Recoil), r.thi_to_kj(theta_p, Ejectile, Recoil)) {
        (Zero, Zero) => {}
        (One(th3), One(k3)) => { println!("{} {}", th3, k3); }
        (Two(th3_1, th3_2), Two(k3_1, k3_2)) => { println!("{} {}", th3_1, k3_1); println!("{} {}", th3_2, k3_2); }
        _ => { eprintln!("th2_to_th3 and th2_to_k3 do not have to same number of results")}
    }
}