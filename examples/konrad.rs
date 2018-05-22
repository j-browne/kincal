extern crate kincal;
use std::env;
use std::process::exit;
use kincal::{ReactionKinematics, Value, Outgoing};

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 5 {
        eprintln!("Usage: kincal <Projectile Name> <Projectile Energy> <Ejectile Theta> <Recoil Excitation Energy>");
        eprintln!();
        eprintln!("This calculates the properties of the recoil from a X(a,p)Y.");
        eprintln!("Allowed projectiles are \"34S\", \"34Cl\", and \"34Ar\".");
        eprintln!();
        eprintln!("Output: <Recoil Angle> <Recoil Energy>");
        exit(1)
    }

    let beam = &args[1];
    let e_in = match args[2].parse::<f64>() {
        Ok(e) => e,
        Err(_) => {
            eprintln!("failed to parse e_in");
            exit(1)
        }
    };
    let theta_p = match args[3].parse::<f64>() {
        Ok(t) => t,
        Err(_) => {
            eprintln!("failed to parse theta_p");
            exit(1)
        }
    };
    let e_ex = match args[4].parse::<f64>() {
        Ok(e) => e,
        Err(_) => {
            eprintln!("failed to parse e_ex");
            exit(1)
        }
    };

    let m = match beam.as_str() {
        "34S" => [33.967867012 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.965902584 * 931.494095 + e_ex],
        "34Cl" => [33.973762491 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.966776314 * 931.494095 + e_ex],
        "34Ar" => [33.970409992 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.962969140 * 931.494095 + e_ex],
        _ => panic!("Unknown beam species"),
    };
    let r = ReactionKinematics::new(m, e_in);

    use Value::{NoVal, OneVal, TwoVal};
    use Outgoing::{Ejectile, Recoil};
    match (r.thi_to_thj(theta_p, Ejectile, Recoil), r.thi_to_kj(theta_p, Ejectile, Recoil)) {
        (NoVal, NoVal) => {}
        (OneVal(th3), OneVal(k3)) => { println!("{} {}", th3, k3); }
        (TwoVal(th3_1, th3_2), TwoVal(k3_1, k3_2)) => { println!("{} {}", th3_1, k3_1); println!("{} {}", th3_2, k3_2); }
        _ => { eprintln!("th2_to_th3 and th2_to_k3 do not have to same number of results")}
    }
}
