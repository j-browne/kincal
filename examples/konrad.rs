extern crate kincal;
use kincal::{
    Outgoing::{Ejectile, Recoil}, ReactionKinematics, Value::{NoVal, OneVal, TwoVal},
};
use std::env;
use std::process::exit;

fn main() {
    let cmd_args: Vec<String> = env::args().collect();

    if cmd_args.len() != 5 {
        eprintln!("Usage: kincal <Projectile Name> <Projectile Energy> <Ejectile Theta> <Recoil Excitation Energy>");
        eprintln!();
        eprintln!("This calculates the properties of the recoil from a X(a,p)Y.");
        eprintln!("Allowed projectiles are \"34S\", \"34Cl\", and \"34Ar\".");
        eprintln!();
        eprintln!("Output: <Recoil Angle> <Recoil Energy>");
        exit(1)
    }

    let beam = &cmd_args[1];
    let e_in = if let Ok(e) = cmd_args[2].parse::<f64>() {
        e
    } else {
        eprintln!("failed to parse e_in");
        exit(1)
    };
    let theta_p = if let Ok(t) = cmd_args[3].parse::<f64>() {
        t
    } else {
        eprintln!("failed to parse theta_p");
        exit(1)
    };
    let e_ex = if let Ok(e) = cmd_args[4].parse::<f64>() {
        e
    } else {
        eprintln!("failed to parse e_ex");
        exit(1)
    };

    let m = match beam.as_str() {
        "34S" => [
            33.967_867_012 * 931.494_095,
            4.001_506_179 * 931.494_095,
            1.007_276_469 * 931.494_095,
            36.965_902_584 * 931.494_095 + e_ex,
        ],
        "34Cl" => [
            33.973_762_491 * 931.494_095,
            4.001_506_179 * 931.494_095,
            1.007_276_469 * 931.494_095,
            36.966_776_314 * 931.494_095 + e_ex,
        ],
        "34Ar" => [
            33.970_409_992 * 931.494_095,
            4.001_506_179 * 931.494_095,
            1.007_276_469 * 931.494_095,
            36.962_969_140 * 931.494_095 + e_ex,
        ],
        _ => panic!("Unknown beam species"),
    };
    let r = ReactionKinematics::new(m, e_in);

    match (
        r.thi_to_thj(theta_p, Ejectile, Recoil),
        r.thi_to_kj(theta_p, Ejectile, Recoil),
    ) {
        (NoVal, NoVal) => {}
        (OneVal(th3), OneVal(k3)) => {
            println!("{} {}", th3, k3);
        }
        (TwoVal(th3_1, th3_2), TwoVal(k3_1, k3_2)) => {
            println!("{} {}", th3_1, k3_1);
            println!("{} {}", th3_2, k3_2);
        }
        _ => eprintln!("th2_to_th3 and th2_to_k3 do not have to same number of results"),
    }
}
