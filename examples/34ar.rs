extern crate kincal;
use kincal::ReactionKinematics;

fn main() {
    let m = [33.970409992 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.962969140 * 931.494095 + 1.37085];
    let r = ReactionKinematics::new(m, 54.1900);

    let ths: Vec<_> = (0 ..= 180).map(|i| i as f64).collect();
    for th in ths.iter() {
        println!("{} {}", th, r.th2_to_k2_plus(*th).unwrap());
    }

    println!();

    let ths: Vec<_> = (0 ..= r.th3_max().unwrap().floor() as usize).map(|i| i as f64).collect();
    for th in ths.iter() {
        println!("{} {}", th, r.th3_to_k3_plus(*th).unwrap());
    }
    for th in ths.iter().rev() {
        println!("{} {}", th, r.th3_to_k3_minus(*th).unwrap());
    }
}