#![feature(try_from)]
extern crate kincal;
use kincal::{Error, ReactionKinematics, Value};
use std::convert::TryInto;

fn main() -> Result<(), Error> {
    let m = [
        33.970_409_992 * 931.494_095,
        4.001_506_179 * 931.494_095,
        1.007_276_469 * 931.494_095,
        36.962_969_140 * 931.494_095 + 1.370_85,
    ];
    let r = ReactionKinematics::new(m, 54.190);

    for recoil_idx in 2..=3 {
        let ths: Vec<_> = (0..=180).map(f64::from).collect();
        let mut k_minus = Vec::with_capacity(ths.len());
        for th in &ths {
            use Value::{NoVal, OneVal, TwoVal};
            match r.th_to_k(*th, recoil_idx.try_into()?) {
                NoVal => {}
                OneVal(k) => {
                    println!("{} {}", th, k);
                }
                TwoVal(k1, k2) => {
                    println!("{} {}", th, k1);
                    k_minus.push((*th, k2));
                }
            }
        }
        for (th, k) in k_minus.iter().rev() {
            println!("{} {}", th, k);
        }
        println!();
    }
    Ok(())
}
