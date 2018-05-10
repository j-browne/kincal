#![feature(try_from)]
extern crate kincal;
use std::convert::TryInto;
use kincal::{ReactionKinematics, Error, Value};

fn main() {
    match main_err() {
        Ok(_) => {}
        Err(e) => { println!("{:?}", e); }
    };
}

fn main_err() -> Result<(), Error> {
    let m = [33.970409992 * 931.494095, 4.001506179 * 931.494095, 1.007276469 * 931.494095, 36.962969140 * 931.494095 + 1.37085];
    let r = ReactionKinematics::new(m, 54.1900);

    for recoil_idx in 2..=3 {
        let ths: Vec<_> = (0 ..= 180).map(|i| i as f64).collect();
        let mut k_minus = Vec::with_capacity(ths.len());
        for th in ths.iter() {
            use Value::{Zero, One, Two};
            match r.th_to_k(*th, recoil_idx.try_into()?) {
                Zero => {}
                One(k) => { println!("{} {}", th, k); }
                Two(k1, k2) => { println!("{} {}", th, k1); k_minus.push((*th, k2)); }
            }
        }
        for (th, k) in k_minus.iter().rev() {
            println!("{} {}", th, k);
        }
        println!();
    }
    Ok(())
}