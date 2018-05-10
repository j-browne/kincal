#![feature(try_from)]
pub mod particle;
mod error;

use particle::{Particle, Outgoing};
pub use error::Error;

#[derive(Debug)]
enum PlusMinus {
    Plus,
    Minus,
}

#[derive(Debug)]
pub enum Value {
    Zero,
    One(f64),
    Two(f64, f64),
}

impl Value {
    pub fn map<F: Fn(f64) -> f64>(&self, f: F) -> Value {
        use Value::{Zero, One, Two};
        match self {
            Zero => Zero,
            One(x) => One(f(*x)),
            Two(x1, x2) => Two(f(*x1), f(*x2)),
        }
    }
}

impl From<f64> for Value {
    fn from(f: f64) -> Self {
        Value::One(f)
    }
}

impl From<(f64, f64)> for Value {
    fn from(f: (f64, f64)) -> Self {
        Value::Two(f.0, f.1)
    }
}

impl From<Option<f64>> for Value {
    fn from(o: Option<f64>) -> Self {
        use Value::{Zero, One};
        match o {
            Some(f) => One(f),
            None => Zero,
        }
    }
}

impl From<(Option<f64>, Option<f64>)> for Value {
    fn from(o: (Option<f64>, Option<f64>)) -> Self {
        use Value::{Zero, One, Two};
        match o {
            (Some(f1), Some(f2)) => Two(f1, f2),
            (Some(f), None) => One(f),
            (None, Some(f)) => One(f),
            (None, None) => Zero,
        }
    }
}

// In `masses`, the indices refer to the following:
// 0: projectile
// 1: target
// 2: ejectile
// 3: recoil
#[derive(Debug)]
pub struct ReactionKinematics {
    masses: [f64; 4],
    energy: f64,
}

impl ReactionKinematics {
    pub fn new(masses: [f64; 4], energy: f64) -> ReactionKinematics {
        ReactionKinematics {
            masses,
            energy
        }
    }
    pub fn masses(&self) -> &[f64; 4] {
        &self.masses
    }

    pub fn mass(&self, part: Particle) -> f64 {
        self.masses[part.index()]
    }

    pub fn set_masses(&mut self, masses: [f64; 4]) {
        self.masses = masses;
    }

    pub fn set_mass(&mut self, part: Particle, mass: f64) {
        self.masses[part.index()] = mass;
    }

    pub fn energy(&self) -> f64 {
        self.energy
    }

    pub fn set_energy(&mut self, energy: f64) {
        self.energy = energy;
    }

    fn s(&self) -> f64 {
        let m0 = self.masses[0];
        let m1 = self.masses[1];
        let e = self.energy;

        (m0 + m1).powi(2) + 2.0 * m1 * e
    }

    fn pcm_i_2(&self) -> f64 {
        let m0_2 = self.masses[0].powi(2);
        let m1_2 = self.masses[1].powi(2);
        let s = self.s();

        ((s - m0_2 - m1_2).powi(2) - 4.0 * m0_2 * m1_2) / (4.0 * s)
    }

    fn pcm_f_2(&self) -> f64 {
        let m2_2 = self.masses[2].powi(2);
        let m3_2 = self.masses[3].powi(2);
        let s = self.s();

        ((s - m2_2 - m3_2).powi(2) - 4.0 * m2_2 * m3_2) / (4.0 * s)
    }

    fn chi(&self) -> f64 {
        let m1 = self.masses[1];
        let m1_2 = m1.powi(2);
        let pcm_i_2 = self.pcm_i_2();
        let pcm_i = pcm_i_2.sqrt();

        ((pcm_i + (m1_2 + pcm_i_2).sqrt())/ m1).ln()
    }

    fn thcm_to_e(&self, thcm: f64, out: Outgoing) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let pcm_f = pcm_f_2.sqrt();
        let m_2 = self.masses[out.index()].powi(2);
        let chi = self.chi();
        let pm = match out { Outgoing::Ejectile => 1.0, Outgoing::Recoil => -1.0 };

        f64::sqrt(pcm_f_2 + m_2) * f64::cosh(chi) + pm * pcm_f * f64::cos(thcm.to_radians()) * f64::sinh(chi)
    }

    pub fn thcm_to_k(&self, thcm: f64, out: Outgoing) -> f64 {
        let e = self.thcm_to_e(thcm, out);
        let m = self.masses[out.index()];

        e - m
    }

    fn th_to_e(&self, th: f64, out: Outgoing) -> Value {
        let p = self.th_to_p(th, out);
        let m = self.masses[out.index()];

        p.map(|p| f64::sqrt(p.powi(2) + m.powi(2)))
    }

    pub fn th_to_k(&self, th: f64, out: Outgoing) -> Value {
        let e = self.th_to_e(th, out);
        let m = self.masses[out.index()];

        e.map(|e| e - m)
    }

    pub fn th_max(&self, out: Outgoing) -> Option<f64> {
        let m = self.masses[out.index()];
        let pcm_f = self.pcm_f_2().sqrt();
        let chi = self.chi();

        if m * chi.sinh() > pcm_f {
            Some((pcm_f / (m * chi.sinh())).asin().to_degrees())
        } else {
            None
        }
    }

    fn th_to_p_pm(&self, th: f64, out: Outgoing, pm: PlusMinus) -> f64 {
        use PlusMinus::{Plus, Minus};
        let pcm_f_2 = self.pcm_f_2();
        let m = self.masses[out.index()];
        let m_2 = m.powi(2);
        let chi = self.chi();
        let pm = match pm { Plus => 1.0, Minus => -1.0 };

        (f64::sqrt(m_2 + pcm_f_2) * f64::cos(th.to_radians()) * f64::sinh(chi) + pm * f64::cosh(chi) * f64::sqrt(pcm_f_2 - m_2 * f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2))) /
        (1.0 + f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2))
    }

    pub fn th_to_p(&self, th: f64, out: Outgoing) -> Value {
        use PlusMinus::{Plus, Minus};
        use Value::{Zero, One, Two};
        let th_max = self.th_max(out);

        match th_max {
            None => Value::One(self.th_to_p_pm(th, out, Plus)),
            Some(th_max) if th == th_max => One(self.th_to_p_pm(th, out, Plus)),
            Some(th_max) if th < th_max => Two(self.th_to_p_pm(th, out, Plus),
                                               self.th_to_p_pm(th, out, Minus)),
            Some(_) => Zero,
        }
    }


    pub fn thcm_to_th(&self, thcm:f64, out: Outgoing) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let m_2 = self.masses[out.index()].powi(2);
        let chi = self.chi();
        let pm_particle = match out { Outgoing::Ejectile => 1.0, Outgoing::Recoil => -1.0 };

        f64::atan2(
            f64::sin(thcm.to_radians()),
            pm_particle * f64::cos(thcm.to_radians()) * f64::cosh(chi) + f64::sqrt(1.0 + m_2 / pcm_f_2) * f64::sinh(chi)
        ).to_degrees()
    }

    fn th_to_thcm_p(&self, th: f64, out: Outgoing, p: f64) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let m_2 = self.masses[out.index()].powi(2);
        let chi = self.chi();
        let pm = match out { Outgoing::Ejectile => 1.0, Outgoing::Recoil => -1.0 };

        f64::atan2(
            p * f64::sin(th.to_radians()) * f64::cosh(chi),
            pm * p * f64::cos(th.to_radians()) - f64::sqrt(pcm_f_2 + m_2)  * f64::sinh(chi)
        ).to_degrees()
    }

    pub fn th_to_thcm(&self, th: f64, out: Outgoing) -> Value {
        let p = self.th_to_p(th, out);

        p.map(|p| self.th_to_thcm_p(th, out, p))
    }


    pub fn thcm_to_th2(&self, thcm: f64) -> Value {
        self.thcm_to_th(thcm, Outgoing::Ejectile).into()
    }

    pub fn thcm_to_th3(&self, thcm: f64) -> Value {
        self.thcm_to_th(thcm, Outgoing::Recoil).into()
    }

    pub fn thcm_to_kcm(&self, _thcm: f64) -> Value {
        unimplemented!();
    }

    pub fn thcm_to_k2(&self, thcm: f64) -> Value {
        self.thcm_to_k(thcm, Outgoing::Ejectile).into()
    }

    pub fn thcm_to_k3(&self, thcm: f64) -> Value {
        self.thcm_to_k(thcm, Outgoing::Recoil).into()
    }

    pub fn th2_to_thcm(&self, th: f64) -> Value {
        self.th_to_thcm(th, Outgoing::Ejectile).into()
    }

    pub fn th2_to_th3(&self, th: f64) -> Value {
        let thcm = self.th_to_thcm(th, Outgoing::Ejectile);

        thcm.map(|thcm| self.thcm_to_th(thcm, Outgoing::Recoil))
    }

    pub fn th2_to_kcm(&self, _th: f64) -> Value {
        unimplemented!();
    }

    pub fn th2_to_k2(&self, th: f64) -> Value {
        self.th_to_k(th, Outgoing::Ejectile)
    }

    pub fn th2_to_k3(&self, th: f64) -> Value {
        let thcm = self.th_to_thcm(th, Outgoing::Ejectile);

        thcm.map(|thcm| self.thcm_to_k(thcm, Outgoing::Recoil))
    }

    pub fn th3_to_thcm(&self, th: f64) -> Value {
        self.th_to_thcm(th, Outgoing::Recoil).into()
    }

    pub fn th3_to_th2(&self, th: f64) -> Value {
        let thcm = self.th_to_thcm(th, Outgoing::Recoil);

        thcm.map(|thcm| self.thcm_to_th(thcm, Outgoing::Ejectile))
    }

    pub fn th3_to_kcm(&self, _th: f64) -> Value {
        unimplemented!();
    }

    pub fn th3_to_k2(&self, th: f64) -> Value {
        let thcm = self.th_to_thcm(th, Outgoing::Recoil);

        thcm.map(|thcm| self.thcm_to_k(thcm, Outgoing::Ejectile))
    }

    pub fn th3_to_k3(&self, th: f64) -> Value {
        self.th_to_k(th, Outgoing::Recoil)
    }
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}