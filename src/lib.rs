mod error;
mod particle;

pub use error::Error;
pub use particle::{Incoming, Outgoing, Particle};

#[derive(Debug, Clone, Copy)]
enum PlusMinus {
    Plus,
    Minus,
}

#[derive(Debug)]
pub enum Value {
    NoVal,
    OneVal(f64),
    TwoVal(f64, f64),
}

impl Value {
    pub fn map<F: Fn(f64) -> f64>(&self, f: F) -> Self {
        use Value::{NoVal, OneVal, TwoVal};
        match self {
            NoVal => NoVal,
            OneVal(x) => OneVal(f(*x)),
            TwoVal(x1, x2) => TwoVal(f(*x1), f(*x2)),
        }
    }
}

impl From<f64> for Value {
    fn from(f: f64) -> Self {
        Value::OneVal(f)
    }
}

impl From<(f64, f64)> for Value {
    fn from(f: (f64, f64)) -> Self {
        Value::TwoVal(f.0, f.1)
    }
}

impl From<Option<f64>> for Value {
    fn from(o: Option<f64>) -> Self {
        use Value::{NoVal, OneVal};
        match o {
            Some(f) => OneVal(f),
            None => NoVal,
        }
    }
}

impl From<(Option<f64>, Option<f64>)> for Value {
    fn from(o: (Option<f64>, Option<f64>)) -> Self {
        use Value::{NoVal, OneVal, TwoVal};
        match o {
            (Some(f1), Some(f2)) => TwoVal(f1, f2),
            (Some(f), None) | (None, Some(f)) => OneVal(f),
            (None, None) => NoVal,
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
    pub fn new(masses: [f64; 4], energy: f64) -> Self {
        Self { masses, energy }
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

        ((pcm_i + (m1_2 + pcm_i_2).sqrt()) / m1).ln()
    }

    pub fn th_max(&self, part: Outgoing) -> Option<f64> {
        let m = self.masses[part.index()];
        let pcm_f = self.pcm_f_2().sqrt();
        let chi = self.chi();

        if m * chi.sinh() > pcm_f {
            Some((pcm_f / (m * chi.sinh())).asin().to_degrees())
        } else {
            None
        }
    }

    pub fn thcm_to_th(&self, thcm: f64, part: Outgoing) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let m_2 = self.masses[part.index()].powi(2);
        let chi = self.chi();
        let pm_particle = match part {
            Outgoing::Ejectile => 1.0,
            Outgoing::Recoil => -1.0,
        };

        f64::atan2(
            f64::sin(thcm.to_radians()),
            pm_particle * f64::cos(thcm.to_radians()) * f64::cosh(chi)
                + f64::sqrt(1.0 + m_2 / pcm_f_2) * f64::sinh(chi),
        ).to_degrees()
    }

    fn thcm_to_e(&self, thcm: f64, part: Outgoing) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let pcm_f = pcm_f_2.sqrt();
        let m_2 = self.masses[part.index()].powi(2);
        let chi = self.chi();
        let pm = match part {
            Outgoing::Ejectile => 1.0,
            Outgoing::Recoil => -1.0,
        };

        f64::sqrt(pcm_f_2 + m_2) * f64::cosh(chi)
            + pm * pcm_f * f64::cos(thcm.to_radians()) * f64::sinh(chi)
    }

    pub fn thcm_to_k(&self, thcm: f64, part: Outgoing) -> f64 {
        let e = self.thcm_to_e(thcm, part);
        let m = self.masses[part.index()];

        e - m
    }

    fn th_to_p_pm(&self, th: f64, part: Outgoing, pm: PlusMinus) -> f64 {
        use PlusMinus::{Minus, Plus};
        let pcm_f_2 = self.pcm_f_2();
        let m = self.masses[part.index()];
        let m_2 = m.powi(2);
        let chi = self.chi();
        let pm = match pm {
            Plus => 1.0,
            Minus => -1.0,
        };

        (f64::sqrt(m_2 + pcm_f_2) * f64::cos(th.to_radians()) * f64::sinh(chi)
            + pm * f64::cosh(chi)
                * f64::sqrt(
                    pcm_f_2 - m_2 * f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2),
                )) / (1.0 + f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2))
    }

    fn th_to_p(&self, th: f64, part: Outgoing) -> Value {
        use PlusMinus::{Minus, Plus};
        use Value::{NoVal, OneVal, TwoVal};
        let th_max = self.th_max(part);

        #[cfg_attr(feature = "cargo-clippy", allow(float_cmp))]
        match th_max {
            None => OneVal(self.th_to_p_pm(th, part, Plus)),
            Some(th_max) if th == th_max => OneVal(self.th_to_p_pm(th, part, Plus)),
            Some(th_max) if th < th_max => TwoVal(
                self.th_to_p_pm(th, part, Plus),
                self.th_to_p_pm(th, part, Minus),
            ),
            Some(_) => NoVal,
        }
    }

    fn th_to_thcm_p(&self, th: f64, part: Outgoing, p: f64) -> f64 {
        let pcm_f_2 = self.pcm_f_2();
        let m_2 = self.masses[part.index()].powi(2);
        let chi = self.chi();
        let pm = match part {
            Outgoing::Ejectile => 1.0,
            Outgoing::Recoil => -1.0,
        };

        f64::atan2(
            p * f64::sin(th.to_radians()) * f64::cosh(chi),
            pm * p * f64::cos(th.to_radians()) - f64::sqrt(pcm_f_2 + m_2) * f64::sinh(chi),
        ).to_degrees()
    }

    pub fn th_to_thcm(&self, th: f64, part: Outgoing) -> Value {
        let p = self.th_to_p(th, part);

        p.map(|p| self.th_to_thcm_p(th, part, p))
    }

    fn th_to_e(&self, th: f64, part: Outgoing) -> Value {
        let p = self.th_to_p(th, part);
        let m = self.masses[part.index()];

        p.map(|p| f64::sqrt(p.powi(2) + m.powi(2)))
    }

    pub fn th_to_k(&self, th: f64, part: Outgoing) -> Value {
        let e = self.th_to_e(th, part);
        let m = self.masses[part.index()];

        e.map(|e| e - m)
    }

    pub fn thi_to_thj(&self, th_i: f64, part_i: Outgoing, part_j: Outgoing) -> Value {
        let thcm = self.th_to_thcm(th_i, part_i);

        thcm.map(|thcm| self.thcm_to_th(thcm, part_j))
    }

    pub fn thi_to_kj(&self, th_i: f64, part_i: Outgoing, part_j: Outgoing) -> Value {
        let thcm = self.th_to_thcm(th_i, part_i);

        thcm.map(|thcm| self.thcm_to_k(thcm, part_j))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
