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

    pub fn set_masses(&mut self, masses: [f64; 4]) {
        self.masses = masses;
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

    fn th_max(&self, idx: usize) -> Option<f64> {
        assert!(idx == 2 || idx == 3);

        let m = self.masses[idx];
        let pcm_f = self.pcm_f_2().sqrt();
        let chi = self.chi();

        if m * chi.sinh() > pcm_f {
            Some((pcm_f / (m * chi.sinh())).asin().to_degrees())
        } else {
            None
        }
    }

    pub fn th2_max(&self) -> Option<f64> {
        self.th_max(2)
    }

    pub fn th3_max(&self) -> Option<f64> {
        self.th_max(3)
    }

    fn thcm_to_th(&self, thcm:f64, idx: usize) -> f64 {
        assert!(idx == 2 || idx == 3);

        let pcm_f_2 = self.pcm_f_2();
        let m_2 = self.masses[idx].powi(2);
        let chi = self.chi();
        let pm = match idx { 2 => 1.0, 3 => -1.0, _ => unreachable!()};

        f64::atan2(
            f64::sin(thcm.to_radians()),
            pm * f64::cos(thcm.to_radians()) * f64::cosh(chi) + f64::sqrt(1.0 + m_2 / pcm_f_2) * f64::sinh(chi)
        ).to_degrees()
    }

    pub fn thcm_to_th2(&self, thcm: f64) -> f64 {
        self.thcm_to_th(thcm, 2)
    }

    pub fn thcm_to_th3(&self, thcm: f64) -> f64 {
        self.thcm_to_th(thcm, 3)
    }

    fn thcm_to_e(&self, thcm: f64, idx: usize) -> f64 {
        assert!(idx == 2 || idx == 3);

        let pcm_f_2 = self.pcm_f_2();
        let pcm_f = pcm_f_2.sqrt();
        let m_2 = self.masses[idx].powi(2);
        let chi = self.chi();
        let pm = match idx { 2 => 1.0, 3 => -1.0, _ => unreachable!()};

        f64::sqrt(pcm_f_2 + m_2) * f64::cosh(chi) + pm * pcm_f * f64::cos(thcm.to_radians()) * f64::sinh(chi)
    }

    fn thcm_to_k(&self, thcm: f64, idx: usize) -> f64 {
        assert!(idx == 2 || idx == 3);

        let e = self.thcm_to_e(thcm, idx);
        let m = self.masses[idx];

        e - m
    }

    pub fn thcm_to_k2(&self, thcm: f64) -> f64 {
        self.thcm_to_k(thcm, 2)
    }

    pub fn thcm_to_k3(&self, thcm: f64) -> f64 {
        self.thcm_to_k(thcm, 3)
    }

    fn th_to_p_pm(&self, th: f64, idx: usize, pm: f64) -> Option<f64> {
        assert!(idx == 2 || idx == 3);
        assert_eq!(pm.abs(), 1.0);

        let pcm_f_2 = self.pcm_f_2();
        let m = self.masses[idx];
        let m_2 = m.powi(2);
        let chi = self.chi();
        let th_max = self.th_max(idx);

        match th_max {
            Some(th_max) if th > th_max => None,
            _ => Some(
                    (f64::sqrt(m_2 + pcm_f_2) * f64::cos(th.to_radians()) * f64::sinh(chi) + pm * f64::cosh(chi) * f64::sqrt(pcm_f_2 - m_2 * f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2))) /
                    (1.0 + f64::sin(th.to_radians()).powi(2) * f64::sinh(chi).powi(2))
                ),
        }
    }

    fn th_to_e_pm(&self, th: f64, idx: usize, pm: f64) -> Option<f64> {
        let p = self.th_to_p_pm(th, idx, pm);
        let m = self.masses[idx];

        p.map(|p| f64::sqrt(p.powi(2) + m.powi(2)))
    }

    fn th_to_k_pm(&self, th: f64, idx: usize, pm: f64) -> Option<f64> {
        let e = self.th_to_e_pm(th, idx, pm);
        let m = self.masses[idx];

        e.map(|e| e - m)
    }

    pub fn th2_to_k2_plus(&self, th: f64) -> Option<f64> {
        self.th_to_k_pm(th, 2, 1.0)
    }

    pub fn th2_to_k2_minus(&self, th: f64) -> Option<f64> {
        self.th_to_k_pm(th, 2, -1.0)
    }

    pub fn th3_to_k3_plus(&self, th: f64) -> Option<f64> {
        self.th_to_k_pm(th, 3, 1.0)
    }

    pub fn th3_to_k3_minus(&self, th: f64) -> Option<f64> {
        self.th_to_k_pm(th, 3, -1.0)
    }

    fn th_to_k(&self, th: f64, idx: usize) -> (Option<f64>, Option<f64>) {
        let th_max = self.th_max(idx);

        match th_max {
            None => (self.th_to_k_pm(th, idx, 1.0), None),
            Some(th_max) if th <= th_max => (self.th_to_k_pm(th, idx, 1.0), self.th_to_k_pm(th, idx, -1.0)),
            _ => (None, None),
        }
    }

    pub fn th2_to_k2(&self, th: f64) -> (Option<f64>, Option<f64>) {
        self.th_to_k(th, 2)
    }

    pub fn th3_to_k3(&self, th: f64) -> (Option<f64>, Option<f64>) {
        self.th_to_k(th, 3)
    }
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
