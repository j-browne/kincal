use error::Error;
use std::convert::TryFrom;

#[derive(Debug, Clone, Copy)]
pub enum Particle {
    In(Incoming),
    Out(Outgoing),
}

impl Particle {
    pub fn index(&self) -> usize {
        use Particle::{In, Out};
        match self {
            In(p) => p.index(),
            Out(p) => p.index(),
        }
    }
}

impl From<Incoming> for Particle {
    fn from(p: Incoming) -> Self {
        Particle::In(p)
    }
}

impl From<Outgoing> for Particle {
    fn from(p: Outgoing) -> Self {
        Particle::Out(p)
    }
}

impl TryFrom<usize> for Particle {
    type Error = Error;
    fn try_from(idx: usize) -> Result<Self, Error> {
        use self::Incoming::{Projectile, Target};
        use self::Outgoing::{Ejectile, Recoil};
        use Particle::{In, Out};
        match idx {
            0 => Ok(In(Projectile)),
            1 => Ok(In(Target)),
            2 => Ok(Out(Ejectile)),
            3 => Ok(Out(Recoil)),
            _ => Err(Error::InvalidIndex),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Incoming {
    Projectile,
    Target,
}

impl Incoming {
    pub fn index(&self) -> usize {
        match self {
            Incoming::Projectile => 0,
            Incoming::Target => 1,
        }
    }
}

impl TryFrom<usize> for Incoming {
    type Error = Error;
    fn try_from(u: usize) -> Result<Self, Error> {
        match Particle::try_from(u)? {
            Particle::In(p) => Ok(p),
            Particle::Out(_) => Err(Error::UnexpectedOutgoing),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Outgoing {
    Ejectile,
    Recoil,
}

impl Outgoing {
    pub fn index(&self) -> usize {
        match self {
            Outgoing::Ejectile => 2,
            Outgoing::Recoil => 3,
        }
    }
}

impl TryFrom<usize> for Outgoing {
    type Error = Error;
    fn try_from(u: usize) -> Result<Self, Error> {
        match Particle::try_from(u)? {
            Particle::Out(p) => Ok(p),
            Particle::In(_) => Err(Error::UnexpectedIncoming),
        }
    }
}
