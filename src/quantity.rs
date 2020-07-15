// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use std::fmt;
use std::fmt::{Display, Formatter};
use std::ops::{Div, Mul};
use serde::Deserialize;

#[derive(Deserialize, Copy, Clone)]
pub struct Units {
    pub f: f32,
    pub l: f32,
    pub t: f32,
}

impl PartialEq for Units {
    fn eq(&self, other: &Self) -> bool {
        self.f == other.f && self.l == other.l && self.t == other.t
    }

    fn ne(&self, other: &Self) -> bool {
        self.f != other.f || self.l != other.l || self.t != other.t
    }
}

impl Units {
    fn pow(&self, exp: f32) -> Self {
        Units {
            f: exp * self.f,
            l: exp * self.l,
            t: exp * self.t,
        }
    }

    fn inv(&self) -> Self {
        Units {
            f: -1.0 * self.f,
            l: -1.0 * self.l,
            t: -1.0 * self.t,
        }
    }

    #[inline]
    fn unitless() -> Units {
        Units {
            f: 0.0,
            l: 0.0,
            t: 0.0,
        }
    }

    #[inline]
    fn force() -> Units {
        Units {
            f: 1.0,
            l: 0.0,
            t: 0.0,
        }
    }

    #[inline]
    fn length() -> Units {
        Units {
            f: 0.0,
            l: 1.0,
            t: 0.0,
        }
    }

    #[inline]
    fn time() -> Units {
        Units {
            f: 0.0,
            l: 0.0,
            t: 1.0,
        }
    }

    #[inline]
    fn diffusion() -> Units {
        Units::length().pow(2.0) / Units::time()
    }

    #[inline]
    fn stress() -> Units {
        Units::force() / Units::length().pow(2.0)
    }

    #[inline]
    fn tinv() -> Units {
        Units::time().inv()
    }

    #[inline]
    fn viscosity() -> Units {
        Units::force() / (Units::length() / Units::time())
    }
}

impl Mul for Units {
    type Output = Units;

    fn mul(self, rhs: Self) -> Self::Output {
        Units {
            f: self.f + rhs.f,
            l: self.l + rhs.l,
            t: self.t + rhs.t,
        }
    }
}

impl Display for Units {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "F^{} L^{} T^{}", self.f, self.l, self.t)
    }
}

impl Div for Units {
    type Output = Units;

    fn div(self, rhs: Self) -> Self::Output {
        Units {
            f: self.f - rhs.f,
            l: self.l - rhs.l,
            t: self.t - rhs.t,
        }
    }
}

pub trait Quantity {
    fn mulf(&self, other: f32) -> Self;

    fn kilo(&self) -> Self;

    fn micro(&self) -> Self;

    fn nano(&self) -> Self;

    fn value(&self) -> f32;

    fn units(&self) -> Units;

    /// Convert to `General` quantity.
    fn g(&self) -> General;

    fn pow(&self, exp: f32) -> General;
}

#[derive(Deserialize, Copy, Clone)]
pub struct General {
    v: f32,
    u: Units,
}

impl General {
    pub fn to_force(&self) -> Result<Force, String> {
        if self.units() == Units::force() {
            Ok(Force(self.value()))
        } else {
            Err(String::from("Quantity does not have units of force."))
        }
    }

    pub fn to_diffusion(&self) -> Result<Diffusion, String> {
        if self.units() == Units::diffusion() {
            Ok(Diffusion(self.value()))
        } else {
            Err(String::from("Quantity does not have units of diffusion."))
        }
    }

    pub fn to_viscosity(&self) -> Result<Viscosity, String> {
        if self.units() == Units::viscosity() {
            Ok(Viscosity(self.value()))
        } else {
            Err(format!(
                "Quantity ({}) does not have units of viscosity ({}).",
                self,
                Units::viscosity()
            ))
        }
    }
}

impl Quantity for General {
    fn mulf(&self, other: f32) -> Self {
        General {
            v: self.v * other,
            u: self.u.clone(),
        }
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.v
    }

    fn units(&self) -> Units {
        self.u.clone()
    }

    fn g(&self) -> General {
        self.clone()
    }

    fn pow(&self, exp: f32) -> General {
        General {
            v: self.v.powf(exp),
            u: self.u.pow(exp),
        }
    }
}

impl Mul for General {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        General {
            v: self.v * rhs.v,
            u: self.u * rhs.u,
        }
    }
}

impl Div for General {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        General {
            v: self.v / rhs.v,
            u: self.u / rhs.u,
        }
    }
}

impl Display for General {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.v, self.u)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Force(pub f32);

impl Quantity for Force {
    fn mulf(&self, other: f32) -> Self {
        Force(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::force()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Clone, Deserialize)]
pub struct Length(pub f32);

impl Quantity for Length {
    fn mulf(&self, other: f32) -> Self {
        Length(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::length()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Time(pub f32);

impl Quantity for Time {
    fn mulf(&self, other: f32) -> Self {
        Time(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::time()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Tinv(pub f32);

impl Quantity for Tinv {
    fn mulf(&self, other: f32) -> Self {
        Tinv(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::tinv()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Diffusion(pub f32);

impl Quantity for Diffusion {
    fn mulf(&self, other: f32) -> Self {
        Diffusion(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::diffusion()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Stress(pub f32);

impl Quantity for Stress {
    fn mulf(&self, other: f32) -> Self {
        Stress(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::stress()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Viscosity(pub f32);

impl Quantity for Viscosity {
    fn mulf(&self, other: f32) -> Self {
        Viscosity(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::viscosity()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Copy, Clone)]
pub struct Unitless(pub f32);

impl Quantity for Unitless {
    fn mulf(&self, other: f32) -> Self {
        Unitless(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mulf(1e3)
    }

    fn micro(&self) -> Self {
        self.mulf(1e-6)
    }

    fn nano(&self) -> Self {
        self.mulf(1e-9)
    }

    fn value(&self) -> f32 {
        self.0
    }

    fn units(&self) -> Units {
        Units::unitless()
    }

    fn g(&self) -> General {
        General {
            v: self.value(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f32) -> General {
        self.g().pow(exp)
    }
}
