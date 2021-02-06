// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::{Display, Formatter};
use std::ops::{Div, Mul};

/// A unit is a quantity written as (F^f)(L^l)(T^t), where
/// F, L and T are units of force, length and time respectively.
///
/// For instance, the units of velocity are: L^1/T^1, while
/// the units of acceleration are: L^1/T^2, which is better
/// written as L^1 T^(-2).
///
/// The reason why we do not use mass as a basic unit, and
/// instead jump directly to force is because wherever we
/// need to use some concept of mass, it comes in the form
/// of a force. That is, recall that force has units
/// M^1 L^1 T^-2 (mass * acceleration).
#[derive(Deserialize, Serialize, Clone, Copy, Default, Debug)]
pub struct Units {
    /// Exponent for force.
    pub f: f64,
    /// Exponent for length.
    pub l: f64,
    /// Exponent for time.
    pub t: f64,
}

impl PartialEq for Units {
    fn eq(&self, other: &Self) -> bool {
        self.f == other.f && self.l == other.l && self.t == other.t
    }
}

impl Units {
    /// Suppose we had units `u` of the form `F^1 T^-3 L^0`,
    /// and a quantity `q` with units `u`.
    ///
    /// Suppose we wanted to raise `q` to the power `n`.
    /// What units should `q` have? By the laws of
    /// exponentiation, it should have the units:
    /// `u^n` = `(F^1 T^-3 L^0)^n` = `F^(n) T^(-3n) L^(0n)`.
    ///
    /// In terms of code, `u` would be defined:
    /// ```
    /// let u = Units {
    ///     f: 1.0,
    ///     t: -3.0,
    ///     l: 0.0,
    /// };
    /// ```
    ///
    /// and `u.pow(n)` will give us the `u` to the nth
    /// power.
    /// ```
    /// let u = Units {
    ///     f: 1.0 * n,
    ///     t: -3.0 * n,
    ///     l: 0.0 * n,
    /// };
    /// ```
    fn pow(&self, exp: f64) -> Units {
        Units {
            f: exp * self.f,
            l: exp * self.l,
            t: exp * self.t,
        }
    }

    /// Take the inverse of units.
    fn inv(&self) -> Units {
        Units {
            f: -1.0 * self.f,
            l: -1.0 * self.l,
            t: -1.0 * self.t,
        }
    }

    /// Units for a unitless quantity.
    #[inline]
    fn unitless() -> Units {
        Units {
            f: 0.0,
            l: 0.0,
            t: 0.0,
        }
    }

    /// Units for a force quantity.
    #[inline]
    fn force() -> Units {
        Units {
            f: 1.0,
            l: 0.0,
            t: 0.0,
        }
    }

    /// Units for a length quantity.
    #[inline]
    fn length() -> Units {
        Units {
            f: 0.0,
            l: 1.0,
            t: 0.0,
        }
    }

    /// Units for a time quantity.
    #[inline]
    fn time() -> Units {
        Units {
            f: 0.0,
            l: 0.0,
            t: 1.0,
        }
    }

    /// Units for a diffusion quantity.
    #[inline]
    fn diffusion() -> Units {
        Units::length().pow(2.0) / Units::time()
    }

    /// Units for a stress quantity.
    #[inline]
    fn stress() -> Units {
        Units::force() / Units::length().pow(2.0)
    }

    /// Units for an inverse time quantity.
    #[inline]
    fn tinv() -> Units {
        Units::time().inv()
    }

    /// Units for a viscosity quantity.
    #[inline]
    fn viscosity() -> Units {
        Units::force() / (Units::length() / Units::time())
    }
}

/// Implement `Mul` trait for Units.
impl Mul for Units {
    type Output = Units;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Self::Output {
        Units {
            f: self.f + rhs.f,
            l: self.l + rhs.l,
            t: self.t + rhs.t,
        }
    }
}

/// Implement `Display` for Units.
impl Display for Units {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "F^{} L^{} T^{}", self.f, self.l, self.t)
    }
}

/// Implement `Div` for units.
impl Div for Units {
    type Output = Units;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Units) -> Self::Output {
        Units {
            f: self.f - rhs.f,
            l: self.l - rhs.l,
            t: self.t - rhs.t,
        }
    }
}

/// Implement a specification (trait) for a `Quantity` object.
pub trait Quantity {
    /// Return a quantity that whose number part is a multiple
    /// of this quantity.
    fn mul_number(&self, multiple: f64) -> Self;

    /// Return a quantity that is `10^-3` times the original.
    fn kilo(&self) -> Self;

    /// Return a quantity that is `10^-6` times the original.
    fn micro(&self) -> Self;

    /// Return a quantity that is `10^-9` times the original.
    fn nano(&self) -> Self;

    /// Return the number (`f64`) part of the quantity (not the units).
    fn number(&self) -> f64;

    /// Return the `Units` of the `Quantity`
    fn units(&self) -> Units;

    /// Convert to `General` quantity.
    fn g(&self) -> General;

    /// Return a quantity that is the `exp`th power of
    /// this quantity.
    fn pow(&self, exp: f64) -> General;
}

/// A general quantity.
#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct General {
    /// Numerical value of this quantity.
    n: f64,
    /// Units of this quantity.
    u: Units,
}

impl General {
    /// Convert the general quantity to a force, if possible.
    pub fn to_force(&self) -> Result<Force, String> {
        if self.units() == Units::force() {
            Ok(Force(self.number()))
        } else {
            Err(String::from(
                "Quantity does not have units of force.",
            ))
        }
    }

    /// Convert the general quantity to a diffusion, if possible.
    pub fn to_diffusion(&self) -> Result<Diffusion, String> {
        if self.units() == Units::diffusion() {
            Ok(Diffusion(self.number()))
        } else {
            Err(String::from(
                "Quantity does not have units of diffusion.",
            ))
        }
    }

    // pub fn to_viscosity(&self) -> Result<Viscosity, String> {
    //     if self.units() == Units::viscosity() {
    //         Ok(Viscosity(self.value()))
    //     } else {
    //         Err(format!(
    //             "Quantity ({}) does not have units of viscosity ({}).",
    //             self,
    //             Units::viscosity()
    //         ))
    //     }
    // }
}

impl Quantity for General {
    fn mul_number(&self, multiplier: f64) -> Self {
        General {
            n: self.n * multiplier,
            u: self.u,
        }
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.n
    }

    fn units(&self) -> Units {
        self.u
    }

    fn g(&self) -> General {
        *self
    }

    fn pow(&self, exp: f64) -> General {
        General {
            n: self.n.powf(exp),
            u: self.u.pow(exp),
        }
    }
}

impl Mul for General {
    type Output = General;

    fn mul(self, rhs: General) -> Self::Output {
        General {
            n: self.n * rhs.n,
            u: self.u * rhs.u,
        }
    }
}

impl Div for General {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        General {
            n: self.n / rhs.n,
            u: self.u / rhs.u,
        }
    }
}

impl Display for General {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.n, self.u)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Force(pub f64);

impl Quantity for Force {
    fn mul_number(&self, other: f64) -> Self {
        Force(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::force()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Length(pub f64);

impl Quantity for Length {
    fn mul_number(&self, other: f64) -> Self {
        Length(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::length()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Time(pub f64);

impl Quantity for Time {
    fn mul_number(&self, other: f64) -> Self {
        Time(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::time()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Tinv(pub f64);

impl Quantity for Tinv {
    fn mul_number(&self, other: f64) -> Self {
        Tinv(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::tinv()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Diffusion(pub f64);

impl Quantity for Diffusion {
    fn mul_number(&self, other: f64) -> Self {
        Diffusion(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::diffusion()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Stress(pub f64);

impl Quantity for Stress {
    fn mul_number(&self, other: f64) -> Self {
        Stress(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::stress()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct Viscosity(pub f64);

impl Quantity for Viscosity {
    fn mul_number(&self, other: f64) -> Self {
        Viscosity(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Viscosity {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::viscosity()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}

#[derive(Deserialize, Clone, Copy, PartialEq, Default, Debug)]
pub struct Unitless(pub f64);

impl Quantity for Unitless {
    fn mul_number(&self, other: f64) -> Self {
        Unitless(self.0 * other)
    }

    fn kilo(&self) -> Self {
        self.mul_number(1e3)
    }

    fn micro(&self) -> Self {
        self.mul_number(1e-6)
    }

    fn nano(&self) -> Self {
        self.mul_number(1e-9)
    }

    fn number(&self) -> f64 {
        self.0
    }

    fn units(&self) -> Units {
        Units::unitless()
    }

    fn g(&self) -> General {
        General {
            n: self.number(),
            u: self.units(),
        }
    }

    fn pow(&self, exp: f64) -> General {
        self.g().pow(exp)
    }
}
