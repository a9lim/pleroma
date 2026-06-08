//! Piecewise-linear rational functions used by thermography.

use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;

/// Exact rational division (panics only on a zero divisor, which the thermograph
/// recursion never produces).
pub(crate) fn rdiv(a: &Rational, b: &Rational) -> Rational {
    a.mul(&b.inv().expect("thermograph: division by zero"))
}

pub(crate) fn req(a: &Rational, b: &Rational) -> bool {
    a.cmp(b) == Ordering::Equal
}

/// A continuous piecewise-linear function of temperature `t >= 0`, stored as
/// ascending breakpoints `(t, value)`. It is constant (`= last value`) beyond
/// the final breakpoint, and the first breakpoint is always at `t = 0`.
#[derive(Clone, Debug)]
pub struct Pl {
    pub(crate) pts: Vec<(Rational, Rational)>,
}

impl Pl {
    /// The constant function `t -> v`.
    pub fn constant(v: Rational) -> Pl {
        Pl {
            pts: vec![(Rational::zero(), v)],
        }
    }

    /// The breakpoints `(t, value)`, ascending in `t`.
    pub fn points(&self) -> &[(Rational, Rational)] {
        &self.pts
    }

    /// Evaluate the wall at temperature `t` (constant-extended past the ends).
    pub fn value_at(&self, t: &Rational) -> Rational {
        let n = self.pts.len();
        if t.cmp(&self.pts[0].0) != Ordering::Greater {
            return self.pts[0].1.clone();
        }
        if t.cmp(&self.pts[n - 1].0) != Ordering::Less {
            return self.pts[n - 1].1.clone();
        }
        for w in self.pts.windows(2) {
            let (ta, va) = (&w[0].0, &w[0].1);
            let (tb, vb) = (&w[1].0, &w[1].1);
            if t.cmp(ta) != Ordering::Less && t.cmp(tb) != Ordering::Greater {
                let frac = rdiv(&t.sub(ta), &tb.sub(ta));
                return va.add(&vb.sub(va).mul(&frac));
            }
        }
        unreachable!("value_at fell through its segments")
    }

    /// Drop duplicate-`t` and collinear interior breakpoints so deep recursions
    /// (sums) don't accrete redundant nodes.
    pub(crate) fn cleaned(mut self) -> Pl {
        let mut dedup: Vec<(Rational, Rational)> = Vec::with_capacity(self.pts.len());
        for p in self.pts.drain(..) {
            match dedup.last() {
                Some(last) if req(&last.0, &p.0) => {}
                _ => dedup.push(p),
            }
        }

        let mut out: Vec<(Rational, Rational)> = Vec::with_capacity(dedup.len());
        for p in dedup {
            while out.len() >= 2 {
                let a = &out[out.len() - 2];
                let b = &out[out.len() - 1];
                let lhs = b.1.sub(&a.1).mul(&p.0.sub(&b.0));
                let rhs = p.1.sub(&b.1).mul(&b.0.sub(&a.0));
                if req(&lhs, &rhs) {
                    out.pop();
                } else {
                    break;
                }
            }
            out.push(p);
        }
        Pl { pts: out }
    }
}

fn merge_ts(f: &Pl, g: &Pl) -> Vec<Rational> {
    let mut ts: Vec<Rational> = Vec::new();
    for (t, _) in f.pts.iter().chain(g.pts.iter()) {
        ts.push(t.clone());
    }
    sort_dedup(&mut ts);
    ts
}

fn sort_dedup(ts: &mut Vec<Rational>) {
    ts.sort_by(|a, b| a.cmp(b));
    ts.dedup_by(|a, b| req(a, b));
}

fn cross(f: &Pl, g: &Pl, a: &Rational, b: &Rational) -> Option<Rational> {
    let (fa, fb) = (f.value_at(a), f.value_at(b));
    let (ga, gb) = (g.value_at(a), g.value_at(b));
    let denom = fb.sub(&fa).sub(&gb.sub(&ga));
    if denom.sign() == Ordering::Equal {
        return None;
    }
    let s = rdiv(&ga.sub(&fa), &denom);
    if s.sign() != Ordering::Greater || s.cmp(&Rational::one()) != Ordering::Less {
        return None;
    }
    Some(a.add(&b.sub(a).mul(&s)))
}

/// Pointwise `max` (or `min`) of two walls, exact (crossings inserted).
pub(crate) fn combine(f: &Pl, g: &Pl, take_max: bool) -> Pl {
    let mut ts = merge_ts(f, g);
    let mut extra = Vec::new();
    for w in ts.windows(2) {
        if let Some(x) = cross(f, g, &w[0], &w[1]) {
            extra.push(x);
        }
    }
    ts.extend(extra);
    sort_dedup(&mut ts);
    let pts = ts
        .iter()
        .map(|t| {
            let (fv, gv) = (f.value_at(t), g.value_at(t));
            let v = if (fv.cmp(&gv) == Ordering::Greater) == take_max {
                fv
            } else {
                gv
            };
            (t.clone(), v)
        })
        .collect();
    Pl { pts }.cleaned()
}

/// Pointwise difference `f - g` (linear on each merged segment, no crossings).
pub(crate) fn sub_pl(f: &Pl, g: &Pl) -> Pl {
    let pts = merge_ts(f, g)
        .into_iter()
        .map(|t| {
            let v = f.value_at(&t).sub(&g.value_at(&t));
            (t, v)
        })
        .collect();
    Pl { pts }.cleaned()
}

/// Pointwise sum `f + g` (linear on each merged segment, no crossings) — the
/// additive twin of [`sub_pl`]. It names the tropical `⊗` on walls (tropical
/// multiplication is ordinary addition of values); see
/// [`crate::games::tropical_thermography`].
pub(crate) fn add_pl(f: &Pl, g: &Pl) -> Pl {
    let pts = merge_ts(f, g)
        .into_iter()
        .map(|t| {
            let v = f.value_at(&t).add(&g.value_at(&t));
            (t, v)
        })
        .collect();
    Pl { pts }.cleaned()
}
