use super::ast::{BinaryOp, Expr, RelOp, StarLiteral, Statement, UnaryOp};
use super::error::{OghamError, OghamErrorKind, OghamResult, Span};
use super::parse::parse_statement;
use super::unparse::unparse_statement;
use crate::clifford::{CliffordAlgebra, Metric, Multivector};
use crate::scalar::{
    checked_factorial_i128, factorial_in_scalar, nim_trace, ExactFieldScalar, FiniteField, Fp, Fpn,
    Integer, IntegerDivExactError, Nimber, Omnific, Ordinal, Poly, Rational, RationalFunction,
    Scalar, Surreal,
};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Display;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EvalLine {
    pub canonical: String,
    pub value: Option<String>,
}

pub fn eval_to_string(world: &str, src: &str) -> OghamResult<String> {
    let mut session = OghamSession::new(world)?;
    let mut out = Vec::new();
    for line in src.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix(":world ") {
            session.set_world(rest)?;
            continue;
        }
        if let Some(value) = session.eval_line(trimmed)?.value {
            out.push(value);
        }
    }
    Ok(out.join("\n"))
}

pub struct OghamSession {
    world: World,
}

impl OghamSession {
    pub fn new(world_decl: &str) -> OghamResult<Self> {
        Ok(OghamSession {
            world: World::from_decl(world_decl)?,
        })
    }

    pub fn set_world(&mut self, world_decl: &str) -> OghamResult<()> {
        self.world = World::from_decl(world_decl)?;
        Ok(())
    }

    pub fn eval_line(&mut self, src: &str) -> OghamResult<EvalLine> {
        let stmt = parse_statement(src)?;
        let canonical = unparse_statement(&stmt);
        let value = self.world.eval_statement(&stmt)?;
        Ok(EvalLine { canonical, value })
    }

    pub fn world_summary(&self) -> String {
        self.world.summary()
    }

    pub fn env_summary(&self) -> Vec<String> {
        self.world.env_summary()
    }
}

enum World {
    Nimber(Runtime<Nimber>),
    Ordinal(Runtime<Ordinal>),
    Surreal(Runtime<Surreal>),
    Omnific(Runtime<Omnific>),
    Integer(Runtime<Integer>),
    Fp2(Runtime<Fp<2>>),
    Fp3(Runtime<Fp<3>>),
    Fp5(Runtime<Fp<5>>),
    Fp7(Runtime<Fp<7>>),
    F4(Runtime<Fpn<2, 2>>),
    F8(Runtime<Fpn<2, 3>>),
    F16(Runtime<Fpn<2, 4>>),
    F9(Runtime<Fpn<3, 2>>),
    F27(Runtime<Fpn<3, 3>>),
    F25(Runtime<Fpn<5, 2>>),
    PolyInt(PolyRuntime<Integer>),
    Poly2(PolyRuntime<Fp<2>>),
    Poly3(PolyRuntime<Fp<3>>),
    Poly5(PolyRuntime<Fp<5>>),
    Poly7(PolyRuntime<Fp<7>>),
    RatFunc2(RatFuncRuntime<Fp<2>>),
    RatFunc3(RatFuncRuntime<Fp<3>>),
    RatFunc5(RatFuncRuntime<Fp<5>>),
    RatFunc7(RatFuncRuntime<Fp<7>>),
}

impl World {
    fn from_decl(decl: &str) -> OghamResult<Self> {
        let decl = decl.trim().strip_prefix(":world ").unwrap_or(decl.trim());
        let mut parts = decl.split_whitespace();
        let name = parts
            .next()
            .ok_or_else(|| parse_error("missing world name"))?;
        let tail: Vec<&str> = parts.collect();
        macro_rules! build_poly {
            ($variant:ident, $ty:ty, $label:expr) => {{
                ensure_function_world_decl(name, &tail)?;
                return Ok(World::$variant(PolyRuntime::<$ty>::new($label)));
            }};
        }
        macro_rules! build_ratfunc {
            ($variant:ident, $ty:ty, $label:expr) => {{
                ensure_function_world_decl(name, &tail)?;
                return Ok(World::$variant(RatFuncRuntime::<$ty>::new($label)));
            }};
        }
        match name {
            "polyint" => build_poly!(PolyInt, Integer, "polyint"),
            "poly2" => build_poly!(Poly2, Fp<2>, "poly2"),
            "poly3" => build_poly!(Poly3, Fp<3>, "poly3"),
            "poly5" => build_poly!(Poly5, Fp<5>, "poly5"),
            "poly7" => build_poly!(Poly7, Fp<7>, "poly7"),
            "ratfunc2" => build_ratfunc!(RatFunc2, Fp<2>, "ratfunc2"),
            "ratfunc3" => build_ratfunc!(RatFunc3, Fp<3>, "ratfunc3"),
            "ratfunc5" => build_ratfunc!(RatFunc5, Fp<5>, "ratfunc5"),
            "ratfunc7" => build_ratfunc!(RatFunc7, Fp<7>, "ratfunc7"),
            _ => {}
        }
        let second = tail
            .first()
            .copied()
            .ok_or_else(|| parse_error("missing world dimension or constructor"))?;
        if name == "nimber" && second.starts_with("gold(") {
            let metric = parse_gold_metric(second)?;
            return Ok(World::Nimber(Runtime::from_metric("nimber", metric)));
        }
        let dim = second
            .parse::<usize>()
            .map_err(|_| parse_error("world dimension must be a usize"))?;
        let rest = decl.split_once(second).map_or("", |(_, tail)| tail).trim();
        macro_rules! build {
            ($variant:ident, $ty:ty, $label:expr) => {
                Ok(World::$variant(build_runtime::<$ty>($label, dim, rest)?))
            };
        }
        match name {
            "nimber" => build!(Nimber, Nimber, "nimber"),
            "ordinal" => build!(Ordinal, Ordinal, "ordinal"),
            "surreal" => build!(Surreal, Surreal, "surreal"),
            "omnific" => build!(Omnific, Omnific, "omnific"),
            "integer" => build!(Integer, Integer, "integer"),
            "fp2" => build!(Fp2, Fp<2>, "fp2"),
            "fp3" => build!(Fp3, Fp<3>, "fp3"),
            "fp5" => build!(Fp5, Fp<5>, "fp5"),
            "fp7" => build!(Fp7, Fp<7>, "fp7"),
            "f4" => build!(F4, Fpn<2, 2>, "f4"),
            "f8" => build!(F8, Fpn<2, 3>, "f8"),
            "f16" => build!(F16, Fpn<2, 4>, "f16"),
            "f9" => build!(F9, Fpn<3, 2>, "f9"),
            "f27" => build!(F27, Fpn<3, 3>, "f27"),
            "f25" => build!(F25, Fpn<5, 2>, "f25"),
            _ => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                format!("unknown world `{name}`"),
            )),
        }
    }

    fn eval_statement(&mut self, stmt: &Statement) -> OghamResult<Option<String>> {
        macro_rules! dispatch {
            ($rt:expr) => {
                $rt.eval_statement(stmt)
            };
        }
        match self {
            World::Nimber(rt) => dispatch!(rt),
            World::Ordinal(rt) => dispatch!(rt),
            World::Surreal(rt) => dispatch!(rt),
            World::Omnific(rt) => dispatch!(rt),
            World::Integer(rt) => dispatch!(rt),
            World::Fp2(rt) => dispatch!(rt),
            World::Fp3(rt) => dispatch!(rt),
            World::Fp5(rt) => dispatch!(rt),
            World::Fp7(rt) => dispatch!(rt),
            World::F4(rt) => dispatch!(rt),
            World::F8(rt) => dispatch!(rt),
            World::F16(rt) => dispatch!(rt),
            World::F9(rt) => dispatch!(rt),
            World::F27(rt) => dispatch!(rt),
            World::F25(rt) => dispatch!(rt),
            World::PolyInt(rt) => dispatch!(rt),
            World::Poly2(rt) => dispatch!(rt),
            World::Poly3(rt) => dispatch!(rt),
            World::Poly5(rt) => dispatch!(rt),
            World::Poly7(rt) => dispatch!(rt),
            World::RatFunc2(rt) => dispatch!(rt),
            World::RatFunc3(rt) => dispatch!(rt),
            World::RatFunc5(rt) => dispatch!(rt),
            World::RatFunc7(rt) => dispatch!(rt),
        }
    }

    fn summary(&self) -> String {
        macro_rules! dispatch {
            ($rt:expr) => {
                $rt.summary()
            };
        }
        match self {
            World::Nimber(rt) => dispatch!(rt),
            World::Ordinal(rt) => dispatch!(rt),
            World::Surreal(rt) => dispatch!(rt),
            World::Omnific(rt) => dispatch!(rt),
            World::Integer(rt) => dispatch!(rt),
            World::Fp2(rt) => dispatch!(rt),
            World::Fp3(rt) => dispatch!(rt),
            World::Fp5(rt) => dispatch!(rt),
            World::Fp7(rt) => dispatch!(rt),
            World::F4(rt) => dispatch!(rt),
            World::F8(rt) => dispatch!(rt),
            World::F16(rt) => dispatch!(rt),
            World::F9(rt) => dispatch!(rt),
            World::F27(rt) => dispatch!(rt),
            World::F25(rt) => dispatch!(rt),
            World::PolyInt(rt) => dispatch!(rt),
            World::Poly2(rt) => dispatch!(rt),
            World::Poly3(rt) => dispatch!(rt),
            World::Poly5(rt) => dispatch!(rt),
            World::Poly7(rt) => dispatch!(rt),
            World::RatFunc2(rt) => dispatch!(rt),
            World::RatFunc3(rt) => dispatch!(rt),
            World::RatFunc5(rt) => dispatch!(rt),
            World::RatFunc7(rt) => dispatch!(rt),
        }
    }

    fn env_summary(&self) -> Vec<String> {
        macro_rules! dispatch {
            ($rt:expr) => {
                $rt.env_summary()
            };
        }
        match self {
            World::Nimber(rt) => dispatch!(rt),
            World::Ordinal(rt) => dispatch!(rt),
            World::Surreal(rt) => dispatch!(rt),
            World::Omnific(rt) => dispatch!(rt),
            World::Integer(rt) => dispatch!(rt),
            World::Fp2(rt) => dispatch!(rt),
            World::Fp3(rt) => dispatch!(rt),
            World::Fp5(rt) => dispatch!(rt),
            World::Fp7(rt) => dispatch!(rt),
            World::F4(rt) => dispatch!(rt),
            World::F8(rt) => dispatch!(rt),
            World::F16(rt) => dispatch!(rt),
            World::F9(rt) => dispatch!(rt),
            World::F27(rt) => dispatch!(rt),
            World::F25(rt) => dispatch!(rt),
            World::PolyInt(rt) => dispatch!(rt),
            World::Poly2(rt) => dispatch!(rt),
            World::Poly3(rt) => dispatch!(rt),
            World::Poly5(rt) => dispatch!(rt),
            World::Poly7(rt) => dispatch!(rt),
            World::RatFunc2(rt) => dispatch!(rt),
            World::RatFunc3(rt) => dispatch!(rt),
            World::RatFunc5(rt) => dispatch!(rt),
            World::RatFunc7(rt) => dispatch!(rt),
        }
    }
}

fn ensure_function_world_decl(name: &str, tail: &[&str]) -> OghamResult<()> {
    if tail.is_empty() || tail == ["0"] {
        Ok(())
    } else {
        Err(parse_error(format!(
            "`{name}` is a function-shaped scalar world; it takes no metric declaration"
        )))
    }
}

struct PolyRuntime<S: PolyWorldCoeff> {
    name: &'static str,
    env: BTreeMap<String, Poly<S>>,
}

impl<S: PolyWorldCoeff> PolyRuntime<S> {
    fn new(name: &'static str) -> Self {
        PolyRuntime {
            name,
            env: BTreeMap::new(),
        }
    }

    fn eval_statement(&mut self, stmt: &Statement) -> OghamResult<Option<String>> {
        match stmt {
            Statement::Binding { name, expr } => {
                if name == "t" {
                    return Err(OghamError::new(
                        OghamErrorKind::Reserved,
                        Span::point(0),
                        format!("`t` is reserved in the `{}` world", self.name),
                    ));
                }
                let value = self.eval_element(expr)?;
                self.env.insert(name.clone(), value);
                Ok(None)
            }
            Statement::Expr(expr) => match expr {
                Expr::Relation { op, lhs, rhs } => {
                    let value = self.eval_relation(*op, lhs, rhs)?;
                    Ok(Some(value.to_string()))
                }
                _ if expression_is_index(expr) => Ok(Some(self.eval_index(expr)?.to_string())),
                _ => Ok(Some(self.eval_element(expr)?.to_string())),
            },
        }
    }

    fn summary(&self) -> String {
        self.name.to_string()
    }

    fn env_summary(&self) -> Vec<String> {
        self.env
            .iter()
            .map(|(name, value)| format!("{name} := {value}"))
            .collect()
    }

    fn eval_relation(&mut self, op: RelOp, lhs: &Expr, rhs: &Expr) -> OghamResult<bool> {
        if expression_is_index(lhs) || expression_is_index(rhs) {
            let lhs = self.eval_index(lhs)?;
            let rhs = self.eval_index(rhs)?;
            return ordered_relation(op, lhs.cmp(&rhs));
        }
        let lhs = self.eval_element(lhs)?;
        let rhs = self.eval_element(rhs)?;
        if op == RelOp::Eq {
            Ok(lhs == rhs)
        } else {
            Err(no_order_error())
        }
    }

    fn eval_element(&mut self, expr: &Expr) -> OghamResult<Poly<S>> {
        match expr {
            Expr::Int(n) => Ok(Poly::constant(S::bare_int(*n, Span::point(0))?)),
            Expr::Star(star) => Ok(Poly::constant(S::star(star, Span::point(0))?)),
            Expr::Omega => Ok(Poly::constant(S::omega(Span::point(0))?)),
            Expr::Blade(_) | Expr::Vector(_) => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "function-shaped worlds do not have Clifford blades or vectors",
            )),
            Expr::Ident(name) => {
                if name == "t" {
                    Ok(Poly::x())
                } else if let Some(value) = self.env.get(name) {
                    Ok(value.clone())
                } else {
                    Err(OghamError::new(
                        OghamErrorKind::Unbound,
                        Span::point(0),
                        format!("unbound identifier `{name}`"),
                    )
                    .with_hint(format!("did you mean `{name} := ...`?")))
                }
            }
            Expr::Call { name, args } => self.eval_call(name, args),
            Expr::Factorial(expr) => {
                let n = self.eval_index(expr)?;
                Ok(Poly::constant(S::factorial(n, Span::point(0))?))
            }
            Expr::Unary { op, expr } => {
                let value = self.eval_element(expr)?;
                match op {
                    UnaryOp::Neg => Ok(value.neg()),
                    UnaryOp::Inv => self.inverse_element(&value),
                }
            }
            Expr::Binary { op, lhs, rhs } => self.eval_binary(*op, lhs, rhs),
            Expr::Relation { .. } => Err(OghamError::new(
                OghamErrorKind::Parse,
                Span::point(0),
                "relations only appear as top-level statements",
            )),
        }
    }

    fn eval_binary(&mut self, op: BinaryOp, lhs: &Expr, rhs: &Expr) -> OghamResult<Poly<S>> {
        if op == BinaryOp::Pow {
            return self.eval_power(lhs, rhs);
        }
        let lhs_v = self.eval_element(lhs)?;
        let rhs_v = self.eval_element(rhs)?;
        match op {
            BinaryOp::Add => Ok(lhs_v.add(&rhs_v)),
            BinaryOp::Sub => Ok(lhs_v.sub(&rhs_v)),
            BinaryOp::Mul => Ok(lhs_v.mul(&rhs_v)),
            BinaryOp::Div => poly_exact_div::<S>(&lhs_v, &rhs_v, Span::point(0)),
            BinaryOp::Rem => poly_rem::<S>(&lhs_v, &rhs_v, Span::point(0)),
            BinaryOp::At => Ok(lhs_v.compose(&rhs_v)),
            BinaryOp::Wedge => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "wedge product belongs to Clifford worlds",
            )),
            BinaryOp::Pow => unreachable!(),
        }
    }

    fn eval_power(&mut self, lhs: &Expr, rhs: &Expr) -> OghamResult<Poly<S>> {
        let base = self.eval_element(lhs)?;
        let exp = self.eval_index(rhs).map_err(|err| {
            if err.kind == OghamErrorKind::IndexSort {
                exp_sort_error()
            } else {
                err
            }
        })?;
        if exp < 0 {
            let inv = self.inverse_element(&base)?;
            let k = exp
                .checked_neg()
                .and_then(|v| u128::try_from(v).ok())
                .ok_or_else(|| overflow("negative exponent magnitude exceeds u128"))?;
            Ok(pow_poly(&inv, k))
        } else {
            let k = u128::try_from(exp).map_err(|_| overflow("exponent exceeds u128"))?;
            Ok(pow_poly(&base, k))
        }
    }

    fn eval_call(&mut self, name: &str, args: &[Expr]) -> OghamResult<Poly<S>> {
        match name {
            "deg" => Err(index_sort_error().with_hint("`deg` returns an Index")),
            "gcd" => {
                expect_arity(name, args, 2)?;
                let lhs = self.eval_element(&args[0])?;
                let rhs = self.eval_element(&args[1])?;
                S::gcd_poly(&lhs, &rhs, Span::point(0))
            }
            _ => Err(OghamError::new(
                OghamErrorKind::UnknownFn,
                Span::point(0),
                format!("unknown function `{name}`"),
            )),
        }
    }

    fn eval_index(&mut self, expr: &Expr) -> OghamResult<i128> {
        match expr {
            Expr::Int(n) => u128_to_i128(*n),
            Expr::Call { name, args } if name == "deg" => {
                expect_arity(name, args, 1)?;
                let value = self.eval_element(&args[0])?;
                let degree = value.degree().ok_or_else(|| {
                    OghamError::new(
                        OghamErrorKind::Domain,
                        Span::point(0),
                        "degree of the zero polynomial is undefined",
                    )
                })?;
                i128::try_from(degree).map_err(|_| overflow("polynomial degree exceeds i128"))
            }
            Expr::Unary {
                op: UnaryOp::Neg,
                expr,
            } => self
                .eval_index(expr)?
                .checked_neg()
                .ok_or_else(|| overflow("index negation overflowed i128")),
            Expr::Binary { op, lhs, rhs } => {
                let lhs = self.eval_index(lhs)?;
                let rhs = self.eval_index(rhs)?;
                eval_index_binary(*op, lhs, rhs)
            }
            _ => Err(index_sort_error()),
        }
    }

    fn inverse_element(&self, value: &Poly<S>) -> OghamResult<Poly<S>> {
        if value.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                Span::point(0),
                "division by zero",
            ));
        }
        value.inv().ok_or_else(|| {
            OghamError::new(
                OghamErrorKind::NotInvertible,
                Span::point(0),
                "polynomial is not a unit",
            )
        })
    }
}

struct RatFuncRuntime<S: OghamScalar + ExactFieldScalar> {
    name: &'static str,
    env: BTreeMap<String, RationalFunction<S>>,
}

impl<S: OghamScalar + ExactFieldScalar> RatFuncRuntime<S> {
    fn new(name: &'static str) -> Self {
        RatFuncRuntime {
            name,
            env: BTreeMap::new(),
        }
    }

    fn eval_statement(&mut self, stmt: &Statement) -> OghamResult<Option<String>> {
        match stmt {
            Statement::Binding { name, expr } => {
                if name == "t" {
                    return Err(OghamError::new(
                        OghamErrorKind::Reserved,
                        Span::point(0),
                        format!("`t` is reserved in the `{}` world", self.name),
                    ));
                }
                let value = self.eval_element(expr)?;
                self.env.insert(name.clone(), value);
                Ok(None)
            }
            Statement::Expr(expr) => match expr {
                Expr::Relation { op, lhs, rhs } => {
                    let value = self.eval_relation(*op, lhs, rhs)?;
                    Ok(Some(value.to_string()))
                }
                _ if expression_is_index(expr) => Ok(Some(self.eval_index(expr)?.to_string())),
                _ => Ok(Some(self.eval_element(expr)?.to_string())),
            },
        }
    }

    fn summary(&self) -> String {
        self.name.to_string()
    }

    fn env_summary(&self) -> Vec<String> {
        self.env
            .iter()
            .map(|(name, value)| format!("{name} := {value}"))
            .collect()
    }

    fn eval_relation(&mut self, op: RelOp, lhs: &Expr, rhs: &Expr) -> OghamResult<bool> {
        let lhs = self.eval_element(lhs)?;
        let rhs = self.eval_element(rhs)?;
        if op == RelOp::Eq {
            Ok(lhs == rhs)
        } else {
            Err(no_order_error())
        }
    }

    fn eval_element(&mut self, expr: &Expr) -> OghamResult<RationalFunction<S>> {
        match expr {
            Expr::Int(n) => Ok(RationalFunction::from_base(S::bare_int(
                *n,
                Span::point(0),
            )?)),
            Expr::Star(star) => Ok(RationalFunction::from_base(S::star(star, Span::point(0))?)),
            Expr::Omega => Ok(RationalFunction::from_base(S::omega(Span::point(0))?)),
            Expr::Blade(_) | Expr::Vector(_) => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "function-shaped worlds do not have Clifford blades or vectors",
            )),
            Expr::Ident(name) => {
                if name == "t" {
                    Ok(RationalFunction::t())
                } else if let Some(value) = self.env.get(name) {
                    Ok(value.clone())
                } else {
                    Err(OghamError::new(
                        OghamErrorKind::Unbound,
                        Span::point(0),
                        format!("unbound identifier `{name}`"),
                    )
                    .with_hint(format!("did you mean `{name} := ...`?")))
                }
            }
            Expr::Call { name, args } => self.eval_call(name, args),
            Expr::Factorial(expr) => {
                let n = self.eval_index(expr)?;
                Ok(RationalFunction::from_base(S::factorial(
                    n,
                    Span::point(0),
                )?))
            }
            Expr::Unary { op, expr } => {
                let value = self.eval_element(expr)?;
                match op {
                    UnaryOp::Neg => Ok(value.neg()),
                    UnaryOp::Inv => self.inverse_element(&value),
                }
            }
            Expr::Binary { op, lhs, rhs } => self.eval_binary(*op, lhs, rhs),
            Expr::Relation { .. } => Err(OghamError::new(
                OghamErrorKind::Parse,
                Span::point(0),
                "relations only appear as top-level statements",
            )),
        }
    }

    fn eval_binary(
        &mut self,
        op: BinaryOp,
        lhs: &Expr,
        rhs: &Expr,
    ) -> OghamResult<RationalFunction<S>> {
        if op == BinaryOp::Pow {
            return self.eval_power(lhs, rhs);
        }
        let lhs_v = self.eval_element(lhs)?;
        let rhs_v = self.eval_element(rhs)?;
        match op {
            BinaryOp::Add => Ok(lhs_v.add(&rhs_v)),
            BinaryOp::Sub => Ok(lhs_v.sub(&rhs_v)),
            BinaryOp::Mul => Ok(lhs_v.mul(&rhs_v)),
            BinaryOp::Div => {
                if rhs_v.is_zero() {
                    Err(OghamError::new(
                        OghamErrorKind::DivisionByZero,
                        Span::point(0),
                        "division by zero",
                    ))
                } else {
                    Ok(lhs_v.mul(&rhs_v.inv().expect("checked nonzero rational function")))
                }
            }
            BinaryOp::Rem => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "function-field worlds are fields; `%` is only active in polynomial worlds",
            )),
            BinaryOp::At => substitute_rational_function(&lhs_v, &rhs_v, Span::point(0)),
            BinaryOp::Wedge => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "wedge product belongs to Clifford worlds",
            )),
            BinaryOp::Pow => unreachable!(),
        }
    }

    fn eval_power(&mut self, lhs: &Expr, rhs: &Expr) -> OghamResult<RationalFunction<S>> {
        let base = self.eval_element(lhs)?;
        let exp = self.eval_index(rhs).map_err(|err| {
            if err.kind == OghamErrorKind::IndexSort {
                exp_sort_error()
            } else {
                err
            }
        })?;
        if exp < 0 {
            let inv = self.inverse_element(&base)?;
            let k = exp
                .checked_neg()
                .and_then(|v| u128::try_from(v).ok())
                .ok_or_else(|| overflow("negative exponent magnitude exceeds u128"))?;
            Ok(pow_rational_function(&inv, k))
        } else {
            let k = u128::try_from(exp).map_err(|_| overflow("exponent exceeds u128"))?;
            Ok(pow_rational_function(&base, k))
        }
    }

    fn eval_call(&mut self, name: &str, _args: &[Expr]) -> OghamResult<RationalFunction<S>> {
        match name {
            "deg" | "gcd" => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                format!("`{name}` is a polynomial-world function, not a ratfunc function"),
            )),
            _ => Err(OghamError::new(
                OghamErrorKind::UnknownFn,
                Span::point(0),
                format!("unknown function `{name}`"),
            )),
        }
    }

    fn eval_index(&mut self, expr: &Expr) -> OghamResult<i128> {
        match expr {
            Expr::Int(n) => u128_to_i128(*n),
            Expr::Call { name, .. } if name == "deg" => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "`deg` is a polynomial-world function, not a ratfunc function",
            )),
            Expr::Unary {
                op: UnaryOp::Neg,
                expr,
            } => self
                .eval_index(expr)?
                .checked_neg()
                .ok_or_else(|| overflow("index negation overflowed i128")),
            Expr::Binary { op, lhs, rhs } => {
                let lhs = self.eval_index(lhs)?;
                let rhs = self.eval_index(rhs)?;
                eval_index_binary(*op, lhs, rhs)
            }
            _ => Err(index_sort_error()),
        }
    }

    fn inverse_element(&self, value: &RationalFunction<S>) -> OghamResult<RationalFunction<S>> {
        if value.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                Span::point(0),
                "division by zero",
            ));
        }
        Ok(value.inv().expect("checked nonzero rational function"))
    }
}

struct Runtime<S: OghamScalar> {
    name: &'static str,
    alg: CliffordAlgebra<S>,
    env: BTreeMap<String, Multivector<S>>,
}

impl<S: OghamScalar> Runtime<S> {
    fn from_metric(name: &'static str, metric: Metric<S>) -> Self {
        Runtime {
            name,
            alg: CliffordAlgebra::new(metric.dim(), metric),
            env: BTreeMap::new(),
        }
    }

    fn eval_statement(&mut self, stmt: &Statement) -> OghamResult<Option<String>> {
        match stmt {
            Statement::Binding { name, expr } => {
                if S::reserved_ident(name) {
                    return Err(OghamError::new(
                        OghamErrorKind::Reserved,
                        Span::point(0),
                        format!("`{name}` is reserved in the `{}` world", self.name),
                    ));
                }
                let value = self.eval_expr(expr)?;
                self.env.insert(name.clone(), value);
                Ok(None)
            }
            Statement::Expr(expr) => match expr {
                Expr::Relation { op, lhs, rhs } => {
                    let value = self.eval_relation(*op, lhs, rhs)?;
                    Ok(Some(value.to_string()))
                }
                _ => Ok(Some(self.eval_expr(expr)?.to_string())),
            },
        }
    }

    fn summary(&self) -> String {
        format!("{} dim {}", self.name, self.alg.dim())
    }

    fn env_summary(&self) -> Vec<String> {
        self.env
            .iter()
            .map(|(name, value)| format!("{name} := {value}"))
            .collect()
    }

    fn eval_relation(&mut self, op: RelOp, lhs: &Expr, rhs: &Expr) -> OghamResult<bool> {
        let lhs = self.eval_expr(lhs)?;
        let rhs = self.eval_expr(rhs)?;
        if op == RelOp::Eq {
            return Ok(lhs == rhs);
        }
        let Some(lhs) = scalar_part(&lhs) else {
            return Err(grade0_error(Span::point(0)));
        };
        let Some(rhs) = scalar_part(&rhs) else {
            return Err(grade0_error(Span::point(0)));
        };
        S::relation(op, &lhs, &rhs, Span::point(0))
    }

    fn eval_expr(&mut self, expr: &Expr) -> OghamResult<Multivector<S>> {
        match expr {
            Expr::Int(n) => Ok(self.alg.scalar(S::bare_int(*n, Span::point(0))?)),
            Expr::Star(star) => Ok(self.alg.scalar(S::star(star, Span::point(0))?)),
            Expr::Omega => Ok(self.alg.scalar(S::omega(Span::point(0))?)),
            Expr::Blade(i) => {
                if *i >= self.alg.dim() {
                    Err(OghamError::new(
                        OghamErrorKind::BladeIndex,
                        Span::point(0),
                        format!("blade e{i} is outside dimension {}", self.alg.dim()),
                    ))
                } else {
                    Ok(self.alg.e(*i))
                }
            }
            Expr::Vector(items) => self.eval_vector(items),
            Expr::Ident(name) => {
                if let Some(value) = self.env.get(name) {
                    Ok(value.clone())
                } else if let Some(x) = S::named_element(name, Span::point(0))? {
                    Ok(self.alg.scalar(x))
                } else {
                    Err(OghamError::new(
                        OghamErrorKind::Unbound,
                        Span::point(0),
                        format!("unbound identifier `{name}`"),
                    )
                    .with_hint(format!("did you mean `{name} := ...`?")))
                }
            }
            Expr::Call { name, args } => self.eval_call(name, args),
            Expr::Factorial(expr) => {
                let n = self.eval_index(expr)?;
                Ok(self.alg.scalar(S::factorial(n, Span::point(0))?))
            }
            Expr::Unary { op, expr } => {
                let value = self.eval_expr(expr)?;
                match op {
                    UnaryOp::Neg => Ok(-value),
                    UnaryOp::Inv => self.inverse_mv(&value),
                }
            }
            Expr::Binary { op, lhs, rhs } => self.eval_binary(*op, lhs, rhs),
            Expr::Relation { .. } => Err(OghamError::new(
                OghamErrorKind::Parse,
                Span::point(0),
                "relations only appear as top-level statements",
            )),
        }
    }

    fn eval_binary(&mut self, op: BinaryOp, lhs: &Expr, rhs: &Expr) -> OghamResult<Multivector<S>> {
        if op == BinaryOp::Pow {
            return self.eval_power(lhs, rhs);
        }
        if op == BinaryOp::At {
            return Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                Span::point(0),
                "evaluation lives in the function-shaped worlds — poly/ratfunc, v1.1",
            ));
        }
        let lhs_v = self.eval_expr(lhs)?;
        let rhs_v = self.eval_expr(rhs)?;
        match op {
            BinaryOp::Add => Ok(lhs_v + rhs_v),
            BinaryOp::Sub => Ok(lhs_v - rhs_v),
            BinaryOp::Mul => self.mul_mv(&lhs_v, &rhs_v),
            BinaryOp::Div => self.div_mv(&lhs_v, &rhs_v),
            BinaryOp::Rem => {
                let Some(lhs_s) = scalar_part(&lhs_v) else {
                    return Err(grade0_error(Span::point(0)));
                };
                let Some(rhs_s) = scalar_part(&rhs_v) else {
                    return Err(grade0_error(Span::point(0)));
                };
                Ok(self.alg.scalar(S::rem(&lhs_s, &rhs_s, Span::point(0))?))
            }
            BinaryOp::Wedge => Ok(self.alg.wedge(&lhs_v, &rhs_v)),
            BinaryOp::Pow | BinaryOp::At => unreachable!(),
        }
    }

    fn eval_power(&mut self, lhs: &Expr, rhs: &Expr) -> OghamResult<Multivector<S>> {
        if lhs.is_omega_atom() {
            if let Err(index_err) = self.eval_index(rhs) {
                if index_err.kind == OghamErrorKind::IndexSort {
                    let exp = self.eval_expr(rhs)?;
                    let Some(exp) = scalar_part(&exp) else {
                        return Err(exp_sort_error());
                    };
                    return Ok(self.alg.scalar(S::omega_pow(exp, Span::point(0))?));
                }
                return Err(index_err);
            }
        }
        let base = self.eval_expr(lhs)?;
        let exp = self.eval_index(rhs).map_err(|err| {
            if err.kind == OghamErrorKind::IndexSort {
                exp_sort_error()
            } else {
                err
            }
        })?;
        if exp < 0 {
            let inv = self.inverse_mv(&base)?;
            let k = exp
                .checked_neg()
                .and_then(|v| u128::try_from(v).ok())
                .ok_or_else(|| overflow("negative exponent magnitude exceeds u128"))?;
            self.pow_mv(&inv, k)
        } else {
            let k = u128::try_from(exp).map_err(|_| overflow("exponent exceeds u128"))?;
            self.pow_mv(&base, k)
        }
    }

    fn eval_vector(&mut self, items: &[Expr]) -> OghamResult<Multivector<S>> {
        if self.alg.dim() == 0 || items.len() != self.alg.dim() {
            return Err(OghamError::new(
                OghamErrorKind::DimMismatch,
                Span::point(0),
                format!(
                    "vector length {} does not match world dimension {}",
                    items.len(),
                    self.alg.dim()
                ),
            ));
        }
        let mut out = self.alg.zero();
        for (i, expr) in items.iter().enumerate() {
            let value = self.eval_expr(expr)?;
            let Some(coeff) = scalar_part(&value) else {
                return Err(grade0_error(Span::point(0)));
            };
            out = self
                .alg
                .add(&out, &self.alg.scalar_mul(&coeff, &self.alg.e(i)));
        }
        Ok(out)
    }

    fn eval_call(&mut self, name: &str, args: &[Expr]) -> OghamResult<Multivector<S>> {
        match name {
            "rev" => {
                expect_arity(name, args, 1)?;
                if self.alg.metric().has_upper() {
                    return Err(OghamError::new(
                        OghamErrorKind::GeneralMetric,
                        Span::point(0),
                        "reverse is undefined for the Chevalley construction",
                    ));
                }
                let x = self.eval_expr(&args[0])?;
                Ok(self.alg.reverse(&x))
            }
            "grade" => {
                expect_arity(name, args, 2)?;
                let x = self.eval_expr(&args[0])?;
                let k = self.eval_index(&args[1])?;
                if k < 0 {
                    return Err(OghamError::new(
                        OghamErrorKind::Domain,
                        Span::point(0),
                        "grade index must be non-negative",
                    ));
                }
                Ok(self.alg.grade_part(&x, k as usize))
            }
            "even" => {
                expect_arity(name, args, 1)?;
                let x = self.eval_expr(&args[0])?;
                Ok(self.alg.even_part(&x))
            }
            "dual" => {
                expect_arity(name, args, 1)?;
                if self.alg.metric().has_upper() {
                    return Err(OghamError::new(
                        OghamErrorKind::GeneralMetric,
                        Span::point(0),
                        "dual is undefined for general-bilinear metrics",
                    ));
                }
                let x = self.eval_expr(&args[0])?;
                self.alg.dual(&x).ok_or_else(|| {
                    OghamError::new(
                        OghamErrorKind::NotInvertible,
                        Span::point(0),
                        "pseudoscalar is not invertible",
                    )
                })
            }
            "frob" => {
                expect_arity(name, args, 1)?;
                let x = self.eval_grade0(&args[0])?;
                Ok(self.alg.scalar(S::frob(&x, Span::point(0))?))
            }
            "tr" => {
                if args.is_empty() || args.len() > 2 {
                    return Err(OghamError::new(
                        OghamErrorKind::Arity,
                        Span::point(0),
                        "`tr` expects one or two arguments",
                    ));
                }
                let x = self.eval_grade0(&args[0])?;
                let m = if args.len() == 2 {
                    Some(self.eval_index(&args[1])?)
                } else {
                    None
                };
                Ok(self.alg.scalar(S::trace(&x, m, Span::point(0))?))
            }
            _ => Err(OghamError::new(
                OghamErrorKind::UnknownFn,
                Span::point(0),
                format!("unknown function `{name}`"),
            )),
        }
    }

    fn eval_grade0(&mut self, expr: &Expr) -> OghamResult<S> {
        let value = self.eval_expr(expr)?;
        scalar_part(&value).ok_or_else(|| grade0_error(Span::point(0)))
    }

    fn eval_index(&mut self, expr: &Expr) -> OghamResult<i128> {
        match expr {
            Expr::Int(n) => u128_to_i128(*n),
            Expr::Unary {
                op: UnaryOp::Neg,
                expr,
            } => self
                .eval_index(expr)?
                .checked_neg()
                .ok_or_else(|| overflow("index negation overflowed i128")),
            Expr::Binary { op, lhs, rhs } => {
                let lhs = self.eval_index(lhs)?;
                let rhs = self.eval_index(rhs)?;
                match op {
                    BinaryOp::Add => lhs
                        .checked_add(rhs)
                        .ok_or_else(|| overflow("index addition overflowed i128")),
                    BinaryOp::Sub => lhs
                        .checked_sub(rhs)
                        .ok_or_else(|| overflow("index subtraction overflowed i128")),
                    BinaryOp::Mul => lhs
                        .checked_mul(rhs)
                        .ok_or_else(|| overflow("index multiplication overflowed i128")),
                    BinaryOp::Pow => {
                        if rhs < 0 {
                            return Err(OghamError::new(
                                OghamErrorKind::Domain,
                                Span::point(0),
                                "index exponent must be non-negative",
                            ));
                        }
                        checked_i128_pow(lhs, rhs as u128)
                    }
                    _ => Err(index_sort_error()),
                }
            }
            _ => Err(index_sort_error()),
        }
    }

    fn inverse_mv(&self, value: &Multivector<S>) -> OghamResult<Multivector<S>> {
        if let Some(s) = scalar_part(value) {
            if s.is_zero() {
                return Err(OghamError::new(
                    OghamErrorKind::DivisionByZero,
                    Span::point(0),
                    "division by zero",
                ));
            }
            return Ok(self.alg.scalar(S::inv_scalar(&s, Span::point(0))?));
        }
        self.alg.multivector_inverse(value).ok_or_else(|| {
            OghamError::new(
                OghamErrorKind::NotInvertible,
                Span::point(0),
                "multivector is not invertible",
            )
        })
    }

    fn div_mv(&self, lhs: &Multivector<S>, rhs: &Multivector<S>) -> OghamResult<Multivector<S>> {
        if rhs.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                Span::point(0),
                "division by zero",
            ));
        }
        if let (Some(a), Some(b)) = (scalar_part(lhs), scalar_part(rhs)) {
            if let Some(out) = S::exact_div(&a, &b, Span::point(0)) {
                return Ok(self.alg.scalar(out?));
            }
        }
        let inv = self.inverse_mv(rhs)?;
        self.mul_mv(lhs, &inv)
    }

    fn mul_mv(&self, lhs: &Multivector<S>, rhs: &Multivector<S>) -> OghamResult<Multivector<S>> {
        if let (Some(a), Some(b)) = (scalar_part(lhs), scalar_part(rhs)) {
            return Ok(self.alg.scalar(S::mul_checked(&a, &b, Span::point(0))?));
        }
        S::mv_mul(&self.alg, lhs, rhs, Span::point(0))
    }

    fn pow_mv(&self, value: &Multivector<S>, k: u128) -> OghamResult<Multivector<S>> {
        if let Some(s) = scalar_part(value) {
            return Ok(self.alg.scalar(S::pow_checked(&s, k, Span::point(0))?));
        }
        S::mv_pow(&self.alg, value, k, Span::point(0))
    }
}

trait PolyWorldCoeff: OghamScalar {
    fn divrem_poly(
        lhs: &Poly<Self>,
        divisor: &Poly<Self>,
        span: Span,
    ) -> OghamResult<(Poly<Self>, Poly<Self>)>;
    fn gcd_poly(lhs: &Poly<Self>, rhs: &Poly<Self>, span: Span) -> OghamResult<Poly<Self>>;
}

impl<const P: u128> PolyWorldCoeff for Fp<P>
where
    Fp<P>: OghamScalar,
{
    fn divrem_poly(
        lhs: &Poly<Self>,
        divisor: &Poly<Self>,
        span: Span,
    ) -> OghamResult<(Poly<Self>, Poly<Self>)> {
        if divisor.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "polynomial division by zero",
            ));
        }
        Ok(lhs.divrem(divisor))
    }

    fn gcd_poly(lhs: &Poly<Self>, rhs: &Poly<Self>, _span: Span) -> OghamResult<Poly<Self>> {
        Ok(lhs.gcd(rhs))
    }
}

impl PolyWorldCoeff for Integer {
    fn divrem_poly(
        lhs: &Poly<Self>,
        divisor: &Poly<Self>,
        span: Span,
    ) -> OghamResult<(Poly<Self>, Poly<Self>)> {
        if divisor.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "polynomial division by zero",
            ));
        }
        if !matches!(divisor.leading(), Some(c) if *c == Integer::one()) {
            return Err(OghamError::new(
                OghamErrorKind::Modulus,
                span,
                "polyint divisors must be monic",
            ));
        }
        Ok(lhs.divrem(divisor))
    }

    fn gcd_poly(lhs: &Poly<Self>, rhs: &Poly<Self>, span: Span) -> OghamResult<Poly<Self>> {
        integer_poly_gcd(lhs, rhs, span)
    }
}

fn poly_rem<S: PolyWorldCoeff>(lhs: &Poly<S>, rhs: &Poly<S>, span: Span) -> OghamResult<Poly<S>> {
    let (_, r) = S::divrem_poly(lhs, rhs, span)?;
    Ok(r)
}

fn poly_exact_div<S: PolyWorldCoeff>(
    lhs: &Poly<S>,
    rhs: &Poly<S>,
    span: Span,
) -> OghamResult<Poly<S>> {
    let (q, r) = S::divrem_poly(lhs, rhs, span)?;
    if r.is_zero() {
        Ok(q)
    } else {
        Err(OghamError::new(
            OghamErrorKind::NotInvertible,
            span,
            format!("polynomial exact division failed with remainder {r}"),
        ))
    }
}

fn pow_poly<S: Scalar>(base: &Poly<S>, mut k: u128) -> Poly<S> {
    if k == 0 {
        return Poly::one();
    }
    let mut acc = Poly::one();
    let mut x = base.clone();
    loop {
        if k & 1 == 1 {
            acc = acc.mul(&x);
        }
        k >>= 1;
        if k == 0 {
            break;
        }
        x = x.mul(&x);
    }
    acc
}

fn pow_rational_function<S: ExactFieldScalar>(
    base: &RationalFunction<S>,
    mut k: u128,
) -> RationalFunction<S> {
    if k == 0 {
        return RationalFunction::one();
    }
    let mut acc = RationalFunction::one();
    let mut x = base.clone();
    loop {
        if k & 1 == 1 {
            acc = acc.mul(&x);
        }
        k >>= 1;
        if k == 0 {
            break;
        }
        x = x.mul(&x);
    }
    acc
}

fn substitute_rational_function<S: OghamScalar + ExactFieldScalar>(
    f: &RationalFunction<S>,
    arg: &RationalFunction<S>,
    span: Span,
) -> OghamResult<RationalFunction<S>> {
    let num = eval_poly_at_rational_function(f.num(), arg);
    let den = eval_poly_at_rational_function(f.den(), arg);
    if den.is_zero() {
        return Err(OghamError::new(
            OghamErrorKind::DivisionByZero,
            span,
            "rational-function evaluation hit a pole",
        ));
    }
    Ok(num.mul(&den.inv().expect("checked nonzero rational function")))
}

fn eval_poly_at_rational_function<S: ExactFieldScalar>(
    poly: &Poly<S>,
    x: &RationalFunction<S>,
) -> RationalFunction<S> {
    let mut acc = RationalFunction::zero();
    for c in poly.coeffs().iter().rev() {
        acc = acc.mul(x).add(&RationalFunction::from_base(c.clone()));
    }
    acc
}

fn expression_is_index(expr: &Expr) -> bool {
    match expr {
        Expr::Call { name, .. } if name == "deg" => true,
        Expr::Unary { expr, .. } => expression_is_index(expr),
        Expr::Binary {
            op: BinaryOp::Add | BinaryOp::Sub | BinaryOp::Mul,
            lhs,
            rhs,
        } => expression_is_index(lhs) || expression_is_index(rhs),
        Expr::Binary {
            op: BinaryOp::Pow,
            lhs,
            rhs,
        } => expression_is_index(lhs) || (plain_index_expr(lhs) && expression_is_index(rhs)),
        _ => false,
    }
}

fn plain_index_expr(expr: &Expr) -> bool {
    match expr {
        Expr::Int(_) => true,
        Expr::Call { name, .. } if name == "deg" => true,
        Expr::Unary {
            op: UnaryOp::Neg,
            expr,
        } => plain_index_expr(expr),
        Expr::Binary {
            op: BinaryOp::Add | BinaryOp::Sub | BinaryOp::Mul | BinaryOp::Pow,
            lhs,
            rhs,
        } => plain_index_expr(lhs) && plain_index_expr(rhs),
        _ => false,
    }
}

fn eval_index_binary(op: BinaryOp, lhs: i128, rhs: i128) -> OghamResult<i128> {
    match op {
        BinaryOp::Add => lhs
            .checked_add(rhs)
            .ok_or_else(|| overflow("index addition overflowed i128")),
        BinaryOp::Sub => lhs
            .checked_sub(rhs)
            .ok_or_else(|| overflow("index subtraction overflowed i128")),
        BinaryOp::Mul => lhs
            .checked_mul(rhs)
            .ok_or_else(|| overflow("index multiplication overflowed i128")),
        BinaryOp::Pow => {
            if rhs < 0 {
                return Err(OghamError::new(
                    OghamErrorKind::Domain,
                    Span::point(0),
                    "index exponent must be non-negative",
                ));
            }
            checked_i128_pow(lhs, rhs as u128)
        }
        _ => Err(index_sort_error()),
    }
}

fn no_order_error() -> OghamError {
    OghamError::new(
        OghamErrorKind::WrongWorld,
        Span::point(0),
        "this world has no canonical order",
    )
}

fn integer_poly_gcd(
    lhs: &Poly<Integer>,
    rhs: &Poly<Integer>,
    span: Span,
) -> OghamResult<Poly<Integer>> {
    let lhs = integer_poly_to_rational(lhs);
    let rhs = integer_poly_to_rational(rhs);
    primitive_integer_poly_from_rational(&lhs.gcd(&rhs), span)
}

fn integer_poly_to_rational(p: &Poly<Integer>) -> Poly<Rational> {
    Poly::new(p.coeffs().iter().map(|c| Rational::from_int(c.0)).collect())
}

fn primitive_integer_poly_from_rational(
    p: &Poly<Rational>,
    span: Span,
) -> OghamResult<Poly<Integer>> {
    if p.is_zero() {
        return Ok(Poly::zero());
    }
    let mut scale = 1i128;
    for c in p.coeffs() {
        scale = lcm_positive_i128(scale, c.denom(), span)?;
    }
    let mut coeffs = Vec::with_capacity(p.coeffs().len());
    for c in p.coeffs() {
        let factor = scale / c.denom();
        coeffs.push(
            c.numer()
                .checked_mul(factor)
                .ok_or_else(|| overflow("integer polynomial gcd coefficient overflowed i128"))?,
        );
    }
    let content = gcd_i128_slice(&coeffs, span)?;
    if content > 1 {
        for c in &mut coeffs {
            *c /= content;
        }
    }
    if coeffs.last().is_some_and(|c| *c < 0) {
        for c in &mut coeffs {
            *c = c.checked_neg().ok_or_else(|| {
                overflow("integer polynomial gcd sign normalization overflowed i128")
            })?;
        }
    }
    Ok(Poly::new(coeffs.into_iter().map(Integer).collect()))
}

fn gcd_i128_slice(values: &[i128], span: Span) -> OghamResult<i128> {
    let mut g = 0u128;
    for value in values {
        g = gcd_u128_local(g, value.unsigned_abs());
    }
    i128::try_from(g).map_err(|_| {
        OghamError::new(
            OghamErrorKind::Overflow,
            span,
            "integer polynomial gcd content exceeds i128",
        )
    })
}

fn lcm_positive_i128(lhs: i128, rhs: i128, span: Span) -> OghamResult<i128> {
    debug_assert!(lhs > 0 && rhs > 0);
    let gcd = gcd_u128_local(lhs as u128, rhs as u128);
    let gcd = i128::try_from(gcd).map_err(|_| {
        OghamError::new(
            OghamErrorKind::Overflow,
            span,
            "integer polynomial denominator gcd exceeds i128",
        )
    })?;
    lhs.checked_div(gcd)
        .and_then(|x| x.checked_mul(rhs))
        .ok_or_else(|| overflow("integer polynomial denominator lcm overflowed i128"))
}

fn gcd_u128_local(mut lhs: u128, mut rhs: u128) -> u128 {
    while rhs != 0 {
        let next = lhs % rhs;
        lhs = rhs;
        rhs = next;
    }
    lhs
}

trait OghamScalar: Scalar + Sized + Display + 'static {
    fn bare_int(n: u128, span: Span) -> OghamResult<Self>;
    fn star(lit: &StarLiteral, span: Span) -> OghamResult<Self>;
    fn omega(span: Span) -> OghamResult<Self>;
    fn omega_pow(_exp: Self, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::ExpSort,
            span,
            "`ω↑s` is only an element-level monomial constructor in surreal-family worlds",
        ))
    }
    fn named_element(_name: &str, _span: Span) -> OghamResult<Option<Self>> {
        Ok(None)
    }
    fn reserved_ident(_name: &str) -> bool {
        false
    }
    fn factorial(n: i128, span: Span) -> OghamResult<Self>;
    fn inv_scalar(value: &Self, span: Span) -> OghamResult<Self> {
        value
            .inv()
            .ok_or_else(|| OghamError::new(OghamErrorKind::NotInvertible, span, "not invertible"))
    }
    fn exact_div(_lhs: &Self, _rhs: &Self, _span: Span) -> Option<OghamResult<Self>> {
        None
    }
    fn rem(_lhs: &Self, _rhs: &Self, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "field worlds have no informative remainder operator",
        ))
    }
    fn relation(_op: RelOp, _lhs: &Self, _rhs: &Self, span: Span) -> OghamResult<bool> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "this world has no canonical order",
        ))
    }
    fn frob(_value: &Self, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`frob` is only available in finite-field worlds",
        ))
    }
    fn trace(_value: &Self, _m: Option<i128>, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`tr` is only available in finite-field worlds",
        ))
    }
    fn mul_checked(lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<Self> {
        Ok(lhs.mul(rhs))
    }
    fn pow_checked(base: &Self, mut k: u128, span: Span) -> OghamResult<Self> {
        if k == 0 {
            return Ok(Self::one());
        }
        let mut acc = Self::one();
        let mut x = base.clone();
        loop {
            if k & 1 == 1 {
                acc = Self::mul_checked(&acc, &x, span)?;
            }
            k >>= 1;
            if k == 0 {
                break;
            }
            x = Self::mul_checked(&x, &x, span)?;
        }
        Ok(acc)
    }
    fn mv_mul(
        alg: &CliffordAlgebra<Self>,
        lhs: &Multivector<Self>,
        rhs: &Multivector<Self>,
        _span: Span,
    ) -> OghamResult<Multivector<Self>> {
        Ok(alg.mul(lhs, rhs))
    }
    fn mv_pow(
        alg: &CliffordAlgebra<Self>,
        value: &Multivector<Self>,
        k: u128,
        _span: Span,
    ) -> OghamResult<Multivector<Self>> {
        Ok(alg.pow(value, k))
    }
}

impl OghamScalar for Nimber {
    fn bare_int(n: u128, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::BareInt,
            span,
            format!("bare integer `{n}` is not a nimber literal"),
        )
        .with_hint(format!("did you mean `*{n}`?")))
    }

    fn star(lit: &StarLiteral, span: Span) -> OghamResult<Self> {
        match lit {
            StarLiteral::Finite(n) => Ok(Nimber(*n)),
            StarLiteral::Cnf(_) => Err(OghamError::new(
                OghamErrorKind::WrongWorld,
                span,
                "transfinite star-literals belong to the `ordinal` world",
            )),
        }
    }

    fn omega(span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`ω` is not a finite nimber literal",
        ))
    }

    fn factorial(n: i128, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::BareInt,
            span,
            format!("`!{n}` would land through a bare integer in a nim-world"),
        ))
    }

    fn relation(op: RelOp, lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<bool> {
        Ok(match op {
            RelOp::Lt | RelOp::Gt => false,
            RelOp::Fuzzy => lhs.fuzzy(rhs),
            RelOp::Eq => lhs == rhs,
        })
    }

    fn frob(value: &Self, _span: Span) -> OghamResult<Self> {
        Ok(value.frobenius())
    }

    fn trace(value: &Self, m: Option<i128>, span: Span) -> OghamResult<Self> {
        let Some(m) = m else {
            return Err(OghamError::new(
                OghamErrorKind::Arity,
                span,
                "`tr` in the nimber world expects `tr(x, m)`",
            ));
        };
        if m <= 0 {
            return Err(OghamError::new(
                OghamErrorKind::Domain,
                span,
                "nimber trace degree must be positive",
            ));
        }
        Ok(Nimber(nim_trace(value.0, m as u128)))
    }
}

impl OghamScalar for Ordinal {
    fn bare_int(n: u128, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::BareInt,
            span,
            format!("bare integer `{n}` is not an ordinal-nimber value"),
        )
        .with_hint(format!("did you mean `*{n}`?")))
    }

    fn star(lit: &StarLiteral, _span: Span) -> OghamResult<Self> {
        Ok(match lit {
            StarLiteral::Finite(n) => Ordinal::from_u128(*n),
            StarLiteral::Cnf(cnf) => cnf.clone(),
        })
    }

    fn omega(span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::BareOrdinal,
            span,
            "bare `ω` is an ordinal address, not a value",
        )
        .with_hint("values are starred here: `*ω`"))
    }

    fn factorial(n: i128, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::BareInt,
            span,
            format!("`!{n}` would land through a bare integer in a nim-world"),
        ))
    }

    fn inv_scalar(value: &Self, span: Span) -> OghamResult<Self> {
        if value.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "division by zero",
            ));
        }
        value.checked_inv().ok_or_else(|| kummer_escape(span))
    }

    fn relation(op: RelOp, lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<bool> {
        Ok(match op {
            RelOp::Lt | RelOp::Gt => false,
            RelOp::Fuzzy => lhs.fuzzy(rhs),
            RelOp::Eq => lhs == rhs,
        })
    }

    fn mul_checked(lhs: &Self, rhs: &Self, span: Span) -> OghamResult<Self> {
        lhs.nim_mul(rhs).ok_or_else(|| kummer_escape(span))
    }

    fn pow_checked(base: &Self, k: u128, span: Span) -> OghamResult<Self> {
        base.nim_pow(k).ok_or_else(|| kummer_escape(span))
    }

    fn mv_mul(
        alg: &CliffordAlgebra<Self>,
        lhs: &Multivector<Self>,
        rhs: &Multivector<Self>,
        span: Span,
    ) -> OghamResult<Multivector<Self>> {
        catch_unwind(AssertUnwindSafe(|| alg.mul(lhs, rhs))).map_err(|_| kummer_escape(span))
    }

    fn mv_pow(
        alg: &CliffordAlgebra<Self>,
        value: &Multivector<Self>,
        k: u128,
        span: Span,
    ) -> OghamResult<Multivector<Self>> {
        catch_unwind(AssertUnwindSafe(|| alg.pow(value, k))).map_err(|_| kummer_escape(span))
    }
}

impl OghamScalar for Surreal {
    fn bare_int(n: u128, _span: Span) -> OghamResult<Self> {
        Ok(Surreal::from_int(u128_to_i128(n)?))
    }

    fn star(_lit: &StarLiteral, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`*3` is a nimber; this is the `surreal` world",
        ))
    }

    fn omega(_span: Span) -> OghamResult<Self> {
        Ok(Surreal::omega())
    }

    fn omega_pow(exp: Self, _span: Span) -> OghamResult<Self> {
        Ok(Surreal::omega_pow(exp))
    }

    fn factorial(n: i128, _span: Span) -> OghamResult<Self> {
        if n < 0 {
            return Err(domain("factorial is only defined for n >= 0"));
        }
        let n = checked_factorial_i128(n).ok_or_else(|| overflow("factorial exceeds i128"))?;
        Ok(Surreal::from_int(n))
    }

    fn inv_scalar(value: &Self, span: Span) -> OghamResult<Self> {
        if value.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "division by zero",
            ));
        }
        value.inv().ok_or_else(|| {
            OghamError::new(
                OghamErrorKind::NotInvertible,
                span,
                "only CNF monomials invert exactly; 1/(ω+1) is an infinite Hahn series",
            )
        })
    }

    fn rem(lhs: &Self, rhs: &Self, span: Span) -> OghamResult<Self> {
        if rhs.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "division by zero",
            ));
        }
        lhs.rem(rhs).ok_or_else(|| modulus_error(span))
    }

    fn relation(op: RelOp, lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<bool> {
        ordered_relation(op, lhs.cmp(rhs))
    }
}

impl OghamScalar for Omnific {
    fn bare_int(n: u128, _span: Span) -> OghamResult<Self> {
        Ok(Omnific::from_int(u128_to_i128(n)?))
    }

    fn star(_lit: &StarLiteral, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`*3` is a nimber; this is the `omnific` world",
        ))
    }

    fn omega(_span: Span) -> OghamResult<Self> {
        Ok(Omnific::omega())
    }

    fn omega_pow(exp: Self, span: Span) -> OghamResult<Self> {
        Omnific::from_surreal(Surreal::omega_pow(exp.inner().clone())).ok_or_else(|| {
            OghamError::new(
                OghamErrorKind::Domain,
                span,
                "omega-power exponent does not produce an omnific integer",
            )
        })
    }

    fn factorial(n: i128, _span: Span) -> OghamResult<Self> {
        if n < 0 {
            return Err(domain("factorial is only defined for n >= 0"));
        }
        let n = checked_factorial_i128(n).ok_or_else(|| overflow("factorial exceeds i128"))?;
        Ok(Omnific::from_int(n))
    }

    fn rem(lhs: &Self, rhs: &Self, span: Span) -> OghamResult<Self> {
        if rhs.is_zero() {
            return Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "division by zero",
            ));
        }
        lhs.rem(rhs).ok_or_else(|| modulus_error(span))
    }

    fn relation(op: RelOp, lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<bool> {
        ordered_relation(op, lhs.cmp(rhs))
    }
}

impl OghamScalar for Integer {
    fn bare_int(n: u128, _span: Span) -> OghamResult<Self> {
        Ok(Integer(u128_to_i128(n)?))
    }

    fn star(_lit: &StarLiteral, span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`*3` is a nimber; this is the `integer` world",
        ))
    }

    fn omega(span: Span) -> OghamResult<Self> {
        Err(OghamError::new(
            OghamErrorKind::WrongWorld,
            span,
            "`ω` belongs to the surreal-family worlds",
        ))
    }

    fn factorial(n: i128, _span: Span) -> OghamResult<Self> {
        if n < 0 {
            return Err(domain("factorial is only defined for n >= 0"));
        }
        let n = checked_factorial_i128(n).ok_or_else(|| overflow("factorial exceeds i128"))?;
        Ok(Integer(n))
    }

    fn exact_div(lhs: &Self, rhs: &Self, span: Span) -> Option<OghamResult<Self>> {
        Some(match lhs.div_exact(rhs) {
            Ok(q) => Ok(q),
            Err(IntegerDivExactError::DivisionByZero) => Err(OghamError::new(
                OghamErrorKind::DivisionByZero,
                span,
                "division by zero",
            )),
            Err(IntegerDivExactError::Remainder(r)) => Err(OghamError::new(
                OghamErrorKind::NotInvertible,
                span,
                format!("integer exact division failed with remainder {r}"),
            )),
        })
    }

    fn rem(lhs: &Self, rhs: &Self, span: Span) -> OghamResult<Self> {
        lhs.rem(rhs).ok_or_else(|| {
            OghamError::new(OghamErrorKind::DivisionByZero, span, "division by zero")
        })
    }

    fn relation(op: RelOp, lhs: &Self, rhs: &Self, _span: Span) -> OghamResult<bool> {
        ordered_relation(op, lhs.cmp(rhs))
    }
}

macro_rules! impl_fp_ogham {
    ($($p:literal),* $(,)?) => {
        $(
            impl OghamScalar for Fp<$p> {
                fn bare_int(n: u128, _span: Span) -> OghamResult<Self> {
                    Ok(Fp::<$p>::from_u128(n))
                }
                fn star(_lit: &StarLiteral, span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "`*3` is a nimber; this is a prime-field world",
                    ))
                }
                fn omega(span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "`ω` belongs to the surreal-family worlds",
                    ))
                }
                fn factorial(n: i128, _span: Span) -> OghamResult<Self> {
                    factorial_in_scalar::<Self>(n).ok_or_else(|| domain("factorial is only defined for n >= 0"))
                }
                fn rem(_lhs: &Self, _rhs: &Self, span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "field worlds have no informative remainder operator",
                    ))
                }
                fn frob(value: &Self, _span: Span) -> OghamResult<Self> {
                    Ok(*value)
                }
                fn trace(value: &Self, m: Option<i128>, span: Span) -> OghamResult<Self> {
                    if m.is_some() {
                        return Err(OghamError::new(
                            OghamErrorKind::Arity,
                            span,
                            "`tr` in prime fields expects one argument",
                        ));
                    }
                    Ok(*value)
                }
            }
        )*
    };
}

macro_rules! impl_fpn_ogham {
    ($(($p:literal, $n:literal)),* $(,)?) => {
        $(
            impl OghamScalar for Fpn<$p, $n> {
                fn bare_int(n: u128, _span: Span) -> OghamResult<Self> {
                    Ok(Fpn::<$p, $n>::constant(n))
                }
                fn star(_lit: &StarLiteral, span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "`*3` is a nimber; this is an extension-field world",
                    ))
                }
                fn omega(span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "`ω` belongs to the surreal-family worlds",
                    ))
                }
                fn named_element(name: &str, _span: Span) -> OghamResult<Option<Self>> {
                    Ok((name == "x").then(Fpn::<$p, $n>::generator))
                }
                fn reserved_ident(name: &str) -> bool {
                    name == "x"
                }
                fn factorial(n: i128, _span: Span) -> OghamResult<Self> {
                    factorial_in_scalar::<Self>(n).ok_or_else(|| domain("factorial is only defined for n >= 0"))
                }
                fn rem(_lhs: &Self, _rhs: &Self, span: Span) -> OghamResult<Self> {
                    Err(OghamError::new(
                        OghamErrorKind::WrongWorld,
                        span,
                        "field worlds have no informative remainder operator",
                    ))
                }
                fn frob(value: &Self, _span: Span) -> OghamResult<Self> {
                    Ok(value.frobenius())
                }
                fn trace(value: &Self, m: Option<i128>, span: Span) -> OghamResult<Self> {
                    if m.is_some() {
                        return Err(OghamError::new(
                            OghamErrorKind::Arity,
                            span,
                            "`tr` in extension fields expects one argument",
                        ));
                    }
                    Ok(value.relative_trace(1))
                }
            }
        )*
    };
}

impl_fp_ogham!(2, 3, 5, 7);
impl_fpn_ogham!((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (5, 2));

fn build_runtime<S: OghamScalar>(
    name: &'static str,
    dim: usize,
    rest: &str,
) -> OghamResult<Runtime<S>> {
    let metric = if rest.trim().is_empty() {
        if dim == 0 {
            Metric::diagonal(Vec::new())
        } else {
            return Err(parse_error(
                "positive-dimensional worlds need `q=[...]` or `grassmann`",
            ));
        }
    } else if rest.contains("grassmann") {
        Metric::grassmann(dim)
    } else {
        let q_src = extract_bracket(rest, "q")?;
        let q = parse_scalar_list::<S>(&q_src)?;
        if q.len() != dim {
            return Err(OghamError::new(
                OghamErrorKind::DimMismatch,
                Span::point(0),
                format!("q length {} does not match dimension {dim}", q.len()),
            ));
        }
        let b = if let Some(b_src) = extract_bracket_opt(rest, "b")? {
            parse_sparse_pairs::<S>(&b_src)?
        } else {
            BTreeMap::new()
        };
        let a = if let Some(a_src) = extract_bracket_opt(rest, "a")? {
            parse_sparse_pairs::<S>(&a_src)?
        } else {
            BTreeMap::new()
        };
        Metric::general(q, b, a)
    };
    Ok(Runtime::from_metric(name, metric))
}

fn parse_gold_metric(src: &str) -> OghamResult<Metric<Nimber>> {
    let inner = src
        .strip_prefix("gold(")
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| parse_error("expected `gold(m,a)`"))?;
    let mut parts = inner.split(',');
    let m = parts
        .next()
        .ok_or_else(|| parse_error("missing gold m"))?
        .trim()
        .parse::<usize>()
        .map_err(|_| parse_error("gold m must be a usize"))?;
    let a = parts
        .next()
        .ok_or_else(|| parse_error("missing gold a"))?
        .trim()
        .parse::<usize>()
        .map_err(|_| parse_error("gold a must be a usize"))?;
    if parts.next().is_some() {
        return Err(parse_error("gold expects exactly two arguments"));
    }
    Ok(crate::forms::gold_form(m, a))
}

fn parse_scalar_list<S: OghamScalar>(src: &str) -> OghamResult<Vec<S>> {
    if src.trim().is_empty() {
        return Ok(Vec::new());
    }
    split_top_level(src, ',')
        .into_iter()
        .map(|part| parse_metric_scalar::<S>(&part))
        .collect()
}

fn parse_sparse_pairs<S: OghamScalar>(src: &str) -> OghamResult<BTreeMap<(usize, usize), S>> {
    let mut out = BTreeMap::new();
    if src.trim().is_empty() {
        return Ok(out);
    }
    for entry in split_top_level(src, ',') {
        let (ij, value) = entry
            .split_once(':')
            .ok_or_else(|| parse_error("sparse metric entries need `(i,j):value`"))?;
        let ij = ij.trim();
        let ij = ij
            .strip_prefix('(')
            .and_then(|s| s.strip_suffix(')'))
            .ok_or_else(|| parse_error("sparse metric key needs `(i,j)`"))?;
        let (i, j) = ij
            .split_once(',')
            .ok_or_else(|| parse_error("sparse metric key needs two indices"))?;
        let i = i
            .trim()
            .parse::<usize>()
            .map_err(|_| parse_error("metric index must be a usize"))?;
        let j = j
            .trim()
            .parse::<usize>()
            .map_err(|_| parse_error("metric index must be a usize"))?;
        out.insert((i, j), parse_metric_scalar::<S>(value)?);
    }
    Ok(out)
}

fn parse_metric_scalar<S: OghamScalar>(src: &str) -> OghamResult<S> {
    let mut rt = Runtime::<S>::from_metric("metric", Metric::diagonal(Vec::new()));
    let stmt = parse_statement(src)?;
    let Statement::Expr(expr) = stmt else {
        return Err(parse_error("metric scalar must be an expression"));
    };
    let value = rt.eval_expr(&expr)?;
    scalar_part(&value).ok_or_else(|| grade0_error(Span::point(0)))
}

fn extract_bracket(rest: &str, key: &str) -> OghamResult<String> {
    extract_bracket_opt(rest, key)?.ok_or_else(|| parse_error(format!("missing `{key}=[...]`")))
}

fn extract_bracket_opt(rest: &str, key: &str) -> OghamResult<Option<String>> {
    let needle = format!("{key}=");
    let Some(start) = rest.find(&needle) else {
        return Ok(None);
    };
    let after = &rest[start + needle.len()..];
    let Some(open) = after.find('[') else {
        return Err(parse_error(format!("`{key}` needs `[...]`")));
    };
    let mut depth = 0i32;
    let mut begin = None;
    for (idx, ch) in after[open..].char_indices() {
        match ch {
            '[' => {
                if depth == 0 {
                    begin = Some(open + idx + ch.len_utf8());
                }
                depth += 1;
            }
            ']' => {
                depth -= 1;
                if depth == 0 {
                    let begin = begin.expect("set at opening bracket");
                    return Ok(Some(after[begin..open + idx].to_string()));
                }
            }
            _ => {}
        }
    }
    Err(parse_error(format!("unterminated `{key}` bracket list")))
}

fn split_top_level(src: &str, delim: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut start = 0usize;
    let mut parens = 0i32;
    let mut brackets = 0i32;
    for (idx, ch) in src.char_indices() {
        match ch {
            '(' => parens += 1,
            ')' => parens -= 1,
            '[' => brackets += 1,
            ']' => brackets -= 1,
            c if c == delim && parens == 0 && brackets == 0 => {
                out.push(src[start..idx].trim().to_string());
                start = idx + ch.len_utf8();
            }
            _ => {}
        }
    }
    out.push(src[start..].trim().to_string());
    out
}

fn scalar_part<S: Scalar>(value: &Multivector<S>) -> Option<S> {
    match value.terms() {
        terms if terms.is_empty() => Some(S::zero()),
        terms if terms.len() == 1 => terms.get(&0).cloned(),
        _ => None,
    }
}

fn expect_arity(name: &str, args: &[Expr], expected: usize) -> OghamResult<()> {
    if args.len() == expected {
        Ok(())
    } else {
        Err(OghamError::new(
            OghamErrorKind::Arity,
            Span::point(0),
            format!("`{name}` expects {expected} argument(s)"),
        ))
    }
}

fn ordered_relation(op: RelOp, cmp: Ordering) -> OghamResult<bool> {
    Ok(match op {
        RelOp::Eq => cmp == Ordering::Equal,
        RelOp::Lt => cmp == Ordering::Less,
        RelOp::Gt => cmp == Ordering::Greater,
        RelOp::Fuzzy => false,
    })
}

fn checked_i128_pow(base: i128, mut exp: u128) -> OghamResult<i128> {
    if exp == 0 {
        return Ok(1);
    }
    let mut acc = 1i128;
    let mut x = base;
    loop {
        if exp & 1 == 1 {
            acc = acc
                .checked_mul(x)
                .ok_or_else(|| overflow("index power overflowed i128"))?;
        }
        exp >>= 1;
        if exp == 0 {
            break;
        }
        x = x
            .checked_mul(x)
            .ok_or_else(|| overflow("index power overflowed i128"))?;
    }
    Ok(acc)
}

fn u128_to_i128(n: u128) -> OghamResult<i128> {
    i128::try_from(n).map_err(|_| overflow("integer literal exceeds i128 in this world"))
}

fn parse_error(message: impl Into<String>) -> OghamError {
    OghamError::new(OghamErrorKind::Parse, Span::point(0), message)
}

fn index_sort_error() -> OghamError {
    OghamError::new(
        OghamErrorKind::IndexSort,
        Span::point(0),
        "expected an Index expression",
    )
}

fn exp_sort_error() -> OghamError {
    OghamError::new(
        OghamErrorKind::ExpSort,
        Span::point(0),
        "exponent must be an Index",
    )
    .with_hint("`↑`/`^` is power; the wedge product is `∧`/`&`")
}

fn grade0_error(span: Span) -> OghamError {
    OghamError::new(
        OghamErrorKind::Grade0,
        span,
        "operation requires a grade-0 element",
    )
}

fn modulus_error(span: Span) -> OghamError {
    OghamError::new(
        OghamErrorKind::Modulus,
        span,
        "moduli here are monic omega-powers: `% ω↑2` truncates the CNF below it",
    )
}

fn kummer_escape(span: Span) -> OghamError {
    OghamError::new(
        OghamErrorKind::KummerEscape,
        span,
        "ordinal nim-product escaped beyond the source-verified tower below ω^(ω^ω)",
    )
    .with_hint("below ω^(ω^ω), primes <= 47 — see OPEN.md")
}

fn overflow(message: impl Into<String>) -> OghamError {
    OghamError::new(OghamErrorKind::Overflow, Span::point(0), message)
}

fn domain(message: impl Into<String>) -> OghamError {
    OghamError::new(OghamErrorKind::Domain, Span::point(0), message)
}
