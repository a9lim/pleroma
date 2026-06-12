use crate::scalar::Ordinal;

#[derive(Clone, Debug, PartialEq)]
pub enum Statement {
    Binding { name: String, expr: Expr },
    Expr(Expr),
}

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Int(u128),
    Bool(bool),
    Star(StarLiteral),
    Omega,
    Blade(usize),
    Vector(Vec<Expr>),
    Tuple(Vec<Expr>),
    Ident(String),
    Lambda {
        binders: Vec<String>,
        body: Box<Expr>,
    },
    Call {
        name: String,
        args: Vec<Expr>,
    },
    Factorial(Box<Expr>),
    Unary {
        op: UnaryOp,
        expr: Box<Expr>,
    },
    Binary {
        op: BinaryOp,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    Ternary {
        cond: Box<Expr>,
        then_expr: Box<Expr>,
        else_expr: Box<Expr>,
    },
    Relation {
        op: RelOp,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
}

impl Expr {
    pub fn is_omega_atom(&self) -> bool {
        matches!(self, Expr::Omega)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum StarLiteral {
    Finite(u128),
    Cnf(Ordinal),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnaryOp {
    Neg,
    Inv,
    Not,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinaryOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
    Wedge,
    Pow,
    At,
    And,
    Or,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RelOp {
    Eq,
    Lt,
    Gt,
    Fuzzy,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sort {
    Element,
    Index,
    Bool,
}
