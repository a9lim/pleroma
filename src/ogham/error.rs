use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Span { start, end }
    }

    pub fn point(pos: usize) -> Self {
        Span {
            start: pos,
            end: pos,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OghamErrorKind {
    Parse,
    Reserved,
    ExpSort,
    IndexSort,
    BoolSort,
    FnSort,
    Shadow,
    SeqValue,
    BareInt,
    BareOrdinal,
    WrongWorld,
    CnfOrder,
    KummerEscape,
    NotInvertible,
    DivisionByZero,
    BladeIndex,
    DimMismatch,
    GeneralMetric,
    Unbound,
    Arity,
    UnknownFn,
    Grade0,
    Modulus,
    Overflow,
    Domain,
}

impl OghamErrorKind {
    pub fn code(self) -> &'static str {
        match self {
            OghamErrorKind::Parse => "E_Parse",
            OghamErrorKind::Reserved => "E_Reserved",
            OghamErrorKind::ExpSort => "E_ExpSort",
            OghamErrorKind::IndexSort => "E_IndexSort",
            OghamErrorKind::BoolSort => "E_BoolSort",
            OghamErrorKind::FnSort => "E_FnSort",
            OghamErrorKind::Shadow => "E_Shadow",
            OghamErrorKind::SeqValue => "E_SeqValue",
            OghamErrorKind::BareInt => "E_BareInt",
            OghamErrorKind::BareOrdinal => "E_BareOrdinal",
            OghamErrorKind::WrongWorld => "E_WrongWorld",
            OghamErrorKind::CnfOrder => "E_CnfOrder",
            OghamErrorKind::KummerEscape => "E_KummerEscape",
            OghamErrorKind::NotInvertible => "E_NotInvertible",
            OghamErrorKind::DivisionByZero => "E_DivisionByZero",
            OghamErrorKind::BladeIndex => "E_BladeIndex",
            OghamErrorKind::DimMismatch => "E_DimMismatch",
            OghamErrorKind::GeneralMetric => "E_GeneralMetric",
            OghamErrorKind::Unbound => "E_Unbound",
            OghamErrorKind::Arity => "E_Arity",
            OghamErrorKind::UnknownFn => "E_UnknownFn",
            OghamErrorKind::Grade0 => "E_Grade0",
            OghamErrorKind::Modulus => "E_Modulus",
            OghamErrorKind::Overflow => "E_Overflow",
            OghamErrorKind::Domain => "E_Domain",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OghamError {
    pub kind: OghamErrorKind,
    pub span: Span,
    pub message: String,
    pub hint: Option<String>,
}

impl OghamError {
    pub fn new(kind: OghamErrorKind, span: Span, message: impl Into<String>) -> Self {
        OghamError {
            kind,
            span,
            message: message.into(),
            hint: None,
        }
    }

    pub fn with_hint(mut self, hint: impl Into<String>) -> Self {
        self.hint = Some(hint.into());
        self
    }
}

impl fmt::Display for OghamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.kind.code(), self.message)?;
        if let Some(hint) = &self.hint {
            write!(f, " ({hint})")?;
        }
        Ok(())
    }
}

impl std::error::Error for OghamError {}

pub type OghamResult<T> = Result<T, OghamError>;
