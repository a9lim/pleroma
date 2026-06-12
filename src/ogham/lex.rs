use super::error::{OghamError, OghamErrorKind, OghamResult, Span};

#[derive(Clone, Debug, PartialEq)]
pub struct Token {
    pub kind: TokenKind,
    pub span: Span,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TokenKind {
    Int(u128),
    Ident(String),
    Blade(usize),
    Omega,
    Star,
    Up,
    Wedge,
    Dot,
    Slash,
    Percent,
    At,
    Bang,
    Question,
    Colon,
    Arrow,
    And,
    Or,
    Not,
    True,
    False,
    Semicolon,
    Eq,
    Less,
    Greater,
    Pipe,
    Assign,
    Plus,
    Minus,
    LParen,
    RParen,
    LBracket,
    RBracket,
    Comma,
}

pub fn lex(src: &str) -> OghamResult<Vec<Token>> {
    let src = src.split_once('#').map_or(src, |(head, _)| head);
    let chars: Vec<(usize, char)> = src.char_indices().collect();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < chars.len() {
        let (pos, ch) = chars[i];
        if ch.is_whitespace() {
            i += 1;
            continue;
        }
        if ch.is_ascii_digit() {
            let start = pos;
            let mut end = pos + ch.len_utf8();
            let mut value = 0u128;
            while i < chars.len() {
                let (p, c) = chars[i];
                if !c.is_ascii_digit() {
                    break;
                }
                let digit = u128::from(c as u8 - b'0');
                value = value
                    .checked_mul(10)
                    .and_then(|v| v.checked_add(digit))
                    .ok_or_else(|| {
                        OghamError::new(
                            OghamErrorKind::Overflow,
                            Span::new(start, p + c.len_utf8()),
                            "integer literal exceeds u128",
                        )
                    })?;
                end = p + c.len_utf8();
                i += 1;
            }
            out.push(Token {
                kind: TokenKind::Int(value),
                span: Span::new(start, end),
            });
            continue;
        }
        if ch == 'e' && i + 1 < chars.len() && chars[i + 1].1.is_ascii_digit() {
            let start = pos;
            i += 1;
            let mut end = start + 1;
            let mut value = 0usize;
            while i < chars.len() {
                let (p, c) = chars[i];
                if !c.is_ascii_digit() {
                    break;
                }
                let digit = usize::from(c as u8 - b'0');
                value = value
                    .checked_mul(10)
                    .and_then(|v| v.checked_add(digit))
                    .ok_or_else(|| {
                        OghamError::new(
                            OghamErrorKind::Overflow,
                            Span::new(start, p + c.len_utf8()),
                            "blade index exceeds usize",
                        )
                    })?;
                end = p + c.len_utf8();
                i += 1;
            }
            out.push(Token {
                kind: TokenKind::Blade(value),
                span: Span::new(start, end),
            });
            continue;
        }
        if ch == 'O' && i + 1 < chars.len() && chars[i + 1].1 == '(' {
            return Err(reserved(Span::new(pos, chars[i + 1].0 + 1)));
        }
        if ch.is_ascii_lowercase() {
            let start = pos;
            let mut s = String::new();
            let mut end = pos + ch.len_utf8();
            while i < chars.len() {
                let (p, c) = chars[i];
                if !(c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_') {
                    break;
                }
                s.push(c);
                end = p + c.len_utf8();
                i += 1;
            }
            let kind = if s == "w" {
                TokenKind::Omega
            } else if s == "and" {
                TokenKind::And
            } else if s == "or" {
                TokenKind::Or
            } else if s == "not" {
                TokenKind::Not
            } else if s == "true" {
                TokenKind::True
            } else if s == "false" {
                TokenKind::False
            } else if s == "e" {
                return Err(OghamError::new(
                    OghamErrorKind::Parse,
                    Span::new(start, end),
                    "`e` needs a blade index, e.g. `e0`",
                ));
            } else {
                TokenKind::Ident(s)
            };
            out.push(Token {
                kind,
                span: Span::new(start, end),
            });
            continue;
        }
        let span = Span::new(pos, pos + ch.len_utf8());
        let kind = match ch {
            'ω' => TokenKind::Omega,
            '*' => TokenKind::Star,
            '^' | '↑' => {
                if i + 1 < chars.len() && matches!(chars[i + 1].1, '^' | '↑') {
                    return Err(reserved(Span::new(
                        pos,
                        chars[i + 1].0 + chars[i + 1].1.len_utf8(),
                    )));
                }
                TokenKind::Up
            }
            '&' | '∧' => TokenKind::Wedge,
            '.' | '⋅' | '·' => TokenKind::Dot,
            '/' => TokenKind::Slash,
            '%' => TokenKind::Percent,
            '@' => TokenKind::At,
            '!' => TokenKind::Bang,
            '?' => TokenKind::Question,
            '=' => {
                if i + 1 < chars.len() && chars[i + 1].1 == '=' {
                    i += 1;
                    out.push(Token {
                        kind: TokenKind::Eq,
                        span: Span::new(pos, chars[i].0 + 1),
                    });
                    i += 1;
                    continue;
                }
                TokenKind::Eq
            }
            '<' => TokenKind::Less,
            '>' => TokenKind::Greater,
            '|' => TokenKind::Pipe,
            ':' => {
                if i + 1 < chars.len() && chars[i + 1].1 == '=' {
                    i += 1;
                    out.push(Token {
                        kind: TokenKind::Assign,
                        span: Span::new(pos, chars[i].0 + 1),
                    });
                    i += 1;
                    continue;
                }
                TokenKind::Colon
            }
            '+' => TokenKind::Plus,
            '-' => TokenKind::Minus,
            ';' => TokenKind::Semicolon,
            '(' => TokenKind::LParen,
            ')' => TokenKind::RParen,
            '[' => TokenKind::LBracket,
            ']' => TokenKind::RBracket,
            ',' => TokenKind::Comma,
            '↦' | '~' => TokenKind::Arrow,
            '{' | '}' => return Err(reserved(span)),
            _ => {
                return Err(OghamError::new(
                    OghamErrorKind::Parse,
                    span,
                    format!("unexpected character `{ch}`"),
                ));
            }
        };
        out.push(Token { kind, span });
        i += 1;
    }
    Ok(out)
}

fn reserved(span: Span) -> OghamError {
    OghamError::new(OghamErrorKind::Reserved, span, "reserved syntax")
        .with_hint("reserved for future games/precision/function syntax")
}
