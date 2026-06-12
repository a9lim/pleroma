use super::ast::{BinaryOp, Expr, RelOp, StarLiteral, Statement, UnaryOp};
use super::error::{OghamError, OghamErrorKind, OghamResult, Span};
use super::lex::{lex, Token, TokenKind};
use crate::scalar::Ordinal;

pub fn parse_statement(src: &str) -> OghamResult<Statement> {
    let tokens = lex(src)?;
    if tokens
        .iter()
        .any(|tok| matches!(tok.kind, TokenKind::Semicolon))
    {
        return Err(OghamError::new(
            OghamErrorKind::SeqValue,
            Span::point(0),
            "sequencing is reserved for value-discarding program statements",
        ));
    }
    let mut parser = Parser { tokens, pos: 0 };
    if parser.tokens.is_empty() {
        return Err(OghamError::new(
            OghamErrorKind::Parse,
            Span::point(0),
            "empty statement",
        ));
    }
    if parser.is_reserved_word_binding() {
        return Err(OghamError::new(
            OghamErrorKind::Reserved,
            parser.span(),
            "reserved word cannot be rebound",
        ));
    }
    let stmt = if let (Some(TokenKind::Ident(name)), Some(TokenKind::Assign)) =
        (parser.peek_kind(), parser.peek_kind_at(1))
    {
        let name = name.clone();
        parser.bump();
        parser.bump();
        let expr = parser.parse_lambda_or_expression()?;
        Statement::Binding { name, expr }
    } else {
        Statement::Expr(parser.parse_lambda_or_expression()?)
    };
    parser.expect_end()?;
    Ok(stmt)
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn peek_kind(&self) -> Option<&TokenKind> {
        self.peek().map(|t| &t.kind)
    }

    fn peek_kind_at(&self, offset: usize) -> Option<&TokenKind> {
        self.tokens.get(self.pos + offset).map(|t| &t.kind)
    }

    fn span(&self) -> Span {
        self.peek().map_or(Span::point(0), |t| t.span)
    }

    fn bump(&mut self) -> Option<Token> {
        let tok = self.tokens.get(self.pos).cloned();
        if tok.is_some() {
            self.pos += 1;
        }
        tok
    }

    fn expect_end(&self) -> OghamResult<()> {
        if let Some(tok) = self.peek() {
            Err(OghamError::new(
                OghamErrorKind::Parse,
                tok.span,
                "unexpected trailing token",
            ))
        } else {
            Ok(())
        }
    }

    fn eat(&mut self, pred: impl FnOnce(&TokenKind) -> bool) -> Option<Token> {
        if self.peek_kind().is_some_and(pred) {
            self.bump()
        } else {
            None
        }
    }

    fn expect(&mut self, pred: impl FnOnce(&TokenKind) -> bool, what: &str) -> OghamResult<Token> {
        self.eat(pred).ok_or_else(|| {
            OghamError::new(
                OghamErrorKind::Parse,
                self.span(),
                format!("expected {what}"),
            )
        })
    }

    fn is_reserved_word_binding(&self) -> bool {
        matches!(
            self.peek_kind(),
            Some(
                TokenKind::And
                    | TokenKind::Or
                    | TokenKind::Not
                    | TokenKind::True
                    | TokenKind::False
            )
        ) && matches!(self.peek_kind_at(1), Some(TokenKind::Assign))
    }

    fn parse_lambda_or_expression(&mut self) -> OghamResult<Expr> {
        if let Some(binders) = self.try_parse_binders()? {
            self.expect(|k| matches!(k, TokenKind::Arrow), "`↦`")?;
            let body = self.parse_lambda_or_expression()?;
            return Ok(Expr::Lambda {
                binders,
                body: Box::new(body),
            });
        }
        self.parse_expression()
    }

    fn try_parse_binders(&mut self) -> OghamResult<Option<Vec<String>>> {
        let save = self.pos;
        let out = match self.peek_kind() {
            Some(TokenKind::Ident(name))
                if matches!(self.peek_kind_at(1), Some(TokenKind::Arrow)) =>
            {
                let name = name.clone();
                self.bump();
                Some(vec![name])
            }
            Some(TokenKind::LParen) => {
                self.bump();
                let mut binders = Vec::new();
                loop {
                    match self.bump() {
                        Some(Token {
                            kind: TokenKind::Ident(name),
                            ..
                        }) => binders.push(name),
                        _ => {
                            self.pos = save;
                            return Ok(None);
                        }
                    }
                    if !matches!(self.peek_kind(), Some(TokenKind::Comma)) {
                        break;
                    }
                    self.bump();
                }
                if !matches!(self.peek_kind(), Some(TokenKind::RParen)) {
                    self.pos = save;
                    return Ok(None);
                }
                self.bump();
                if matches!(self.peek_kind(), Some(TokenKind::Arrow)) {
                    Some(binders)
                } else {
                    self.pos = save;
                    return Ok(None);
                }
            }
            _ => None,
        };
        if out.is_none() {
            self.pos = save;
        }
        Ok(out)
    }

    fn parse_expression(&mut self) -> OghamResult<Expr> {
        let expr = self.parse_or()?;
        if !matches!(self.peek_kind(), Some(TokenKind::Question)) {
            return Ok(expr);
        }
        self.bump();
        let then_expr = self.parse_additive()?;
        self.expect(|k| matches!(k, TokenKind::Colon), "`:`")?;
        let else_expr = self.parse_additive()?;
        Ok(Expr::Ternary {
            cond: Box::new(expr),
            then_expr: Box::new(then_expr),
            else_expr: Box::new(else_expr),
        })
    }

    fn parse_or(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_and()?;
        while matches!(self.peek_kind(), Some(TokenKind::Or)) {
            self.bump();
            let rhs = self.parse_and()?;
            expr = Expr::Binary {
                op: BinaryOp::Or,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_and(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_not()?;
        while matches!(self.peek_kind(), Some(TokenKind::And)) {
            self.bump();
            let rhs = self.parse_not()?;
            expr = Expr::Binary {
                op: BinaryOp::And,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_not(&mut self) -> OghamResult<Expr> {
        if matches!(self.peek_kind(), Some(TokenKind::Not)) {
            self.bump();
            let expr = self.parse_not()?;
            return Ok(Expr::Unary {
                op: UnaryOp::Not,
                expr: Box::new(expr),
            });
        }
        self.parse_relation()
    }

    fn parse_relation(&mut self) -> OghamResult<Expr> {
        let lhs = self.parse_additive()?;
        let Some(op) = self.parse_relop() else {
            return Ok(lhs);
        };
        let rhs = self.parse_additive()?;
        if self.parse_relop().is_some() {
            return Err(OghamError::new(
                OghamErrorKind::Parse,
                self.span(),
                "relations are top-level and non-associative",
            ));
        }
        Ok(Expr::Relation {
            op,
            lhs: Box::new(lhs),
            rhs: Box::new(rhs),
        })
    }

    fn parse_relop(&mut self) -> Option<RelOp> {
        match self.peek_kind()? {
            TokenKind::Eq => {
                self.bump();
                Some(RelOp::Eq)
            }
            TokenKind::Less => {
                self.bump();
                Some(RelOp::Lt)
            }
            TokenKind::Greater => {
                self.bump();
                Some(RelOp::Gt)
            }
            TokenKind::Pipe => {
                self.bump();
                Some(RelOp::Fuzzy)
            }
            _ => None,
        }
    }

    fn parse_additive(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_mulexpr()?;
        loop {
            let op = match self.peek_kind() {
                Some(TokenKind::Plus) => BinaryOp::Add,
                Some(TokenKind::Minus) => BinaryOp::Sub,
                _ => break,
            };
            self.bump();
            let rhs = self.parse_mulexpr()?;
            expr = Expr::Binary {
                op,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_mulexpr(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_wedge()?;
        loop {
            let op = match self.peek_kind() {
                Some(TokenKind::Dot) => BinaryOp::Mul,
                Some(TokenKind::Slash) => BinaryOp::Div,
                Some(TokenKind::Percent) => BinaryOp::Rem,
                _ => break,
            };
            self.bump();
            let rhs = self.parse_wedge()?;
            expr = Expr::Binary {
                op,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_wedge(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_unary()?;
        while matches!(self.peek_kind(), Some(TokenKind::Wedge)) {
            self.bump();
            let rhs = self.parse_unary()?;
            expr = Expr::Binary {
                op: BinaryOp::Wedge,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_unary(&mut self) -> OghamResult<Expr> {
        let mut ops = Vec::new();
        loop {
            match self.peek_kind() {
                Some(TokenKind::Minus) => {
                    self.bump();
                    ops.push(UnaryOp::Neg);
                }
                Some(TokenKind::Slash) => {
                    self.bump();
                    ops.push(UnaryOp::Inv);
                }
                _ => break,
            }
        }
        let mut expr = self.parse_power()?;
        for op in ops.into_iter().rev() {
            expr = Expr::Unary {
                op,
                expr: Box::new(expr),
            };
        }
        Ok(expr)
    }

    fn parse_power(&mut self) -> OghamResult<Expr> {
        let base = self.parse_appl()?;
        if !matches!(self.peek_kind(), Some(TokenKind::Up)) {
            return Ok(base);
        }
        self.bump();
        let rhs = if matches!(self.peek_kind(), Some(TokenKind::Minus))
            && matches!(self.peek_kind_at(1), Some(TokenKind::Int(_)))
        {
            self.bump();
            let tok = self.bump().expect("peeked int");
            let TokenKind::Int(n) = tok.kind else {
                unreachable!()
            };
            Expr::Unary {
                op: UnaryOp::Neg,
                expr: Box::new(Expr::Int(n)),
            }
        } else {
            self.parse_power()?
        };
        Ok(Expr::Binary {
            op: BinaryOp::Pow,
            lhs: Box::new(base),
            rhs: Box::new(rhs),
        })
    }

    fn parse_appl(&mut self) -> OghamResult<Expr> {
        let mut expr = self.parse_atom()?;
        while matches!(self.peek_kind(), Some(TokenKind::At)) {
            self.bump();
            let rhs = self.parse_appl_arg()?;
            expr = Expr::Binary {
                op: BinaryOp::At,
                lhs: Box::new(expr),
                rhs: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_appl_arg(&mut self) -> OghamResult<Expr> {
        if !matches!(self.peek_kind(), Some(TokenKind::LParen)) {
            return self.parse_atom();
        }
        self.bump();
        let first = self.parse_lambda_or_expression()?;
        if !matches!(self.peek_kind(), Some(TokenKind::Comma)) {
            self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
            return Ok(first);
        }
        let mut items = vec![first];
        while matches!(self.peek_kind(), Some(TokenKind::Comma)) {
            self.bump();
            items.push(self.parse_lambda_or_expression()?);
        }
        self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
        Ok(Expr::Tuple(items))
    }

    fn parse_atom(&mut self) -> OghamResult<Expr> {
        let tok = self.bump().ok_or_else(|| {
            OghamError::new(OghamErrorKind::Parse, Span::point(0), "expected atom")
        })?;
        match tok.kind {
            TokenKind::Int(n) => Ok(Expr::Int(n)),
            TokenKind::True => Ok(Expr::Bool(true)),
            TokenKind::False => Ok(Expr::Bool(false)),
            TokenKind::Star => self.parse_star(),
            TokenKind::Omega => Ok(Expr::Omega),
            TokenKind::Blade(i) => Ok(Expr::Blade(i)),
            TokenKind::Ident(name) => {
                if matches!(self.peek_kind(), Some(TokenKind::LParen)) {
                    self.bump();
                    let mut args = Vec::new();
                    if !matches!(self.peek_kind(), Some(TokenKind::RParen)) {
                        loop {
                            args.push(self.parse_expression()?);
                            if !matches!(self.peek_kind(), Some(TokenKind::Comma)) {
                                break;
                            }
                            self.bump();
                        }
                    }
                    self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
                    Ok(Expr::Call { name, args })
                } else {
                    Ok(Expr::Ident(name))
                }
            }
            TokenKind::Bang => {
                let expr = if matches!(self.peek_kind(), Some(TokenKind::LParen)) {
                    self.bump();
                    let expr = self.parse_expression()?;
                    self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
                    expr
                } else {
                    self.parse_atom()?
                };
                Ok(Expr::Factorial(Box::new(expr)))
            }
            TokenKind::LParen => {
                let expr = self.parse_lambda_or_expression()?;
                self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
                Ok(expr)
            }
            TokenKind::LBracket => {
                let mut items = Vec::new();
                if !matches!(self.peek_kind(), Some(TokenKind::RBracket)) {
                    loop {
                        items.push(self.parse_lambda_or_expression()?);
                        if !matches!(self.peek_kind(), Some(TokenKind::Comma)) {
                            break;
                        }
                        self.bump();
                    }
                }
                self.expect(|k| matches!(k, TokenKind::RBracket), "`]`")?;
                Ok(Expr::Vector(items))
            }
            _ => Err(OghamError::new(
                OghamErrorKind::Parse,
                tok.span,
                "expected atom",
            )),
        }
    }

    fn parse_star(&mut self) -> OghamResult<Expr> {
        match self.peek_kind() {
            Some(TokenKind::Int(n)) => {
                let n = *n;
                self.bump();
                Ok(Expr::Star(StarLiteral::Finite(n)))
            }
            Some(TokenKind::Omega) => {
                self.bump();
                Ok(Expr::Star(StarLiteral::Cnf(Ordinal::omega())))
            }
            Some(TokenKind::LParen) => {
                self.bump();
                let cnf = self.parse_cnf()?;
                self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
                Ok(Expr::Star(StarLiteral::Cnf(cnf)))
            }
            _ => Ok(Expr::Star(StarLiteral::Finite(1))),
        }
    }

    fn parse_cnf(&mut self) -> OghamResult<Ordinal> {
        let mut terms = Vec::<(Ordinal, u128)>::new();
        loop {
            terms.push(self.parse_cnf_term()?);
            if !matches!(self.peek_kind(), Some(TokenKind::Plus)) {
                break;
            }
            self.bump();
        }
        for pair in terms.windows(2) {
            if pair[0].0.cmp(&pair[1].0) != std::cmp::Ordering::Greater {
                return Err(OghamError::new(
                    OghamErrorKind::CnfOrder,
                    self.span(),
                    "CNF exponents must be strictly descending",
                )
                .with_hint("CNF indices are structural: write `*(ω + 1)`, not `*(1 + ω)`"));
            }
        }
        let mut out = Ordinal::from_u128(0);
        for (exp, coeff) in terms {
            let term = if exp.is_zero() {
                Ordinal::from_u128(coeff)
            } else {
                Ordinal::monomial(exp, coeff)
            };
            out = out.nim_add(&term);
        }
        Ok(out)
    }

    fn parse_cnf_term(&mut self) -> OghamResult<(Ordinal, u128)> {
        match self.bump() {
            Some(Token {
                kind: TokenKind::Int(n),
                ..
            }) => Ok((Ordinal::from_u128(0), n)),
            Some(Token {
                kind: TokenKind::Omega,
                ..
            }) => {
                let exp = if matches!(self.peek_kind(), Some(TokenKind::Up)) {
                    self.bump();
                    self.parse_cnf_exp()?
                } else {
                    Ordinal::from_u128(1)
                };
                let coeff = if matches!(self.peek_kind(), Some(TokenKind::Dot)) {
                    self.bump();
                    match self.bump() {
                        Some(Token {
                            kind: TokenKind::Int(n),
                            ..
                        }) => n,
                        Some(tok) => {
                            return Err(OghamError::new(
                                OghamErrorKind::Parse,
                                tok.span,
                                "expected finite CNF coefficient",
                            ));
                        }
                        None => {
                            return Err(OghamError::new(
                                OghamErrorKind::Parse,
                                Span::point(0),
                                "expected finite CNF coefficient",
                            ));
                        }
                    }
                } else {
                    1
                };
                Ok((exp, coeff))
            }
            Some(tok) => Err(OghamError::new(
                OghamErrorKind::Parse,
                tok.span,
                "expected CNF term",
            )),
            None => Err(OghamError::new(
                OghamErrorKind::Parse,
                Span::point(0),
                "expected CNF term",
            )),
        }
    }

    fn parse_cnf_exp(&mut self) -> OghamResult<Ordinal> {
        match self.bump() {
            Some(Token {
                kind: TokenKind::Int(n),
                ..
            }) => Ok(Ordinal::from_u128(n)),
            Some(Token {
                kind: TokenKind::Omega,
                ..
            }) => {
                if matches!(self.peek_kind(), Some(TokenKind::Up)) {
                    self.bump();
                    let exp = self.parse_cnf_exp()?;
                    Ok(Ordinal::omega_pow(exp))
                } else {
                    Ok(Ordinal::omega())
                }
            }
            Some(Token {
                kind: TokenKind::LParen,
                ..
            }) => {
                let cnf = self.parse_cnf()?;
                self.expect(|k| matches!(k, TokenKind::RParen), "`)`")?;
                Ok(cnf)
            }
            Some(tok) => Err(OghamError::new(
                OghamErrorKind::Parse,
                tok.span,
                "expected CNF exponent",
            )),
            None => Err(OghamError::new(
                OghamErrorKind::Parse,
                Span::point(0),
                "expected CNF exponent",
            )),
        }
    }
}
