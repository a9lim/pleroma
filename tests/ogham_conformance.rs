use ogdoad::ogham::{EvalLine, OghamError, OghamErrorKind, OghamSession};

#[derive(Debug)]
enum Outcome {
    Ok(EvalLine),
    Err(OghamError),
}

#[test]
fn ogham_conformance_corpus() {
    let corpus = include_str!("../docs/ogham/conformance.txt");
    let mut session: Option<OghamSession> = None;
    let mut pending: Option<(usize, String, Outcome)> = None;
    let lines = corpus.lines().collect::<Vec<_>>();
    let mut idx = 0usize;
    while idx < lines.len() {
        let raw = lines[idx];
        let line_no = idx + 1;
        let line = raw.trim();
        idx += 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(decl) = line.strip_prefix("@world ") {
            finish_pending(&mut pending);
            session = Some(
                OghamSession::new(decl)
                    .unwrap_or_else(|err| panic!("line {line_no}: world failed: {err}")),
            );
            continue;
        }
        if let Some(input) = line.strip_prefix("> ") {
            finish_pending(&mut pending);
            let mut input = input.to_string();
            while idx < lines.len() {
                let cont = lines[idx].trim();
                if let Some(rest) = cont.strip_prefix(">> ") {
                    input.push('\n');
                    input.push_str(rest);
                    idx += 1;
                } else {
                    break;
                }
            }
            let sess = session
                .as_mut()
                .unwrap_or_else(|| panic!("line {line_no}: statement before @world"));
            let outcome = match sess.eval_line(&input) {
                Ok(line) => Outcome::Ok(line),
                Err(err) => Outcome::Err(err),
            };
            pending = Some((line_no, input, outcome));
            continue;
        }
        if line.starts_with(">> ") {
            panic!("line {line_no}: continuation without input");
        }
        if let Some(expected) = line.strip_prefix("~ ") {
            let Some((input_line, input, outcome)) = pending.as_ref() else {
                panic!("line {line_no}: canonical expectation without input");
            };
            match outcome {
                Outcome::Ok(out) => assert_eq!(
                    out.canonical, expected,
                    "line {input_line}: canonical echo for `{input}`"
                ),
                Outcome::Err(err) => {
                    panic!("line {input_line}: expected canonical echo but got {err}")
                }
            }
            continue;
        }
        if let Some(expected) = line.strip_prefix("= ") {
            let Some((input_line, input, outcome)) = pending.take() else {
                panic!("line {line_no}: value expectation without input");
            };
            match outcome {
                Outcome::Ok(out) => assert_eq!(
                    out.value.as_deref(),
                    Some(expected),
                    "line {input_line}: value for `{input}`"
                ),
                Outcome::Err(err) => panic!("line {input_line}: expected value but got {err}"),
            }
            continue;
        }
        if let Some(expected) = line.strip_prefix("! ") {
            let Some((input_line, input, outcome)) = pending.take() else {
                panic!("line {line_no}: error expectation without input");
            };
            let (kind, needle) = expected
                .split_once(':')
                .map_or((expected, ""), |(kind, needle)| {
                    (kind.trim(), needle.trim())
                });
            match outcome {
                Outcome::Err(err) => {
                    assert_eq!(
                        err.kind.code(),
                        kind,
                        "line {input_line}: error kind for `{input}`"
                    );
                    let haystack = format!("{err}");
                    assert!(
                        needle.is_empty() || haystack.contains(needle),
                        "line {input_line}: expected error substring `{needle}` in `{haystack}`"
                    );
                }
                Outcome::Ok(out) => panic!("line {input_line}: expected error but got {out:?}"),
            }
            continue;
        }
        panic!("line {line_no}: unknown corpus directive `{line}`");
    }
    finish_pending(&mut pending);
}

fn finish_pending(pending: &mut Option<(usize, String, Outcome)>) {
    let Some((line_no, input, outcome)) = pending.take() else {
        return;
    };
    match outcome {
        Outcome::Ok(out) => {
            assert!(
                out.value.is_none(),
                "line {line_no}: `{input}` produced unexpected value {:?}",
                out.value
            );
        }
        Outcome::Err(err) => panic!("line {line_no}: `{input}` failed unexpectedly: {err}"),
    }
}

#[test]
fn error_kind_codes_are_stable() {
    assert_eq!(OghamErrorKind::BareInt.code(), "E_BareInt");
    assert_eq!(OghamErrorKind::KummerEscape.code(), "E_KummerEscape");
}
