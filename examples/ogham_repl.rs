use ogdoad::ogham::{needs_continuation, OghamSession};
use std::io::{self, Write};

fn main() {
    let mut session = OghamSession::new("integer 0").expect("default ogham world");
    println!("ogham — {}", session.world_summary());
    let stdin = io::stdin();
    let mut pending = String::new();
    loop {
        if pending.is_empty() {
            print!("og> ");
        } else {
            print!(">> ");
        }
        io::stdout().flush().expect("flush prompt");
        let mut line = String::new();
        if stdin.read_line(&mut line).expect("read line") == 0 {
            break;
        }
        let line = line.trim();
        if pending.is_empty() && line.is_empty() {
            continue;
        }
        if pending.is_empty() {
            match line {
                ":quit" | ":q" => break,
                ":help" => {
                    println!(":world <decl>  change world");
                    println!(":env           show bindings");
                    println!(":quit          exit");
                    continue;
                }
                ":env" => {
                    println!("{}", session.world_summary());
                    for binding in session.env_summary() {
                        println!("{binding}");
                    }
                    continue;
                }
                _ => {}
            }
        }
        if pending.is_empty() {
            if let Some(rest) = line.strip_prefix(":world ") {
                match session.set_world(rest) {
                    Ok(()) => println!("{}", session.world_summary()),
                    Err(err) => eprintln!("{err}"),
                }
                continue;
            }
        }
        if !pending.is_empty() {
            pending.push('\n');
        }
        pending.push_str(line);
        match needs_continuation(&pending) {
            Ok(true) => continue,
            Ok(false) => {}
            Err(err) => {
                eprintln!("{err}");
                pending.clear();
                continue;
            }
        }
        match session.eval_line(&pending) {
            Ok(out) => {
                if out.canonical != pending {
                    println!("{}", out.canonical);
                }
                if let Some(value) = out.value {
                    println!("{value}");
                }
            }
            Err(err) => eprintln!("{err}"),
        }
        pending.clear();
    }
}
