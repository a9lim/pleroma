use proptest::test_runner::Config as ProptestConfig;

pub fn proptest_config(default_cases: u32) -> ProptestConfig {
    let cases = std::env::var("OGDOAD_PROPTEST_CASES")
        .or_else(|_| std::env::var("PROPTEST_CASES"))
        .ok()
        .and_then(|raw| raw.parse::<u32>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default_cases);
    let mut config = ProptestConfig::with_cases(cases);
    config.failure_persistence = None;
    config
}
