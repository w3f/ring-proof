use std::path::Path;

fn main() {
    let contracts_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("contracts");
    let config = foundry_config::Config::load_with_root(&contracts_dir)
        .unwrap()
        .canonic();
    let project = config.project().unwrap();
    let _output = project.compile().unwrap();
}
