use std::fs;
use std::process::Command;

fn pseudo_sequence(length: usize, mut state: u32) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            BASES[(state >> 30) as usize]
        })
        .collect()
}

#[test]
fn cli_maps_fastq_through_the_ordered_worker_pool() {
    let root = std::env::temp_dir().join(format!(
        "rs-lra-cli-smoke-{}-{}",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    fs::create_dir_all(&root).unwrap();

    let reference_sequence = pseudo_sequence(5_000, 91);
    let read_sequence = &reference_sequence[100..2_700];
    let reference_path = root.join("reference.fa");
    let reads_path = root.join("reads.fq");
    let output_path = root.join("output.sam");
    let reference_text = String::from_utf8(reference_sequence.clone()).unwrap();
    fs::write(&reference_path, format!(">chr0\n{reference_text}\n")).unwrap();
    fs::write(
        &reads_path,
        format!(
            "@read0\n{}\n+\n{}\n",
            String::from_utf8(read_sequence.to_vec()).unwrap(),
            "!".repeat(read_sequence.len())
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_rs-lra"))
        .args([
            "--reference",
            reference_path.to_str().unwrap(),
            "--reads",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            "--workers",
            "2",
            "--chunk-size",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let sam = fs::read_to_string(&output_path).unwrap();
    let record = sam
        .lines()
        .find(|line| line.starts_with("read0\t"))
        .unwrap();
    let fields: Vec<_> = record.split('\t').collect();
    assert_eq!(fields[0], "read0");
    assert_eq!(fields[2], "chr0");
    assert_eq!(fields[3], "101");
    assert_eq!(fields[5], "2600M");
    assert_eq!(fields[10], "!".repeat(read_sequence.len()));

    fs::remove_dir_all(root).unwrap();
}
