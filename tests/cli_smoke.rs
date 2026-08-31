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

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
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

#[test]
fn cli_emits_reverse_complemented_sequence_and_reversed_quality() {
    let root = std::env::temp_dir().join(format!(
        "rs-lra-cli-reverse-{}-{}",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    fs::create_dir_all(&root).unwrap();

    let reference_sequence = pseudo_sequence(5_000, 193);
    let forward_slice = &reference_sequence[100..2_700];
    let read_sequence = reverse_complement(forward_slice);
    let qualities: String = (0..read_sequence.len())
        .map(|index| char::from(b'!' + (index % 40) as u8))
        .collect();
    let reference_path = root.join("reference.fa");
    let reads_path = root.join("reads.fq");
    let output_path = root.join("output.sam");
    fs::write(
        &reference_path,
        format!(
            ">chr0\n{}\n",
            String::from_utf8(reference_sequence.clone()).unwrap()
        ),
    )
    .unwrap();
    fs::write(
        &reads_path,
        format!(
            "@reverse\n{}\n+\n{}\n",
            String::from_utf8(read_sequence).unwrap(),
            qualities
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
            "1",
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
        .find(|line| line.starts_with("reverse\t"))
        .unwrap();
    let fields: Vec<_> = record.split('\t').collect();
    assert_eq!(fields[1], "16");
    assert_eq!(fields[3], "101");
    assert_eq!(
        fields[9],
        String::from_utf8(forward_slice.to_vec()).unwrap()
    );
    let expected_quality: String = qualities.chars().rev().collect();
    assert_eq!(fields[10], expected_quality);

    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_marks_reads_shorter_than_the_anchor_k_as_unmapped() {
    let root = std::env::temp_dir().join(format!(
        "rs-lra-cli-short-{}-{}",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    fs::create_dir_all(&root).unwrap();

    let reference_path = root.join("reference.fa");
    let reads_path = root.join("reads.fq");
    let output_path = root.join("output.sam");
    fs::write(&reference_path, ">chr0\nACGTACGTACGTACGT\n").unwrap();
    fs::write(&reads_path, "@short\nACG\n+\n!\"#\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_rs-lra"))
        .args([
            "--reference",
            reference_path.to_str().unwrap(),
            "--reads",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            "--workers",
            "1",
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
        .find(|line| line.starts_with("short\t"))
        .unwrap();
    let fields: Vec<_> = record.split('\t').collect();
    assert_eq!(fields[1], "4");
    assert_eq!(fields[2], "*");
    assert_eq!(fields[9], "ACG");
    assert_eq!(fields[10], "!\"#");

    fs::remove_dir_all(root).unwrap();
}
