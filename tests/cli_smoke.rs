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
    assert!(fields.contains(&"AS:i:5200"));
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

#[test]
fn dual_affine_cli_preserves_long_indels_on_both_strands() {
    let root = std::env::temp_dir().join(format!("rs-lra-dual-indels-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let reference = pseudo_sequence(5_000, 717);
    let mut deletion = reference[100..2700].to_vec();
    deletion.drain(1200..1240);
    let mut insertion = reference[100..2700].to_vec();
    insertion.splice(1200..1200, pseudo_sequence(25, 881));
    let reads = [
        ("del", deletion.clone(), 40, "40D"),
        ("del_rev", reverse_complement(&deletion), 40, "40D"),
        ("ins", insertion.clone(), 25, "25I"),
        ("ins_rev", reverse_complement(&insertion), 25, "25I"),
    ];
    let reference_path = root.join("ref.fa");
    let reads_path = root.join("reads.fq");
    let output_path = root.join("out.sam");
    fs::write(
        &reference_path,
        format!(">chr0\n{}\n", String::from_utf8(reference).unwrap()),
    )
    .unwrap();
    let fastq: String = reads
        .iter()
        .map(|(name, q, _, _)| {
            format!(
                "@{name}\n{}\n+\n{}\n",
                String::from_utf8(q.clone()).unwrap(),
                "I".repeat(q.len())
            )
        })
        .collect();
    fs::write(&reads_path, fastq).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rs-lra"))
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
            "--dual-affine",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let sam = fs::read_to_string(&output_path).unwrap();
    for (name, q, nm, gap) in reads {
        let fields: Vec<_> = sam
            .lines()
            .find(|line| line.starts_with(&format!("{name}\t")))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(fields[2], "chr0");
        assert_eq!(fields[3], "101");
        assert_eq!(fields[1], if name.ends_with("_rev") { "16" } else { "0" });
        assert!(fields[5].contains(gap), "{name}: {}", fields[5]);
        assert!(
            fields.contains(&format!("NM:i:{nm}").as_str()),
            "{name}: {fields:?}"
        );
        assert_eq!(fields[9].len(), q.len());
        // These fixtures contain only one indel and otherwise exact matches.
        let matches = q.len() as i32 - if gap.ends_with('I') { nm } else { 0 };
        let expected_score = 2 * matches - (6 + 2 * nm).min(24 + nm);
        assert!(
            fields.contains(&format!("AS:i:{expected_score}").as_str()),
            "{name}: {fields:?}"
        );
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn secondary_records_preserve_primary_and_calibration_caps_unique_mapq() {
    let root = std::env::temp_dir().join(format!("rs-lra-secondary-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let repeated = String::from_utf8(pseudo_sequence(8000, 937)).unwrap();
    let unique = String::from_utf8(pseudo_sequence(8000, 117)).unwrap();
    fs::write(
        root.join("ref.fa"),
        format!(">a\n{repeated}\n>b\n{repeated}\n>unique\n{unique}\n"),
    )
    .unwrap();
    fs::write(
        root.join("reads.fa"),
        format!(
            ">repeat\n{}\n>unique\n{}\n",
            &repeated[100..5100],
            &unique[100..5100]
        ),
    )
    .unwrap();
    fs::write(
        root.join("caps.tsv"),
        (0..=60)
            .map(|q| format!("{q}\t{}\n", q.min(12)))
            .collect::<String>(),
    )
    .unwrap();
    let run = |secondary: usize, calibrated: bool| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rs-lra"));
        command.args([
            "--reference",
            root.join("ref.fa").to_str().unwrap(),
            "--reads",
            root.join("reads.fa").to_str().unwrap(),
            "--output",
            root.join("out.sam").to_str().unwrap(),
            "--workers",
            "1",
            "--secondary",
            &secondary.to_string(),
        ]);
        if calibrated {
            command.arg("--mapq-calibration").arg(root.join("caps.tsv"));
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        fs::read_to_string(root.join("out.sam"))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let baseline = run(0, false);
    let enabled = run(1, false);
    let primary = |records: &Vec<String>| {
        records
            .iter()
            .filter(|line| line.split('\t').nth(1).unwrap().parse::<u16>().unwrap() & 0x900 == 0)
            .cloned()
            .collect::<Vec<_>>()
    };
    assert_eq!(primary(&baseline), primary(&enabled));
    let alternatives: Vec<_> = enabled
        .iter()
        .filter(|line| line.split('\t').nth(1).unwrap().parse::<u16>().unwrap() & 0x100 != 0)
        .collect();
    assert_eq!(alternatives.len(), 1, "{enabled:?}");
    let fields: Vec<_> = alternatives[0].split('\t').collect();
    assert_eq!(fields[0], "repeat");
    assert_eq!(fields[4], "0");
    assert!(!alternatives[0].contains("SA:Z:"));
    assert!(fields.contains(&"NM:i:0"));
    let calibrated = run(1, true);
    let unique_record = calibrated
        .iter()
        .find(|r| r.starts_with("unique\t"))
        .unwrap();
    assert_eq!(unique_record.split('\t').nth(4), Some("12"));
    fs::remove_dir_all(root).unwrap();
}
