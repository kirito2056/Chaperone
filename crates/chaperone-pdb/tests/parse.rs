use chaperone_pdb::{parse, ParseError};

#[allow(clippy::too_many_arguments)]
fn atom(
    record: &str,
    serial: i32,
    name: &str,
    alt: char,
    res: &str,
    chain: char,
    seq: i32,
    x: f64,
    y: f64,
    z: f64,
    element: &str,
) -> String {
    format!(
        "{record:<6}{serial:>5} {name:<4}{alt}{res:>3} {chain}{seq:>4}    \
         {x:8.3}{y:8.3}{z:8.3}  1.00  0.00          {element:>2}"
    )
}

const DX: f64 = 3.7;
const DY: f64 = 0.866;

fn zigzag(i: usize) -> (f64, f64, f64) {
    (
        DX * i as f64,
        if i.is_multiple_of(2) { 0.0 } else { DY },
        0.0,
    )
}

fn chain_lines(n: usize) -> Vec<String> {
    let mut lines = Vec::new();
    for i in 0..n {
        let (x, y, z) = zigzag(i);
        let seq = i as i32 + 1;
        lines.push(atom(
            "ATOM",
            2 * seq - 1,
            " N  ",
            ' ',
            "ALA",
            'A',
            seq,
            x,
            y + 1.0,
            z,
            "N",
        ));
        lines.push(atom(
            "ATOM",
            2 * seq,
            " CA ",
            ' ',
            "ALA",
            'A',
            seq,
            x,
            y,
            z,
            "C",
        ));
    }
    lines
}

fn chain_text(n: usize) -> String {
    chain_lines(n).join("\n")
}

#[test]
fn parses_a_minimal_chain() {
    let s = parse(&chain_text(5), None).expect("valid chain");
    assert_eq!(s.chain, b'A');
    assert_eq!(s.len(), 5);
    assert_eq!(s.residues[0].name_str(), "ALA");
    assert_eq!(s.residues[4].seq, 5);
    for d in s.peptide_bond_lengths() {
        assert!((d - 3.8).abs() < 1e-3, "peptide bond {d} is not 3.8");
    }
    for r in &s.residues {
        assert_eq!(r.heavy.len(), 2, "N and CA should both be kept");
    }
}

#[test]
fn record_name_is_read_from_columns_not_substring() {
    let mut lines = vec![
        "REMARK 200  X-RAY GENERATOR MODEL              : ROTATING ANODE".to_string(),
        "REMARK 280 MODEL BUILT WITH O".to_string(),
    ];
    lines.extend(chain_lines(4));

    let s = parse(&lines.join("\n"), None)
        .expect("REMARK text containing MODEL must not terminate parsing");
    assert_eq!(s.len(), 4);
}

#[test]
fn coordinates_that_run_together_are_split_by_column() {
    let mut lines = Vec::new();
    for i in 0..4 {
        let seq = i + 1;
        let x = -100.0 - 3.7 * i as f64;
        let y = -100.0 + if i % 2 == 0 { 0.0 } else { DY };
        lines.push(atom(
            "ATOM", seq, " CA ", ' ', "ALA", 'A', seq, x, y, -100.0, "C",
        ));
    }
    let text = lines.join("\n");
    assert!(
        text.contains("-100.000-100.000"),
        "fixture must actually produce touching coordinate fields"
    );

    let s = parse(&text, None).expect("fixed-column parsing must survive touching fields");
    assert_eq!(s.len(), 4);
    assert_eq!(s.residues[0].ca, [-100.0, -100.0, -100.0]);
    assert!((s.residues[1].ca[0] + 103.7).abs() < 1e-9);
}

#[test]
fn alternate_locations_do_not_duplicate_residues() {
    let mut lines = chain_lines(4);
    let (x, y, z) = zigzag(1);
    lines.insert(
        3,
        atom("ATOM", 100, " CA ", 'A', "ALA", 'A', 2, x, y, z, "C"),
    );
    lines.insert(
        4,
        atom("ATOM", 101, " CA ", 'B', "ALA", 'A', 2, x + 9.0, y, z, "C"),
    );

    let s = parse(&lines.join("\n"), None).expect("altLoc B must be dropped");
    assert_eq!(s.len(), 4, "altLoc must not create an extra residue");
    assert!(
        (s.residues[1].ca[0] - x).abs() < 1e-9,
        "conformer A should win, got {:?}",
        s.residues[1].ca
    );
}

#[test]
fn water_and_ligands_are_excluded() {
    let mut lines = chain_lines(4);
    lines.push(atom(
        "HETATM", 900, " O  ", ' ', "HOH", 'A', 200, 50.0, 50.0, 50.0, "O",
    ));
    lines.push(atom(
        "HETATM", 901, "CA  ", ' ', " CA", 'A', 201, 60.0, 60.0, 60.0, "CA",
    ));

    let s = parse(&lines.join("\n"), None).expect("HETATM must be skipped");
    assert_eq!(
        s.len(),
        4,
        "water and the calcium ion must not become residues"
    );
}

#[test]
fn only_the_first_model_is_read() {
    let mut lines = vec!["MODEL        1".to_string()];
    lines.extend(chain_lines(4));
    lines.push("ENDMDL".to_string());
    lines.push("MODEL        2".to_string());
    lines.extend(chain_lines(4));
    lines.push("ENDMDL".to_string());

    let s = parse(&lines.join("\n"), None).expect("first model only");
    assert_eq!(s.len(), 4, "the second model must not be appended");
}

#[test]
fn hydrogens_are_excluded_by_element_column() {
    let mut lines = chain_lines(3);
    let (x, y, z) = zigzag(0);
    lines.push(atom(
        "ATOM",
        500,
        " HA ",
        ' ',
        "ALA",
        'A',
        1,
        x,
        y,
        z + 1.0,
        "H",
    ));
    lines.push(atom(
        "ATOM",
        501,
        "1HG2",
        ' ',
        "ALA",
        'A',
        1,
        x,
        y,
        z + 2.0,
        "H",
    ));

    let s = parse(&lines.join("\n"), None).expect("hydrogens skipped");
    assert_eq!(s.residues[0].heavy.len(), 2, "no hydrogen may reach heavy");
}

#[test]
fn hydrogens_are_inferred_when_the_element_column_is_absent() {
    let mut lines: Vec<String> = chain_lines(3).iter().map(|l| l[..54].to_string()).collect();
    let (x, y, z) = zigzag(0);
    let mut h = atom("ATOM", 500, " HA ", ' ', "ALA", 'A', 1, x, y, z + 1.0, "H");
    h.truncate(54);
    lines.push(h);

    let s = parse(&lines.join("\n"), None).expect("short lines must still parse");
    assert_eq!(
        s.residues[0].heavy.len(),
        2,
        "hydrogen must be inferred from the atom name when column 77-78 is missing"
    );
}

#[test]
fn a_chain_break_is_rejected() {
    let mut lines = chain_lines(3);
    lines.extend(vec![atom(
        "ATOM", 99, " CA ", ' ', "ALA", 'A', 40, 40.0, 0.0, 0.0, "C",
    )]);

    match parse(&lines.join("\n"), None) {
        Err(ParseError::ChainBreak { distance, .. }) => {
            assert!(
                distance > 4.5,
                "distance {distance} should exceed the cutoff"
            )
        }
        other => panic!("expected a chain break error, got {other:?}"),
    }
}

#[test]
fn implausibly_close_residues_are_rejected() {
    let mut lines = chain_lines(2);
    lines.push(atom(
        "ATOM",
        99,
        " CA ",
        ' ',
        "ALA",
        'A',
        3,
        DX,
        DY + 0.5,
        0.0,
        "C",
    ));

    match parse(&lines.join("\n"), None) {
        Err(ParseError::ImplausiblyClose { distance, .. }) => {
            assert!(
                distance < 2.5,
                "distance {distance} should be below the floor"
            )
        }
        other => panic!("expected an implausibly-close error, got {other:?}"),
    }
}

#[test]
fn a_residue_without_ca_is_rejected() {
    let mut lines = chain_lines(3);
    lines.retain(|l| &l[12..16] != " CA " || &l[22..26] != "   2");

    match parse(&lines.join("\n"), None) {
        Err(ParseError::MissingCa { seq, .. }) => assert_eq!(seq, 2),
        other => panic!("expected a missing CA error, got {other:?}"),
    }
}

#[test]
fn a_second_chain_is_ignored_unless_requested() {
    let mut lines = chain_lines(4);
    for i in 0..4 {
        let (x, y, z) = zigzag(i);
        let seq = i as i32 + 1;
        lines.push(atom(
            "ATOM",
            800 + seq,
            " CA ",
            ' ',
            "GLY",
            'B',
            seq,
            x + 30.0,
            y,
            z,
            "C",
        ));
    }
    let text = lines.join("\n");

    let a = parse(&text, None).expect("first chain by default");
    assert_eq!(a.chain, b'A');
    assert_eq!(a.len(), 4);

    let b = parse(&text, Some(b'B')).expect("chain B on request");
    assert_eq!(b.chain, b'B');
    assert_eq!(b.len(), 4);
    assert_eq!(b.residues[0].name_str(), "GLY");
}

#[test]
fn an_empty_file_is_rejected() {
    assert!(matches!(parse("", None), Err(ParseError::NoResidues)));
    assert!(matches!(
        parse("HEADER    NOTHING HERE\nEND", None),
        Err(ParseError::NoResidues)
    ));
}

#[test]
fn contact_order_counts_sequence_separation() {
    let contacts = [(0usize, 3usize), (1, 5)];
    let co = chaperone_pdb::relative_contact_order(&contacts, 10).unwrap();
    assert!((co - (3.0 + 4.0) / (2.0 * 10.0)).abs() < 1e-12);
    assert!(chaperone_pdb::relative_contact_order(&[], 10).is_none());
}
