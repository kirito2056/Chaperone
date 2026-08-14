use std::fmt;

pub const PEPTIDE_MIN: f64 = 2.5;
pub const PEPTIDE_MAX: f64 = 4.5;

#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    NoResidues,
    BadNumber {
        line: usize,
        field: &'static str,
    },
    ShortLine {
        line: usize,
        len: usize,
    },
    MissingCa {
        index: usize,
        seq: i32,
    },
    ChainBreak {
        index: usize,
        seq_a: i32,
        seq_b: i32,
        distance: f64,
    },
    ImplausiblyClose {
        index: usize,
        seq_a: i32,
        seq_b: i32,
        distance: f64,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::NoResidues => write!(f, "no protein residues found"),
            ParseError::BadNumber { line, field } => {
                write!(f, "line {line}: could not parse {field}")
            }
            ParseError::ShortLine { line, len } => {
                write!(f, "line {line}: truncated ATOM record ({len} bytes)")
            }
            ParseError::MissingCa { index, seq } => {
                write!(f, "residue {seq} (index {index}) has no CA atom")
            }
            ParseError::ChainBreak {
                index,
                seq_a,
                seq_b,
                distance,
            } => write!(
                f,
                "chain break between residues {seq_a} and {seq_b} (index {index}): \
                 CA-CA distance {distance:.3} A exceeds {PEPTIDE_MAX} A"
            ),
            ParseError::ImplausiblyClose {
                index,
                seq_a,
                seq_b,
                distance,
            } => write!(
                f,
                "residues {seq_a} and {seq_b} (index {index}) are {distance:.3} A apart, \
                 below the {PEPTIDE_MIN} A floor; the coordinates look misparsed"
            ),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Residue {
    pub name: [u8; 3],
    pub seq: i32,
    pub icode: u8,
    pub ca: [f64; 3],
    pub heavy: Vec<[f64; 3]>,
}

impl Residue {
    pub fn name_str(&self) -> &str {
        std::str::from_utf8(&self.name).unwrap_or("???")
    }
}

#[derive(Debug, Clone)]
pub struct Structure {
    pub chain: u8,
    pub residues: Vec<Residue>,
}

fn field(bytes: &[u8], from: usize, to: usize) -> &[u8] {
    let end = to.min(bytes.len());
    if from >= end {
        return &[];
    }
    let slice = &bytes[from..end];
    let start = slice.iter().position(|b| !b.is_ascii_whitespace());
    let stop = slice.iter().rposition(|b| !b.is_ascii_whitespace());
    match (start, stop) {
        (Some(a), Some(b)) => &slice[a..=b],
        _ => &[],
    }
}

fn parse_f64(bytes: &[u8], line: usize, name: &'static str) -> Result<f64, ParseError> {
    std::str::from_utf8(bytes)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or(ParseError::BadNumber { line, field: name })
}

fn parse_i32(bytes: &[u8], line: usize, name: &'static str) -> Result<i32, ParseError> {
    std::str::from_utf8(bytes)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or(ParseError::BadNumber { line, field: name })
}

fn is_hydrogen(raw_name: &[u8], element: &[u8]) -> bool {
    if !element.is_empty() {
        return element == b"H" || element == b"D";
    }
    match raw_name.first() {
        Some(b'H') => true,
        Some(b' ') | Some(b'0'..=b'9') => raw_name.get(1) == Some(&b'H'),
        _ => false,
    }
}

pub fn parse(text: &str, chain: Option<u8>) -> Result<Structure, ParseError> {
    let mut residues: Vec<Residue> = Vec::new();
    let mut selected_chain: Option<u8> = chain;
    let mut current_key: Option<(i32, u8)> = None;

    for (number, line) in text.lines().enumerate() {
        let bytes = line.as_bytes();
        let record = field(bytes, 0, 6);

        if record == b"ENDMDL" {
            break;
        }
        if record != b"ATOM" {
            continue;
        }
        if bytes.len() < 54 {
            return Err(ParseError::ShortLine {
                line: number + 1,
                len: bytes.len(),
            });
        }

        let alt_loc = bytes[16];
        if alt_loc != b' ' && alt_loc != b'A' {
            continue;
        }

        let chain_id = bytes[21];
        match selected_chain {
            None => selected_chain = Some(chain_id),
            Some(want) if want != chain_id => continue,
            Some(_) => {}
        }

        let raw_name = &bytes[12..16];
        let element = field(bytes, 76, 78);
        if is_hydrogen(raw_name, element) {
            continue;
        }

        let seq = parse_i32(field(bytes, 22, 26), number + 1, "resSeq")?;
        let icode = bytes[26];
        let pos = [
            parse_f64(field(bytes, 30, 38), number + 1, "x")?,
            parse_f64(field(bytes, 38, 46), number + 1, "y")?,
            parse_f64(field(bytes, 46, 54), number + 1, "z")?,
        ];

        if current_key != Some((seq, icode)) {
            let mut name = [b' '; 3];
            for (slot, byte) in name.iter_mut().zip(field(bytes, 17, 20)) {
                *slot = *byte;
            }
            residues.push(Residue {
                name,
                seq,
                icode,
                ca: [f64::NAN; 3],
                heavy: Vec::new(),
            });
            current_key = Some((seq, icode));
        }

        let residue = residues.last_mut().expect("residue pushed above");
        residue.heavy.push(pos);
        if raw_name == b" CA " {
            residue.ca = pos;
        }
    }

    if residues.is_empty() {
        return Err(ParseError::NoResidues);
    }

    for (index, residue) in residues.iter().enumerate() {
        if residue.ca.iter().any(|c| c.is_nan()) {
            return Err(ParseError::MissingCa {
                index,
                seq: residue.seq,
            });
        }
    }

    for index in 1..residues.len() {
        let d = distance(residues[index - 1].ca, residues[index].ca);
        if d > PEPTIDE_MAX {
            return Err(ParseError::ChainBreak {
                index,
                seq_a: residues[index - 1].seq,
                seq_b: residues[index].seq,
                distance: d,
            });
        }
        if d < PEPTIDE_MIN {
            return Err(ParseError::ImplausiblyClose {
                index,
                seq_a: residues[index - 1].seq,
                seq_b: residues[index].seq,
                distance: d,
            });
        }
    }

    Ok(Structure {
        chain: selected_chain.expect("residues implies a chain"),
        residues,
    })
}

pub fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

impl Structure {
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    pub fn ca_positions(&self) -> Vec<[f64; 3]> {
        self.residues.iter().map(|r| r.ca).collect()
    }

    pub fn peptide_bond_lengths(&self) -> Vec<f64> {
        (1..self.len())
            .map(|i| distance(self.residues[i - 1].ca, self.residues[i].ca))
            .collect()
    }

    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.len() as f64;
        let mut centre = [0.0; 3];
        for r in &self.residues {
            for (c, v) in centre.iter_mut().zip(r.ca) {
                *c += v;
            }
        }
        for c in centre.iter_mut() {
            *c /= n;
        }

        let mut sum = 0.0;
        for r in &self.residues {
            let d = distance(r.ca, centre);
            sum += d * d;
        }
        (sum / n).sqrt()
    }

    pub fn heavy_atom_contacts(&self, cutoff: f64, min_sep: usize) -> Vec<(usize, usize)> {
        let cutoff_sq = cutoff * cutoff;
        let mut contacts = Vec::new();

        for i in 0..self.len() {
            for j in (i + min_sep.max(1))..self.len() {
                let touching = self.residues[i].heavy.iter().any(|a| {
                    self.residues[j].heavy.iter().any(|b| {
                        let dx = a[0] - b[0];
                        let dy = a[1] - b[1];
                        let dz = a[2] - b[2];
                        dx * dx + dy * dy + dz * dz < cutoff_sq
                    })
                });
                if touching {
                    contacts.push((i, j));
                }
            }
        }

        contacts
    }

    pub fn ca_contacts(&self, cutoff: f64, min_sep: usize) -> Vec<(usize, usize)> {
        let cutoff_sq = cutoff * cutoff;
        let mut contacts = Vec::new();

        for i in 0..self.len() {
            for j in (i + min_sep.max(1))..self.len() {
                let d = distance(self.residues[i].ca, self.residues[j].ca);
                if d * d < cutoff_sq {
                    contacts.push((i, j));
                }
            }
        }

        contacts
    }
}

pub fn relative_contact_order(contacts: &[(usize, usize)], n: usize) -> Option<f64> {
    if contacts.is_empty() || n == 0 {
        return None;
    }
    let total: usize = contacts.iter().map(|(i, j)| j - i).sum();
    Some(total as f64 / (contacts.len() as f64 * n as f64))
}
