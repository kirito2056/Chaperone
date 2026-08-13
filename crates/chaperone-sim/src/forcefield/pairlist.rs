pub struct PairList {
    pub i: Vec<u32>,
    pub j: Vec<u32>,
}

impl PairList {
    pub fn new() -> Self {
        PairList {
            i: Vec::new(),
            j: Vec::new(),
        }
    }

    pub fn push(&mut self, i: u32, j: u32) {
        self.i.push(i);
        self.j.push(j);
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn all_pairs(n: usize, min_sequence_separation: usize) -> Self {
        let mut pairs = PairList::new();
        for i in 0..n {
            for j in (i + min_sequence_separation)..n {
                pairs.push(i as u32, j as u32);
            }
        }
        pairs
    }
}

impl Default for PairList {
    fn default() -> Self {
        Self::new()
    }
}
