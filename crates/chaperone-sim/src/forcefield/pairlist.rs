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
        assert!(i != j, "self-pair ({i}, {j}) is not a valid interaction");
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
        let separation = min_sequence_separation.max(1);
        let mut pairs = PairList::new();
        for i in 0..n {
            for j in (i + separation)..n {
                pairs.push(i as u32, j as u32);
            }
        }
        pairs
    }

    pub fn exclude(&self, excluded_i: &[u32], excluded_j: &[u32]) -> Self {
        assert_eq!(excluded_i.len(), excluded_j.len());

        let key = |i: u32, j: u32| ((i as u64) << 32) | (j as u64);
        let excluded: std::collections::HashSet<u64> = excluded_i
            .iter()
            .zip(excluded_j)
            .map(|(&i, &j)| key(i.min(j), i.max(j)))
            .collect();

        let mut kept = PairList::new();
        for p in 0..self.len() {
            let (i, j) = (self.i[p], self.j[p]);
            if !excluded.contains(&key(i.min(j), i.max(j))) {
                kept.push(i, j);
            }
        }
        kept
    }

    pub fn validate(&self, n: usize) {
        for p in 0..self.len() {
            let (i, j) = (self.i[p] as usize, self.j[p] as usize);
            assert!(
                i < n && j < n,
                "pair ({i}, {j}) out of range for system of {n} atoms"
            );
        }
    }
}

impl Default for PairList {
    fn default() -> Self {
        Self::new()
    }
}
