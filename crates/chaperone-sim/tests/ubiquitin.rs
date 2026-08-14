use chaperone_sim::forcefield::{angle, bond, dihedral, native};
use chaperone_sim::model;
use chaperone_sim::scenario;
use chaperone_sim::system::Real;

const PATH: &str = "../../data/pdb/1UBQ.pdb";

fn load() -> chaperone_pdb::Structure {
    let text = std::fs::read_to_string(PATH).unwrap_or_else(|e| {
        panic!("{PATH}: {e}\nfetch it with: curl -o data/pdb/1UBQ.pdb https://files.rcsb.org/download/1UBQ.pdb")
    });
    chaperone_pdb::parse(&text, None).expect("1UBQ must parse cleanly")
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn ubiquitin_structure_matches_reference_measurements() {
    let s = load();

    assert_eq!(s.chain, b'A');
    assert_eq!(s.len(), 76, "ubiquitin has 76 residues");
    assert_eq!(s.residues[0].seq, 1);
    assert_eq!(s.residues[75].seq, 76);

    let bonds = s.peptide_bond_lengths();
    let lo = bonds.iter().cloned().fold(Real::INFINITY, Real::min);
    let hi = bonds.iter().cloned().fold(Real::NEG_INFINITY, Real::max);
    assert!(
        (3.75..=3.86).contains(&lo) && (3.75..=3.86).contains(&hi),
        "peptide CA-CA range [{lo:.3}, {hi:.3}] left [3.75, 3.86]"
    );

    let rg = s.radius_of_gyration();
    assert!(
        (rg - 11.49).abs() < 0.02,
        "radius of gyration {rg:.3}, expected 11.49"
    );
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn ubiquitin_contact_maps_match_reference_counts() {
    let s = load();

    let heavy = s.heavy_atom_contacts(model::CONTACT_CUTOFF, model::MIN_SEQUENCE_SEPARATION);
    let ca = s.ca_contacts(7.5, model::MIN_SEQUENCE_SEPARATION);

    assert_eq!(
        heavy.len(),
        180,
        "heavy-atom 4.5 A with sep >= 3; the ~110-130 figure in version_02 came from \
         occlusion-filtered maps (CSU / shadow), not a plain distance cutoff"
    );
    assert_eq!(ca.len(), 155, "CA-CA 7.5 A with sep >= 3");

    let heavy_set: std::collections::HashSet<_> = heavy.iter().copied().collect();
    let ca_set: std::collections::HashSet<_> = ca.iter().copied().collect();
    let shared = heavy_set.intersection(&ca_set).count();

    assert_eq!(
        shared, 120,
        "the two criteria are different maps, not a rescaling"
    );
    assert_eq!(
        heavy.len() - shared,
        60,
        "contacts only the heavy-atom map sees"
    );
    assert_eq!(ca.len() - shared, 35, "contacts only the CA map sees");

    let model_co = chaperone_pdb::relative_contact_order(&heavy, s.len()).unwrap();
    assert!(
        (model_co - 0.281).abs() < 5e-3,
        "model contact order {model_co:.4}, expected 0.281"
    );

    let plaxco = s.heavy_atom_contacts(model::PLAXCO_CUTOFF, 1);
    let plaxco_co = chaperone_pdb::relative_contact_order(&plaxco, s.len()).unwrap();
    assert!(
        (plaxco_co - 0.196).abs() < 5e-3,
        "Plaxco contact order {plaxco_co:.4}, expected 0.196 (literature ~0.19)"
    );
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn ubiquitin_go_model_is_force_free_at_the_native_state() {
    let s = load();
    let (mut sys, ff) = model::go_model(
        &s,
        scenario::BOND_K,
        scenario::ANGLE_K,
        scenario::K_PHI1,
        scenario::K_PHI3,
        scenario::EPS,
        scenario::SIGMA,
    );

    assert_eq!(ff.bonds.len(), 75);
    assert_eq!(ff.angles.len(), 74);
    assert_eq!(ff.dihedrals.len(), 73);
    assert_eq!(ff.native.len(), 180);
    assert_eq!(ff.repulsion_pairs.len(), 2701 - 180);

    let lo = ff
        .native
        .sigma
        .iter()
        .cloned()
        .fold(Real::INFINITY, Real::min);
    let hi = ff
        .native
        .sigma
        .iter()
        .cloned()
        .fold(Real::NEG_INFINITY, Real::max);
    let mean = ff.native.sigma.iter().sum::<Real>() / ff.native.len() as Real;
    assert!(
        (lo - 3.93).abs() < 0.02 && (hi - 11.12).abs() < 0.02 && (mean - 6.82).abs() < 0.02,
        "sigma range [{lo:.3}, {hi:.3}] mean {mean:.3}, expected [3.93, 11.12] mean 6.82"
    );

    sys.clear_forces();
    bond::accumulate(&mut sys, &ff.bonds, ff.bond_k);
    assert!(sys.max_force() < 1e-9, "bond force at the native state");

    sys.clear_forces();
    angle::accumulate(&mut sys, &ff.angles, ff.angle_k);
    assert!(sys.max_force() < 1e-9, "angle force at the native state");

    sys.clear_forces();
    dihedral::accumulate(&mut sys, &ff.dihedrals, ff.k_phi1, ff.k_phi3);
    assert!(sys.max_force() < 1e-9, "dihedral force at the native state");

    sys.clear_forces();
    native::accumulate(&mut sys, &ff.native, ff.eps);
    assert!(
        sys.max_force() < 1e-9,
        "native contact force at the native state; sigma_ij must be the CA-CA distance"
    );
}
