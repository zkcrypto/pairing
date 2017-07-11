use super::*;
use ::*;

fn test_vectors<G: CurveProjective, E: EncodedPoint<Affine=G::Affine>>(expected: &[u8])
{
    let mut e = G::one();

    let mut v = vec![];
    {
        let mut expected = expected;
        for _ in 0..1000 {
            let e_affine = e.to_affine();
            let encoded = E::from_affine(e_affine).unwrap();
            v.extend_from_slice(encoded.as_ref());

            let mut decoded = E::empty();
            decoded.as_mut().copy_from_slice(&expected[0..E::size()]);
            expected = &expected[E::size()..];
            let decoded = decoded.into_affine().unwrap();
            assert_eq!(e_affine, decoded);

            e.add_assign(&G::one());
        }
    }

    assert_eq!(&v[..], expected);
}

#[test]
fn test_g1_uncompressed_vectors() {
    test_vectors::<G1, G1Uncompressed>(include_bytes!("g1_uncompressed_test_vectors.dat"));
}

#[test]
fn test_g1_compressed_vectors() {
    test_vectors::<G1, G1Compressed>(include_bytes!("g1_compressed_test_vectors.dat"));
}

#[test]
fn test_g2_uncompressed_vectors() {
    test_vectors::<G2, G2Uncompressed>(include_bytes!("g2_uncompressed_test_vectors.dat"));
}

#[test]
fn test_g2_compressed_vectors() {
    test_vectors::<G2, G2Compressed>(include_bytes!("g2_compressed_test_vectors.dat"));
}

