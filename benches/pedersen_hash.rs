#![feature(test)]

extern crate fil_sapling_crypto;
extern crate paired;
extern crate rand;
extern crate test;

use fil_sapling_crypto::jubjub::JubjubBls12;
use fil_sapling_crypto::pedersen_hash::{pedersen_hash, Personalization};
use paired::bls12_381::Bls12;
use rand::{thread_rng, Rng};

fn bench_pedersen_hash_aux(b: &mut test::Bencher, params: JubjubBls12) {
    let rng = &mut thread_rng();
    let bits = (0..510).map(|_| rng.gen()).collect::<Vec<bool>>();
    let personalization = Personalization::MerkleTree(31);

    b.iter(|| pedersen_hash::<Bls12, _>(personalization, bits.clone(), &params));
}

#[bench]
fn bench_pedersen_hash(b: &mut test::Bencher) {
    let params = JubjubBls12::new();
    bench_pedersen_hash_aux(b, params);
}

#[bench]
fn bench_pedersen_hash_16(b: &mut test::Bencher) {
    let params = JubjubBls12::new_with_window_size(16);
    bench_pedersen_hash_aux(b, params);
}

#[bench]
fn bench_pedersen_hash_17(b: &mut test::Bencher) {
    let params = JubjubBls12::new_with_window_size(17);
    bench_pedersen_hash_aux(b, params);
}

#[bench]
fn bench_pedersen_hash_18(b: &mut test::Bencher) {
    let params = JubjubBls12::new_with_window_size(18);
    bench_pedersen_hash_aux(b, params);
}

#[bench]
fn bench_pedersen_hash_19(b: &mut test::Bencher) {
    let params = JubjubBls12::new_with_window_size(19);
    bench_pedersen_hash_aux(b, params);
}
