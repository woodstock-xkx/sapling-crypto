use ff::{Field, PrimeField, PrimeFieldRepr};
use jubjub::*;

#[derive(Copy, Clone)]
pub enum Personalization {
    NoteCommitment,
    MerkleTree(usize),
}

impl Personalization {
    pub fn get_bits(&self) -> Vec<bool> {
        match *self {
            Personalization::NoteCommitment => vec![true, true, true, true, true, true],
            Personalization::MerkleTree(num) => {
                assert!(num < 63);

                (0..6).map(|i| (num >> i) & 1 == 1).collect()
            }
        }
    }
}

pub fn pedersen_hash<E, I>(
    personalization: Personalization,
    bits: I,
    params: &E::Params,
) -> edwards::Point<E, PrimeOrder>
where
    I: IntoIterator<Item = bool>,
    E: JubjubEngine,
{
    let mut bits = personalization
        .get_bits()
        .into_iter()
        .chain(bits.into_iter());

    let mut result = edwards::Point::zero();
    let mut generators = params.pedersen_hash_exp_table().iter();

    let n_groups = JubjubBls12::pedersen_scalar_n();
    let bits_per_iteration = n_groups * 3;

    loop {
        let mut simple_scalar_table = params.pedersen_hash_scalar_table().iter();
        let mut scalar_table = params.pedersen_hash_scalar_n_table().iter();
        let mut acc = E::Fs::zero();
        let mut chunks_remaining = params.pedersen_hash_chunks_per_generator();
        let mut encountered_bits = false;
        let mut stashed_acc = E::Fs::zero();
        let mut stashed_bits = Vec::new(); // FIXME: Reuse allocation?

        let stashed_chunks_remaining = chunks_remaining;
        let mut incomplete_final_bits = false;
        'outer: while let Some(a) = bits.next() {
            stashed_bits.push(a);

            let table = scalar_table.next().expect("not enough scalar chunks");
            encountered_bits = true;

            let mut index = 0;
            if a {
                index += 1
            };
            let mut x = 2;

            let mut bit_count = 1;
            for _ in 0..bits_per_iteration - 1 {
                let unwrapped_bit = bits.next();

                match unwrapped_bit {
                    Some(bit) => {
                        bit_count += 1;
                        if bit_count % 3 == 0 {
                            chunks_remaining -= 1;
                        }
                        if chunks_remaining == 0 {
                            stashed_bits.push(bit);
                            if (bits_per_iteration - bit_count) >= 3 {
                                incomplete_final_bits = true;
                            }
                            break 'outer;
                        }
                    }
                    None => {
                        if (bits_per_iteration - bit_count) >= 3 {
                            incomplete_final_bits = true;
                            break 'outer;
                        }
                    }
                }
                let bit = unwrapped_bit.unwrap_or(false);
                stashed_bits.push(bit);
                if bit {
                    index += x;
                }
                x <<= 1;
            }

            let scalar_for_bits = &table[index];
            acc.add_assign(scalar_for_bits);

            if chunks_remaining == 0 {
                stashed_acc = acc.clone();
                break;
            }
        }

        if incomplete_final_bits {
            acc = stashed_acc;
            chunks_remaining = stashed_chunks_remaining;
            let mut bits = stashed_bits.into_iter();

            while let Some(a) = bits.next() {
                let table = simple_scalar_table
                    .next()
                    .expect("not enough scalar chunks");

                encountered_bits = true;

                let b = bits.next().unwrap_or(false);
                let c = bits.next().unwrap_or(false);

                {
                    let mut index = 0;
                    if a {
                        index += 1
                    };
                    if b {
                        index += 2
                    };
                    if c {
                        index += 4
                    };

                    let tmp = table[index];
                    acc.add_assign(&tmp);
                }

                chunks_remaining -= 1;

                if chunks_remaining == 0 {
                    break;
                }
            }
        }

        if !encountered_bits {
            break;
        }

        let mut table: &[Vec<edwards::Point<E, _>>] =
            &generators.next().expect("we don't have enough generators");
        let window = JubjubBls12::pedersen_hash_exp_window_size();
        let window_mask = (1 << window) - 1;

        let mut acc = acc.into_repr();

        let mut tmp = edwards::Point::zero();

        while !acc.is_zero() {
            let i = (acc.as_ref()[0] & window_mask) as usize;

            tmp = tmp.add(&table[0][i], params);

            acc.shr(window);
            table = &table[1..];
        }

        result = result.add(&tmp, params);
    }

    result
}
