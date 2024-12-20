use bitvec::vec::BitVec;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{BitAnd, BitOr, BitOrAssign, Not, Shl, ShlAssign, Shr, ShrAssign};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct BitPackedKmer<const K: usize, Integer> {
    kmer: Integer,
}

#[derive(Debug, Clone, Eq, PartialEq, PartialOrd, Ord)]
pub struct BitPackedVectorKmer {
    kmer: BitVec,
}

pub trait Kmer: Ord + Sized + Clone {
    fn reverse_complement(&self) -> Self;

    fn canonical(&self) -> Self {
        let reverse_complement = self.reverse_complement();
        if &reverse_complement < self {
            reverse_complement
        } else {
            self.clone()
        }
    }

    fn predecessor(&self, character: u8) -> Self;

    fn successor(&self, character: u8) -> Self;
}

impl<
        const K: usize,
        Integer: Default + Shl<i32, Output = Integer> + BitOr<Integer, Output = Integer> + From<u8>,
    > FromIterator<u8> for BitPackedKmer<K, Integer>
{
    fn from_iter<Iter: IntoIterator<Item = u8>>(iter: Iter) -> Self {
        assert!(2 * K <= 8 * std::mem::size_of::<Integer>());

        let iter = iter.into_iter();
        let size = iter.size_hint();
        assert_eq!(Some(size.0), size.1);
        assert_eq!(K, size.0);

        BitPackedKmer {
            kmer: iter.fold(Integer::default(), |result, character| {
                let bits = match character {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    other => panic!("Not a DNA character: {other}"),
                }
                .into();

                let result = result << 2;
                result | bits
            }),
        }
    }
}

impl FromIterator<u8> for BitPackedVectorKmer {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        let iter = iter.into_iter();
        let kmer = BitVec::with_capacity(iter.size_hint().0 * 2);
        BitPackedVectorKmer {
            kmer: iter.fold(kmer, |mut result, character| {
                let bits = match character {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    other => panic!("Not a DNA character: {other}"),
                };

                result.push(bits & 2 != 0);
                result.push(bits & 1 != 0);
                result
            }),
        }
    }
}

impl<
        const K: usize,
        Integer: BitAnd<Integer, Output = Integer>
            + ShlAssign<i32>
            + Shl<usize, Output = Integer>
            + Shr<usize, Output = Integer>
            + From<u8>
            + Copy,
    > Display for BitPackedKmer<K, Integer>
where
    usize: TryFrom<Integer>,
    <usize as TryFrom<Integer>>::Error: Debug,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        static CHARACTERS: [char; 4] = ['A', 'C', 'G', 'T'];

        let mut current = self.kmer << (std::mem::size_of::<Integer>() * 8 - 2 * K);
        let mask_shift = std::mem::size_of::<Integer>() * 8 - 2;
        let mask = Integer::from(3u8) << mask_shift;
        for _ in 0..K {
            let bits = (current & mask) >> mask_shift;
            current <<= 2;
            let character = CHARACTERS[usize::try_from(bits).unwrap()];
            write!(f, "{character}",)?;
        }

        Ok(())
    }
}

impl Display for BitPackedVectorKmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        assert_eq!(self.kmer.len() % 2, 0);
        for bits in self.kmer.chunks(2) {
            write!(
                f,
                "{}",
                match (bits[0], bits[1]) {
                    (false, false) => 'A',
                    (false, true) => 'C',
                    (true, false) => 'G',
                    (true, true) => 'T',
                }
            )?;
        }

        Ok(())
    }
}

impl<
        const K: usize,
        Integer: BitAnd<Integer, Output = Integer>
            + BitOrAssign<Integer>
            + Not<Output = Integer>
            + ShlAssign<i32>
            + ShrAssign<i32>
            + From<u8>
            + Copy
            + Ord,
    > Kmer for BitPackedKmer<K, Integer>
{
    fn reverse_complement(&self) -> Self {
        let mut source = !self.kmer;
        let mut result = 0.into();
        for _ in 0..K {
            result <<= 2;
            result |= source & 3.into();
            source >>= 2;
        }

        BitPackedKmer { kmer: result }
    }

    fn predecessor(&self, character: u8) -> Self {
        let mut character_bits = Integer::from(match character {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => panic!("Not a DNA character: {other}"),
        });
        character_bits <<= (i32::try_from(K).unwrap() - 1) * 2;

        let mut kmer = self.kmer;
        kmer >>= 2;
        kmer |= character_bits;

        Self { kmer }
    }

    fn successor(&self, character: u8) -> Self {
        let character_bits = match character {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => panic!("Not a DNA character: {other}"),
        };

        let mut kmer = self.kmer;
        kmer <<= 2;
        kmer |= character_bits.into();

        // Clear high bits.
        let mut mask = Integer::from(3);
        mask <<= i32::try_from(K).unwrap() * 2;
        mask = !mask;
        kmer = kmer & mask;

        Self { kmer }
    }
}

impl Kmer for BitPackedVectorKmer {
    fn reverse_complement(&self) -> Self {
        assert_eq!(self.kmer.len() % 2, 0);
        Self {
            kmer: self
                .kmer
                .chunks(2)
                .rev()
                .flat_map(|bits| [!bits[0], !bits[1]])
                .collect(),
        }
    }

    fn predecessor(&self, character: u8) -> Self {
        let bits = match character {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => panic!("Not a DNA character: {other}"),
        };

        let mut kmer = self.kmer.clone();
        kmer.shift_right(2);
        kmer.set(1, bits & 1 != 0);
        kmer.set(0, bits & 2 != 0);

        Self { kmer }
    }

    fn successor(&self, character: u8) -> Self {
        let bits = match character {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => panic!("Not a DNA character: {other}"),
        };

        let mut kmer = self.kmer.clone();
        kmer.shift_left(2);
        let kmer_len = kmer.len();
        kmer.set(kmer_len - 1, bits & 1 != 0);
        kmer.set(kmer_len - 2, bits & 2 != 0);

        Self { kmer }
    }
}

#[cfg(test)]
mod tests {
    use crate::kmer::{BitPackedVectorKmer, Kmer};
    use crate::BitPackedKmer;

    #[test]
    fn test_k31_display() {
        let kmer: String = "A".repeat(31);
        let bit_packed_kmer = BitPackedKmer::<31, u64>::from_iter(kmer.bytes());
        assert_eq!(format!("{bit_packed_kmer}"), kmer);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            BitPackedKmer::<3, u8>::from_iter("AAA".as_bytes().iter().copied())
                .reverse_complement(),
            BitPackedKmer::<3, u8>::from_iter("TTT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedKmer::<3, u8>::from_iter("ACA".as_bytes().iter().copied())
                .reverse_complement(),
            BitPackedKmer::<3, u8>::from_iter("TGT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedKmer::<3, u8>::from_iter("ACC".as_bytes().iter().copied())
                .reverse_complement(),
            BitPackedKmer::<3, u8>::from_iter("GGT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedKmer::<51, u128>::from_iter(
                "ACAACAACAACAACAACAACAACAACAACAACAACAACAACAACATTTTTT"
                    .as_bytes()
                    .iter()
                    .copied()
            )
            .reverse_complement(),
            BitPackedKmer::<51, u128>::from_iter(
                "AAAAAATGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGT"
                    .as_bytes()
                    .iter()
                    .copied()
            )
        );
        assert_eq!(
            BitPackedKmer::<4, u8>::from_iter("ACAA".as_bytes().iter().copied())
                .reverse_complement(),
            BitPackedKmer::<4, u8>::from_iter("TTGT".as_bytes().iter().copied())
        );

        assert_eq!(
            BitPackedVectorKmer::from_iter("AAA".as_bytes().iter().copied()).reverse_complement(),
            BitPackedVectorKmer::from_iter("TTT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("ACA".as_bytes().iter().copied()).reverse_complement(),
            BitPackedVectorKmer::from_iter("TGT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("ACC".as_bytes().iter().copied()).reverse_complement(),
            BitPackedVectorKmer::from_iter("GGT".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter(
                "ACAACAACAACAACAACAACAACAACAACAACAACAACAACAACATTTTTT"
                    .as_bytes()
                    .iter()
                    .copied()
            )
            .reverse_complement(),
            BitPackedVectorKmer::from_iter(
                "AAAAAATGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGT"
                    .as_bytes()
                    .iter()
                    .copied()
            )
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("ACAA".as_bytes().iter().copied()).reverse_complement(),
            BitPackedVectorKmer::from_iter("TTGT".as_bytes().iter().copied())
        );
    }

    #[test]
    fn test_predecessor_successor() {
        assert_eq!(
            BitPackedKmer::<3, u8>::from_iter("GGG".as_bytes().iter().copied()).successor(b'C'),
            BitPackedKmer::<3, u8>::from_iter("GGC".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedKmer::<3, u8>::from_iter("GGG".as_bytes().iter().copied()).predecessor(b'C'),
            BitPackedKmer::<3, u8>::from_iter("CGG".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied()).successor(b'C'),
            BitPackedVectorKmer::from_iter("GGC".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied()).predecessor(b'C'),
            BitPackedVectorKmer::from_iter("CGG".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied()).successor(b'A'),
            BitPackedVectorKmer::from_iter("GGA".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied()).predecessor(b'A'),
            BitPackedVectorKmer::from_iter("AGG".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied())
                .successor(b'A')
                .predecessor(b'A'),
            BitPackedVectorKmer::from_iter("AGG".as_bytes().iter().copied())
        );
        assert_eq!(
            BitPackedVectorKmer::from_iter("GGG".as_bytes().iter().copied())
                .predecessor(b'A')
                .successor(b'A'),
            BitPackedVectorKmer::from_iter("GGA".as_bytes().iter().copied())
        );
    }
}
