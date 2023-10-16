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

pub trait Kmer {
    fn reverse_complement(&self) -> Self;
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
            + Copy,
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
}
