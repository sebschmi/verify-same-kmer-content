use std::fmt::{Display, Formatter};
use std::ops::{BitAnd, BitOr, BitOrAssign, Not, Shl, ShlAssign, ShrAssign};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct BitPackedKmer<const K: usize, Integer> {
    kmer: Integer,
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

impl<
        const K: usize,
        Integer: BitAnd<Integer, Output = Integer>
            + ShlAssign<i32>
            + Shl<usize, Output = Integer>
            + From<u8>
            + Into<u128>
            + Copy,
    > Display for BitPackedKmer<K, Integer>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut current = self.kmer << (std::mem::size_of::<Integer>() * 8 - 2 * K);
        let mask = (3 << (std::mem::size_of::<Integer>() * 8 - 2)).into();
        for _ in 0..K {
            let bits = (current & mask).into() >> (std::mem::size_of::<Integer>() * 8 - 2);
            current <<= 2;
            write!(
                f,
                "{}",
                match bits {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!("Masking with 3 cannot result in any other number"),
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

pub trait Kmer {
    fn reverse_complement(&self) -> Self;
}

#[cfg(test)]
mod tests {
    use crate::kmer::Kmer;
    use crate::BitPackedKmer;

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
    }
}
