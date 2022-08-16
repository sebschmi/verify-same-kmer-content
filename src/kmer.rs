use std::fmt::{Display, Formatter};
use std::ops::{BitAnd, BitOr, Shl, ShrAssign};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct Kmer<const K: usize, Integer> {
    kmer: Integer,
}

impl<
        const K: usize,
        Integer: Default + Shl<i32, Output = Integer> + BitOr<Integer, Output = Integer> + From<u8>,
    > FromIterator<u8> for Kmer<K, Integer>
{
    fn from_iter<Iter: IntoIterator<Item = u8>>(iter: Iter) -> Self {
        assert!(2 * K <= 8 * std::mem::size_of::<Integer>());

        let iter = iter.into_iter();
        let size = iter.size_hint();
        assert_eq!(Some(size.0), size.1);
        assert_eq!(K, size.0);

        Kmer {
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
        Integer: BitAnd<Integer, Output = Integer> + ShrAssign<i32> + From<u8> + Into<u128> + Copy,
    > Display for Kmer<K, Integer>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut current = self.kmer;
        for _ in 0..K {
            let bits = (current & 3.into()).into();
            current >>= 2;
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
