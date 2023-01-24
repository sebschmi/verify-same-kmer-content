use log::warn;
use std::collections::VecDeque;
use std::io::{BufReader, Read};
use std::marker::PhantomData;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum State {
    None,
    GfaS,
    GfaSequence,
    FaId,
    FaSequence,
    Eof,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum Format {
    None,
    Gfa,
    Fa,
}

pub struct KmerIterator<InputReader: Read, KmerType> {
    input: BufReader<InputReader>,
    k: usize,
    state: State,
    format: Format,
    buffer: VecDeque<u8>,
    character_buffer: [u8; 1],
    sequence_count: usize,
    character_count: usize,
    panic_on_parse_error: bool,
    kmer_type: PhantomData<KmerType>,
}

impl<InputReader: Read, KmerType> KmerIterator<InputReader, KmerType> {
    pub fn new(input: InputReader, k: usize, panic_on_parse_error: bool) -> Self {
        Self {
            input: BufReader::with_capacity(16 * 1024 * 1024, input),
            k,
            state: State::None,
            format: Format::None,
            buffer: Default::default(),
            character_buffer: Default::default(),
            sequence_count: 0,
            character_count: 0,
            panic_on_parse_error,
            kmer_type: Default::default(),
        }
    }

    fn read_char(&mut self) -> Option<u8> {
        let read = self.input.read(&mut self.character_buffer).unwrap();
        if read == 1 {
            Some(self.character_buffer[0])
        } else {
            None
        }
    }

    pub fn sequence_count(&self) -> usize {
        self.sequence_count
    }

    pub fn character_count(&self) -> usize {
        self.character_count
    }
}

impl<InputReader: Read, KmerType: FromIterator<u8>> Iterator
    for KmerIterator<InputReader, KmerType>
{
    type Item = KmerType;

    fn next(&mut self) -> Option<KmerType> {
        while self.state != State::Eof {
            match self.state {
                State::None => loop {
                    let character = self.read_char();
                    if character == Some(b'S') {
                        if self.format == Format::Fa {
                            if self.panic_on_parse_error {
                                panic!("Found GFA within fasta");
                            } else {
                                warn!("Found GFA within fasta");
                            }
                        } else {
                            self.format = Format::Gfa;
                        }

                        self.state = State::GfaS;
                        break;
                    } else if character == Some(b'>') {
                        if self.format == Format::Gfa {
                            if self.panic_on_parse_error {
                                panic!("Found fasta within GFA");
                            } else {
                                warn!("Found fasta within GFA");
                            }
                        } else {
                            self.format = Format::Fa;
                        }

                        self.state = State::FaId;
                        break;
                    } else if character.is_none() {
                        self.state = State::Eof;
                        break;
                    }
                },
                State::GfaS => {
                    let character = self.read_char();
                    if character == Some(b'\t') {
                        loop {
                            let character = self.read_char();
                            if character == Some(b'\t') {
                                self.sequence_count += 1;
                                self.state = State::GfaSequence;
                                break;
                            } else if character.is_none() {
                                self.state = State::Eof;
                                break;
                            }
                        }
                    }
                }
                State::GfaSequence => {
                    while self.state == State::GfaSequence {
                        let character = self.read_char();
                        if let Some(character) = character {
                            let character = character.to_ascii_uppercase();
                            match character {
                                b'A' | b'C' | b'G' | b'T' => {
                                    self.buffer.push_back(character);
                                }
                                _ => {
                                    self.state = State::None;
                                }
                            }
                        } else {
                            self.state = State::Eof;
                        }

                        assert!(self.buffer.len() <= self.k);
                        if self.buffer.len() == self.k {
                            let kmer = self.buffer.iter().copied().collect();
                            self.character_count += 1;
                            self.buffer.pop_front();
                            return Some(kmer);
                        }
                    }

                    self.character_count += self.buffer.len();
                    self.buffer.clear();
                }
                State::FaId => loop {
                    let character = self.read_char();
                    if character == Some(b'\n') {
                        self.sequence_count += 1;
                        self.state = State::FaSequence;
                        break;
                    } else if character.is_none() {
                        self.state = State::Eof;
                        break;
                    }
                },
                State::FaSequence => {
                    while self.state == State::FaSequence {
                        let character = self.read_char();
                        if let Some(character) = character {
                            let character = character.to_ascii_uppercase();
                            match character {
                                b'A' | b'C' | b'G' | b'T' => {
                                    self.buffer.push_back(character);
                                }
                                b'\n' => { /* ignore newlines */ }
                                b'>' => {
                                    self.state = State::FaId;
                                }
                                _ => {
                                    self.state = State::None;
                                }
                            }
                        } else {
                            self.state = State::Eof;
                        }

                        assert!(self.buffer.len() <= self.k);
                        if self.buffer.len() == self.k {
                            let kmer = self.buffer.iter().copied().collect();
                            self.character_count += 1;
                            self.buffer.pop_front();
                            return Some(kmer);
                        }
                    }

                    self.character_count += self.buffer.len();
                    self.buffer.clear();
                }
                State::Eof => unreachable!("Loop is not entered when self.state == State::Eof"),
            }
        }

        assert_eq!(self.state, State::Eof);
        if self.format == Format::None {
            if self.panic_on_parse_error {
                panic!("Found no kmers");
            } else {
                warn!("Found no kmers");
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::{initialise_logging, BitPackedKmer, KmerIterator};
    use log::LevelFilter;

    #[test]
    fn test_simple_fa() {
        initialise_logging(LevelFilter::Debug);
        let tigs = ">b\nAAAC\n>\nCAGT\n>a\nCCC";
        let mut iterator = KmerIterator::<_, BitPackedKmer<3, u8>>::new(tigs.as_bytes(), 3, true);
        let kmers: Vec<_> = iterator.by_ref().collect();
        assert_eq!(
            kmers,
            vec![
                BitPackedKmer::from_iter("AAA".as_bytes().iter().copied()),
                BitPackedKmer::from_iter("AAC".as_bytes().iter().copied()),
                BitPackedKmer::from_iter("CAG".as_bytes().iter().copied()),
                BitPackedKmer::from_iter("AGT".as_bytes().iter().copied()),
                BitPackedKmer::from_iter("CCC".as_bytes().iter().copied()),
            ]
        );
        assert_eq!(iterator.sequence_count(), 3);
        assert_eq!(iterator.character_count(), 11);
    }
}
