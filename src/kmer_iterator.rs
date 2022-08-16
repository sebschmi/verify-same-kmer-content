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
}

impl<InputReader: Read, KmerType: FromIterator<u8>> Iterator
    for KmerIterator<InputReader, KmerType>
{
    type Item = KmerType;

    fn next(&mut self) -> Option<KmerType> {
        loop {
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

                        self.sequence_count += 1;
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

                        self.sequence_count += 1;
                        self.state = State::FaId;
                        break;
                    } else if character == None {
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
                                self.state = State::GfaSequence;
                                break;
                            } else if character == None {
                                self.state = State::Eof;
                                break;
                            }
                        }
                    }
                }
                State::GfaSequence => loop {
                    let character = self.read_char();
                    if let Some(character) = character {
                        let character = character.to_ascii_uppercase();
                        match character {
                            b'A' | b'C' | b'G' | b'T' => {
                                self.buffer.push_back(character);
                            }
                            _ => {
                                self.state = State::None;
                                self.buffer.clear();
                                break;
                            }
                        }
                        assert!(self.buffer.len() <= self.k);
                        if self.buffer.len() == self.k {
                            let kmer = self.buffer.iter().copied().collect();
                            self.buffer.pop_front();
                            return Some(kmer);
                        }
                    } else {
                        self.state = State::Eof;
                        self.buffer.clear();
                        break;
                    }
                },
                State::FaId => loop {
                    let character = self.read_char();
                    if character == Some(b'\n') {
                        self.state = State::FaSequence;
                        break;
                    } else if character == None {
                        self.state = State::Eof;
                        break;
                    }
                },
                State::FaSequence => {
                    loop {
                        let character = self.read_char();
                        if let Some(character) = character {
                            let character = character.to_ascii_uppercase();
                            match character {
                                b'A' | b'C' | b'G' | b'T' => {
                                    self.buffer.push_back(character);
                                }
                                b'\n' => { /* ignore newlines */ }
                                _ => {
                                    self.state = State::None;
                                    self.buffer.clear();
                                    break;
                                }
                            }
                            assert!(self.buffer.len() <= self.k);
                            if self.buffer.len() == self.k {
                                let kmer = self.buffer.iter().copied().collect();
                                self.buffer.pop_front();
                                return Some(kmer);
                            }
                        } else {
                            self.state = State::Eof;
                            self.buffer.clear();
                            break;
                        }
                    }
                }
                State::Eof => break,
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
