use log::warn;
use std::collections::VecDeque;
use std::io::{BufReader, Read};
use std::marker::PhantomData;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
#[allow(non_camel_case_types)]
enum State {
    None,
    GFA_S,
    GFA_Sequence,
    FA_ID,
    FA_Sequence,
    EOF,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum Format {
    None,
    GFA,
    FA,
}

pub struct KmerIterator<InputReader: Read, KmerType> {
    input: BufReader<InputReader>,
    k: usize,
    state: State,
    format: Format,
    buffer: VecDeque<u8>,
    character_buffer: [u8; 1],
    kmer_type: PhantomData<KmerType>,
}

impl<InputReader: Read, KmerType> KmerIterator<InputReader, KmerType> {
    pub fn new(input: InputReader, k: usize) -> Self {
        Self {
            input: BufReader::with_capacity(16 * 1024 * 1024, input),
            k,
            state: State::None,
            format: Format::None,
            buffer: Default::default(),
            character_buffer: Default::default(),
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
                        if self.format == Format::FA {
                            warn!("Found GFA within fasta");
                        } else {
                            self.format = Format::GFA;
                        }

                        self.state = State::GFA_S;
                        break;
                    } else if character == Some(b'>') {
                        if self.format == Format::GFA {
                            warn!("Found fasta within GFA");
                        } else {
                            self.format = Format::FA;
                        }

                        self.state = State::FA_ID;
                        break;
                    } else if character == None {
                        self.state = State::EOF;
                        break;
                    }
                },
                State::GFA_S => {
                    let character = self.read_char();
                    if character == Some(b'\t') {
                        loop {
                            let character = self.read_char();
                            if character == Some(b'\t') {
                                self.state = State::GFA_Sequence;
                                break;
                            } else if character == None {
                                self.state = State::EOF;
                                break;
                            }
                        }
                    }
                }
                State::GFA_Sequence => loop {
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
                        self.state = State::EOF;
                        self.buffer.clear();
                        break;
                    }
                },
                State::FA_ID => loop {
                    let character = self.read_char();
                    if character == Some(b'\n') {
                        self.state = State::FA_Sequence;
                        break;
                    } else if character == None {
                        self.state = State::EOF;
                        break;
                    }
                },
                State::FA_Sequence => {
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
                            self.state = State::EOF;
                            self.buffer.clear();
                            break;
                        }
                    }
                }
                State::EOF => break,
            }
        }

        assert_eq!(self.state, State::EOF);
        if self.format == Format::None {
            warn!("Found no kmers");
        }
        None
    }
}
