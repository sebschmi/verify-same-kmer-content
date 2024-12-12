use crate::kmer::{BitPackedKmer, BitPackedVectorKmer, Kmer};
use crate::kmer_iterator::KmerIterator;
use clap::Parser;
use log::{debug, error, info, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, TermLogger, TerminalMode};
use std::cmp::Ordering;
use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Mutex;

mod kmer;
mod kmer_iterator;

static LOGGING_INITIALISED: Mutex<bool> = Mutex::new(false);

pub fn initialise_logging(log_level: LevelFilter) {
    let mut logging_initialised = LOGGING_INITIALISED.lock().unwrap();

    if !*logging_initialised {
        CombinedLogger::init(vec![TermLogger::new(
            log_level,
            Default::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        )])
        .unwrap();

        info!("Logging initialised successfully");
        *logging_initialised = true;
    }
}

/// Verify that an SPSS contains the same kmer content as a set of unitigs.
#[derive(Parser, Debug)]
pub struct Config {
    /// The desired log level.
    #[clap(short, long, default_value = "Info")]
    log_level: LevelFilter,

    /// The kmer size.
    #[clap(short)]
    k: usize,

    /// Skip the actual verification, and only compute statistics.
    #[clap(long)]
    do_not_verify: bool,

    /// Do not print warnings during parsing, but instead abort if there is any warning.
    #[clap(long)]
    panic_on_parse_error: bool,

    /// Do not treat k-mers in the test tigs as missing if they are missing due to cuttlefish2's error.
    ///
    /// This allows k-mers to be missing if they are not part of any k+1-mer.
    /// See [this github issue][1] for details.
    ///
    /// [1]: https://github.com/COMBINE-lab/cuttlefish/issues/36
    #[clap(long)]
    allow_cuttlefish2_errors: bool,

    /// A file containing the ground truth kmer set as unitigs.
    #[clap(index = 1)]
    unitigs: PathBuf,

    /// A file containing the test kmer set as any set of strings.
    #[clap(index = 2)]
    test_tigs: PathBuf,
}

#[derive(Debug)]
enum Error {
    Mismatch,
    IllegalKmerSize {
        #[allow(dead_code)]
        kmer_size: usize,
    },
}

fn compare_kmer_sets<KmerType: FromIterator<u8> + Ord + Clone + Display + Kmer>(
    unitigs: impl Read,
    test_tigs: impl Read,
    config: Config,
) -> Result<(), Error> {
    let mut kmer_iter_unitigs =
        KmerIterator::<_, KmerType>::new(unitigs, config.k, config.panic_on_parse_error);
    let mut kmer_iter_test_tigs =
        KmerIterator::<_, KmerType>::new(test_tigs, config.k, config.panic_on_parse_error);

    let (has_superfluous_kmers_unitigs, has_superfluous_kmers_test_tigs) = if !config.do_not_verify
    {
        info!("Reading first input file");
        let mut kmers_unitigs: Vec<_> = kmer_iter_unitigs
            .by_ref()
            .map(|kmer| Kmer::canonical(&kmer))
            .collect();
        let input_unitig_kmer_amount = kmers_unitigs.len();
        info!("Sorting kmers in first input file");
        kmers_unitigs.sort_unstable();

        info!("Removing duplicates from first input file");
        let mut previous_kmer = None;
        kmers_unitigs.retain(|kmer| {
            if let Some(previous_kmer) = previous_kmer.as_mut() {
                let result = kmer != previous_kmer;
                *previous_kmer = kmer.clone();
                result
            } else {
                previous_kmer = Some(kmer.clone());
                true
            }
        });

        let kmers_unitigs = kmers_unitigs;
        let duplicate_unitig_kmer_amount = input_unitig_kmer_amount - kmers_unitigs.len();
        debug!(
            "Duplicate kmers: {duplicate_unitig_kmer_amount}/{input_unitig_kmer_amount} ({:.0}%)",
            duplicate_unitig_kmer_amount as f64 / input_unitig_kmer_amount as f64
        );

        assert_eq!(
            kmers_unitigs.len() + duplicate_unitig_kmer_amount,
            kmer_iter_unitigs.character_count()
                - kmer_iter_unitigs.sequence_count() * (config.k - 1),
            "unitigs: character_count: {}; sequence_count: {}; k: {}",
            kmer_iter_unitigs.character_count(),
            kmer_iter_unitigs.sequence_count(),
            config.k
        );

        let unitig_kmers_without_superstrings = if config.allow_cuttlefish2_errors {
            info!("Collecting kmers without superstrings");
            kmers_unitigs
                .iter()
                .filter(|&kmer| !has_superstring(kmer, &kmers_unitigs))
                .cloned()
                .collect()
        } else {
            Vec::new()
        };
        debug_assert!(unitig_kmers_without_superstrings.is_sorted());
        for kmer in &unitig_kmers_without_superstrings {
            debug!("Unitig kmer without superstrings: {kmer}");
        }

        info!("Reading second input file");
        let mut kmers_test_tigs: Vec<_> = kmer_iter_test_tigs
            .by_ref()
            .map(|kmer| Kmer::canonical(&kmer))
            .collect();
        let input_test_tig_kmer_amount = kmers_test_tigs.len();
        info!("Sorting kmers in second input file");
        kmers_test_tigs.sort_unstable();

        info!("Removing duplicates from second input file");
        let mut previous_kmer = None;
        kmers_test_tigs.retain(|kmer| {
            if let Some(previous_kmer) = previous_kmer.as_mut() {
                let result = kmer != previous_kmer;
                *previous_kmer = kmer.clone();
                result
            } else {
                previous_kmer = Some(kmer.clone());
                true
            }
        });

        let kmers_test_tigs = kmers_test_tigs;
        let duplicate_test_tig_kmer_amount = input_test_tig_kmer_amount - kmers_test_tigs.len();
        debug!(
            "Duplicate kmers: {duplicate_test_tig_kmer_amount}/{input_test_tig_kmer_amount} ({:.0}%)",
            duplicate_test_tig_kmer_amount as f64 / input_test_tig_kmer_amount as f64
        );

        assert_eq!(
            kmers_test_tigs.len() + duplicate_test_tig_kmer_amount,
            kmer_iter_test_tigs.character_count()
                - kmer_iter_test_tigs.sequence_count() * (config.k - 1),
            "unitigs: character_count: {}; sequence_count: {}; k: {}",
            kmer_iter_test_tigs.character_count(),
            kmer_iter_test_tigs.sequence_count(),
            config.k
        );

        info!("Comparing kmer content");
        let mut unitig_kmer_iterator = kmers_unitigs.iter().peekable();
        let mut test_tig_kmer_iterator = kmers_test_tigs.iter().peekable();
        let mut superfluous_unitig_kmer_count = 0usize;
        let mut superfluous_test_tig_kmer_count = 0usize;

        while let (Some(unitig_kmer), Some(test_tig_kmer)) =
            (unitig_kmer_iterator.peek(), test_tig_kmer_iterator.peek())
        {
            match unitig_kmer.cmp(test_tig_kmer) {
                Ordering::Less => {
                    if unitig_kmers_without_superstrings
                        .binary_search(unitig_kmer)
                        .is_err()
                    {
                        superfluous_unitig_kmer_count += 1;
                        debug!("Unitigs contain kmer that is missing in test tigs: {unitig_kmer}");
                    }
                    unitig_kmer_iterator.next().unwrap();
                }
                Ordering::Equal => {
                    unitig_kmer_iterator.next().unwrap();
                    test_tig_kmer_iterator.next().unwrap();
                }
                Ordering::Greater => {
                    superfluous_test_tig_kmer_count += 1;
                    debug!("Test tigs contains kmer that is missing in unitigs: {test_tig_kmer}");
                    test_tig_kmer_iterator.next().unwrap();
                }
            }
        }

        if superfluous_unitig_kmer_count != 0 {
            info!(
                "Test tigs miss {superfluous_unitig_kmer_count} kmers that are present in unitigs"
            );
        }
        if superfluous_test_tig_kmer_count != 0 {
            info!("Test tigs contain {superfluous_test_tig_kmer_count} kmers that are not present in unitigs");
        }

        (
            superfluous_unitig_kmer_count != 0,
            superfluous_test_tig_kmer_count != 0,
        )
    } else {
        info!("Reading first input file");
        assert!(kmer_iter_unitigs.by_ref().all(|_| true));

        info!("Reading second input file");
        assert!(kmer_iter_test_tigs.by_ref().all(|_| true));
        (false, false)
    };

    let unitigs_sequence_size = kmer_iter_unitigs.character_count();
    let test_tigs_sequence_size = kmer_iter_test_tigs.character_count();
    let unitigs_string_count = kmer_iter_unitigs.sequence_count();
    let test_tigs_string_count = kmer_iter_test_tigs.sequence_count();
    let compression_rate = test_tigs_sequence_size as f64 / unitigs_sequence_size as f64;
    let string_count_rate = test_tigs_string_count as f64 / unitigs_string_count as f64;
    let unique_kmer_count = unitigs_sequence_size - unitigs_string_count * (config.k - 1);
    let test_tigs_kmer_count = test_tigs_sequence_size - test_tigs_string_count * (config.k - 1);

    std::io::stdout().flush().unwrap();
    std::io::stderr().flush().unwrap();
    println!("ground_truth_size:   {unitigs_sequence_size}");
    println!("test_size: {test_tigs_sequence_size}");

    println!("ground_truth_str_cnt:   {unitigs_string_count}");
    println!("test_str_cnt: {test_tigs_string_count}");

    println!("compression_rate: {compression_rate}");
    println!("str_cnt_rate: {string_count_rate}");

    println!("unique_kmer_count: {unique_kmer_count}");
    std::io::stdout().flush().unwrap();
    std::io::stderr().flush().unwrap();

    if !has_superfluous_kmers_unitigs && !has_superfluous_kmers_test_tigs {
        match unique_kmer_count.cmp(&test_tigs_kmer_count) {
            Ordering::Greater => {
                debug!("Unitig kmer count: {unique_kmer_count}");
                debug!("Test tigs kmer count: {test_tigs_kmer_count}");
                if config.allow_cuttlefish2_errors {
                    debug!("Missing kmers in test tigs are ignored because cuttlefish2 errors are allowed.");
                    info!("Success!");
                    Ok(())
                } else {
                    error!("Test tigs are missing kmers. Note that the test tigs are assumed to contain no duplicate kmers.");
                    Err(Error::Mismatch)
                }
            }
            Ordering::Equal => {
                info!("Success!");
                Ok(())
            }
            Ordering::Less => {
                debug!("Unitig kmer count: {unique_kmer_count}");
                debug!("Test tigs kmer count: {test_tigs_kmer_count}");
                info!("Test tigs contain more kmers than unitigs. This may happen if they contain duplicates.");
                info!("Success!");
                Ok(())
            }
        }
    } else if !has_superfluous_kmers_unitigs {
        error!("Test tigs contain kmers that are missing in unitigs");
        Err(Error::Mismatch)
    } else if !has_superfluous_kmers_test_tigs {
        error!("Test tigs miss kmers that are present in unitigs");
        Err(Error::Mismatch)
    } else {
        error!("Test tigs both miss kmers and contain kmers that are not present in unitigs");
        Err(Error::Mismatch)
    }
}

fn has_superstring<KmerType: FromIterator<u8> + Ord + Clone + Display + Kmer>(
    kmer: &KmerType,
    all_kmers: &[KmerType],
) -> bool {
    debug_assert!(all_kmers.is_sorted());

    for &character in b"ACGT" {
        let predecessor = kmer.predecessor(character);
        let successor = kmer.successor(character);
        let predecessor_rc = predecessor.reverse_complement();
        let successor_rc = successor.reverse_complement();

        if all_kmers.binary_search(&predecessor).is_ok()
            || all_kmers.binary_search(&successor).is_ok()
            || all_kmers.binary_search(&predecessor_rc).is_ok()
            || all_kmers.binary_search(&successor_rc).is_ok()
        {
            return true;
        }
    }

    false
}

fn main() -> Result<(), Error> {
    let config = Config::parse();
    initialise_logging(config.log_level);
    debug!("{config:?}");

    let unitigs_file = File::open(&config.unitigs)
        .unwrap_or_else(|_| panic!("--unitigs points to a file: {:?}", &config.unitigs));
    let test_tigs_file = File::open(&config.test_tigs)
        .unwrap_or_else(|_| panic!("--test-tigs points to a file: {:?}", &config.test_tigs));

    // This is not the most clever way to handle different kmer sizes in the type system, but it gets the job done.
    // It results in larger binary sizes, but therefore we can have e.g. a Display implementation for Kmer.
    match config.k {
        0 => {
            error!("Kmer size cannot be zero");
            Err(Error::IllegalKmerSize {
                kmer_size: config.k,
            })
        }
        1 => compare_kmer_sets::<BitPackedKmer<1, u8>>(unitigs_file, test_tigs_file, config),
        2 => compare_kmer_sets::<BitPackedKmer<2, u8>>(unitigs_file, test_tigs_file, config),
        3 => compare_kmer_sets::<BitPackedKmer<3, u8>>(unitigs_file, test_tigs_file, config),
        4 => compare_kmer_sets::<BitPackedKmer<4, u8>>(unitigs_file, test_tigs_file, config),
        5 => compare_kmer_sets::<BitPackedKmer<5, u16>>(unitigs_file, test_tigs_file, config),
        6 => compare_kmer_sets::<BitPackedKmer<6, u16>>(unitigs_file, test_tigs_file, config),
        7 => compare_kmer_sets::<BitPackedKmer<7, u16>>(unitigs_file, test_tigs_file, config),
        8 => compare_kmer_sets::<BitPackedKmer<8, u16>>(unitigs_file, test_tigs_file, config),
        9 => compare_kmer_sets::<BitPackedKmer<9, u32>>(unitigs_file, test_tigs_file, config),
        10 => compare_kmer_sets::<BitPackedKmer<10, u32>>(unitigs_file, test_tigs_file, config),
        11 => compare_kmer_sets::<BitPackedKmer<11, u32>>(unitigs_file, test_tigs_file, config),
        12 => compare_kmer_sets::<BitPackedKmer<12, u32>>(unitigs_file, test_tigs_file, config),
        13 => compare_kmer_sets::<BitPackedKmer<13, u32>>(unitigs_file, test_tigs_file, config),
        14 => compare_kmer_sets::<BitPackedKmer<14, u32>>(unitigs_file, test_tigs_file, config),
        15 => compare_kmer_sets::<BitPackedKmer<15, u32>>(unitigs_file, test_tigs_file, config),
        16 => compare_kmer_sets::<BitPackedKmer<16, u32>>(unitigs_file, test_tigs_file, config),
        17 => compare_kmer_sets::<BitPackedKmer<17, u64>>(unitigs_file, test_tigs_file, config),
        18 => compare_kmer_sets::<BitPackedKmer<18, u64>>(unitigs_file, test_tigs_file, config),
        19 => compare_kmer_sets::<BitPackedKmer<19, u64>>(unitigs_file, test_tigs_file, config),
        20 => compare_kmer_sets::<BitPackedKmer<20, u64>>(unitigs_file, test_tigs_file, config),
        21 => compare_kmer_sets::<BitPackedKmer<21, u64>>(unitigs_file, test_tigs_file, config),
        22 => compare_kmer_sets::<BitPackedKmer<22, u64>>(unitigs_file, test_tigs_file, config),
        23 => compare_kmer_sets::<BitPackedKmer<23, u64>>(unitigs_file, test_tigs_file, config),
        24 => compare_kmer_sets::<BitPackedKmer<24, u64>>(unitigs_file, test_tigs_file, config),
        25 => compare_kmer_sets::<BitPackedKmer<25, u64>>(unitigs_file, test_tigs_file, config),
        26 => compare_kmer_sets::<BitPackedKmer<26, u64>>(unitigs_file, test_tigs_file, config),
        27 => compare_kmer_sets::<BitPackedKmer<27, u64>>(unitigs_file, test_tigs_file, config),
        28 => compare_kmer_sets::<BitPackedKmer<28, u64>>(unitigs_file, test_tigs_file, config),
        29 => compare_kmer_sets::<BitPackedKmer<29, u64>>(unitigs_file, test_tigs_file, config),
        30 => compare_kmer_sets::<BitPackedKmer<30, u64>>(unitigs_file, test_tigs_file, config),
        31 => compare_kmer_sets::<BitPackedKmer<31, u64>>(unitigs_file, test_tigs_file, config),
        32 => compare_kmer_sets::<BitPackedKmer<32, u64>>(unitigs_file, test_tigs_file, config),
        33 => compare_kmer_sets::<BitPackedKmer<33, u128>>(unitigs_file, test_tigs_file, config),
        34 => compare_kmer_sets::<BitPackedKmer<34, u128>>(unitigs_file, test_tigs_file, config),
        35 => compare_kmer_sets::<BitPackedKmer<35, u128>>(unitigs_file, test_tigs_file, config),
        36 => compare_kmer_sets::<BitPackedKmer<36, u128>>(unitigs_file, test_tigs_file, config),
        37 => compare_kmer_sets::<BitPackedKmer<37, u128>>(unitigs_file, test_tigs_file, config),
        38 => compare_kmer_sets::<BitPackedKmer<38, u128>>(unitigs_file, test_tigs_file, config),
        39 => compare_kmer_sets::<BitPackedKmer<39, u128>>(unitigs_file, test_tigs_file, config),
        40 => compare_kmer_sets::<BitPackedKmer<40, u128>>(unitigs_file, test_tigs_file, config),
        41 => compare_kmer_sets::<BitPackedKmer<41, u128>>(unitigs_file, test_tigs_file, config),
        42 => compare_kmer_sets::<BitPackedKmer<42, u128>>(unitigs_file, test_tigs_file, config),
        43 => compare_kmer_sets::<BitPackedKmer<43, u128>>(unitigs_file, test_tigs_file, config),
        44 => compare_kmer_sets::<BitPackedKmer<44, u128>>(unitigs_file, test_tigs_file, config),
        45 => compare_kmer_sets::<BitPackedKmer<45, u128>>(unitigs_file, test_tigs_file, config),
        46 => compare_kmer_sets::<BitPackedKmer<46, u128>>(unitigs_file, test_tigs_file, config),
        47 => compare_kmer_sets::<BitPackedKmer<47, u128>>(unitigs_file, test_tigs_file, config),
        48 => compare_kmer_sets::<BitPackedKmer<48, u128>>(unitigs_file, test_tigs_file, config),
        49 => compare_kmer_sets::<BitPackedKmer<49, u128>>(unitigs_file, test_tigs_file, config),
        50 => compare_kmer_sets::<BitPackedKmer<50, u128>>(unitigs_file, test_tigs_file, config),
        51 => compare_kmer_sets::<BitPackedKmer<51, u128>>(unitigs_file, test_tigs_file, config),
        52 => compare_kmer_sets::<BitPackedKmer<52, u128>>(unitigs_file, test_tigs_file, config),
        53 => compare_kmer_sets::<BitPackedKmer<53, u128>>(unitigs_file, test_tigs_file, config),
        54 => compare_kmer_sets::<BitPackedKmer<54, u128>>(unitigs_file, test_tigs_file, config),
        55 => compare_kmer_sets::<BitPackedKmer<55, u128>>(unitigs_file, test_tigs_file, config),
        56 => compare_kmer_sets::<BitPackedKmer<56, u128>>(unitigs_file, test_tigs_file, config),
        57 => compare_kmer_sets::<BitPackedKmer<57, u128>>(unitigs_file, test_tigs_file, config),
        58 => compare_kmer_sets::<BitPackedKmer<58, u128>>(unitigs_file, test_tigs_file, config),
        59 => compare_kmer_sets::<BitPackedKmer<59, u128>>(unitigs_file, test_tigs_file, config),
        60 => compare_kmer_sets::<BitPackedKmer<60, u128>>(unitigs_file, test_tigs_file, config),
        61 => compare_kmer_sets::<BitPackedKmer<61, u128>>(unitigs_file, test_tigs_file, config),
        62 => compare_kmer_sets::<BitPackedKmer<62, u128>>(unitigs_file, test_tigs_file, config),
        63 => compare_kmer_sets::<BitPackedKmer<63, u128>>(unitigs_file, test_tigs_file, config),
        64 => compare_kmer_sets::<BitPackedKmer<64, u128>>(unitigs_file, test_tigs_file, config),
        _ => compare_kmer_sets::<BitPackedVectorKmer>(unitigs_file, test_tigs_file, config),
    }
}

#[cfg(test)]
mod tests {
    use crate::{compare_kmer_sets, initialise_logging, BitPackedKmer, Config};
    use log::LevelFilter;

    #[test]
    fn test_simple() {
        initialise_logging(LevelFilter::Debug);
        let unitigs = ">a\nTAAACTG";
        let test_tigs = ">\nTAAAC\n>\nCAGT\n";

        let result = compare_kmer_sets::<BitPackedKmer<3, u8>>(
            unitigs.as_bytes(),
            test_tigs.as_bytes(),
            Config {
                log_level: LevelFilter::Debug,
                k: 3,
                do_not_verify: false,
                panic_on_parse_error: true,
                allow_cuttlefish2_errors: false,
                unitigs: Default::default(),
                test_tigs: Default::default(),
            },
        );

        assert!(result.is_ok(), "Expected ok result, but got {result:?}");
    }

    #[test]
    fn test_self_complemental_node() {
        initialise_logging(LevelFilter::Debug);
        let unitigs = ">a\nTAATTACTG";
        let test_tigs = ">\nTAATTA\n>\nCAGTAA\n";
        assert!(compare_kmer_sets::<BitPackedKmer<4, u8>>(
            unitigs.as_bytes(),
            test_tigs.as_bytes(),
            Config {
                log_level: LevelFilter::Debug,
                k: 4,
                do_not_verify: false,
                panic_on_parse_error: true,
                allow_cuttlefish2_errors: false,
                unitigs: Default::default(),
                test_tigs: Default::default(),
            }
        )
        .is_ok());
    }
}
