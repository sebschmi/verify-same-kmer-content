use crate::kmer::{BitPackedKmer, Kmer};
use crate::kmer_iterator::KmerIterator;
use clap::Parser;
use log::{debug, error, info, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, TermLogger, TerminalMode};
use std::collections::BTreeSet;
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

fn compare_kmer_sets<KmerType: FromIterator<u8> + Ord + Copy + Display + Kmer>(
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
        let mut kmers_unitigs: Vec<_> = kmer_iter_unitigs.by_ref().collect();
        kmers_unitigs.sort_unstable();
        let kmers_unitigs = kmers_unitigs;
        let mut kmers_unitigs_visited = vec![false; kmers_unitigs.len()];

        info!("Reading second input file");
        let mut superfluous_kmers_test_tigs = BTreeSet::new();
        for kmer in kmer_iter_test_tigs.by_ref() {
            let mut found = false;
            if let Ok(index) = kmers_unitigs.binary_search(&kmer) {
                kmers_unitigs_visited[index] = true;
                found = true;
            }
            if let Ok(index) = kmers_unitigs.binary_search(&kmer.reverse_complement()) {
                kmers_unitigs_visited[index] = true;
                found = true;
            }

            if !found {
                superfluous_kmers_test_tigs.insert(kmer);
            }
        }

        let superfluous_kmers_test_tigs = superfluous_kmers_test_tigs;
        let superfluous_kmers_unitigs: Vec<_> = kmers_unitigs_visited
            .into_iter()
            .enumerate()
            .filter_map(|(index, visited)| {
                if visited {
                    None
                } else {
                    Some(kmers_unitigs[index])
                }
            })
            .collect();

        for kmer in &superfluous_kmers_unitigs {
            debug!("Unitigs contain kmer that is missing in test tigs: {kmer}");
        }
        for kmer in &superfluous_kmers_test_tigs {
            debug!("Test tigs contains kmer that is missing in unitigs: {kmer}");
        }

        assert_eq!(
            kmers_unitigs.len(),
            kmer_iter_unitigs.character_count()
                - kmer_iter_unitigs.sequence_count() * (config.k - 1),
            "character_count: {}; sequence_count: {}; k: {}",
            kmer_iter_unitigs.character_count(),
            kmer_iter_unitigs.sequence_count(),
            config.k
        );
        (
            !superfluous_kmers_unitigs.is_empty(),
            !superfluous_kmers_test_tigs.is_empty(),
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
        info!("Success!");
        Ok(())
    } else if !has_superfluous_kmers_unitigs {
        error!("Unitigs contain kmers that are missing in test tigs");
        Err(Error::Mismatch)
    } else if !has_superfluous_kmers_test_tigs {
        error!("Test tigs contains kmers that are missing in unitigs");
        Err(Error::Mismatch)
    } else {
        error!("Unitigs and test tigs contain kmers that are missing in each other");
        Err(Error::Mismatch)
    }
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
        other => {
            error!("Unsupported kmer size: {} > 64", other);
            Err(Error::IllegalKmerSize {
                kmer_size: config.k,
            })
        }
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
        assert!(compare_kmer_sets::<BitPackedKmer<3, u8>>(
            unitigs.as_bytes(),
            test_tigs.as_bytes(),
            Config {
                log_level: LevelFilter::Debug,
                k: 3,
                do_not_verify: false,
                panic_on_parse_error: true,
                unitigs: Default::default(),
                test_tigs: Default::default(),
            }
        )
        .is_ok());
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
                unitigs: Default::default(),
                test_tigs: Default::default(),
            }
        )
        .is_ok());
    }
}
