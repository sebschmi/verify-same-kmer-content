use crate::kmer::Kmer;
use crate::kmer_iterator::KmerIterator;
use clap::Parser;
use log::{debug, error, info, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, TermLogger, TerminalMode};
use std::collections::BTreeSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;

mod kmer;
mod kmer_iterator;

pub fn initialise_logging(log_level: LevelFilter) {
    CombinedLogger::init(vec![TermLogger::new(
        log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )])
    .unwrap();

    info!("Logging initialised successfully");
}

#[derive(Parser, Debug)]
pub struct Config {
    #[clap(short, long, default_value = "Info")]
    log_level: LevelFilter,

    #[clap(short)]
    k: usize,

    #[clap(long)]
    panic_on_parse_error: bool,

    #[clap(index = 0)]
    file1: PathBuf,

    #[clap(index = 1)]
    file2: PathBuf,
}

fn compare_kmer_sets<KmerType: FromIterator<u8> + Ord + Copy + Display>(
    input1: impl Read,
    input2: impl Read,
    config: Config,
) {
    let mut kmer_iter1 =
        KmerIterator::<_, KmerType>::new(input1, config.k, config.panic_on_parse_error);
    let mut kmer_iter2 =
        KmerIterator::<_, KmerType>::new(input2, config.k, config.panic_on_parse_error);

    info!("Reading first input file");
    let kmers1: BTreeSet<_> = kmer_iter1.by_ref().collect();
    let kmers1: Vec<_> = kmers1.into_iter().collect();
    let mut kmers1_visited = vec![false; kmers1.len()];

    // assert that kmers1 is sorted
    debug_assert!(kmers1.windows(2).all(|w| w[0] <= w[1]));

    info!("Reading second input file");
    let mut superfluous_kmers2 = BTreeSet::new();
    let mut kmers2_len = 0;
    for kmer in kmer_iter2.by_ref() {
        kmers2_len += 1;
        if let Ok(index) = kmers1.binary_search(&kmer) {
            kmers1_visited[index] = true;
        } else {
            superfluous_kmers2.insert(kmer);
        }
    }
    let superfluous_kmers2 = superfluous_kmers2;
    let superfluous_kmers1: Vec<_> = kmers1_visited
        .into_iter()
        .enumerate()
        .filter_map(|(index, visited)| if visited { None } else { Some(kmers1[index]) })
        .collect();

    for kmer in &superfluous_kmers1 {
        debug!("File1 contains kmer that is missing in file2: {kmer}");
    }
    for kmer in &superfluous_kmers2 {
        debug!("File2 contains kmer that is missing in file1: {kmer}");
    }

    let file1_sequence_size = kmers1.len() + kmer_iter1.sequence_count() * (config.k - 1);
    let file2_sequence_size = kmers2_len + kmer_iter2.sequence_count() * (config.k - 1);
    let file1_string_count = kmer_iter1.sequence_count();
    let file2_string_count = kmer_iter2.sequence_count();
    let compression_rate = file2_sequence_size as f64 / file1_sequence_size as f64;
    let string_count_rate = file2_string_count as f64 / file1_string_count as f64;
    let unique_kmer_count = kmers1.len();

    std::io::stdout().flush().unwrap();
    std::io::stderr().flush().unwrap();
    println!("ground_truth_size:   {file1_sequence_size}");
    println!("test_size: {file2_sequence_size}");

    println!("ground_truth_str_cnt:   {file1_string_count}");
    println!("test_str_cnt: {file2_string_count}");

    println!("compression_rate: {compression_rate}");
    println!("str_cnt_rate: {string_count_rate}");

    println!("unique_kmer_count: {unique_kmer_count}");
    std::io::stdout().flush().unwrap();
    std::io::stderr().flush().unwrap();

    if superfluous_kmers1.is_empty() && superfluous_kmers2.is_empty() {
        info!("Success!");
    } else if superfluous_kmers1.is_empty() {
        error!("File1 contains kmers that are missing in file2");
    } else if superfluous_kmers2.is_empty() {
        error!("File2 contains kmers that are missing in file1");
    } else {
        error!("File1 and file2 contain kmers that are missing in each other");
    }
}

fn main() {
    let config = Config::parse();
    initialise_logging(config.log_level);

    let file1 = File::open(&config.file1)
        .unwrap_or_else(|_| panic!("--file1 points to a file: {:?}", &config.file1));
    let file2 = File::open(&config.file2)
        .unwrap_or_else(|_| panic!("--file2 points to a file: {:?}", &config.file2));

    match config.k {
        0 => error!("Kmer size cannot be zero"),
        1 => compare_kmer_sets::<Kmer<1, u32>>(file1, file2, config),
        2 => compare_kmer_sets::<Kmer<2, u32>>(file1, file2, config),
        3 => compare_kmer_sets::<Kmer<3, u32>>(file1, file2, config),
        4 => compare_kmer_sets::<Kmer<4, u32>>(file1, file2, config),
        5 => compare_kmer_sets::<Kmer<5, u32>>(file1, file2, config),
        6 => compare_kmer_sets::<Kmer<6, u32>>(file1, file2, config),
        7 => compare_kmer_sets::<Kmer<7, u32>>(file1, file2, config),
        8 => compare_kmer_sets::<Kmer<8, u32>>(file1, file2, config),
        9 => compare_kmer_sets::<Kmer<9, u32>>(file1, file2, config),
        10 => compare_kmer_sets::<Kmer<10, u32>>(file1, file2, config),
        11 => compare_kmer_sets::<Kmer<11, u32>>(file1, file2, config),
        12 => compare_kmer_sets::<Kmer<12, u32>>(file1, file2, config),
        13 => compare_kmer_sets::<Kmer<13, u32>>(file1, file2, config),
        14 => compare_kmer_sets::<Kmer<14, u32>>(file1, file2, config),
        15 => compare_kmer_sets::<Kmer<15, u32>>(file1, file2, config),
        16 => compare_kmer_sets::<Kmer<16, u32>>(file1, file2, config),
        17 => compare_kmer_sets::<Kmer<17, u64>>(file1, file2, config),
        18 => compare_kmer_sets::<Kmer<18, u64>>(file1, file2, config),
        19 => compare_kmer_sets::<Kmer<19, u64>>(file1, file2, config),
        20 => compare_kmer_sets::<Kmer<20, u64>>(file1, file2, config),
        21 => compare_kmer_sets::<Kmer<21, u64>>(file1, file2, config),
        22 => compare_kmer_sets::<Kmer<22, u64>>(file1, file2, config),
        23 => compare_kmer_sets::<Kmer<23, u64>>(file1, file2, config),
        24 => compare_kmer_sets::<Kmer<24, u64>>(file1, file2, config),
        25 => compare_kmer_sets::<Kmer<25, u64>>(file1, file2, config),
        26 => compare_kmer_sets::<Kmer<26, u64>>(file1, file2, config),
        27 => compare_kmer_sets::<Kmer<27, u64>>(file1, file2, config),
        28 => compare_kmer_sets::<Kmer<28, u64>>(file1, file2, config),
        29 => compare_kmer_sets::<Kmer<29, u64>>(file1, file2, config),
        30 => compare_kmer_sets::<Kmer<30, u64>>(file1, file2, config),
        31 => compare_kmer_sets::<Kmer<31, u64>>(file1, file2, config),
        32 => compare_kmer_sets::<Kmer<32, u64>>(file1, file2, config),
        33 => compare_kmer_sets::<Kmer<33, u128>>(file1, file2, config),
        34 => compare_kmer_sets::<Kmer<34, u128>>(file1, file2, config),
        35 => compare_kmer_sets::<Kmer<35, u128>>(file1, file2, config),
        36 => compare_kmer_sets::<Kmer<36, u128>>(file1, file2, config),
        37 => compare_kmer_sets::<Kmer<37, u128>>(file1, file2, config),
        38 => compare_kmer_sets::<Kmer<38, u128>>(file1, file2, config),
        39 => compare_kmer_sets::<Kmer<39, u128>>(file1, file2, config),
        40 => compare_kmer_sets::<Kmer<40, u128>>(file1, file2, config),
        41 => compare_kmer_sets::<Kmer<41, u128>>(file1, file2, config),
        42 => compare_kmer_sets::<Kmer<42, u128>>(file1, file2, config),
        43 => compare_kmer_sets::<Kmer<43, u128>>(file1, file2, config),
        44 => compare_kmer_sets::<Kmer<44, u128>>(file1, file2, config),
        45 => compare_kmer_sets::<Kmer<45, u128>>(file1, file2, config),
        46 => compare_kmer_sets::<Kmer<46, u128>>(file1, file2, config),
        47 => compare_kmer_sets::<Kmer<47, u128>>(file1, file2, config),
        48 => compare_kmer_sets::<Kmer<48, u128>>(file1, file2, config),
        49 => compare_kmer_sets::<Kmer<49, u128>>(file1, file2, config),
        50 => compare_kmer_sets::<Kmer<50, u128>>(file1, file2, config),
        51 => compare_kmer_sets::<Kmer<51, u128>>(file1, file2, config),
        52 => compare_kmer_sets::<Kmer<52, u128>>(file1, file2, config),
        53 => compare_kmer_sets::<Kmer<53, u128>>(file1, file2, config),
        54 => compare_kmer_sets::<Kmer<54, u128>>(file1, file2, config),
        55 => compare_kmer_sets::<Kmer<55, u128>>(file1, file2, config),
        56 => compare_kmer_sets::<Kmer<56, u128>>(file1, file2, config),
        57 => compare_kmer_sets::<Kmer<57, u128>>(file1, file2, config),
        58 => compare_kmer_sets::<Kmer<58, u128>>(file1, file2, config),
        59 => compare_kmer_sets::<Kmer<59, u128>>(file1, file2, config),
        60 => compare_kmer_sets::<Kmer<60, u128>>(file1, file2, config),
        61 => compare_kmer_sets::<Kmer<61, u128>>(file1, file2, config),
        62 => compare_kmer_sets::<Kmer<62, u128>>(file1, file2, config),
        63 => compare_kmer_sets::<Kmer<63, u128>>(file1, file2, config),
        64 => compare_kmer_sets::<Kmer<64, u128>>(file1, file2, config),
        other => error!("Unsupported kmer size: {} > 64", other),
    }
}
