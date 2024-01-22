use clap::{Parser, Subcommand};
use bio::io::fastq;
use bio::pattern_matching::myers::Myers;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Split { 
        /*
        #[arg(short, long, value_name = "prefix")]
        prefix: String,
        #[arg(short, long, value_name = "prefix")]
        suffix: String,
        */
        fastq: String
        },
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Split { fastq } => {
            split_reads(fastq);
        }
    }
}

fn filter_overlaps(occ: Vec<(usize, usize, u8)>) -> Vec<(usize, usize, u8)> {

    let mut out = Vec::<(usize, usize, u8)>::new();
    // lazy naive method n^2, discard any hits that overlap a hit with a lower edit distance
    for a in occ {

        let mut add = true;
        for b in out.iter_mut() {
            // check overlap
            let overlaps = a.0 <= b.1 && b.0 <= a.1;
            if overlaps {
                if a.2 <= b.2 { 
                    // better hit, replace b
                    *b = a;
                }

                add = false; 
                break;
            }
        }

        if add { out.push(a); } // no overlaps, add to collection
    }

    return out;
}

fn display_matches(patterns: &mut [(Myers<u128>, Myers<u128>); 2], strand_seq: &[[&str; 2]; 2], seq: &[u8], max_ed: u8) -> () {

    let mut match_strings = Vec::new();
    for si in 0..=1 {
        let mut strand_matches = vec![' '; seq.len()];

        let occ_prefix: Vec<_> = filter_overlaps(patterns[si].0.find_all(seq, max_ed).collect());
        let occ_suffix: Vec<_> = filter_overlaps(patterns[si].1.find_all(seq, max_ed).collect());

        for hit in &occ_prefix {

            let mut e = hit.0 + strand_seq[si][0].len(); 
            e = if e > strand_matches.len() { strand_matches.len() } else { e };
            strand_matches.splice(hit.0 .. e, strand_seq[si][0].chars());
        }

        for hit in &occ_suffix {
            let mut e = hit.0 + strand_seq[si][1].len(); 
            e = if e > strand_matches.len() { strand_matches.len() } else { e };
            strand_matches.splice(hit.0 .. e, strand_seq[si][1].chars());
        }

        let sm: String = strand_matches.into_iter().collect();
        println!("POC: {:?}", occ_prefix);
        println!("SOC: {:?}", occ_suffix);
        match_strings.push(sm);
    }

    let pp_seq = String::from_utf8_lossy(seq);
    let l = pp_seq.len();
    let step = 80;
    for i in (0..l).step_by(step) {
        let m = if i + step >= l { l - 1 } else { i + step };
        println!("S1: {}", &match_strings[0][i..m]);
        println!("FQ: {}", &pp_seq[i..m]);
        println!("S2: {}", &match_strings[1][i..m]);
        println!("");
    }
}

fn split_reads(input_fastq: &str) {
    let reader = fastq::Reader::from_file(input_fastq).expect("Could not open input fastq");
    let max_ed = 14;
    let display = false;


    let strand_seq = [ [ "GAATCCTCGGATTCCATGATCGTTACATGATTTTCTGTTGGTGCTGATATTGC", 
                         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAACTTTCTTGTCAC" ],
                       [ "GTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG", 
                         "GCAATATCAGCACCAACAGAAAATCATGTAACGATCATGGAATCCGAGGATTC" ] ];

    let mut strand_pattern = [ (Myers::<u128>::new(strand_seq[0][0].bytes()), Myers::<u128>::new(strand_seq[0][1].bytes())),
                               (Myers::<u128>::new(strand_seq[1][0].bytes()), Myers::<u128>::new(strand_seq[1][1].bytes())) ];
    
    for r in reader.records() {
        if let Ok(record) = r {
            let seq = record.seq();
            let qual = record.qual();

            if display {
                display_matches(&mut strand_pattern, &strand_seq, &seq, max_ed);
            } else {
                let mut subs = Vec::new();

                for si in 0..=1 {

                    let occ_prefix: Vec<_> = filter_overlaps(strand_pattern[si].0.find_all(seq, max_ed).collect());
                    let occ_suffix: Vec<_> = filter_overlaps(strand_pattern[si].1.find_all(seq, max_ed).collect());

                    for idx in 0..occ_prefix.len() {
                        let start = occ_prefix[idx];
                        let end = if idx < occ_suffix.len() { Some(occ_suffix[idx]) } else { None };
                        subs.push( (start, end, si) );
                    }
                }

                subs.sort_by(|a, b| a.0.0.cmp( & b.0.0 ));
                let mut count = 0;
                for (start, end, strand) in subs {

                    let s = start.1 + 1; // endpoint of front match
                    let e = if let Some(t) = end { t.0 } else { seq.len() };
                    if e < s { continue; }
                    println!("@{}_{} strand={} start={:?} end={:?}", record.id(), count, strand, start, end);
                    count += 1;
                    let subseq = &seq[s..e];
                    let subqual = &qual[s..e];
                    println!("{}", std::str::from_utf8(&subseq).unwrap());
                    println!("+");
                    println!("{}", std::str::from_utf8(&subqual).unwrap());
                }
            }
        }
    }
}
