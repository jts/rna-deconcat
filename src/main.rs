use clap::{Parser, Subcommand};
use bio::io::fasta;
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
        barcode: String,
        fastq: String
        },
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Split { barcode, fastq } => {
            split_reads(barcode, fastq);
        }
    }
}

fn get_best_match(occ: Vec<(usize, usize, u8)>) -> Option<(usize, usize, u8)> {

    let mut out: Option<(usize, usize, u8)> = None;

    for a in occ {
        if out.is_none() || a.2 < out.unwrap().2 {
            out = Some(a);
        }
    }

    return out;
}

fn split_reads(input_barcode: &str, input_fastq: &str) {

    let max_ed = 4;
    let barcode_reader = fasta::Reader::from_file(input_barcode).expect("Could not open input fasta");
    let mut patterns = Vec::new();
    for r in barcode_reader.records() {
        if let Ok(record) = r {
            patterns.push( (record.id().to_owned(), Myers::<u64>::new(record.seq()) ));
        }
    }

    let reader = fastq::Reader::from_file(input_fastq).expect("Could not open input fastq");
    
    for r in reader.records() {
        if let Ok(record) = r {

            let seq = record.seq();
            let mut hits = Vec::new();
            
            // find lowest distance hit between each barcode and read sequence
            for p in patterns.iter_mut() {
                hits.push( get_best_match(p.1.find_all(seq, max_ed).collect()) );
            }

            // parse vector of hits to find best (lowest edit distance)
            let mut best_barcode: Option<usize> = None;
            let mut best_hit: Option<(usize, usize, u8)>= None;
            for (idx, h) in hits.iter().enumerate() {
                if h.is_some() {
                    let t = h.unwrap();
                    let ed = t.2;
                    if best_hit.is_none() || ed < best_hit.unwrap().2 {
                        best_barcode = Some(idx);
                        best_hit = Some(t);
                    }
                }
            }

            if best_hit.is_some() {
                let b = best_hit.unwrap();
                println!("{}\t{}\t{}\t{}\t{}", record.id(), patterns[best_barcode.unwrap()].0, b.2, b.0, b.1); 
            } else {
                println!("{}\t-\t-\t-\t-", record.id());
            }
        }
    }
}
