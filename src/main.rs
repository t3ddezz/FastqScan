use clap::{Arg, Command};
use std::fs::File;
use flate2::read::GzDecoder;
use std::io::Read;

fn main() {
    let matches = Command::new("FastQScan")
        .author("Leeroy Jenkins")
        .version("-")
        .about("Fast and safe Q/C for FASTQ files")
        .arg(Arg::new("Path_to_file_1").short('1').long("read1").required(true).takes_value(true))
        .arg(Arg::new("Path_to_file_2").short('2').long("read2").required(true).takes_value(true))
        .get_matches();
    //value_of: verwendet, um den Wert des CLI-Arguments Path_to_file_1 aus den matches zu holen
    //unwrap() wird verwendet, um aus Option<&str> direkt ein &str zu machen
    let r1_path = matches.value_of("Path_to_file_1").unwrap();
    let reader_r1 = open_fasta_file(r1_path);
    let avg_qualities_r1 = average_base_quality(reader_r1);
    println!("{:?}", avg_qualities_r1);
    
}

fn average_base_quality<R: Read>(mut reader: R) -> Vec<f64> {
    let mut buffer = String::new();
    reader.read_to_string(&mut buffer).expect("Fehler beim Lesen der Datei");

    let mut total_quality: Vec<u64> = Vec::new();
    let mut read_count: u64 = 0;
    //buffer.lines() gibt jede Zeile als Iterator zurück
    //enumerate() fügt einen Index i hinzu, der zählt, welche Zeile gerade gelesen wird
    //i % 4 == 3 prüft, ob die Zeile eine Qualitätszeile ist (weil FASTQ 4-Zeilen-Blöcke hat)
    for (i, line) in buffer.lines().enumerate() {
        if i % 4 == 3 { // Qualitätsscore-Zeilen
            read_count += 1;
            for (pos, ch) in line.chars().enumerate() {
                let phred_score = (ch as u8) - 33;
                if pos >= total_quality.len() {
                    total_quality.push(phred_score as u64);
                } else {
                    total_quality[pos] += phred_score as u64;
                }
            }
        }
    }
    total_quality as f64
}

fn open_fasta_file(path: &str)-> GzDecoder<File>{
    let file = File::open(path).expect("Fehler beim Öffnen der Datei");
    GzDecoder::new(file)
}