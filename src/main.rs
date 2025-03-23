use clap::{Arg, Command}; 
use std::fs::File; 
use std::io::Read; 
use flate2::read::GzDecoder; 

fn main() {
    let matches = Command::new("FastQScan") 
        .author("Leeroy Jenkins") 
        .version("-") 
        .about("Fast and safe Q/C for FASTQ files") 
        .arg(Arg::new("read1").short('1').long("read1").required(true).value_name("FILE")) 
        .arg(Arg::new("read2").short('2').long("read2").required(false).value_name("FILE")) 
        .get_matches(); 
    
    // value_of: verwendet, um den Wert des CLI-Arguments "read1" aus den matches zu holen 
    // unwrap() wird verwendet, um aus Option<&str> direkt ein &str zu machen 
    let r1_path = matches.get_one::<String>("read1").expect("Missing read1 file"); 
    let reader_r1 = open_fasta_file(r1_path); 
    
    let (avg_qualities_r1, avg_quality_r1) = compute_quality_metrics(reader_r1);

    println!("Average per position: {:?}", avg_qualities_r1);
    println!("Average quality of all reads: {:.2}", avg_quality_r1);
}

// Diese Funktion öffnet eine gzip-komprimierte FASTQ-Datei 
fn open_fasta_file(path: &str) -> GzDecoder<File> { 
    let file = File::open(path).expect("Fehler beim Öffnen der Datei"); 
    GzDecoder::new(file) // Entpackt die gz-Datei direkt beim Lesen 
}

fn average_base_quality<R: Read>(mut reader: R) -> (Vec<f64>, f64) { 
    let mut buffer = String::new(); 
    reader.read_to_string(&mut buffer).expect("Fehler beim Lesen der Datei"); 
    
    let mut total_quality: Vec<u64> = Vec::new(); 
    let mut read_count: u64 = 0; 
    let mut total_quality_sum: u64 = 0;
    let mut total_base_count: u64 = 0;

    // buffer.lines() gibt jede Zeile als Iterator zurück 
    // enumerate() fügt einen Index i hinzu, der zählt, welche Zeile gerade gelesen wird 
    // i % 4 == 3 prüft, ob die Zeile eine Qualitätszeile ist (weil FASTQ 4-Zeilen-Blöcke hat) 
    for (i, line) in buffer.lines().enumerate() { 
        if i % 4 == 3 { // Qualitätsscore-Zeilen 
            read_count += 1; 
            for (pos, ch) in line.chars().enumerate() { 
                let phred_score = (ch as u8) - 33; // ASCII-Offset für Phred-Qualität 
                if pos >= total_quality.len() { 
                    total_quality.push(phred_score as u64); 
                } else { 
                    total_quality[pos] += phred_score as u64; 
                }
                // Gesamtqualität berechnen
                total_quality_sum += phred_score as u64;
                total_base_count += 1;
            } 
        } 
    } 
    
    // Durchschnitt pro Position berechnen
    let avg_per_position = total_quality_per_pos.iter().map(|&sum| sum as f64 / read_count as f64).collect();

    // Durchschnitt über alle Reads berechnen
    let avg_quality_all_reads = total_quality_sum as f64 / total_base_count as f64;

    (avg_per_position, avg_quality_all_reads)
}
