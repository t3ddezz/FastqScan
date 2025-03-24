use clap::{Arg, Command}; 
use std::fs::File; 
use std::io::Read; 
use flate2::read::GzDecoder; 
use std::io::{BufRead, BufReader};
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
    let reader_r1_for_qual = open_fasta_file(r1_path);
    let reader_r1_for_base = open_fasta_file(r1_path);

    let (avg_qualities_r1, avg_quality_r1) = average_base_quality(reader_r1_for_qual);

    let mut buf_reader = BufReader::new(reader_r1_for_base);
    let base_counts_r1 = count_base_composition(&mut buf_reader);


    println!("Average per position: {:?}", avg_qualities_r1);
    println!("Average quality of all reads: {:.2}", avg_quality_r1);
}

// Diese Funktion öffnet eine gzip-komprimierte FASTQ-Datei 
fn open_fasta_file(path: &str) -> GzDecoder<File> { 
    let file = File::open(path).expect("Fehler beim Öffnen der Datei"); 
    GzDecoder::new(file) // Entpackt die gz-Datei direkt beim Lesen 
}


fn average_base_quality<R: Read>(reader: R) -> (Vec<f64>, f64) {
    let buffer = BufReader::new(reader);

    let mut total_quality: Vec<u64> = Vec::new();
    let mut read_count: u64 = 0;
    let mut total_quality_sum: u64 = 0;
    let mut total_base_count: u64 = 0;

    for (i, line_result) in buffer.lines().enumerate() {
        let line = line_result.expect("Fehler beim Lesen der Zeile");

        if i % 4 == 3 {
            read_count += 1;

            for (pos, ch) in line.chars().enumerate() {
                let phred_score = (ch as u8) - 33;

                if pos >= total_quality.len() {
                    total_quality.push(phred_score as u64);
                } else {
                    total_quality[pos] += phred_score as u64;
                }

                total_quality_sum += phred_score as u64;
                total_base_count += 1;
            }
        }
    }

    let avg_per_position = total_quality
        .iter()
        .map(|&sum| sum as f64 / read_count as f64)
        .collect();

    let avg_quality_all_reads = total_quality_sum as f64 / total_base_count as f64;

    (avg_per_position, avg_quality_all_reads)
}

// Funktion zur Zählung der Basenzusammensetzung (A, C, G, T, N) an jeder Position in Sequenzen aus einem FASTQ-Format
fn count_base_composition<R: BufRead>(reader: &mut R) -> Vec<[u64; 5]> {
    // Initialisiert einen Vektor für die Zählungen: 
    // Für jede Position in den Sequenzen wird ein Array [A, C, G, T, N] gehalten
    let mut base_counts: Vec<[u64; 5]> = Vec::new();
    let mut line = String::new();
    let mut line_index = 0;

    // Solange noch Zeilen gelesen werden können
    // unwrap_or(0) damit kein panic 
    while reader.read_line(&mut line).unwrap_or(0) > 0 {
        // Entfernt Zeilenumbruch am Ende der Zeile
        let trimmed = line.trim_end();
        
        if line_index % 4 == 1 {
            // Iteration über die Zeichen der Sequenz und deren Position
            for (pos, ch) in trimmed.chars().enumerate() {
                // Bestimmen, welche Base es ist (Index: 0=A, 1=C, 2=G, 3=T, 4=N)
                let base_index = match ch {
                    'A' | 'a' => 0,
                    'C' | 'c' => 1,
                    'G' | 'g' => 2,
                    'T' | 't' => 3,
                    'N' | 'n' => 4,
                    _ => continue, // Ignoriert andere Zeichen
                };

                // Wenn die aktuelle Position noch nicht im Vektor existiert, erweitere ihn
                if pos >= base_counts.len() {
                    base_counts.push([0; 5]); // Initialisiert Zählarray für neue Position
                }

                // Erhöhe den Zähler für die jeweilige Base an der aktuellen Position
                base_counts[pos][base_index] += 1;
            }
        }

        // Zeilenpuffer leeren für nächste Zeile
        line.clear();
        line_index += 1;
    }

    base_counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_base_composition_simple() {
        let fastq_data = b"@read1\nACGT\n+\n!!!!\n@read2\nTGCA\n+\n####\n";

        
        let mut reader = Cursor::new(fastq_data);

        let result = count_base_composition(&mut reader);

        
        assert_eq!(result[0], [1, 0, 0, 1, 0]);
        
        assert_eq!(result[1], [0, 1, 1, 0, 0]);
        
        assert_eq!(result[2], [0, 1, 1, 0, 0]);
        
        assert_eq!(result[3], [1, 0, 0, 1, 0]);
    }
}
