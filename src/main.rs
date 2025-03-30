mod runners; // Hier fügen wir das runners Modul ein

use clap::{Arg, Command};
use std::fs::File;
use std::io::{BufReader};
use flate2::read::GzDecoder; // Importiere den GzDecoder für das Entpacken
use runners::{WorkflowRunner, BaseQualityPosStatistic, ReadQualityStatistic, BaseCompositionStatistic, GcContentStatistic, LengthDistributionStatistic};
use serde_json::to_string_pretty; // Für die hübsche JSON-Ausgabe

fn main() {
    // Command-line Argument Parsing
    let matches = Command::new("FastQScan") 
        .author("Leeroy Jenkins") 
        .version("-") 
        .about("Fast and safe Q/C for FASTQ files") 
        .arg(Arg::new("read1").short('1').long("read1").required(true).value_name("FILE"))
        .arg(Arg::new("read2").short('2').long("read2").required(true).value_name("FILE"))
        .get_matches(); 
    
    let r1_path = matches.get_one::<String>("read1").expect("Missing read1 file");
    let r2_path = matches.get_one::<String>("read2").expect("Missing read2 file");

    println!("Verarbeite Datei 1: {}", r1_path);
    println!("Verarbeite Datei 2: {}", r2_path);

    // Öffnen der GZIP-komprimierten FASTQ-Dateien
    let r1_file = File::open(r1_path).expect("Fehler beim Öffnen der Datei 1");
    let r2_file = File::open(r2_path).expect("Fehler beim Öffnen der Datei 2");
    
    let r1_decoder = GzDecoder::new(r1_file); // GZIP-Decoder für Datei 1
    let r2_decoder = GzDecoder::new(r2_file); // GZIP-Decoder für Datei 2
    
    let r1_reader = BufReader::new(r1_decoder); // Entpackten Reader für Datei 1
    let r2_reader = BufReader::new(r2_decoder); // Entpackten Reader für Datei 2

    // Debug-Ausgabe, um sicherzustellen, dass das Entpacken funktioniert
    println!("Dateien erfolgreich entpackt, beginne mit dem Lesen der Dateien...");

    // Erstelle den WorkflowRunner mit den gewünschten Statistiken
    let mut runner = WorkflowRunner {
        statistics: vec![
            Box::new(BaseQualityPosStatistic::new()), // Aufruf des neuen Konstruktors
            Box::new(ReadQualityStatistic::new()),
            Box::new(BaseCompositionStatistic::new()),
            Box::new(GcContentStatistic::new()),
            Box::new(LengthDistributionStatistic::new()),
        ],
    };

    // Verarbeite die erste Datei mit dem Runner
    println!("Verarbeitung der Datei 1 gestartet...");
    runner.process(r1_reader);
    
    // Verarbeite die zweite Datei mit dem Runner
    println!("Verarbeitung der Datei 2 gestartet...");
    runner.process(r2_reader);

    // Ergebnisse abrufen und als JSON ausgeben
    let results = runner.finalize();

    // Alle Statistiken im JSON-Format ausgeben
    let mut json_output = Vec::new();
    for statistic in results {
        let report = statistic.report();
        json_output.push(report);
    }

    // Die JSON-Daten ausgeben
    match to_string_pretty(&json_output) {
        Ok(json) => println!("{}", json),
        Err(e) => eprintln!("Fehler beim Serialisieren der JSON-Daten: {}", e),
    }

    println!("Verarbeitung abgeschlossen.");
}
