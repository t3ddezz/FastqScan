use std::io::{self, BufRead};

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub trait Statistic {
    fn process(&mut self, record: &FastqRecord);
    fn report(&self) -> serde_json::Value;
}


pub struct BaseQualityPosStatistic {
    avg_per_position: Vec<f64>, 
}

impl BaseQualityPosStatistic {
    pub fn new() -> Self {
        BaseQualityPosStatistic {
            avg_per_position: Vec::new(),
        }
    }
}

impl Statistic for BaseQualityPosStatistic {
    fn process(&mut self, record: &FastqRecord) {
        let total_quality: Vec<u64> = record.qual.iter().map(|&qual| qual as u64).collect();
        self.avg_per_position = total_quality.iter()
            .map(|&sum| sum as f64 / record.seq.len() as f64)
            .collect(); 
    }

    fn report(&self) -> serde_json::Value {
        let mut result = Vec::new();
        for (pos, avg_quality) in self.avg_per_position.iter().enumerate() {
            result.push(serde_json::json!({
                "position": pos + 1,
                "average_quality": avg_quality,
            }));
        }
        serde_json::json!({ "base_quality_position": result })
    }
}


pub struct ReadQualityStatistic {
    avg_quality_all_reads: f64, 
}

impl ReadQualityStatistic {
    pub fn new() -> Self {
        ReadQualityStatistic {
            avg_quality_all_reads: 0.0,
        }
    }
}

impl Statistic for ReadQualityStatistic {
    fn process(&mut self, record: &FastqRecord) {
        let total_quality: u64 = record.qual.iter().map(|&qual| qual as u64).sum();
        self.avg_quality_all_reads = total_quality as f64 / record.seq.len() as f64; 
    }

    fn report(&self) -> serde_json::Value {
        serde_json::json!({
            "read_quality": self.avg_quality_all_reads
        })
    }
}


pub struct BaseCompositionStatistic {
    proportions: Vec<[f64; 5]>, 
}

impl BaseCompositionStatistic {
    pub fn new() -> Self {
        BaseCompositionStatistic {
            proportions: Vec::new(),
        }
    }

    fn char_to_base_index(ch: char) -> usize {
        match ch {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            'N' | 'n' => 4,
            _ => return 5,
        }
    }
}

impl Statistic for BaseCompositionStatistic {
    fn process(&mut self, record: &FastqRecord) {
        let mut counts = vec![[0; 5]; record.seq.len()]; 

        for (pos, &ch) in record.seq.iter().enumerate() {
            let base_idx = Self::char_to_base_index(ch as char);
            if base_idx < 5 {
                counts[pos][base_idx] += 1;
            }
        }

        self.proportions = counts.into_iter()
            .map(|count| {
                let total: f64 = count.iter().sum::<i32>() as f64;
                let proportion: Vec<f64> = count.iter().map(|&count| count as f64 / total).collect();
                [proportion[0], proportion[1], proportion[2], proportion[3], proportion[4]]
            }).collect();
    }

    fn report(&self) -> serde_json::Value {
        let result: Vec<serde_json::Value> = self.proportions.iter()
            .map(|&proportion| serde_json::json!({
                "A": proportion[0],
                "C": proportion[1],
                "G": proportion[2],
                "T": proportion[3],
                "N": proportion[4],
            }))
            .collect();
        serde_json::json!({ "base_compositions": result })
    }
}


pub struct GcContentStatistic {
    gc_per_position: Vec<f64>,
    gc_per_read: f64,
}

impl GcContentStatistic {
    pub fn new() -> Self {
        GcContentStatistic {
            gc_per_position: Vec::new(),
            gc_per_read: 0.0,
        }
    }

    fn calculate_gc_content(seq: &Vec<u8>) -> f64 {
        let gc_count = seq.iter().filter(|&&base| base == b'G' || base == b'C').count();
        gc_count as f64 / seq.len() as f64
    }
}

impl Statistic for GcContentStatistic {
    fn process(&mut self, record: &FastqRecord) {
        self.gc_per_read = Self::calculate_gc_content(&record.seq);

        self.gc_per_position = record.seq.iter()
            .map(|&base| {
                if base == b'G' || base == b'C' {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
    }

    fn report(&self) -> serde_json::Value {
        serde_json::json!({
            "gc_content_per_position": self.gc_per_position,
            "gc_content_per_read": self.gc_per_read,
        })
    }
}


pub struct LengthDistributionStatistic {
    lengths: Vec<usize>,
}

impl LengthDistributionStatistic {
    pub fn new() -> Self {
        LengthDistributionStatistic {
            lengths: Vec::new(),
        }
    }
}

impl Statistic for LengthDistributionStatistic {
    fn process(&mut self, record: &FastqRecord) {
        self.lengths.push(record.seq.len());
    }

    fn report(&self) -> serde_json::Value {
        serde_json::json!({
            "length_distribution": self.lengths
        })
    }
}

pub struct WorkflowRunner {
    pub statistics: Vec<Box<dyn Statistic>>,
}

impl WorkflowRunner {
    pub fn process<R>(&mut self, mut read: R)
    where
        R: BufRead,
    {
        let mut record = FastqRecord::default();
        let mut count = 0;
        let mut error_count = 0; 
    
        println!("Starte die Verarbeitung der FASTQ-Datei...");
    
        while let Ok(true) = WorkflowRunner::parse_record(&mut read, &mut record) {
            if record.seq.is_empty() || record.qual.is_empty() {
                println!("Fehler: Leere Sequenz oder Qualitätswerte gefunden.");
                error_count += 1;
                if error_count > 5 { 
                    println!("Zu viele fehlerhafte Datensätze gefunden. Beende die Verarbeitung.");
                    break;
                }
                continue; 
            }
    
            count += 1;
            if count % 100 == 0 {
                println!("Verarbeite den {}. Record", count);
            }
            println!("Verarbeite Record: {:?}", record);
    
            for statistic in self.statistics.iter_mut() {
                statistic.process(&record);
            }
        }
    
        println!("Verarbeitung der Datei abgeschlossen.");
    }
    


    
    pub fn parse_record<R>(read: &mut R, record: &mut FastqRecord) -> io::Result<bool>
    where
        R: BufRead,
    {
        let mut line = String::new();

        
        println!("Beginne mit dem Einlesen der Zeilen...");

        
        read.read_line(&mut line)?;
        line.clear(); 

        
        read.read_line(&mut line)?;
        if line.trim().is_empty() {
            println!("Fehler: Leere Sequenz gefunden.");
            return Ok(false); 
        }
        record.seq = line.trim().bytes().collect(); 
        println!("Sequenz: {:?}", record.seq); 

        line.clear(); 

        
        read.read_line(&mut line)?;
        line.clear(); 

        
        read.read_line(&mut line)?;
        if line.trim().is_empty() {
            println!("Fehler: Leere Qualitätswerte gefunden.");
            return Ok(false); 
        }
        record.qual = line.trim().bytes().collect(); 
        println!("Qualitätswerte: {:?}", record.qual); 

        Ok(true) 
    }


    pub fn finalize(self) -> Vec<Box<dyn Statistic>> {
        self.statistics
    }
}


