mod tests {
    use super::*; 

    #[test]
    fn test_parse_valid_record() {
        let data = b"@seq1\nACGT\n+\n!!!!\n";
        let mut cursor = std::io::Cursor::new(data);
        let mut record = FastqRecord::default();
        let result = WorkflowRunner::parse_record(&mut cursor, &mut record); 
        assert!(result.is_ok(), "Expected parse_record to succeed.");
        assert_eq!(record.seq, b"ACGT".to_vec(), "Expected valid sequence.");
        assert_eq!(record.qual, b"!!!!".to_vec(), "Expected valid quality.");
    }

    #[test]
    fn test_parse_invalid_record() {
        let data = b"@seq2\n\n+\n\n"; 
        let mut cursor = std::io::Cursor::new(data);
        let mut record = FastqRecord::default();
        let result = WorkflowRunner::parse_record(&mut cursor, &mut record);
        
        assert!(result.is_ok(), "Expected parse_record to succeed.");
        assert!(record.seq.is_empty(), "Expected empty sequence.");
        assert!(record.qual.is_empty(), "Expected empty quality.");
    }

    #[test]
    fn test_process_with_valid_data() {
        let data = b"@seq1\nACGT\n+\n!!!!\n@seq2\nCGTA\n+\n!!!!\n";
        let mut cursor = std::io::Cursor::new(data);
        let mut runner = WorkflowRunner {
            statistics: vec![
                Box::new(ReadQualityStatistic::new()), 
            ],
        };
        runner.process(&mut cursor);

        
        assert!(!runner.finalize().is_empty(), "Expected statistics to be generated.");
    }
}
