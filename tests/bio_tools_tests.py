#!/usr/bin/env python3
"""
Unit tests for bio_tool module.
"""

import unittest
import os
import tempfile
import shutil
import logging
import sys
import time
from io import StringIO
from unittest.mock import patch
from Bio import SeqIO

# Import colorama for colored output
try:
    from colorama import init, Fore, Style
    COLORAMA_AVAILABLE = True
    init()  # Initialize colorama
except ImportError:
    COLORAMA_AVAILABLE = False

# Import the module
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    from bio_tool import (
        BiologicalSequence, DNASequence, RNASequence, AminoAcidSequence,
        filter_fastq, run_sequence_operation, setup_logging
    )
except ImportError:
    print("Could not import bio_tool module. Make sure it's in your PYTHONPATH.")
    sys.exit(1)


class TestResult(unittest.TestResult):
    """Custom TestResult with progress indicators."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.timeTaken = 0
        self.current_test = 0
        self.total_tests = 0
    
    def startTest(self, test):
        """Called when a test begins."""
        super().startTest(test)
        self.current_test += 1
        test_name = test.id().split('.')[-1]
        
        # Progress indicator with color
        if COLORAMA_AVAILABLE:
            progress = f"[{self.current_test}/{self.total_tests}] "
            print(f"{Fore.CYAN}{progress}{Style.RESET_ALL}{test_name}... ", end="")
        else:
            progress = f"[{self.current_test}/{self.total_tests}] "
            print(f"{progress}{test_name}... ", end="")
        
        sys.stdout.flush()
    
    def addSuccess(self, test):
        """Called when a test succeeds."""
        super().addSuccess(test)
        if COLORAMA_AVAILABLE:
            print(f"{Fore.GREEN}PASS{Style.RESET_ALL}")
        else:
            print("PASS")
    
    def addError(self, test, err):
        """Called when a test raises an unexpected exception."""
        super().addError(test, err)
        if COLORAMA_AVAILABLE:
            print(f"{Fore.RED}ERROR{Style.RESET_ALL}")
        else:
            print("ERROR")
    
    def addFailure(self, test, err):
        """Called when a test fails."""
        super().addFailure(test, err)
        if COLORAMA_AVAILABLE:
            print(f"{Fore.RED}FAIL{Style.RESET_ALL}")
        else:
            print("FAIL")
    
    def addSkip(self, test, reason):
        """Called when a test is skipped."""
        super().addSkip(test, reason)
        if COLORAMA_AVAILABLE:
            print(f"{Fore.YELLOW}SKIP{Style.RESET_ALL} ({reason})")
        else:
            print(f"SKIP ({reason})")


# Enhanced test runner with progress indicators
class ColorTextTestRunner(unittest.TextTestRunner):
    """Custom test runner with colored output and progress indicators."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.success_count = 0
        self.failure_count = 0
        self.error_count = 0
        self.skipped_count = 0
    
    def _print_colored(self, text, color, bold=False):
        """Print text with specified color if colorama is available."""
        if COLORAMA_AVAILABLE:
            style = Style.BRIGHT if bold else ""
            print(f"{style}{color}{text}{Style.RESET_ALL}")
        else:
            print(text)
    
    def run(self, test):
        """Run the test suite with progress indicators."""
        print("\n" + "=" * 70)
        if COLORAMA_AVAILABLE:
            self._print_colored("BIO TOOLS TEST SUITE", Fore.CYAN, bold=True)
        else:
            print("BIO TOOLS TEST SUITE")
        print("=" * 70)
        
        result = super().run(test)
        
        print("\n" + "-" * 70)
        if COLORAMA_AVAILABLE:
            if result.wasSuccessful():
                self._print_colored("ALL TESTS PASSED", Fore.GREEN, bold=True)
            else:
                self._print_colored("TESTS FAILED", Fore.RED, bold=True)
        else:
            if result.wasSuccessful():
                print("ALL TESTS PASSED")
            else:
                print("TESTS FAILED")
        
        # Print summary
        total = result.testsRun
        success = total - len(result.failures) - len(result.errors) - len(result.skipped)
        
        print(f"\nRan {total} tests in {result.timeTaken:.3f}s")
        
        # Display statistics with colors if available
        if COLORAMA_AVAILABLE:
            success_msg = f"Success: {success}/{total}"
            failure_msg = f"Failures: {len(result.failures)}"
            error_msg = f"Errors: {len(result.errors)}"
            skipped_msg = f"Skipped: {len(result.skipped)}"
            
            self._print_colored(success_msg, Fore.GREEN)
            if result.failures:
                self._print_colored(failure_msg, Fore.RED)
            else:
                print(failure_msg)
            if result.errors:
                self._print_colored(error_msg, Fore.RED)
            else:
                print(error_msg)
            if result.skipped:
                self._print_colored(skipped_msg, Fore.YELLOW)
            else:
                print(skipped_msg)
        else:
            print(f"Success: {success}/{total}")
            print(f"Failures: {len(result.failures)}")
            print(f"Errors: {len(result.errors)}")
            print(f"Skipped: {len(result.skipped)}")
        
        print("-" * 70)
        return result


class TestDNASequence(unittest.TestCase):
    """Tests for DNASequence class."""
    
    def test_valid_dna_sequence(self):
        """Test initialization with valid DNA sequence."""
        dna = DNASequence("ATCG")
        self.assertEqual(str(dna), "ATCG")
        
    def test_invalid_dna_sequence(self):
        """Test initialization with invalid DNA sequence raises error."""
        with self.assertRaises(ValueError):
            DNASequence("AXYZ")
    
    def test_reverse_complement(self):
        """Test reverse complement operation."""
        dna = DNASequence("ATCG")
        rev_comp = dna.reverse_complement()
        self.assertEqual(str(rev_comp), "CGAT")
        
    def test_transcribe(self):
        """Test transcription to RNA."""
        dna = DNASequence("ATCG")
        rna = dna.transcribe()
        self.assertEqual(str(rna), "AUCG")
        self.assertIsInstance(rna, RNASequence)


class TestAminoAcidSequence(unittest.TestCase):
    """Tests for AminoAcidSequence class."""
    
    def test_amino_acid_composition(self):
        """Test amino acid composition calculation."""
        aa_seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
        composition = aa_seq.amino_acid_composition()
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            self.assertEqual(composition.get(aa, 0), 1)


class TestFASTQFiltering(unittest.TestCase):
    """Tests for FASTQ filtering functionality."""
    
    def setUp(self):
        """Set up temporary directory for test files."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a test FASTQ file
        self.fastq_content = """@seq1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@seq2
GCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIII
"""
        self.test_fastq = os.path.join(self.temp_dir, "test.fastq")
        with open(self.test_fastq, "w") as f:
            f.write(self.fastq_content)
        
        # Create a test output path
        self.output_fastq = os.path.join(self.temp_dir, "output.fastq")
    
    def tearDown(self):
        """Remove temporary directory and files."""
        shutil.rmtree(self.temp_dir)
    
    def test_filter_fastq_file(self):
        """Test filtering FASTQ file and writing output."""
        # Test with default parameters (should keep all sequences)
        result = filter_fastq(
            self.test_fastq,
            self.output_fastq,
            gc_bounds=(0.0, 1.0),
            length_bounds=(0, 100),
            quality_threshold=0
        )
        
        # Check that output file exists
        self.assertTrue(os.path.exists(self.output_fastq))
        
        # Check content (should match input with both sequences)
        with open(self.output_fastq, "r") as f:
            content = f.read()
        self.assertIn("@seq1", content)
        self.assertIn("@seq2", content)
        
        # Check with proper file closing
        record_count = 0
        with open(self.output_fastq, "r") as handle:
            record_count = sum(1 for _ in SeqIO.parse(handle, "fastq"))
        self.assertEqual(record_count, 2)
    
    def test_filter_fastq_gc_content(self):
        """Test filtering FASTQ by GC content."""
        # Create test with different GC content sequences
        gc_rich_fastq = """@seq1
ATCGATCGATCGATCG
+
IIIIIIIIIIIIIIII
@seq2
GGGGCCCCGGGGCCCC
+
IIIIIIIIIIIIIIII
"""
        test_file = os.path.join(self.temp_dir, "gc_test.fastq")
        with open(test_file, "w") as f:
            f.write(gc_rich_fastq)
            
        output_file = os.path.join(self.temp_dir, "gc_output.fastq")
        
        # Filter for high GC content (>75%)
        result = filter_fastq(
            test_file,
            output_file,
            gc_bounds=(0.75, 1.0),  # Only sequences with GC content >= 75%
            length_bounds=(0, 100),
            quality_threshold=0
        )
        
        # Check content (should only have seq2 with 100% GC)
        with open(output_file, "r") as f:
            content = f.read()
        self.assertNotIn("@seq1", content)  # seq1 has 50% GC content
        self.assertIn("@seq2", content)     # seq2 has 100% GC content
        
        # Additional check: count records using SeqIO.parse with proper closing
        record_count = 0
        with open(output_file, "r") as handle:
            record_count = sum(1 for _ in SeqIO.parse(handle, "fastq"))
        self.assertEqual(record_count, 1, 
                         "Expected exactly one record in filtered output")
    
    def test_file_exists_error(self):
        """Test that FileExistsError is raised when output file exists."""
        # Create the output file first
        test_output = os.path.join(self.temp_dir, "existing_output.fastq")
        with open(test_output, "w") as f:
            f.write("test")
        
        # Attempt to filter without overwrite flag
        with self.assertRaises(FileExistsError):
            filter_fastq(
                self.test_fastq,
                test_output
            )


class TestLogging(unittest.TestCase):
    """Tests for logging functionality."""
    
    def setUp(self):
        """Set up temporary log file."""
        self.log_file = os.path.join(tempfile.gettempdir(), "test_log.log")
        # Remove the log file if it exists
        if os.path.exists(self.log_file):
            os.remove(self.log_file)
    
    def tearDown(self):
        """Remove temporary log file."""
        if os.path.exists(self.log_file):
            os.remove(self.log_file)
    
    def test_logging_setup(self):
        """Test that logging is set up correctly and writes to file."""
        # Reset any existing logger configuration
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        
        # Set up logging with absolute path
        setup_logging(self.log_file, verbosity=1)
        
        # Log a test message
        logging.info("Test info message")
        logging.error("Test error message")
        
        # Force logger to flush handlers
        for handler in logging.root.handlers:
            handler.flush()
        
        # Check that the log file exists and contains the messages
        self.assertTrue(os.path.exists(self.log_file), 
                       f"Log file not found at {self.log_file}")
        
        with open(self.log_file, "r") as f:
            log_content = f.read()
        
        self.assertIn("Test info message", log_content)
        self.assertIn("Test error message", log_content)
        
        # Additional check: verify log levels are correctly recorded
        self.assertIn("INFO", log_content, 
                   "INFO log level not found in log file")
        self.assertIn("ERROR", log_content, 
                      "ERROR log level not found in log file")
        
        # Make sure log entries are properly formatted with levels
        info_line = next((line for line in log_content.split('\n') 
                          if "Test info message" in line), "")
        error_line = next((line for line in log_content.split('\n') 
                           if "Test error message" in line), "")
        
        self.assertIn("INFO", info_line, 
                   "INFO log level missing from info message line")
        self.assertIn("ERROR", error_line, 
                      "ERROR log level missing from error message line")


if __name__ == "__main__":
    # Create test suite
    test_loader = unittest.TestLoader()
    test_suite = test_loader.loadTestsFromTestCase(TestDNASequence)
    test_suite.addTests(test_loader.loadTestsFromTestCase(TestAminoAcidSequence))
    test_suite.addTests(test_loader.loadTestsFromTestCase(TestFASTQFiltering))
    test_suite.addTests(test_loader.loadTestsFromTestCase(TestLogging))
    
    # Count total tests
    test_count = test_suite.countTestCases()
    
    # Run with custom test runner
    runner = ColorTextTestRunner(verbosity=2)
    result = TestResult()
    result.total_tests = test_count
    
    # Initialize progress bar
    start_time = time.time()
    
    # Run tests
    test_suite.run(result)
    
    # Calculate time taken
    result.timeTaken = time.time() - start_time
    
    # Print final results
    runner.run = lambda x: result
    runner.run(test_suite)