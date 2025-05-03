use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $dna  =  'GGGTCGACGTCCGGTTCTACTTCCGGAACTAGCGGGTCAGGAGGAGGATCAACTGGCTCTGGATCGACCTCAGGTTCAACCTCGGGTACTTCCGGATCGACTAGTAGTTCAACAACTACGGGAGGAAGCACTAGCACGGGCTCAACATCGGGAACGAGCACTACCACAACCTCGGGATCGGGCGGAAGCACAACTGGGTCGGGGGGGTCGACATCCACTAGTGGCGGTAGT';
my $dna2prot =  Dna2prot
             -> new()
             -> find_and_replace_self_matches($dna)
             -> save_dna('test_gst_dna_sequence_self_matches.txt')
             -> save_textbox('test_gst_dna_self_matches_alignment.txt');

