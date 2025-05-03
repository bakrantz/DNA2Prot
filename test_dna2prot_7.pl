use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;

# Primer has three mismatches (8.6% mismatch) from PA_DNA
my $dna2prot = Dna2prot
  -> new()
  -> open_dna('PA_DNA.txt')
  -> calculate_primer_tm('GTCAAAAATAAAACTACTATTCTTTCACCATGGATT')
  -> save_textbox('PA_primer.txt');
