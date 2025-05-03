use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8"; #Make a path to the module
use Dna2prot;

my $dna2prot1 =  Dna2prot
              -> new()
              -> open_dna('PA_DNA_read1.txt')
              -> save_dna('PA_DNA_read1.txt');

my $dna1 = $dna2prot1 -> current_dna_sequence;

my $dna2prot2 =  Dna2prot
              -> new()
              -> open_dna('PA_DNA_read2.txt')
              -> save_dna('PA_DNA_read2.txt');

my $dna2 = $dna2prot2 -> current_dna_sequence;

my $dna2prot3 =  Dna2prot
              -> new()
              -> align_and_join($dna1, $dna2);

my $aligned_joined_dna = $dna2prot3 -> current_dna_sequence;

$dna2prot3    -> find_largest_orf
              -> save_dna('PA_DNA_joined.txt')
              -> find_largest_orf;
$dna2prot3    -> protein_name('PA');
$dna2prot3    -> dna_name('PA_DNA');
$dna2prot3    -> save_record('PA_DNA_record_2.txt');
