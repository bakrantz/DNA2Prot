use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8"; #Make a path to the module
use Dna2prot;

my $dna2prot = Dna2prot
          -> new
          -> align_and_join_files('PA_DNA_read1.txt', 'PA_DNA_read2.txt');

$dna2prot -> find_largest_orf
          -> save_dna('PA_DNA_joined.txt');
$dna2prot -> protein_name('PA');
$dna2prot -> dna_name('PA_DNA');
$dna2prot -> save_record('PA_DNA_record_3.txt');
