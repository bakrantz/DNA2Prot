use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $dna2prot =  Dna2prot
             -> new();
$dna2prot    -> verbose(1);
$dna2prot    -> open_dna('PA_DNA.txt')
             -> find_largest_orf
             -> reverse_translate
             -> save_dna('PA_reverse_translate.txt');
