use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $dna2prot =  Dna2prot
             -> new();
$dna2prot    -> verbose(1);
$dna2prot    -> open_protein('UBA52.txt')
             -> reverse_translate;
my $dna      = $dna2prot -> current_dna_sequence;
$dna2prot    -> dna_synthesis($dna, 'LCR')
             -> save_textbox('UBA52_synthesis_LCR.txt')

