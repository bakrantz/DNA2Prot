use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;

my $dna2prot = Dna2prot
  -> new()
  -> open_dna('PA_DNA.txt')
  -> find_largest_orf();
 $dna2prot -> protein_name('PA'); 
 $dna2prot -> site_directed_mutagenesis(70, 'TYR')
  -> save_textbox('PA_mutant_primer.txt');
