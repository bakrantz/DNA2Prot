use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;

my $dna2prot = Dna2prot
  -> new()
  -> open_dna('PA_DNA.txt')
  -> save_dna('PA_DNA.txt')
  -> find_largest_orf
  -> verbose_translate
  -> save_textbox('PA_verbose_translation.txt');
