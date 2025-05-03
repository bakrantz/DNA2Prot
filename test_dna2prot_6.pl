use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;

my $dna2prot = Dna2prot
  -> new()
  -> open_dna('PA_DNA.txt')
  -> find_restriction_sites
  -> save_textbox('PA_restriction_sites.txt');
