use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $protein  =  'MGIKLPWRTEFGSTZ';
my $dna2prot =  Dna2prot
             -> new()
             -> reverse_translate($protein)
             -> verbose_translate
             -> save_textbox('PA_verbose_translation_reverse_translate.txt');

