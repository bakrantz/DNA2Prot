use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8"; #Make a path to the module
use Dna2prot;

my $dna2prot =  Dna2prot
             -> new()
             -> open_dna('PA_DNA.txt')
             -> save_dna('PA_DNA.txt')
             -> find_largest_orf;
   $dna2prot -> protein_name('PA');
   $dna2prot -> save_record('PA_DNA_record.txt');
