use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $dna2prot =  Dna2prot
             -> new()
             -> open_protein('PA_protein.txt')
             -> save_protein('PA_protein.txt')
             -> calculate_protein_parameters;
my $isoelectric_point = $dna2prot -> isoelectric_point;
print "PA's isoelectric point is $isoelectric_point\n";
   $dna2prot -> open_dna('PA_DNA.txt')
             -> save_record('PA_DNA_record.txt');
