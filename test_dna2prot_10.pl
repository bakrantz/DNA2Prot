use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Dna2prot"; #Make a path to the module
use Dna2prot;
my $dna2prot =  Dna2prot
             -> new()
             -> open_dna('PA_DNA.txt')
             -> find_largest_orf;
$dna2prot    -> amino_acid_numbers_offset(-20);
$dna2prot    -> protein_name('PA');
$dna2prot    -> site_directed_mutagenesis(50, 'PHE')
             -> save_textbox('PA_mutagenesis_primer_G50F.txt')
             -> transfer_mutant_to_current
             -> verbose_translate
             -> save_textbox('PA_verbose_translation_G50F.txt');
