my $file = 'UBA52_synthesis.txt';

my $dna = 'ATGCAGATATTTGTGAAGACATTGACAGGGAAGACTATTACTCTGGAGGTTGAGCCCTCGGACACAATCGAGAATGTCAAGGCAAAAATCCAAGATAAAGAAGGTATCCCTCCTGACCAGCAACGACTTATATTTGCTGGCAAGCAACTAGAAGATGGGCGGACACTTTCGGATTACAATATACAAAAGGAAAGTACGCTCCATCTAGTGCTACGGCTTAGGGGTGGTATAATCGAGCCTTCACTTCGGCAGCTCGCCCAAAAATACAACTGCGATAAAATGATCTGCCGTAAGTGTTATGCTCGTTTACACCCCAGGGCAGTGAACTGCAGAAAAAAGAAATGCGGCCACACGAACAACCTCCGACCAAAGAAAAAAGTTAAG';

open(FH, '<', $file) || die "Cannot open $file. $!";
while (my $line = <FH>) {
  $line = ts($line);
  my ($name, $oligo) = split(/\|/, $line);
  $oligo =~ s/[^ATCG]//gi;
  $oligo = calculate_reverse_complement($oligo);
  print $name . '|' . $oligo . "\n";
};
close(FH);

$dna = calculate_reverse_complement($dna);

print "\n$dna\n";

sub calculate_reverse_complement {
  my $oligo = shift;
  $oligo =~ tr/[A,C,G,T]/[t,g,c,a]/;
  return uc(reverse($oligo))
}

#remove trailing white spaces, i.e. the carriage return and linefeed for inputing PC/DOS compatible file data
sub ts { my $i = shift; $i =~ s/\s+$//; return $i }
