#####################################################################################
#################                     Dna2prot.pm                    ################
#################              Krantz Lab       Dec 2022             ################
#####################################################################################
#################  Module to process DNA and protein sequence data   ################
#####################################################################################
package Dna2prot;
################################# THE CONSTRUCTOR ###################################
sub new {
  my $class                                   = shift;
  my $self                                    = {};                 #the anon hash
     $self->{CODONS}                          = [];
     $self->{AMINO_ACIDS}                     = [];
     $self->{AMINO_ACID_NUMBERS}              = [];
     $self->{MUTANT_AMINO_ACIDS}              = [];
     $self->{MUTANT_CODONS}                   = [];
     $self->{MUTANT_AMINO_ACID_NUMBERS}       = [];
     $self->{OLIGOS}                          = [];
     $self->{GIBSON_SEQUENCES}                = [];
     $self->{PROTEIN_LIBRARY}                 = [];
     $self->{DNA_FILE}                        = undef;
     $self->{RECORD_FILE}                     = undef;
     $self->{PROTEIN_FILE}                    = undef;
     $self->{PROTEIN_LIBRARY_FILE}            = undef;
     $self->{CURRENT_DNA_SEQUENCE}            = undef;
     $self->{CURRENT_PROTEIN_SEQUENCE}        = undef;
     $self->{AMINO_ACID_NUMBERS_OFFSET}       = undef;
     $self->{MUTANT_DNA_SEQUENCE}             = undef;
     $self->{MUTANT_PROTEIN_SEQUENCE}         = undef;
     $self->{MUTANT_DNA_NAME}                 = undef;
     $self->{MUTANT_PROTEIN_NAME}             = undef;
     $self->{DNA_NAME}                        = undef;
     $self->{PROTEIN_NAME}                    = undef;
     $self->{TEXTBOX}                         = undef;
     $self->{PRIMER_TM}                       = undef;
     $self->{PRIMER_PERCENT_GC}               = undef;
     $self->{PRIMER_LENGTH}                   = undef;
     $self->{PRIMER_PERCENT_MISMATCH}         = undef; 
     $self->{PRIMER_FORWARD_NAME}             = undef;
     $self->{PRIMER_REVERSE_NAME}             = undef;
     $self->{PRIMER_FORWARD_SEQUENCE}         = undef;
     $self->{PRIMER_REVERSE_SEQUENCE}         = undef;  
     $self->{PROTEIN_LENGTH}                  = undef;
     $self->{EXTINCTION_COEFFICIENT}          = undef;
     $self->{MOLECULAR_WEIGHT}                = undef;
     $self->{ISOELECTRIC_POINT}               = undef;
     $self->{VERBOSE}                         = undef;
  return bless($self, $class)                                  #bless and return thy self
}

############################### METHOD SUBS ########################################
sub open_dna {
  my ($self, $file) = @_;
  $file = dna_file($self) unless $file;  
  open(FH, '<', dna_file($self, $file)) || die "Cannot open $file in open_dna. $!";
  print "Opening DNA data file $file.\n" if verbose($self);
  my $gt = chr(62); #greater than character denotes beginning of the name of the dna sequence
  my ($raw_dna_seq, $name) = (undef, []);
  while (my $line = <FH>) {
    $line = ts($line); #remove trailing white space(s) from newlines
    if ($line =~ /^$gt\w*/) { @$name = $line =~ m/^$gt(\w*)/ }
    else { $raw_dna_seq .= $line };
  };
  close(FH);
  #Next remove non-dna characters, whitespaces, etc. and store dna sequence in current_dna_sequence
  fasta_dna($self, uc($raw_dna_seq));
  dna_name($self, $name->[0]) if $name->[0]; #won't overwrite name in module if not given in file
  return $self
}

sub open_protein {
  my ($self, $file) = @_;
  $file = protein_file($self) unless $file;  
  open(FH, '<', protein_file($self, $file)) || die "Cannot open $file in open_protein. $!";
  print "Opening protein data file $file.\n" if verbose($self);
  my $gt = chr(62); #greater than character denotes beginning of the name of the protein sequence
  my ($raw_protein_seq, $name) = (undef, []);  
  while (my $line = <FH>) {
    $line = ts($line); #remove trailing white space(s) from newlines
    if ($line =~ /^$gt\w*/) { @$name = $line =~ m/^$gt(\w*)/ }
    else { $raw_protein_seq .= $line };
  };
  close(FH);
  #Next remove non-protein characters, spaces, etc. and store protein sequence in current_protein_sequence
  fasta_protein($self, uc($raw_protein_seq)); 
  protein_name($self, $name->[0]) if $name->[0]; #won't overwrite protein name in module if not given in file 
  return $self
}

sub save_protein {
  my ($self, $file, $protein, $name) = @_;
  $file = protein_file($self) unless $file;
  $protein = current_protein_sequence($self) unless $protein;
  $name = protein_name($self) unless $name;
  my $cr = chr(13) . chr(10);
  my $first_line = '>' . protein_name($self, $name) . $cr;
  open(FH, '>', protein_file($self, $file)) || die "Cannot save $file in save_protein. $!";
  print "Saving protein data file $file.\n" if verbose($self);
  print FH $first_line;
  print FH current_protein_sequence($self, $protein), $cr;
  close(FH);
  return $self
}

sub save_dna {
  my ($self, $file, $dna, $name) = @_;
  $file = dna_file($self) unless $file;
  $dna = current_dna_sequence($self) unless $dna;
  $name = dna_name($self) unless $name;
  my $cr = chr(13) . chr(10);
  my $first_line = '>' . dna_name($self, $name) . $cr;
  open(FH, '>', dna_file($self, $file)) || die "Cannot save $file in save_dna. $!";
  print "Saving dna data file $file.\n" if verbose($self);
  print FH $first_line;
  print FH current_dna_sequence($self, $dna), $cr;
  close(FH);
  return $self
}

sub save_textbox {
  my ($self, $outfile) = @_;
  my $textbox = textbox($self);
  open(FH, '>', $outfile) || die "Cannot save $outfile in save_textbox. $!";
  print "Saving textbox data file $outfile.\n" if verbose($self);
  print FH $textbox;
  close(FH);
  return $self
}

sub save_record {
  my ($self, $outfile) = @_;
  $outfile = record_file($self) unless $outfile;
  open(FH, '>', record_file($self, $outfile)) || die "Cannot save $outfile in save_record. $!";
  print "Saving DNA/protein record file $outfile.\n" if verbose($self);  
  my ($pound, $cr) = (chr(35), chr(13).chr(10)); #pound to signify comments in output file record 
  calculate_protein_parameters($self);
  print FH "$pound Protein name: ", protein_name($self), $cr;
  print FH "$pound Protein length: ", protein_length($self), " residues$cr";
  print FH "$pound Molecular weight: ", molecular_weight($self), " Da$cr";
  print FH "$pound Isoelectric point (pI): ", isoelectric_point($self), $cr;
  print FH "$pound Extinction coefficient: ", extinction_coefficient($self), " M-1 cm-1$cr";
  print FH $cr;
  print FH ">", protein_name($self), $cr;
  print FH current_protein_sequence($self), $cr;
  print FH $cr;
  print FH ">", dna_name($self), $cr;
  print FH current_dna_sequence($self), $cr;
  close(FH);
  return $self
}

sub open_record {
  my ($self, $file) = @_;
  $file = record_file($self) unless $file;
  open(FH, '<', record_file($self, $file)) || die "Cannot open $file in open_record. $!";
  print "Opening DNA/protein record file $file.\n" if verbose($self); 
  my $pound = chr(35); #pound denotes header and comments
  my $gt = chr(62); #greater than character denotes beginning of the name of the sequences
  my @header = ();
  my @names = ();
  my @sequences =();
  while (my $line = <FH>) {
    $line = ts($line);
    if ($line =~ /^$pound/) {
      push @header, $line
    }
    elsif ($line =~ /^$gt/) {
      my @n = $line =~ m/^$gt(\w*)/;
      push @names, $n[0];
      my $seq = <FH>; #next line must be dna or protein sequence string
      push @sequences, ts($seq);
    };
  };
  close(FH);
  protein_name($self, $names[0]); #first name is protein name
  current_protein_sequence($self, $sequences[0]); #first sequence is protein sequence
  dna_name($self, $names[1]); #second name is dna name
  current_dna_sequence($self, $sequences[1]); #second sequence in dna sequence
  return $self
}

sub open_protein_library {
  my ($self, $file) = @_;
  $file = protein_library_file($self) unless $file;
  open(FH, '<', protein_library_file($self, $file)) || die "Cannot open $file in open_protein_library. $!";
  print "Opening protein library file $file.\n" if verbose($self);
  my $pound = chr(35); #pound denotes header and comments
  my $gt = chr(62); #greater than character denotes beginning of the name of the sequences  
  my ($i, $header, $protein_library) = (0, [], []);
  while (my $line = <FH>) {
    $line = ts($line);
    if ($line =~ /^$pound/) {
      push @$header, $line
    }
    elsif ($line =~ /^$gt/) {
      my @name = $line =~ m/^$gt(\w*)/; #extract name after greater than sign
      my $sequence = <FH>; #next line must be protein sequence string
      $protein_library -> [$i]   -> {NAME}     = $name[0];
      $protein_library -> [$i++] -> {SEQUENCE} = ts($sequence);
    };
  };
  close(FH);
  protein_library($self, @$protein_library);
  return $self
}

sub save_protein_library {
  my ($self, $file, $protein_library) = @_;
  $file = protein_library_file($self) unless $file;
  $protein_library = [ protein_library($self) ] unless $protein_library;
  open(FH, '>', protein_library_file($self, $file)) || die "Cannot save $file in save_protein_library. $!";
  print "Saving protein library file $file in save_protein_library.\n" if verbose($self);
  my ($gt, $cr) = (chr(62), chr(13) . chr(10)); #greater than character denotes beginning of the name of the sequences
  foreach my $protein (@$protein_library) { print FH $gt . $protein -> {NAME} . $cr . $protein -> {SEQUENCE} . $cr };
  close(FH);
  return $self
}

#Appends current_protein_sequence with protein_name to the protein_library and saves the library file
sub append_protein_library {
  my ($self, $protein_sequence, $protein_name, $file) = @_;
  $protein_sequence = current_protein_sequence($self) unless $protein_sequence;
  $protein_name = protein_name($self) unless $protein_name;
  $file = protein_library_file($self) unless $file;
  open_protein_library($self, $file);
  my $protein_library = [ protein_library($self) ];
  my $index = @$protein_library;
  $protein_library -> [$index] -> {NAME}     = protein_name($self, $protein_name);
  $protein_library -> [$index] -> {SEQUENCE} = current_protein_sequence($self, $protein_sequence);  
  protein_library($self, @$protein_library);
  save_protein_library($self, $file);
  return $self
}

#Searches for best match of provided protein or current_protein_sequence with those sequences in the library
sub search_protein_library {
  my ($self, $protein, $protein_library_file) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  $protein_library_file = protein_library_file($self) unless $protein_library_file;
  open_protein_library($self, $protein_library_file) if $protein_library_file;
  my $protein_library = [ protein_library($self) ];
  #Find longest common pepitde comparing protein to library
  #Also find length difference between target sequence and each sequence in library
  for (my $i = 0; $i < scalar(@$protein_library); $i++) {
    my $lcs = longest_common_sequence($protein_library -> [$i] -> {SEQUENCE}, $protein);
    $protein_library -> [$i] -> {LONGEST_COMMON_SEQUENCE} = $lcs;
    my $length_difference = length($protein) - length($protein_library -> [$i] -> {SEQUENCE});
    $protein_library -> [$i] -> {LENGTH_DIFFERENCE} = $length_difference;
  };
  #Score by the length of the lcs
  my $high_score = length($protein_library -> [0] -> {LONGEST_COMMON_SEQUENCE});
  $protein_library -> [0] -> {SCORE} = $high_score;
  for (my $i = 1; $i < scalar(@$protein_library); $i++) {  
    my $current_score = length($protein_library -> [$i] -> {LONGEST_COMMON_SEQUENCE});
    $protein_library -> [$i] -> {SCORE} = $current_score;
    if ($current_score >= $high_score) { $high_score = $current_score };
  };
  #Find the candidates' indexes that have the highest score. There could be more than one with the same high score.
  my $candidates = [];
  for (my $i = 0; $i < scalar(@$protein_library); $i++) { if ($protein_library -> [$i] -> {SCORE} == $high_score) { push @$candidates, $i } };
  my $best = undef;
  if (scalar(@$candidates) > 1) {
    #when there is more than one candidate pick the one closest in length to target sequence
    my $best_abs_length_difference = abs($protein_library -> [0] -> {LENGTH_DIFFERENCE});
    $best = 0;
    for (my $i = 1; $i < @$candidates; $i++) {
      my $current_abs_length_difference = abs($protein_library -> [$i] -> {LENGTH_DIFFERENCE});
      if ($current_abs_length_difference < $best_abs_length_difference) {
	$best = $candidate -> [$i];
	$best_abs_length_difference = $current_abs_length_difference;
      };
    };
  }
  elsif (scalar(@$candidates) == 1) { $best = $candidates -> [0] };
  #Align target protein and best candidate from library and determine percentage match
  my $best_library_protein = $protein_library -> [$best] -> {SEQUENCE};
  my $lcs = $protein_library -> [$best] -> {LONGEST_COMMON_SEQUENCE};
  my $protein_index = index($protein, $lcs);
  my $best_index = index($best_library_protein, $lcs);
  #determine the starting index in protein and best protein given the alignment of the lcs
  my ($start_protein_index, $start_best_index, $end_protein_index, $end_best_index) = (undef, undef, undef, undef);
  my ($length_protein, $length_best) = (length($protein), length($best_library_protein));
  if ($protein_index > $best_index) {
    $start_protein_index = $protein_index - $best_index;
    $start_best_index = 0;
    $end_protein_index = $start_protein_index + $length_best - 1;
    if ($end_protein_index < $length_protein - 1) {
      $end_best_index = $length_best - 1;
    }
    elsif ($end_protein_index > $length_protein - 1) {
      $end_best_index = $length_best - 1 - ($end_protein_index - $length_protein + 1);
      $end_protein_index = $length_protein - 1;
    }
    elsif ($end_protein_index == $length_protein - 1) {
      $end_protein_index = $length_protein - 1;
      $end_best_index = $length_best - 1;
    }
  }
  elsif ($protein_index < $best_index) {
    $start_protein_index = 0;
    $start_best_index = $best_index - $protein_index;
    $end_best_index = $start_best_index + $length_protein - 1;
    if ($end_best_index < $length_best - 1) {
      $end_protein_index = $length_protein - 1;
    }
    elsif ($end_best_index > $length_best - 1) {
      $end_protein_index = $length_protein - 1 - ($end_best_index - $length_best + 1);
      $end_best_index = $length_best - 1;
    }
    elsif ($end_best_index == $length_best - 1) {
      $end_protein_index = $length_protein - 1;
      $end_best_index = $length_best - 1;
    }
  }
  elsif ($protein_index == $best_index) {
    $start_protein_index = 0;
    $start_best_index = 0;
    if ($length_protein < $length_best) {
      $end_protein_index = $length_protein - 1;
      $end_best_index = $length_protein - 1;      
    }
    elsif ($length_protein > $length_best) {
      $end_protein_index = $length_best - 1;
      $end_best_index = $length_best - 1; 
    }
    elsif ($length_protein == $length_best) {
      $end_protein_index = $length_protein - 1;
      $end_best_index = $length_best - 1;
    }
  };
  my $point_mutations = [];
  my $j = $start_protein_index;
  for (my $i = $start_best_index; $i <= $end_best_index; $i++) {
    my $best_residue_number = $i + 1; 
    my $best_amino_acid    = substr($best_library_protein, $i, 1);
    my $protein_amino_acid = substr($protein, $j, 1);
    if ($best_amino_acid ne $protein_amino_acid) { push @$point_mutations, $best_amino_acid . $best_residue_number . $protein_amino_acid };
    $j++;
  };
  #percent_mutations is calculated for the aligned region
  my $percent_mutations = 100 * scalar(@$point_mutations) / ($end_best_index - $start_best_index + 1);
  my ($name, $mutations, $mutant_name) = (undef, undef, undef);
  if ($percent_mutations > 80) {
    print "No significant alignment observed with the protein library.\n" if verbose($self);
  }
  else {
    $name = $protein_library -> [$best] -> {NAME};
    if (scalar(@$point_mutations)) { $mutations = join(chr(95), @$point_mutations) } else { $mutations = 'wildtype' };
    $mutant_name = $name . chr(95) . $mutations;
    print "Best alignment to $name. Protein's name is $mutant_name.\n" if verbose($self);
    protein_name($self, $mutant_name);
  };
  return $self
}


#Aligns and joins two DNA sequencing reads from the same construct.
#Two DNA sequence filenames are given as arguments. See also align_and_join below.
sub align_and_join_files {
  my ($self, $file1, $file2) = @_;
  die "Two files need to be given in align_and_join_files.\n" unless ($file1 && $file2);
  open_dna($self, $file1);
  my $dna1 = current_dna_sequence($self);
  open_dna($self, $file2);
  my $dna2 = current_dna_sequence($self);  
  align_and_join($self, $dna1, $dna2);
  return $self
}

#Takes two dna sequences as arguments. Finds largest common dna string and stitches dna sequences together
#searches for best alignment when second dna string is either forward or reverse complemented orientation
sub align_and_join {
  my ($self, $dna1, $dna2) = @_;
  reverse_complement($self, $dna2);
  my $dna2_rc = current_dna_sequence($self);  
  my $lcs     = longest_common_sequence($dna1, $dna2);
  my $lcs_rc  = longest_common_sequence($dna1, $dna2_rc);
  my ($dna1_fragment, $dna2_fragment) = (undef, undef);
  if (length($lcs) >= length($lcs_rc)) {
    print "Sequences show greater alignment when both are forward DNA sequences. Joining the aligned sequences.\n" if verbose($self);
    $dna1_fragment = substr($dna1, 0, index($dna1, $lcs) + length($lcs));
    $dna2_fragment = substr($dna2, index($dna2, $lcs) + length($lcs));	   
  }
  else {
    print "Sequences show greater alignment when second DNA sequence in reverse complemented. Joining the aligned sequences.\n" if verbose($self);
    $dna1_fragment = substr($dna1, 0, index($dna1, $lcs_rc) + length($lcs_rc));
    $dna2_fragment = substr($dna2_rc, index($dna2_rc, $lcs_rc) + length($lcs_rc));    
  };
  current_dna_sequence($self, $dna1_fragment . $dna2_fragment);
  return $self
}

sub longest_common_sequence {
  my ($seq1, $seq2) = @_;
  my $lcs = undef;
  my @best;
  "$seq1\n$seq2" =~ /(.+).*\n.*\1(?{ $best[length $1]{$1}++})(*FAIL)/;
  $lcs = $_ for sort keys %{ $best[-1] };
  return $lcs
}

#finds and replaces self matching DNA sequences with silent mutations to preserve amino acid coding
# $dna is the dna sequence. $rfdna is the reading frame of the amino acid coding sequence.
# $max is the maximum allowable self match
sub find_and_replace_self_matches {
  my $self = shift;
  my ($dna, $rfdna, $max) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  $max = 9 unless $max; #Can't make max much smaller 
  my $old_dna = $dna; #remember unmodified DNA sequence
  my $match = longest_self_match($dna, $max); #find longest self-matching sequence in $dna
  my $counter = 0;
  while (length($match) >= $max && $counter < 10000) {
    my $pos = index($dna, $match); #position in $dna where self-match sequence $match aligns with dna sequence
    if ($pos != -1) {
      my $rfmatch = (3 - (($pos - $rfdna) % 3)) % 3; #$rfmatch is the reading frame of the self-matched sequence, called $match
      my $i = $rfmatch + $pos; #first position of a codon that can be silently mutated
      my @pos_to_mut = ();
      do { push @pos_to_mut, $i; $i += 3 } while $i < ($pos + length($match)); #all codon positions that can be silently mutated
      print "Number of positions to mutate in find_and_replace_self_matches is: @pos_to_mut\n" if verbose($self);
      my $pos_picked = $pos_to_mut[int rand @pos_to_mut]; #pick random position to silent mutate
      my $codon = substr($dna, $pos_picked, 3); #codon located at picked position
      print "Original codon is: $codon at $pos_picked\n" if verbose($self);
      my %codons = (); #a hash of possible codons to replace selected site to silently mutate
      my $aa = +{codon_table()} -> {$codon};
      #the amino acids W and M have only one possible codon and no silent mutation can be selected
      if ($aa ne 'W' && $aa ne 'M') { 
        foreach my $key (keys %{ +{codon_table()} }) { $codons{$key} = $aa if (+{codon_table()} -> {$key} eq $aa)};
        delete $codons{$codon};
      };
      if (scalar(keys %codons)) {
        my $codon_sel = [ keys %codons ] -> [ int rand keys %codons ]; #pick random silent mutation codon
        print "Selected codon is: $codon_sel at $pos_picked\n" if verbose($self);
        $dna = substr($dna, 0, $pos_picked) . $codon_sel . substr($dna, $pos_picked + 3); #replace silent mutation codon in dna
      };
    };
    $match = longest_self_match($dna, $max);
    $counter++;
  };
  my ($comp_out, $textbox, $cr) = (undef, undef, chr(13).chr(10));
  for (0..(length($dna)-1)) { if (substr($old_dna, $_, 1) eq substr($dna, $_, 1)) { $comp_out .= '.' } else { $comp_out .= '*' } };
  #Output to textbox new DNA sequence aligned to previous input DNA sequence 
  $textbox = '>New_DNA_Sequence'.$cr.$dna.$cr.'>Alignment'.$cr.$comp_out.$cr.'>Original_DNA_Sequence'.$cr.$old_dna.$cr;
  textbox($self, $textbox);
  print $textbox if verbose($self);
  current_dna_sequence($self, $dna); #Store silently mutated new DNA sequence in current_dna_sequence
  return $self
}

#Returns longest internal self-match within a DNA string
sub longest_self_match {
  my ($dna, $max) = @_;
  my %hash_matches = ();
  $dna =~ s/(?=(.{$max}))/$hash_matches{$1}++;$1/eg;
  foreach my $key (keys %hash_matches) { delete $hash_matches{$key} if ($hash_matches{$key} == 1) };
  my $longest_self_match = [ sort { length($b) <=> length($a) } keys %hash_matches ] -> [0];
  print "Longest self match is: $longest_self_match\n" if verbose($self);
  return $longest_self_match
}

#DNA/gene synthesis for 500 bp or less sequences
#Generates oligos for either Ligase Chain Reaction (LCR) or Polymerase Chain Reaction (PCR) method.
sub dna_synthesis {
  my ($self, $dna, $method, $target_tm) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  current_dna_sequence($self, $dna);
  $method = 'PCR' unless $method;
  $target_tm = 55 unless $target_tm;
  my ($oligos, $hybridization) = ([], []);
  if ($method eq 'PCR') {
    #Build PCR hybridization units
    my $index_hybrid_units = 0;
    my $five_prime_pos   = 19;
    my $three_prime_pos  = $five_prime_pos + 15;
    do {
      my $oligo_length     = $three_prime_pos - $five_prime_pos + 1;
      my $oligo_sequence   = substr($dna, $five_prime_pos, $oligo_length);
      my $simple_tm        = calculate_simple_tm($oligo_sequence);
      #optimize tm of hybridization unit
      while ($simple_tm < $target_tm && $oligo_length <= 20 && $three_prime_pos < length($dna) - 20) {
        $three_prime_pos++;  
        $oligo_length = $three_prime_pos - $five_prime_pos + 1;
        $oligo_sequence = substr($dna, $five_prime_pos, $oligo_length);
        $simple_tm = calculate_simple_tm($oligo_sequence);
      };
      $hybridization -> [$index_hybrid_units] -> {'FIVE_PRIME_POSITION'}  = $five_prime_pos;
      $hybridization -> [$index_hybrid_units] -> {'THREE_PRIME_POSITION'} = $three_prime_pos;    
      $hybridization -> [$index_hybrid_units] -> {'OLIGO_LENGTH'}         = $oligo_length;
      $hybridization -> [$index_hybrid_units] -> {'OLIGO_SEQUENCE'}       = $oligo_sequence;
      $hybridization -> [$index_hybrid_units] -> {'SIMPLE_TM'}            = $simple_tm;
      $five_prime_pos = $three_prime_pos + 1;
      $three_prime_pos  = $five_prime_pos + 15;
      $index_hybrid_units++;
    } while ($three_prime_pos < length($dna) - 20);
    #Build PCR primer set using hybridization units
    my $index_primers = 0;
    $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = 0;
    $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'};    
    $oligos -> [$index_primers] -> {'OLIGO_LENGTH'}         = $hybridization -> [0] -> {'THREE_PRIME_POSITION'} + 1;
    $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'}       = substr($dna, 0, $oligos -> [$index_primers] -> {'OLIGO_LENGTH'});
    $oligos -> [$index_primers] -> {'DIRECTION'}            = 'FORWARD';    
    $index_primers++;
    my $max_hybrid = @$hybridization;
    my $max_hybrid_minus_one = $max_hybrid - 1;
    for ( my $i = 0; $i < $max_hybrid_minus_one; $i++ ) {
      $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
      $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
      #set directions of the internal oligos
      if ($index_primers <= int($max_hybrid / 2)) { $oligos -> [$index_primers] -> {'DIRECTION'} = 'FORWARD' }
      else { $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE' };      
      $index_primers++;
    };
    $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'} = $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};
    $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = length($dna) - 1;    
    $oligos -> [$index_primers] -> {'OLIGO_LENGTH'} = length($dna) - $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};
    $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'OLIGO_LENGTH'});
    $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE';
    #Reverse complement all REVERSE direction oligos
    for ( my $i = 0; $i < @$oligos; $i++ ) {
      if ($oligos -> [$i] -> {'DIRECTION'} eq 'REVERSE') { $oligos -> [$i] -> {'OLIGO_SEQUENCE'} = calculate_reverse_complement($oligos -> [$i] -> {'OLIGO_SEQUENCE'}); };
    };
    #Name the FORWARD oligo sequences
    for ( my $i = 0; $i < @$oligos; $i++ ) { if ($oligos->[$i]->{'DIRECTION'} eq 'FORWARD') { my $name_index = $i + 1; $oligos->[$i]->{'NAME'} = 'F' . $name_index } };
    #Name the REVERSE oligo sequences
    my $rev_name_index = 1;
    for ( my $i = @$oligos - 1; $i >= 0; $i-- ) { if ($oligos->[$i]->{'DIRECTION'} eq 'REVERSE') { $oligos->[$i]->{'NAME'} = 'R' . $rev_name_index++ } };
  }
  elsif ($method eq 'LCR') {
    my $index_hybrid_units = 0;
    my $five_prime_pos   = 0;
    my $three_prime_pos  = $five_prime_pos + 15;
    do {
      my $oligo_length     = $three_prime_pos - $five_prime_pos + 1;
      my $oligo_sequence   = substr($dna, $five_prime_pos, $oligo_length);
      my $simple_tm        = calculate_simple_tm($oligo_sequence);
      #optimize tm of hybridization unit
      while ($simple_tm < $target_tm && $oligo_length <= 20 && $three_prime_pos < length($dna) - 15) {
        $three_prime_pos++;  
        $oligo_length = $three_prime_pos - $five_prime_pos + 1;
        $oligo_sequence = substr($dna, $five_prime_pos, $oligo_length);
        $simple_tm = calculate_simple_tm($oligo_sequence);
      };
      $hybridization -> [$index_hybrid_units] -> {'FIVE_PRIME_POSITION'}  = $five_prime_pos;
      $hybridization -> [$index_hybrid_units] -> {'THREE_PRIME_POSITION'} = $three_prime_pos;    
      $hybridization -> [$index_hybrid_units] -> {'OLIGO_LENGTH'}         = $oligo_length;
      $hybridization -> [$index_hybrid_units] -> {'OLIGO_SEQUENCE'}       = $oligo_sequence;
      $hybridization -> [$index_hybrid_units] -> {'SIMPLE_TM'}            = $simple_tm;
      $five_prime_pos = $three_prime_pos + 1;
      $three_prime_pos  = $five_prime_pos + 15;
      $index_hybrid_units++;      
    } while ($three_prime_pos < length($dna) - 15);
    #Make last hybridization unit   
    $hybridization -> [$index_hybrid_units] -> {'FIVE_PRIME_POSITION'}  = $five_prime_pos;
    $hybridization -> [$index_hybrid_units] -> {'THREE_PRIME_POSITION'} = length($dna) - 1;
    $hybridization -> [$index_hybrid_units] -> {'OLIGO_LENGTH'}         = length($dna) - $five_prime_pos;
    $hybridization -> [$index_hybrid_units] -> {'OLIGO_SEQUENCE'}       = substr($dna, $five_prime_pos, length($dna) - $five_prime_pos);
    $hybridization -> [$index_hybrid_units] -> {'SIMPLE_TM'}            = calculate_simple_tm($hybridization -> [$index_hybrid_units] -> {'OLIGO_SEQUENCE'});    
    #Generate the FORWARD direction oligo primers for LCR method
    my $max_hybrid = @$hybridization;
    my $max_hybrid_minus_one = $max_hybrid - 1;
    #Forward primer design depends of whether maximum number of hybridization units are ODD or EVEN
    if ($max_hybrid % 2) {
      #Generate forward primers for LCR synthesis with ODD number of hybridization units
      my $index_primers = 0;
      for ($i = 0; $i < $max_hybrid_minus_one; $i += 2) {
        $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
        $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
        $oligos -> [$index_primers] -> {'DIRECTION'} = 'FORWARD'; 	
        $index_primers++;
      };
      $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [-1] -> {'THREE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [-1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'} + 1;
      $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
      $oligos -> [$index_primers] -> {'DIRECTION'} = 'FORWARD';
      #Name the FORWARD oligo sequences
      for ( my $i = 0; $i < @$oligos; $i++ ) { if ($oligos ->[$i]->{'DIRECTION'} eq 'FORWARD') { my $name_index = $i + 1; $oligos->[$i]->{'NAME'} = 'F' . $name_index } };
      #Generate reverse primers for LCR synthesis with ODD number of hybridization units
      $index_primers++;
      $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [0] -> {'FIVE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'} - $hybridization -> [0] -> {'FIVE_PRIME_POSITION'} + 1;
      $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $hybridization -> [0] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
      $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE';
      $index_primers++;
      for ( my $i = 1; $i < $max_hybrid; $i += 2) {
        $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
        $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
        $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE'; 	
        $index_primers++;
      };
      #Name the REVERSE oligo sequences
      my $rev_name_index = 1;
      for ( my $i = @$oligos - 1; $i >= 0; $i-- ) { if ($oligos ->[$i]->{'DIRECTION'} eq 'REVERSE') { $oligos->[$i]->{'NAME'} = 'R' . $rev_name_index++ } };
    }
    else {
      #Generate forward primers for LCR synthesis with EVEN number of hybridization units
      my $index_primers = 0;
      for ( my $i = 0; $i < $max_hybrid; $i += 2) {
        $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
        $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
        $oligos -> [$index_primers] -> {'DIRECTION'} = 'FORWARD'; 	
        $index_primers++;
      };
      #Name the FORWARD oligo sequences
      for ( my $i = 0; $i < @$oligos; $i++ ) { if ($oligos ->[$i]->{'DIRECTION'} eq 'FORWARD') { my $name_index = $i + 1; $oligos->[$i]->{'NAME'} = 'F' . $name_index } };
      #Generate reverse primers for LCR synthesis with EVEN number of hybridization units
      $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [0] -> {'FIVE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'} - $hybridization -> [0] -> {'FIVE_PRIME_POSITION'} + 1;
      $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $hybridization -> [0] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
      $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE';
      $index_primers++;      
      for ( my $i = 1; $i < $max_hybrid_minus_one; $i += 2 ) {
        $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
        $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
        $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
        $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE'; 	
        $index_primers++;
      };
      $oligos -> [$index_primers] -> {'FIVE_PRIME_POSITION'}  = $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'THREE_PRIME_POSITION'} = $hybridization -> [-1] -> {'THREE_PRIME_POSITION'};
      $oligos -> [$index_primers] -> {'LENGTH'} = $hybridization -> [-1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'} + 1;
      $oligos -> [$index_primers] -> {'OLIGO_SEQUENCE'} = substr($dna, $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'}, $oligos -> [$index_primers] -> {'LENGTH'});
      $oligos -> [$index_primers] -> {'DIRECTION'} = 'REVERSE';
      #Name the REVERSE oligo sequences
      my $rev_name_index = 1;
      for ( my $i = @$oligos - 1; $i >= 0; $i-- ) { if ($oligos ->[$i]->{'DIRECTION'} eq 'REVERSE') { $oligos->[$i]->{'NAME'} = 'R' . $rev_name_index++ } };
    };
    #Reverse complement all REVERSE direction oligos in LCR synthesis
    for ( my $i = 0; $i < @$oligos; $i++ ) {
      if ($oligos -> [$i] -> {'DIRECTION'} eq 'REVERSE') { $oligos->[$i]->{'OLIGO_SEQUENCE'} = calculate_reverse_complement($oligos->[$i]->{'OLIGO_SEQUENCE'}); };
    };
  }
  else {
    die "Method $method not recognized. Use either PCR or LCR methods. $!";  
  };
  #Print hybridization unit output if in verbose mode
  print "DNA Synthesis for $method method: the Hybridization Units\n" if verbose($self);
  foreach my $hybrid_unit (@$hybridization) { 
    print '5prime pos '  . $hybrid_unit -> {'FIVE_PRIME_POSITION'} if verbose($self);
    print ' 3prime pos ' . $hybrid_unit -> {'THREE_PRIME_POSITION'} if verbose($self);
    print ' Length '     . $hybrid_unit -> {'OLIGO_LENGTH'} if verbose($self);
    print ' Tm '         . $hybrid_unit -> {'SIMPLE_TM'} if verbose($self);
    print ' Sequence '   . $hybrid_unit -> {'OLIGO_SEQUENCE'} . "\n" if verbose($self);
  };
  #Build textbox of primer sequence output
  print "Primer sequences for DNA Synthesis using $method method\n" if verbose($self);
  my ($textbox, $cr, $prime, $pipe) = (undef, chr(13).chr(10), chr(39), chr(124));
  for ( my $i = 0; $i < @$oligos; $i++ ) { $textbox .= $oligos->[$i]->{'NAME'} . $pipe . '5' . $prime . '-' . $oligos->[$i]->{'OLIGO_SEQUENCE'} . '-3' . $prime . $cr };
  print $textbox if verbose($self);
  textbox($self, $textbox);
  oligos($self, @$oligos);
  return $self
}

sub calculate_simple_tm {
  my $oligo = shift;
  my $g_count = ($oligo =~ s/G/G/g);
  my $c_count = ($oligo =~ s/C/C/g);
  my $a_count = ($oligo =~ s/A/A/g);
  my $t_count = ($oligo =~ s/T/T/g);
  my $simple_tm = ($g_count + $c_count) * 4 + ($a_count + $t_count) * 2;
  return $simple_tm
}

sub calculate_reverse_complement {
  my $oligo = shift;
  $oligo =~ tr/[A,C,G,T]/[t,g,c,a]/;
  return uc(reverse($oligo))
}

#Gibson assembly of dsDNA fragments to make larger genes >500 bp
# Outputs an array of hashes containing the sequences of ~500 bp with overlapping ends to direct the assembly.
# Productive assembly has been achieved for DNA fragments with as little as a 12 bp overlap,
# however, it depends on the GC content of the overlap. We recommend using at least 15 bp overlaps,
# or more, for dsDNA assembly with a Tm ≥ 48°C (AT pair = 2°C and GC pair = 4°C). Increasing the length
# of overlap between fragments also reduces the amount of DNA needed for assembly.
sub gibson_assembly {
  my ($self, $dna, $target_tm) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  current_dna_sequence($self, $dna);
  $target_tm = 55 unless $target_tm; #using simple_tm method
  my ($hybridization, $gibson_sequences) = ([], []);
  my $length_dna = length($dna);
  my $avg_gibson_fragment_length = 500;
  my $number_hybrid = int($length_dna / $avg_gibson_fragment_length);
  if ($length_dna / $avg_gibson_fragment_length - $number_hybrid < 0.1) {
    do { $avg_gibson_fragment_length++ } while ( $length_dna / $avg_gibson_fragment_length - int($length_dna / $avg_gibson_fragment_length) > 0.01);
    $number_hybrid--;
  };
  die "Target DNA sequence is too small. Use dna_synthesis to design primers to synthesize smaller target DNA sequence ~500 bp in length" if ($number_hybrid == 0);
  #Build hybridization units approximately avg_gibson_fragment_length apart
  my $five_prime_pos = $avg_gibson_fragment_length - 15;
  my $three_prime_pos = $avg_gibson_fragment_length;
  for my $index_hybrid (0..($number_hybrid - 1)) {
    my $length_hybrid = $three_prime_pos - $five_prime_pos + 1;
    my $sequence_hybrid = substr($dna, $five_prime_pos, $length_hybrid);
    my $hybrid_tm = calculate_simple_tm($sequence_hybrid);
    while ($hybrid_tm < $target_tm && $length_hybrid < 40) {
      $three_prime_pos++;
      $length_hybrid = $three_prime_pos - $five_prime_pos + 1;
      $sequence_hybrid = substr($dna, $five_prime_pos, $length_hybrid);
      $hybrid_tm = calculate_simple_tm($sequence_hybrid);
    };
    $hybridization -> [$index_hybrid] -> {'FIVE_PRIME_POSITION'}  = $five_prime_pos;
    $hybridization -> [$index_hybrid] -> {'THREE_PRIME_POSITION'} = $three_prime_pos;
    $hybridization -> [$index_hybrid] -> {'LENGTH'}               = $length_hybrid;
    $hybridization -> [$index_hybrid] -> {'SEQUENCE'}             = $sequence_hybrid;
    $hybridization -> [$index_hybrid] -> {'TM'}                   = $hybrid_tm;
    $five_prime_pos = $three_prime_pos + $avg_gibson_fragment_length - 15;
    $three_prime_pos = $five_prime_pos + 15;
  };
  #Build Gibson fragment sequences
  my $gibson_index = 0;
  $gibson_sequences -> [$gibson_index] -> {'FIVE_PRIME_POSITION'}  = 0;
  $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} = $hybridization -> [0] -> {'THREE_PRIME_POSITION'};
  $gibson_sequences -> [$gibson_index] -> {'LENGTH'}               = $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} + 1;
  $gibson_sequences -> [$gibson_index] -> {'SEQUENCE'}             = substr($dna, 0, $gibson_sequences -> [$gibson_index] -> {'LENGTH'});
  my ($dna_from, $dna_to) = (1, $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} + 1);
  $gibson_sequences -> [$gibson_index] -> {'NAME'}                 = dna_name($self) . '_' . $dna_from . '-' . $dna_to;
  $gibson_index++;
  if ($number_hybrid > 1) {
    for my $i (0..($number_hybrid - 2)) {
      $gibson_sequences -> [$gibson_index] -> {'FIVE_PRIME_POSITION'} = $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'};
      $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'};
      $gibson_sequences -> [$gibson_index] -> {'LENGTH'} = $hybridization -> [$i + 1] -> {'THREE_PRIME_POSITION'} - $hybridization -> [$i] -> {'FIVE_PRIME_POSITION'} + 1;
      $gibson_sequences -> [$gibson_index] -> {'SEQUENCE'} = substr($dna, $hybridization->[$i]->{'FIVE_PRIME_POSITION'}, $gibson_sequences->[$gibson_index]->{'LENGTH'});
      ($dna_from, $dna_to) = ($gibson_sequences -> [$gibson_index] -> {'FIVE_PRIME_POSITION'} + 1, $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} + 1);
      $gibson_sequences -> [$gibson_index] -> {'NAME'} = dna_name($self) . '_' . $dna_from . '-' . $dna_to;   
      $gibson_index++;
    };
  };
  $gibson_sequences -> [$gibson_index] -> {'FIVE_PRIME_POSITION'} = $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};  
  $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} = $length_dna - 1;
  $gibson_sequences -> [$gibson_index] -> {'LENGTH'} = $length_dna - $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'};
  $gibson_sequences -> [$gibson_index] -> {'SEQUENCE'} = substr($dna, $hybridization -> [-1] -> {'FIVE_PRIME_POSITION'}, $gibson_sequences -> [$gibson_index] -> {'LENGTH'});
  ($dna_from, $dna_to) = ($gibson_sequences -> [$gibson_index] -> {'FIVE_PRIME_POSITION'} + 1, $gibson_sequences -> [$gibson_index] -> {'THREE_PRIME_POSITION'} + 1);
  $gibson_sequences -> [$gibson_index] -> {'NAME'} = dna_name($self) . '_' . $dna_from . '-' . $dna_to;   
  #Print hybridization unit output if in verbose mode
  print "Gibson Assembly Hybridization Units\n" if verbose($self);
  foreach my $hybrid_unit (@$hybridization) { 
    print '5prime pos '  . $hybrid_unit -> {'FIVE_PRIME_POSITION'} if verbose($self);
    print ' 3prime pos ' . $hybrid_unit -> {'THREE_PRIME_POSITION'} if verbose($self);
    print ' Length '     . $hybrid_unit -> {'LENGTH'} if verbose($self);
    print ' Tm '         . $hybrid_unit -> {'TM'} if verbose($self);
    print ' Sequence '   . $hybrid_unit -> {'SEQUENCE'} . "\n" if verbose($self);
  };
  #Build textbox of Gibson fragment sequences as fasta output
  print "Gibson fragment sequences\n" if verbose($self);
  my ($textbox, $cr, $gt) = (undef, chr(13).chr(10), chr(62));
  for ( my $i = 0; $i < @$gibson_sequences; $i++ ) { $textbox .= $gt . $gibson_sequences->[$i]->{'NAME'} . $cr . $gibson_sequences->[$i]->{'SEQUENCE'} . $cr; };
  print $textbox if verbose($self);
  textbox($self, $textbox);
  gibson_sequences($self, @$gibson_sequences);  
  return $self
}

#Searches for largest open reading frame in three forward and three reverse reading frames
#outputs DNA sequence of open reading frame and the protein orf itself to current sequence properties
sub find_largest_orf {
  my ($self, $dna) = @_;
  $dna = current_dna_sequence($self) unless $dna;  
  reverse_complement($self, $dna);
  my $rc_dna = current_dna_sequence($self);
  my ($translations, $proteins) = ({}, {});
  for (0..2) {
    translate($self, $dna, $_);
    my $fwd_key = 'FORWARD_FRAME_'.$_;
    $translations -> {$fwd_key} -> {PROTEIN_SEQUENCE}   = current_protein_sequence($self);
    $translations -> {$fwd_key} -> {DNA_SEQUENCE}       = $dna;
    $translations -> {$fwd_key} -> {FRAME}              = $_;
    $translations -> {$fwd_key} -> {FORWARD_OR_REVERSE} = 'FORWARD';    
    translate($self, $rc_dna, $_);
    my $rev_key = 'REVERSE_FRAME_'.$_;
    $translations -> {$rev_key} -> {PROTEIN_SEQUENCE}   = current_protein_sequence($self);
    $translations -> {$rev_key} -> {DNA_SEQUENCE}       = $rc_dna;
    $translations -> {$rev_key} -> {FRAME}              = $_;
    $translations -> {$rev_key} -> {FORWARD_OR_REVERSE} = 'REVERSE';   
  };
  foreach my $frame (@{ [ keys %$translations ] }) {
    my $translation = $translations -> {$frame} -> {PROTEIN_SEQUENCE};
    my @orfs        = $translation  =~ m/([M][^Z]*[Z])/gi; #regexp grabs all orf sequences with start Met and ending stop codon
    my @orf_at_end  = $translation  =~ m/([M][^Z]*)$/i; #regexp to grab orf at the end of raw reading frame translation string when stop codon is missing
    push @orfs, @orf_at_end; #add the two orf arrays together
    my $i = 0;
    foreach my $orf (@orfs) {
      my $orf_key = $frame.'_SEQ_'.$i++; 
      $proteins -> {$orf_key} -> {ORF}                = $orf;     
      $proteins -> {$orf_key} -> {PROTEIN_SEQUENCE}   = $translations -> {$frame} -> {PROTEIN_SEQUENCE};
      $proteins -> {$orf_key} -> {DNA_SEQUENCE}       = $translations -> {$frame} -> {DNA_SEQUENCE};
      $proteins -> {$orf_key} -> {FRAME}              = $translations -> {$frame} -> {FRAME};
      $proteins -> {$orf_key} -> {FORWARD_OR_REVERSE} = $translations -> {$frame} -> {FORWARD_OR_REVERSE};
    };
  };
  my @orf_keys = sort { length( $proteins -> {$b} -> {ORF} ) <=> length( $proteins -> {$a} -> {ORF} ) } @{ [ keys %$proteins ] };
  print "Searching for largest orf in " . scalar(@orf_keys) . " reading frames.\n" if verbose($self);
  my $orf_key = $orf_keys[0]; #largest orf hash key is first in the array
  my $largest_orf_protein = $proteins -> {$orf_key} -> {ORF};
  my $translation         = $proteins -> {$orf_key} -> {PROTEIN_SEQUENCE};
  my $dna_seq             = $proteins -> {$orf_key} -> {DNA_SEQUENCE};
  my $frame_position      = $proteins -> {$orf_key} -> {FRAME};
  #Next find position of starting methionine of largest orf in the raw translation of the reading frame
  my $met_position        = index($translation, $largest_orf_protein);
  print "Found largest orf at Met position $met_position. It is approximately ". length($largest_orf_protein) . " residues long.\n" if verbose($self);
  #Next extract the corresponding dna seq to the protein ORF starting of course at Met-1
  my $largest_orf_dna     = substr($dna_seq, $frame_position + $met_position * 3, length($largest_orf_protein) * 3); 
  current_protein_sequence($self, $largest_orf_protein);
  current_dna_sequence($self, $largest_orf_dna);
  translate($self);
  return $self
}

#Outputs to textbox the amino acid numbering, the amino acid residue, and codon
#An amino acid residue numbering offset can be set to account for amino terminal his6 tags or secretion sequences
sub verbose_translate {
  my $self = shift;
  my ($dna, $offset, $position, $chars_per_line) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  $chars_per_line = 80 unless $chars_per_line;
  $offset = amino_acid_numbers_offset($self) unless $offset;
  translate($self, current_dna_sequence($self, $dna), $position, amino_acid_numbers_offset($self, $offset));
  print "Verbose translation with amino acid numbering offset of $offset.\n" if verbose($self);
  my @codons = codons($self);
  my @amino_acids = amino_acids($self);
  my @amino_acid_numbers        = amino_acid_numbers($self);
  my @amino_acids_3_letter      = map { +{ amino_acid_3_letter() } -> {$_} } @amino_acids;
  my ($pipe, $cr, $space, $textbox)    = (chr(124), chr(13).chr(10), chr(32), undef);
  my @amino_acid_number_lengths = map { length } @amino_acid_numbers;
  my @amino_acid_number_pads    = map { if ($_ < 3) { 3 - $_ } else { 0 } } @amino_acid_number_lengths;
  my @amino_acid_pads           = map { if ($_ > 3) { $_ - 3 } else { 0 } } @amino_acid_number_lengths;
  my ($amino_acids_line, $codons_line, $numbers_line) = (undef, undef, undef); 
  my $i = 0;
  while ($i < scalar(@amino_acids)) {
    ($amino_acids_line, $codons_line, $numbers_line) = ($pipe, $pipe, $pipe);
    while ((length($amino_acids_line) <= $chars_per_line) && ($i < scalar(@amino_acids))) {
      $amino_acids_line .= $amino_acids_3_letter[$i] . $space x $amino_acid_pads[$i]        . $pipe;
      $codons_line      .= $codons[$i]               . $space x $amino_acid_pads[$i]        . $pipe;
      $numbers_line     .= $amino_acid_numbers[$i]   . $space x $amino_acid_number_pads[$i] . $pipe;
      $i++;
    };
    $numbers_line     .= $cr;
    $codons_line      .= $cr;
    $amino_acids_line .= $cr;
    $textbox .= $numbers_line . $amino_acids_line . $codons_line . $cr; 
  };
  print $textbox if verbose($self);
  textbox($self, $textbox);
  return $self
}

#Translates dna either passed as a variable or the current_dna_sequence
#starts at the DNA position passed as second argument to the method
#Stores resulting protein sequence in current_protein_sequence property
#stores outputs as arrays in the module for codons, amino_acids, and amino_acid_numbers
sub translate {
  my ($self, $dna, $position, $offset) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  current_dna_sequence($self, $dna);
  my ($protein, $codons, $amino_acids, $number) = (undef, [], [], 0);
  my $codon_table = +{ codon_table() };
  while ($position < length($dna)) { 
    my $codon = substr($dna, $position, 3);
    push @$codons, $codon;
    my $amino_acid = $codon_table -> {$codon};
    push @$amino_acids, $amino_acid;
    $protein .= $amino_acid;
    $position += 3;
  };
  current_protein_sequence($self, $protein);
  codons($self, @$codons);
  amino_acids($self, @$amino_acids);
  #Default residue numbering maps index zero to one so residue one is first amino acid and codon
  my $amino_acid_numbers = [ map { ++$number } @$amino_acids ]; #map the residue numbers indexed to starting value of one
  #Use offset to allow for proper amino acid numbering when a secretion tag or his6 tag
  #or other extraneous amino acids are appended or deleted from the amino terminus of the orf
  $offset = amino_acid_numbers_offset($self) unless $offset;
  amino_acid_numbers_offset($self, $offset);
  $amino_acid_numbers = [ map { $_ + $offset } @$amino_acid_numbers ] if $offset;   
  amino_acid_numbers($self, @$amino_acid_numbers);
  return $self
}

#Method takes a protein sequence and calculates a possible DNA sequence.
#When more than one codon are possible for a given amino acid, method randomly selects a possible codon.
#Stores new DNA sequence in current_dna_sequence property. Also creates and stores codons, amino_acids, and amino_acid_numbers array properties.
sub reverse_translate {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  current_protein_sequence($self, $protein);
  my $codon_table = +{ codon_table() }; #Codon to amino acid lookup table
  my ($number, $i, $dna, $amino_acids, $codons, $amino_acid_numbers) = (0, 0, undef, [], [], []);
  while ($i < length($protein)) {
    my $codon = undef;
    my $amino_acid = substr($protein, $i, 1);
    my @possible_codons = grep { $codon_table->{$_} eq $amino_acid } keys %$codon_table; #array of possible codons that match given amino acid
    my $number_of_codons = scalar(@possible_codons);   
    if ($number_of_codons == 1) { $codon = $possible_codons[0] }
    elsif ($number_of_codons > 1) { $codon = $possible_codons[ int(rand($number_of_codons)) ] };
    $dna .= $codon;
    push @$amino_acids, $codon_table -> {$codon} if $codon;
    push @$codons, $codon if $codon;
    $i++;
  };
  codons($self, @$codons);
  amino_acids($self, @$amino_acids);
  $amino_acid_numbers = [ map { ++$number } @$amino_acids ]; #map the residue numbers indexed to starting value of one
  my $offset = amino_acid_numbers_offset($self);
  $amino_acid_numbers = [ map { $_ + $offset } @$amino_acid_numbers ] if $offset; #if offset is set then re-map shifted amino acid numbers  
  amino_acid_numbers($self, @$amino_acid_numbers);  
  current_dna_sequence($self, $dna);
  return $self
}

#Method performs site-directed mutagenesis generating quikchange primers of optimized Tm
#mutagenesis primer should be 25 to 45 bases long, end and begin in a G or C, and have a target Tm of 78 degrees C
sub site_directed_mutagenesis {
  my ($self, $amino_acid_position_to_mutate, $amino_acid_mutation, $codon_mutation, $target_primer_tm) = @_;
  $amino_acid_mutation = uc($amino_acid_mutation);
  $codon_mutation = uc($codon_mutation);
  print "Starting site-directed mutagenesis at amino acid position $amino_acid_position_to_mutate, replacing with $amino_acid_mutation.\n" if verbose($self); 
  translate($self); #generate wild-type codons, amino acids, and amino acid numbering arrays in module
  my ($amino_acids, $codons, $amino_acid_numbers) = ( [ amino_acids($self) ], [ codons($self) ], [ amino_acid_numbers($self) ] );
  my $wild_type_dna = undef;
  foreach my $codon (@$codons) { $wild_type_dna .= $codon };
  #index in computer array starts at zero but residue positions typically start at one and depend on a specified offset
  #Next step finds index of the given real amino acid number
  my ($index_to_mutate) = grep { $amino_acid_numbers -> [$_] == $amino_acid_position_to_mutate } (0 .. @$amino_acid_numbers - 1);
  my $amino_acid_to_mutate = $amino_acids -> [ $index_to_mutate ];
  my $codon_to_mutate = $codons -> [ $index_to_mutate ];
  my ($mutant_amino_acids, $mutant_codons, $mutant_amino_acid_numbers) = ($amino_acids, $codons, $amino_acid_numbers);
  #Parse whether $amino_acid_mutation is one-letter or three-letter abbreviation
  my @possible_mutant_codons = (); #array of possible codons that match amino acid selected for mutation
  my $codon_table         = +{ codon_table() }; #Codon to amino acid lookup table
  my $amino_acid_3_letter = +{ amino_acid_3_letter() }; #one-letter to three-letter amino acid abbreviation lookup table
  if (length($amino_acid_mutation) == 1 && $amino_acid_mutation =~ /[ACDEFGHIKLMNPQRSTVWYZ]/) {
    @possible_mutant_codons = grep { $codon_table->{$_} eq $amino_acid_mutation } keys %$codon_table;
  }   
  elsif (length($amino_acid_mutation) == 3 && (grep { $_ eq $amino_acid_mutation } values %$amino_acid_3_letter) ) {
    my @one_letter_amino_acid_mutation = grep { $amino_acid_3_letter -> {$_} eq $amino_acid_mutation } keys %$amino_acid_3_letter;
    $amino_acid_mutation = $one_letter_amino_acid_mutation[0]; #grab corresponding one-letter abbreviation for amino_acid_mutation
    @possible_mutant_codons = grep { $codon_table->{$_} eq $amino_acid_mutation } keys %$codon_table;
  }
  else {
    die "Cannot recognize amino acid mutation in site_directed_mutagenesis. Use valid one-letter or three-letter amino acid abbreviation for mutation. Aborting.\n";
  };
  my $is_possible_codon_boolean = undef;	      
  $is_possible_codon_boolean = grep /$codon_mutation/, @possible_mutant_codons if $codon_mutation; #check if codon_mutation is in possible_mutant_codons
  if ($codon_mutation && !$is_possible_codon_boolean) { die "Codon desired for site-directed mutation does not match selected amino acid mutation. Aborting.\n"; }
  elsif (!$codon_mutation) {
    #Select optimal codon_mutation from possible_mutant_codons that with lowest mismatch score
    my $mismatch_score = 4; #start score at one greater than max possible mismatch score for a three-base codon
    foreach my $pos_codon (@possible_mutant_codons) {
      my $current_mismatch_score = ( $codon_to_mutate ^ $pos_codon ) =~ tr/\0//c; #Calculated as bit-wise XOR comparison where nulls count as mismatches
      if ($current_mismatch_score <= $mismatch_score) {
	$codon_mutation = $pos_codon;
        $mismatch_score = $current_mismatch_score;
      };
    };
  };
  #make amino acid and codon substitutions into arrays
  $mutant_amino_acids -> [ $index_to_mutate ] = $amino_acid_mutation;
  $mutant_codons -> [ $index_to_mutate ] = $codon_mutation;
  #Build mutant protein and dna text strings
  my $mutant_protein = undef;
  foreach my $aa (@$mutant_amino_acids) { $mutant_protein .= $aa };
  my $mutant_dna = undef;
  foreach my $codon (@$mutant_codons) { $mutant_dna .= $codon };
  #store mutant properties in module
  mutant_amino_acids($self, @$mutant_amino_acids);
  mutant_codons($self, @$mutant_codons);
  mutant_amino_acid_numbers($self, @$mutant_amino_acid_numbers);
  mutant_protein_sequence($self, $mutant_protein);
  mutant_dna_sequence($self, $mutant_dna);
  #Generate protein name and dna name for mutant
  my $mutant_protein_name = protein_name($self) . '_' . $amino_acid_to_mutate . $amino_acid_position_to_mutate . $amino_acid_mutation;
  my $mutant_dna_name     = dna_name($self)     . '_' . $amino_acid_to_mutate . $amino_acid_position_to_mutate . $amino_acid_mutation;
  mutant_protein_name($self, $mutant_protein_name);
  mutant_dna_name($self, $mutant_dna_name);
  #Design Tm-optimized Quikchange primers with GC clamp at 5' and 3' ends
  $target_primer_tm = 78 unless $target_primer_tm; #recommended default tm of 78 degrees C for quikchange mutagenesis
  my $minus_from_center = 12; #start out -12 bases 5' from center of mutagenized codon
  my $plus_from_center  = 12; #start out +12 bases 3' from center of mutagenized codon
  #get starting primer sequence
  my $primer_sequence = get_primer_sequence($mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
  #Pre-iterations of increasing primer length moving to 5' and 3' end GC clamps
  ($primer_sequence, $minus_from_center) = find_5_prime_gc_clamp($primer_sequence, $mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
  ($primer_sequence, $plus_from_center)  = find_3_prime_gc_clamp($primer_sequence, $mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center); 
  print "Generating Quikchange primer pair for mutation with annealing temperature of $target_primer_tm C.\n" if verbose($self); 
  #calculate primer tm
  calculate_primer_tm($self, $primer_sequence, $wild_type_dna);
  my $tm = primer_tm($self);
  while ($tm < $target_primer_tm) {
    if    ($minus_from_center > $plus_from_center)  {
      $plus_from_center++;
      $primer_sequence = get_primer_sequence($mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
      ($primer_sequence, $plus_from_center)  = find_3_prime_gc_clamp($primer_sequence, $mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
    }
    elsif ($minus_from_center <= $plus_from_center)  {
      $minus_from_center++;
      $primer_sequence = get_primer_sequence($mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
      ($primer_sequence, $minus_from_center) = find_5_prime_gc_clamp($primer_sequence, $mutant_dna, $index_to_mutate, $minus_from_center, $plus_from_center);
    };
   calculate_primer_tm($self, $primer_sequence, $wild_type_dna);
   $tm = primer_tm($self);   
  };
  #output primer information to textbox
  #Update fwd and rev primer names with point mutation appended to the names
  primer_forward_name($self, $mutant_protein_name . '.F');
  primer_reverse_name($self, $mutant_protein_name . '.R');
  #make new textbox of primer pair properties
  my ($textbox, $cr, $prime, $degree, $pipe) = (undef, chr(13).chr(10), chr(39), chr(176), chr(124));
  $textbox .= primer_forward_name($self) . $pipe . '5' . $prime . '-' . primer_forward_sequence($self) . '-3' . $prime . $cr;
  $textbox .= primer_reverse_name($self) . $pipe . '5' . $prime . '-' . primer_reverse_sequence($self) . '-3' . $prime . $cr;
  $textbox .= 'Primer length: ' . primer_length($self) . ' bases'. $cr; 
  $textbox .= 'Melting midpoint (Tm): ' . primer_tm($self) . $degree . 'C' . $cr;
  $textbox .= 'Percent GC: ' . primer_percent_gc($self) . $cr;
  $textbox .= 'Percent mismatch: ' . primer_percent_mismatch($self) . $cr;
  textbox($self, $textbox);
  print $textbox if verbose($self);
  return $self
}  

sub get_primer_sequence {
  my ($dna, $index_to_mutate, $minus_from_center, $plus_from_center) = @_;
  my $dna_center_position  = $index_to_mutate * 3; #find start position of mutagenized codon in dna
  my $dna_5_prime_position = $dna_center_position - $minus_from_center;
  my $primer_length = $minus_from_center + $plus_from_center + 3;
  #return primer sequence as substr of dna
  return substr($dna, $dna_5_prime_position, $primer_length)
}

sub find_5_prime_gc_clamp {
 my ($primer_sequence, $dna, $index_to_mutate, $minus_from_center, $plus_from_center) = @_;
  #walk base by base to a GC clamp on the 5' end of primer sequence
  my $five_prime_base = substr($primer_sequence, 0, 1);
  while ($five_prime_base =~ /[AT]/ ) {
    $minus_from_center++;
    $primer_sequence = get_primer_sequence($dna, $index_to_mutate, $minus_from_center, $plus_from_center);
    $five_prime_base = substr($primer_sequence, 0, 1);
  };
  return $primer_sequence, $minus_from_center
}

sub find_3_prime_gc_clamp {
 my ($primer_sequence, $dna, $index_to_mutate, $minus_from_center, $plus_from_center) = @_;
  #walk base by base to a GC clamp on the 3' end of primer sequence
  my $three_prime_base = substr($primer_sequence, -1, 1);
  while ($three_prime_base =~ /[AT]/ ) {
    $plus_from_center++;
    $primer_sequence = get_primer_sequence($dna, $index_to_mutate, $minus_from_center, $plus_from_center);
    $three_prime_base = substr($primer_sequence, -1, 1);
  };
  return $primer_sequence, $plus_from_center
}

#moves mutant amino acid, codons and numbering and dna and protein sequences to current arrays and scalars
#useful to further mutagenize a mutant sequence or verbose translate the mutant to check if the point mutation was made correctly
sub transfer_mutant_to_current {
  my $self = shift;
  print "Transferring mutant DNA and protein properties to current sequences and arrays.\n" if verbose($self);
  amino_acids($self, mutant_amino_acids($self));
  codons($self, mutant_codons($self));
  amino_acid_numbers($self, mutant_amino_acid_numbers($self));
  current_protein_sequence($self, mutant_protein_sequence($self));
  current_dna_sequence($self, mutant_dna_sequence($self));
  protein_name($self, mutant_protein_name($self));
  dna_name($self, mutant_dna_name($self));
  return $self
}

#Calculates primer melting midpoint for PCR and quikchange site-directed mutagenesis
#Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch
#where N = total number of bases
sub calculate_primer_tm {
  my ($self, $oligo, $dna) = @_; #dna is the target dna sequence. oligo is the primer
  $dna = current_dna_sequence($self) unless $dna;
  current_dna_sequence($self, $dna);
  $oligo = uc($oligo);
  #make reverse complement primer
  my $reverse_oligo = reverse($oligo);
  $reverse_oligo =~ tr/[A,C,G,T]/[t,g,c,a]/; #complement reversed oligo
  $reverse_oligo = uc($reverse_oligo); #Back to uppercase sequence
  #calculate n and percent_gc
  my $n = length($oligo);
  my $g_count = ($oligo =~ s/G/G/g);
  my $c_count = ($oligo =~ s/C/C/g);
  my $gc_count = $g_count + $c_count;
  my $percent_gc = undef; 
  $percent_gc = 100 * $gc_count / $n if $n;
  $percent_gc = sprintf("%.1f", $percent_gc); #round to one decimal place
  #Determine percent mismatch by bitwise string comparison of oligo to target dna
  my $percent_mismatch = undef;
  my $largest_match = longest_common_sequence($oligo, $dna); #find anneal site in target dna
  my $start_position_of_match_in_oligo = index($oligo, $largest_match);
  my $start_position_of_match_in_dna   = index($dna, $largest_match);
  my $five_prime_position_of_dna = $start_position_of_match_in_dna - $start_position_of_match_in_oligo; 
  my $target_dna_fragment = substr($dna, $five_prime_position_of_dna, length($oligo)); #extract target dna fragment same length as primer oligo
  #Bitwise XOR string comparison, where regular expression counts instances of nulls as mismatches
  #note comparison is good for comparing two strings of the same length
  my $mismatch_count = ( $oligo ^ $target_dna_fragment ) =~ tr/\0//c; 
  $percent_mismatch = 100 * $mismatch_count / $n if $n;
  $percent_mismatch = sprintf("%.1f", $percent_mismatch); #round to one decimal place
  #Calculate tm
  my $tm = undef;
  $tm = 81.5 + 0.41 * $percent_gc - 675 / $n - $percent_mismatch if $n;
  $tm = sprintf("%.1f", $tm); #round to one decimal place
  #Store primer properties in module
  primer_forward_sequence($self, $oligo);  
  primer_reverse_sequence($self, $reverse_oligo);
  primer_tm($self, $tm);
  primer_percent_gc($self, $percent_gc);
  primer_length($self, $n);    
  primer_percent_mismatch($self, $percent_mismatch);
  my ($forward_name, $reverse_name) = (dna_name($self).'.F', dna_name($self).'.R');
  primer_forward_name($self, $forward_name);                  
  primer_reverse_name($self, $reverse_name);
  my ($textbox, $cr, $prime, $degree, $pipe) = (undef, chr(13).chr(10), chr(39), chr(176), chr(124));
  $textbox .= $forward_name.$pipe.'5'.$prime.'-'.$oligo.'-3'.$prime.$cr;
  $textbox .= $reverse_name.$pipe.'5'.$prime.'-'.$reverse_oligo.'-3'.$prime.$cr;
  $textbox .= "Primer length: $n$cr"; 
  $textbox .= 'Melting midpoint (Tm): '.$tm.$degree.'C'.$cr;
  $textbox .= "Percent GC: $percent_gc$cr";
  $textbox .= "Percent mismatch: $percent_mismatch$cr";
  textbox($self, $textbox);
  print $textbox if verbose($self);
  return $self
}

sub find_restriction_sites {
  my ($self, $dna)  = @_;
  $dna = current_dna_sequence($self) unless $dna;
  current_dna_sequence($self, $dna);
  print "Searching for resistriction sites.\n" if verbose($self);
  my $rebase = open_rebase(); 
  my $hash = {};
  #search for each enzyme site or sites in restriction enzyme database hash 
  foreach my $enz ( @{ [ keys %$rebase ] } ) {
    my $sites = $rebase -> {$enz} -> {'SITE'};
    my $positions = [];
    foreach my $site (@$sites) {
      #find all instances of $site sequence in $dna sequence
      my $offset   = 0;
      my $position = index($dna, $site, $offset);
      while ($position != -1) {   
        push @$positions, $position + 1; #dna position is one greater than computer position which indexes to zero
	$offset   = $position + 1; #move offset to allow for additional sites to be captured
	$position = index($dna, $site, $offset); #search next position
      };      
    };
    if (scalar(@$positions)) {
      $hash -> {$enz} -> {'NAME'}           = $rebase -> {$enz} -> {'NAME'};
      $hash -> {$enz} -> {'DEGEN_CUT_SITE'} = $rebase -> {$enz} -> {'DEGEN_CUT_SITE'};
      $hash -> {$enz} -> {'POSITION'}       = $positions;
    };
  };
  #Output the restriction enzyme sites to a textbox which can be saved with save_textbox
  my ($textbox, $cr, $space, $pound) = (undef, chr(13).chr(10), chr(32), chr(35));
  $textbox .= "$pound Restriction Enzyme Map$cr";
  foreach my $enz ( @{ [ sort { $a cmp $b } @{ [ keys %$hash ] } ] } ) {
    $textbox .= $hash -> {$enz} -> {'NAME'} . $space
             . $hash -> {$enz} -> {'DEGEN_CUT_SITE'} . $space
             . join($space, @{ [ sort { $a <=> $b } @{ $hash -> {$enz} -> {'POSITION'} } ] } ) . $cr;
  };
  print $textbox if verbose($self);
  textbox($self, $textbox);
  return $self
}

#DNA Restriction Enzyme Utility for Finding Cut Sites
#REBASE file format is DNA Strider format delimited by commas, where the header is commented out with #'s
#also any restriction enzyme commented out with a # is not added to database
sub open_rebase {
    my ($file, $pound, $comma, $rebase) = (shift, chr(35), chr(44), {});
    $file = 'rebase.txt' unless $file;
    print "Opening Rebase database in file $file.\n" if verbose($self);
    open(FH, "<$file") || die "Cannot open Rebase file, $file. $!";
    #EXAMPLE: TseI,g/cwgc,
    while (my $line = <FH>) {
        if (($line =~ /$comma[RYMKSWBDHVNATCG\/]+$comma/i) && ($line !~ /^$pound/)) {
            my ($name, $degen_site) = split(/$comma/, $line);
            $degen_site = uc($degen_site);
            my $degen_site_with_cut_site = $degen_site;
            $degen_site =~ s/\///; #remove slash / char which marks the cut site in REBASE DNA Strider format
#            print "Restriction enzyme: $name\t$degen_site_with_cut_site\t$degen_site\n";
            my $site = calculate_degenerate_sequences($degen_site);
#            print scalar(@$site)."numbdegen_site\n";
            my $key = uc($name);
            $rebase->{$key}->{'NAME'} = $name;  
            $rebase->{$key}->{'DEGEN_SITE'} = $degen_site;
            $rebase->{$key}->{'DEGEN_CUT_SITE'} = $degen_site_with_cut_site;
            $rebase->{$key}->{'SITE'} = [@$site];
        }
    }
    close(FH);
    return $rebase
}

sub calculate_degenerate_sequences {
    my ($nnn, $degen, $seqs) = (shift, [], []);
    my $degen_code = +{ degenerate_code() };
    #do loop until $degen does not have any more degen bases
    if ($nnn =~ /[RYMKSWBDHVN]/) {
        push @$degen, $nnn;
        do {
            my $degen_seq = shift(@$degen);
            my ($degen_base) = ($degen_seq =~ /([RYMKSWBDHVN])/);
            my $pos = index $degen_seq, $degen_base;
            my $poss_bases = $degen_code->{uc($degen_base)};
            foreach ( @{ [ split(//,$poss_bases) ] } ) {
                my $new_seq = uc(substr($degen_seq,0,$pos).$_.substr($degen_seq,$pos+1));
                if ($new_seq =~ /[RYMKSWBDHVN]/) { push @$degen, $new_seq } else { push @$seqs, $new_seq }
            }
        } while scalar(@$degen)
    }
    else { return [$nnn] };
    return [@$seqs]
}

sub calculate_protein_parameters {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  print "Calculating protein parameters.\n" if verbose($self);
  fasta_protein($self, $protein);
  calculate_protein_length($self);
  calculate_molecular_weight($self);
  calculate_extinction_coefficient($self);
  calculate_isoelectric_point($self);
  return $self
}

sub calculate_molecular_weight {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  current_protein_sequence($self, $protein);
  my ($mw, $position) = (0, 0);
  my $amino_acid_molecular_weight = +{ amino_acid_molecular_weight() };
  while ($position < length($protein)) {
    my $residue = substr($protein, $position++, 1);
    #mw data in residue_mw hash are for amino acids not residue mass, hence subtract a water (18 Da)
    my $residue_mw = $amino_acid_molecular_weight -> {$residue} - 18.01528; 
    $mw += $residue_mw if ($residue ne 'Z'); #don't add weight or do anything for the stop codon
  };
  $mw += 18.01528 if ($mw > 0); #add weight of N- and C-termini, which is equiv. to a water
  $mw = sprintf("%.1f", $mw);
  print "Calculated molecular weight is $mw Daltons\n" if verbose($self);
  molecular_weight($self, $mw);
  return $self
}

sub calculate_protein_length {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  current_protein_sequence($self, $protein);  
  $protein =~ s/[Z]//g; #removes stop codon 'Z' prior to counting residues with length
  protein_length($self, length($protein));
  print "Protein is " . length($protein) . " residues.\n" if verbose($self);
  return $self
}

#Extinction coefficients of Tyr and Trp are from https://web.expasy.org/protparam/
sub calculate_extinction_coefficient {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  current_protein_sequence($self, $protein);
  my ($ext_coef, $position, $aa_ext_coefs) = (0, 0, {});
  $aa_ext_coefs -> {W} = 5500;
  $aa_ext_coefs -> {Y} = 1490;
  while ($position < length($protein)) { $ext_coef += $aa_ext_coefs -> { substr($protein, $position++, 1) } };
  print "Extinction coefficient is $ext_coef M-1 cm-1\n" if verbose($self);
  extinction_coefficient($self, $ext_coef);
  return $self
}

sub calculate_isoelectric_point {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  current_protein_sequence($self, $protein);  
  my ($ph, $pi, $min_absolute_charge) = (1, 1, 5);
  while ($ph >= 1 && $ph <= 13) {
    my $alpha_amino_charge    = calculate_charge('ALPHA_AMINO', $ph);   
    my $alpha_carboxyl_charge = calculate_charge('ALPHA_CARBOXYL', $ph);
    my $residues_charge = 0;
    my $position = 0;
    while ($position < length($protein)) { $residues_charge += calculate_charge(substr($protein, $position++, 1), $ph) };
    my $current_absolute_charge = abs( $alpha_amino_charge + $alpha_carboxyl_charge + $residues_charge );    
    if ($current_absolute_charge < $min_absolute_charge) {
      $pi = $ph;
      $min_absolute_charge = $current_absolute_charge;
    };
    $ph += 0.03;
  };
  $pi = sprintf("%.2f", $pi);
  print "Isoelectric point is $pi\n" if verbose($self);
  isoelectric_point($self, $pi);
  return $self
}

sub calculate_charge {
  my ($group, $ph) = @_;
  my $amino_acid_pka = +{ amino_acid_pka() };
  my $max_charge = $amino_acid_pka -> {$group} -> {CHARGE};
  my $pka        = $amino_acid_pka -> {$group} -> {PKA};
  my $charge = 0;
  if    ($max_charge ==  1) { $charge = $max_charge / (1 + 10**($ph-$pka)) }
  elsif ($max_charge == -1) { $charge = $max_charge / (1 + 10**($pka-$ph)) }
  else { $charge = 0 };
  return $charge
}

sub reverse_complement {
  my ($self, $dna) = @_;
  $dna = current_dna_sequence($self) unless $dna;
  complement($self, $dna);
  current_dna_sequence($self, reverse(current_dna_sequence($self)));
  return $self 
}

sub complement {
  my ($self, $dna) = @_;
  $dna =  current_dna_sequence($self) unless $dna;
  $dna =~ tr/[A,C,G,T]/[t,g,c,a]/;
  current_dna_sequence($self, uc($dna));
  return $self
}

sub fasta_dna {
  my ($self, $dna) = @_;
  $dna = current_dna_sequence($self) unless $dna;  
  $dna =~ s/[^ATCG]//gi;
  current_dna_sequence($self, $dna);
  return $self  
}

sub fasta_protein {
  my ($self, $protein) = @_;
  $protein = current_protein_sequence($self) unless $protein;
  $protein =~ s/[^ACDEFGHIKLMNPQRSTVWYZ]//gi;
  current_protein_sequence($self, $protein);  
  return $self
}

#remove trailing white spaces, i.e. the carriage return and linefeed for inputing PC/DOS compatible file data
sub ts { my $i = shift; $i =~ s/\s+$//; return $i }

############################# PROPERTY SETS/GETS ###################################

### ARRAY PROPERTIES

sub codons {
  my $self = shift;
  if (@_) { @{ $self->{CODONS} } = @_ };  
  return @{ $self->{CODONS} }
}

sub amino_acids {
  my $self = shift;
  if (@_) { @{ $self->{AMINO_ACIDS} } = @_ };
  return @{ $self->{AMINO_ACIDS} }
}

sub amino_acid_numbers {
  my $self = shift;
  if (@_) { @{ $self->{AMINO_ACID_NUMBERS} } = @_ };
  return @{ $self->{AMINO_ACID_NUMBERS} }
}

sub mutant_amino_acids {
  my $self = shift;
  if (@_) { @{ $self->{MUTANT_AMINO_ACIDS} } = @_ };
  return @{ $self->{MUTANT_AMINO_ACIDS} }
}

sub mutant_codons {
  my $self = shift;
  if (@_) { @{ $self->{MUTANT_CODONS} } = @_ };
  return @{ $self->{MUTANT_CODONS} }
}

sub mutant_amino_acid_numbers {
  my $self = shift;
  if (@_) { @{ $self->{MUTANT_AMINO_ACID_NUMBERS} } = @_ };
  return @{ $self->{MUTANT_AMINO_ACID_NUMBERS} } 
}

sub oligos {
  my $self = shift;
  if (@_) { @{ $self->{OLIGOS} } = @_ };
  return @{ $self->{OLIGOS} }
}

sub gibson_sequences {
  my $self = shift;
  if (@_) { @{ $self->{GIBSON_SEQUENCES} } = @_ };
  return @{ $self->{GIBSON_SEQUENCES} }
}

sub protein_library {
  my $self = shift;
  if (@_) { @{ $self->{PROTEIN_LIBRARY} } = @_ }; 
  return @{ $self->{PROTEIN_LIBRARY} }
}

### SCALAR PROPERTIES

sub record_file {
  my $self = shift;
  if (@_) { $self->{RECORD_FILE} = shift };
  return $self->{RECORD_FILE}
}

sub dna_file {
  my $self = shift;
  if (@_) { $self->{DNA_FILE} = shift };
  return $self->{DNA_FILE}
}

sub protein_file {
  my $self = shift;
  if (@_) { $self->{PROTEIN_FILE} = shift };
  return $self->{PROTEIN_FILE}
}

sub protein_library_file {
  my $self = shift;
  if (@_) { $self->{PROTEIN_LIBRARY_FILE} = shift };
  return $self->{PROTEIN_LIBRARY_FILE}
}

sub current_dna_sequence {
  my $self = shift;
  if (@_) { $self->{CURRENT_DNA_SEQUENCE} = uc(shift) };
  return uc($self->{CURRENT_DNA_SEQUENCE})
}

sub current_protein_sequence {
  my $self = shift;
  if (@_) { $self->{CURRENT_PROTEIN_SEQUENCE} = uc(shift) };
  return uc($self->{CURRENT_PROTEIN_SEQUENCE})
}

sub mutant_dna_sequence {
  my $self = shift;
  if (@_) { $self->{MUTANT_DNA_SEQUENCE} = shift };
  return $self->{MUTANT_DNA_SEQUENCE}
}

sub mutant_protein_sequence {
  my $self = shift;
  if (@_) { $self->{MUTANT_PROTEIN_SEQUENCE} = shift };
  return $self->{MUTANT_PROTEIN_SEQUENCE}
}

sub dna_name {
  my $self = shift;
  if (@_) { $self->{DNA_NAME} = shift };
  return $self->{DNA_NAME}
}

sub protein_name {
  my $self = shift;
  if (@_) { $self->{PROTEIN_NAME} = shift };
  return $self->{PROTEIN_NAME}
}

sub mutant_dna_name {
  my $self = shift;
  if (@_) { $self->{MUTANT_DNA_NAME} = shift };
  return $self->{MUTANT_DNA_NAME}
}

sub mutant_protein_name {
  my $self = shift;
  if (@_) { $self->{MUTANT_PROTEIN_NAME} = shift };    
  return $self->{MUTANT_PROTEIN_NAME}
}

sub amino_acid_numbers_offset {
  my $self = shift;
  if (@_) { $self->{AMINO_ACID_NUMBERS_OFFSET} = shift };
  return $self->{AMINO_ACID_NUMBERS_OFFSET}
}

sub textbox {
  my $self = shift;
  if (@_) { $self->{TEXTBOX} = shift };
  return $self->{TEXTBOX}
}

sub primer_forward_sequence {
  my $self = shift;
  if (@_) { $self->{PRIMER_FORWARD_SEQUENCE} = shift };
  return $self->{PRIMER_FORWARD_SEQUENCE}   
}
  
sub primer_reverse_sequence {
  my $self = shift;
  if (@_) { $self->{PRIMER_REVERSE_SEQUENCE} = shift };
  return $self->{PRIMER_REVERSE_SEQUENCE}  
}

sub primer_tm {
  my $self = shift;
  if (@_) { $self->{PRIMER_TM} = shift };
  return $self->{PRIMER_TM}
}

sub primer_percent_gc {
  my $self = shift;
  if (@_) { $self->{PRIMER_PERCENT_GC} = shift };
  return $self->{PRIMER_PERCENT_GC}
}

sub primer_length {
  my $self = shift;  
  if (@_) { $self->{PRIMER_LENGTH} = shift };
  return $self->{PRIMER_LENGTH} 
}

sub primer_percent_mismatch {
  my $self = shift;
  if (@_) { $self->{PRIMER_PERCENT_MISMATCH} = shift };
  return $self->{PRIMER_PERCENT_MISMATCH}
}

sub primer_forward_name {
  my $self = shift;
  if (@_) { $self->{PRIMER_FORWARD_NAME} = shift };
  return $self->{PRIMER_FORWARD_NAME}
}
                  
sub primer_reverse_name {
  my $self = shift;
  if (@_) { $self->{PRIMER_REVERSE_NAME} = shift };
  return $self->{PRIMER_REVERSE_NAME}
}

sub protein_length {
  my $self = shift;
  if (@_) { $self->{PROTEIN_LENGTH} = shift };
  return $self->{PROTEIN_LENGTH}
}

sub extinction_coefficient {
  my $self = shift;
  if (@_) { $self->{EXTINCTION_COEFFICIENT} = shift };
  return $self->{EXTINCTION_COEFFICIENT}
}

sub molecular_weight {
  my $self = shift;
  if (@_) { $self->{MOLECULAR_WEIGHT} = shift };
  return $self->{MOLECULAR_WEIGHT}
}

sub isoelectric_point {
  my $self = shift;
  if (@_) { $self->{ISOELECTRIC_POINT} = shift };
  return $self->{ISOELECTRIC_POINT}
}

sub verbose {
  my $self = shift;
  if (@_) { $self->{VERBOSE} = shift };
  return $self->{VERBOSE}
}

### LOOKUP TABLES OF NUCLEIC ACID AND PROTEIN DATA

sub codon_table {
    my %codon_table = (
	"TTT" => "F",
	"TTC" => "F",
	"TTA" => "L",
	"TTG" => "L",
	"TCT" => "S",
	"TCC" => "S",
	"TCA" => "S",
	"TCG" => "S",
	"TAT" => "Y",
	"TAC" => "Y",
	"TAA" => "Z",
	"TAG" => "Z",
	"TGT" => "C",
	"TGC" => "C",
	"TGA" => "Z",
	"TGG" => "W",
	"CTT" => "L",
	"CTC" => "L",
	"CTA" => "L",
	"CTG" => "L",
	"CCT" => "P",
	"CCC" => "P",
	"CCA" => "P",
	"CCG" => "P",
	"CAT" => "H",
	"CAC" => "H",
	"CAA" => "Q",
	"CAG" => "Q",
	"CGT" => "R",
	"CGC" => "R",
	"CGA" => "R",
	"CGG" => "R",
	"ATT" => "I",
	"ATC" => "I",
	"ATA" => "I",
	"ATG" => "M",
	"ACT" => "T",
	"ACC" => "T",
	"ACA" => "T",
	"ACG" => "T",
	"AAT" => "N",
	"AAC" => "N",
	"AAA" => "K",
	"AAG" => "K",
	"AGT" => "S",
	"AGC" => "S",
	"AGA" => "R",
	"AGG" => "R",
	"GTT" => "V",
	"GTC" => "V",
	"GTA" => "V",
	"GTG" => "V",
	"GCT" => "A",
	"GCC" => "A",
	"GCA" => "A",
	"GCG" => "A",
	"GAT" => "D",
	"GAC" => "D",
	"GAA" => "E",
	"GAG" => "E",
	"GGT" => "G",
	"GGC" => "G",
	"GGA" => "G",
	"GGG" => "G",
    );
    return %codon_table
};

#amino acid molecular weights are for the whole amino acid not the residue mass
#subtract 18 to get the residue mass
sub amino_acid_molecular_weight {
  my $amino_acid_molecular_weight = {};
  $amino_acid_molecular_weight -> {A} =	89.09;
  $amino_acid_molecular_weight -> {R} =	174.20;
  $amino_acid_molecular_weight -> {N} =	132.12;
  $amino_acid_molecular_weight -> {D} =	133.10;
  $amino_acid_molecular_weight -> {C} =	121.15;
  $amino_acid_molecular_weight -> {E} =	147.13;
  $amino_acid_molecular_weight -> {Q} =	146.15;
  $amino_acid_molecular_weight -> {G} =	75.07;
  $amino_acid_molecular_weight -> {H} =	155.16;
  $amino_acid_molecular_weight -> {I} =	131.17;
  $amino_acid_molecular_weight -> {L} =	131.17;
  $amino_acid_molecular_weight -> {K} =	146.19;
  $amino_acid_molecular_weight -> {M} =	149.21;
  $amino_acid_molecular_weight -> {F} =	165.19;
  $amino_acid_molecular_weight -> {P} =	115.13;
  $amino_acid_molecular_weight -> {S} =	105.09;
  $amino_acid_molecular_weight -> {T} =	119.12;
  $amino_acid_molecular_weight -> {W} =	204.23;
  $amino_acid_molecular_weight -> {Y} =	181.19;
  $amino_acid_molecular_weight -> {V} =	117.15;
  $amino_acid_molecular_weight -> {Z} =	0;  
  return %$amino_acid_molecular_weight
}

#Three-letter amino acid names
sub amino_acid_3_letter {
  my $amino_acid_3_letter = {};
  $amino_acid_3_letter -> {A} =	'ALA';
  $amino_acid_3_letter -> {R} =	'ARG';
  $amino_acid_3_letter -> {N} =	'ASN';
  $amino_acid_3_letter -> {D} = 'ASP';
  $amino_acid_3_letter -> {C} =	'CYS';
  $amino_acid_3_letter -> {E} =	'GLU';
  $amino_acid_3_letter -> {Q} =	'GLN';
  $amino_acid_3_letter -> {G} = 'GLY';
  $amino_acid_3_letter -> {H} = 'HIS';
  $amino_acid_3_letter -> {I} =	'ILE';
  $amino_acid_3_letter -> {L} =	'LEU';
  $amino_acid_3_letter -> {K} =	'LYS';
  $amino_acid_3_letter -> {M} =	'MET';
  $amino_acid_3_letter -> {F} =	'PHE';
  $amino_acid_3_letter -> {P} =	'PRO';
  $amino_acid_3_letter -> {S} =	'SER';
  $amino_acid_3_letter -> {T} =	'THR';
  $amino_acid_3_letter -> {W} =	'TRP';
  $amino_acid_3_letter -> {Y} =	'TYR';
  $amino_acid_3_letter -> {V} =	'VAL';
  $amino_acid_3_letter -> {Z} =	'STP';  
  return %$amino_acid_3_letter
}

#amino acid pKa values are from Bjellqvist
# Bjellqvist B, Basse B, Olsen E, Celis JE. Reference points for comparisons of two‐dimensional maps of proteins
# from different human cell types defined in a pH scale where isoelectric points correlate with polypeptide compositions.
# Electrophoresis. 1994;15(1):529–39.
sub amino_acid_pka {
  my $amino_acid_pka = {};
  #Alpha amino and carboxyl group pKas for N- and C-termini
  $amino_acid_pka -> {ALPHA_AMINO} -> {PKA}       = 7.50;
  $amino_acid_pka -> {ALPHA_CARBOXYL} -> {PKA}    = 3.55;
  #sidechain pKas
  $amino_acid_pka -> {C} -> {PKA}                 = 9.00;
  $amino_acid_pka -> {D} -> {PKA}                 = 4.05;
  $amino_acid_pka -> {E} -> {PKA}                 = 4.45;
  $amino_acid_pka -> {H} -> {PKA}                 = 5.98;
  $amino_acid_pka -> {K} -> {PKA}                 = 10.00;
  $amino_acid_pka -> {R} -> {PKA}                 = 12.00;
  $amino_acid_pka -> {Y} -> {PKA}                 = 10.00;
  #charges of ionized group
  $amino_acid_pka -> {ALPHA_AMINO} -> {CHARGE}    =  1;
  $amino_acid_pka -> {ALPHA_CARBOXYL} -> {CHARGE} = -1;
  $amino_acid_pka -> {C} -> {CHARGE}              = -1;
  $amino_acid_pka -> {D} -> {CHARGE}              = -1;
  $amino_acid_pka -> {E} -> {CHARGE}              = -1;
  $amino_acid_pka -> {H} -> {CHARGE}              =  1;
  $amino_acid_pka -> {K} -> {CHARGE}              =  1;
  $amino_acid_pka -> {R} -> {CHARGE}              =  1;
  $amino_acid_pka -> {Y} -> {CHARGE}              = -1; 
  return %$amino_acid_pka
}

# These degenerate bases (RYMKSWBDHVN) are weird and translate to: ATGC
# Recognition sequences for restriction enzymes are given using these standard abbreviations
# (Eur. J. Biochem. 150: 1-5, 1985) to represent ambiguity.
#R = G or A 
#Y = C or T
#M = A or C
#K = G or T
#S = G or C
#W = A or T
#B = not A (C or G or T)
#D = not C (A or G or T)
#H = not G (A or C or T)
#V = not T (A or C or G)
#N = A or C or G or T

sub degenerate_code {
  my $degen_code = {};
  $degen_code->{R} = 'GA';
  $degen_code->{Y} = 'CT';
  $degen_code->{M} = 'AC';
  $degen_code->{K} = 'GT';
  $degen_code->{S} = 'GC';
  $degen_code->{W} = 'AT';
  $degen_code->{B} = 'CGT';
  $degen_code->{D} = 'AGT';
  $degen_code->{H} = 'ACT';
  $degen_code->{V} = 'ACG';
  $degen_code->{N} = 'ACGT';
  return %$degen_code
}

1; #always return '1' at end of a Perl module
