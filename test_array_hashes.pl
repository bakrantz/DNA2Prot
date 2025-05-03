my $array_hashes = [];

$array_hashes -> [0] -> {'KEY1'} = 1;
$array_hashes -> [0] -> {'KEY2'} = 2;
$array_hashes -> [0] -> {'KEY3'} = 3;

$array_hashes -> [1] -> {'KEY1'} = 4;
$array_hashes -> [1] -> {'KEY2'} = 5;
$array_hashes -> [1] -> {'KEY3'} = 6;

$array_hashes -> [2] -> {'KEY1'} = 7;
$array_hashes -> [2] -> {'KEY2'} = 8;
$array_hashes -> [2] -> {'KEY3'} = 9;

print "Number of hashes in array is ". scalar(@$array_hashes)  ."\n";

my $index = 0;

for (my $i = 0; $i < @$array_hashes; $i++) {
  my $value = $array_hashes -> [$i] -> {'KEY1'};
  print "index is $index value is $value\n";
  $index++;
};
print "Index is $index\n";

print "Reverse order\n";
for (my $i = @$array_hashes - 1; $i >= 0; $i--) {
  my $value = $array_hashes -> [$i] -> {'KEY1'};
  print "index is $index value is $value\n"; 
};
