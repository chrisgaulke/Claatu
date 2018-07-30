#!/usr/bin/perl


use strict; use warnings; 

<>; 

while (my $line = <>){
	chomp $line; 
	my @array = split "\t", $line; 
	my $seq   = $array[0]; 
	my $tax   = join "; ", ("k__$array[1]", "p__$array[2]", "c__$array[3]", "o__$array[4]", "f__$array[5]", "g__$array[6]", "s__"); 
	print "$seq\t$tax\n"; 
}


