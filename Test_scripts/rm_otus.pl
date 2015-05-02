#!/usr/bin/perl 


##############
#            #
# rm_otus.pl #
#            #
##############

use strict; use warnings; 
use Data::Dumper;

# remove select OTUs from greengene files

# in_rep: input file gg_rep_set 
# in_aln: input file gg_rep_set_aln
# out_rep: output rep_set file clean of duplicates
# out_aln: output aln file clean of duplicates
# names_file: output file that contains identical OTUs 

my $usage = "Usage: rm_otus.pl <in_rep> <in_aln> <out_rep> <out_aln> <names_file>";

die $usage if (@ARGV != 5); 

my $in_rep  = $ARGV[0];
my $in_aln  = $ARGV[1];
my $out_rep = $ARGV[2];
my $out_aln = $ARGV[3];
my $names   = $ARGV[4];

open(my $aln, '<', $in_aln)       or die "$! ... Can't open the $in_aln file\n";

print "prepping hashes\n"; 
 
#my %aln_hash;
#my %names; 
my $aln_hash;
my $names_hash; 

while(my $header = <$aln>){
	chomp $header;
	my $seq = <$aln>;
	chomp $seq; 
	if ( exists $aln_hash->{$seq} ){
		push(@{$names_hash->{$aln_hash->{$seq}}}, $header);
		next;
	}else{
		$aln_hash->{$seq} = $header;
		$names_hash->{$header} = [($header)];
	}
}
close($aln);

#print Dumper($aln_hash);
#print Dumper($names_hash);


#print out the cleaned alignment
print "writing $out_aln\n"; 
open(my $aln_out, '>', $out_aln)  or die "$! ... Can't open the $out_aln file\n";

while( my($key, $value) = each %$aln_hash){
	print $aln_out "$value\n";
	print $aln_out "$key\n";
}
close($aln_out);



#Write out the names file
print "writing $names\n"; 
open(my $names_file, '>', $names) or die "$! ... Can't open the $names file\n";

while( my($key, $value) = each %$names_hash){
	print $names_file "$key\t@{$names_hash->{$key}}\n";
	#print $names_file "$key\n";
}
close($names_file);

print "Finding duplicate sequences\n"; 
#Get duplicate seqs names
my %drops; 
#while( my $key = %$names_hash ){
foreach my $key (keys %$names_hash){
	#print "$names_hash->{$key}\n";
	my $number = @{$names_hash->{$key}}; 
	if ($number == 2){
		my $keep = pop @{ $names_hash->{ $key } }; 
		$drops{$keep} =1 ;
		#push (@drops, $keep); 
	}elsif($number > 2){
		my @arr = @{ $names_hash->{ $key } };
		shift @arr;
		foreach my $elem(@arr){
			$drops{$elem} = 1;
			#push (@drops, $elem);
		} 
	}else{
		next;
	}
}
#foreach my $key (keys %drops){
#	print "$key...$drops{$key\n}";
#}
#print Dumper(%drops);


print "Cleaning the rep_set file\n";
#remove duplicates from the rep set
open(my $rep, '<', $in_rep)  or die "$! ... Can't open the $in_rep file\n";
open(my $rep_out, '>', $out_rep)  or die "$! ... Can't open the $out_rep file\n";


while( my $header = <$rep> ){
	chomp $header;
	my $seq = <$rep>; 	
	chomp $seq;
	if (exists $drops{$header}){
		next;
	}else{
		print $rep_out "$header\n$seq\n";
	}
}


close($rep);
close($rep_out);
print "Done!!!\n"; 

__END__



