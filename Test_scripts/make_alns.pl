#!/usr/bin/perl

################
#              #
# make_alns.pl #
#              #
################

use strict; use warnings;

#This is a wrapper that will make infernal alignments for each fasta file in the input directory

#usage: make_alns.pl <full_path_to_input_directory> <full_path_to_profile> <full_path_to_ref_align> extension

#ARGS 
#<full_path_to_input_directory>: Full path to directory containing the input fasta files
#<full_path_to_profile>        : Full path to profile alignment (usually .cm file) for input in to cmalign
#<full_path_to_ref_align>      : Full path to alignment of reference seqs (those used in the tree that will be used downstream with pplacer)
#extension                     : File extension used for fasta files (i.e., fna, fa, fasta, etc. ) only one extension supported do not include '.'

die "wrong number of arguments" if(@ARGV != 4);  
 
my $fnadir  = $ARGV[0];
my $profile = $ARGV[1];
my $ref		= $ARGV[2];
my $ext		= $ARGV[3];
opendir(DIR, $fnadir) or die $!;

my @files; 

while ( my $file = readdir(DIR) ){
	
	if ( $file =~ m/\.$ext/ ) {
		chomp $file;
		push ( @files, $file );
	}else{ 
		next;	
	}
}

close(DIR);
my $dir_name = "Query_alns"; 
my $ndir = system("mkdir $dir_name");

print "Making query alignments\n"; 

foreach my $entry( @files ){
	chomp $entry;
	my @fname = split('\.', $entry);
	my $fname = $fname[0]; 
	#print "@fname\n\n\n\n"; 

	my $aln   = system("cmalign --hbanded --sub --dnaout --outformat pfam -g --notrunc -o $dir_name/$fname\_aln.sto $profile $entry 2> $dir_name/$fname.log"); 			
}

my $merge_dir_name = "Merged_alns"; 
my $merge_aln_dir  = system("mkdir $merge_dir_name");


opendir(NDIR, $dir_name) or die $!;

my @alns; 
while ( my $aln_file = readdir(NDIR) ){
	
	if ( $aln_file =~ m/\.sto/ ) {
		chomp $aln_file;
		push ( @alns, $aln_file );
	}else{ 
		next;
	}
}


print "Making merged alignments\n"; 

foreach my $aln_f ( @alns ){
	chomp $aln_f;
	my @fname = split('\.', $aln_f);
	my $fname = $fname[0]; 
	my $scall = system("esl-alimerge --dna -o $merge_dir_name/merged_$fname\_ref.sto $ref $dir_name/$aln_f 2> $merge_dir_name/$fname.log"); 			
}

close(NDIR);

__END__


