#!/bin/perl
#Author: ALin
#Purpose: To parsed a parsed blast result and get the entries of identical query sequence. Criteria for a known miRNA are as following: 1. plus strand; 2. hit cover >= 0.9; 3. No indel; 4. the first 2 - 7 nt identical; 5. number of mismatch <= 2.
#Usage: perl parse_parsed_blast_known.pl <in parsed blast file> <output>

use strict;

if(@ARGV != 2){
        print "#Usage: perl parse_parsed_blast_known.pl <in parsed blast file> <output>\n";
        exit;
}

my $in = shift @ARGV;
my $out = shift @ARGV;
my $new_line_flag = 0;
my %ids = ();

my $hit_cov_cutoff = 0.9;
my $max_num_mismatch = 2;

open(IN, "<", $in) or die"Cannot open $in!\n";
open(OUT, ">", $out) or die"Cannot create $out!\n";

while(<IN>){
	chomp;
	my $strand_flag = 0;
	my $hit_cov_flag = 0;
	my $indel_flag = 0;
	my $seed_flag = 0;
	my $mismatch_flag = 0;

	my $line = $_;
	my @line = split('\t', $line);
	
	
	#Criterion 1: strand
	if($line[15] == 1){
		$strand_flag = 1;
	}
	else{
		next;
	}

	#Criterion 2: hit cover

	if($line[19] >= $hit_cov_cutoff){
		$hit_cov_flag = 1;
	}
	else{
		next;
	}
	
	#Criterion 3: indel
	if(($line[22] =~ /-/) || ($line[23] =~ /-/)){
		next;
	}
	else{
		$indel_flag = 1;
	}
	
	#Criterion 4: seed match
	if($line[13] <= 2){
		my $seed_start = 2 - $line[13];
		my $query_seed = substr($line[22], $seed_start, 6);
		my $hit_seed = substr($line[23], $seed_start, 6);
		
		if($query_seed eq $hit_seed){
			$seed_flag = 1;
		}
		else{
			next;
		}
	}
	else{
		next;
	}
	
	#Criterion 5: mismatch
	my $length_aligned = $line[14] - $line[13] + 1;
	my $num_mismatch = $line[2] - $length_aligned;
	for(my $i = 0; $i < $length_aligned; $i++){
		my $temp_query_nt = substr($line[22], $i, 1);
		my $temp_hit_nt = substr($line[23], $i, 1);
		if($temp_query_nt ne $temp_hit_nt){
		$num_mismatch++;
		}
	}
	if($num_mismatch <= $max_num_mismatch){
		$mismatch_flag = 1;
	}
	else{
		next;
	}

	if(($strand_flag == 1) && ($hit_cov_flag  == 1) && ($indel_flag == 1) && ($seed_flag == 1) && ($mismatch_flag == 1)){
		if($new_line_flag == 0){
			$new_line_flag = 1;
		}
		else{
			print OUT "\n";
		}
		print OUT $line;
	}
}

#close IN;

#foreach my $id (keys %ids){
#	if($new_line_flag == 0){
#		$new_line_flag = 1;
#	}
#	else{
#		print OUT "\n";
#	}
#	print OUT $id;
#}


close OUT;








