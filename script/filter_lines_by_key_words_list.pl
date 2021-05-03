#!/bin/perl
#Author: ALin
#Purpose: To filter lines in any files by looking at a list of words in specific column (0-based).
#Usage: filter_lines_by_key_words_list.pl <in file> <output> <list of key words> <column A> [column B]...

use strict;

if(@ARGV < 4){
	print "Usage: perl filter_lines_by_key_words_list.pl <in file> <output> <list of key words> <column A (0-based)> [column B]...\n";
	exit;
}

my $in = shift @ARGV;
my $out = shift @ARGV;
my $list = shift @ARGV;
my %key_words = ();
my @columns = ();

while(@ARGV > 0){
	my $temp_column = shift @ARGV;
	print "Column $temp_column (+1) will be checked...\n";
        push(@columns, $temp_column);
}

open(IN, "<", $in) or die"Cannot open $in!\n";
open(LIST, "<", $list) or die"Cannot open $list!\n";
open(OUT, ">", $out) or die"Cannot create $out!\n";

while(<LIST>){
	chomp;
	my $key_word = $_;
	$key_word =~ s/\r//;
	$key_words{$key_word} = 1;
}
close LIST;

my $num = (keys %key_words);
my $num_c = @columns;
print "Total entries in the list: $num\n$num_c column(s)\n";

my $new_line_flag = 0;

LINE:while(<IN>){
	chomp;
	my $line = $_;
	$line =~ s/\r//;
	my @line = split('\t', $line);
	foreach my $i (@columns){
		if(exists $key_words{$line[$i]}){
			#print "$key_words{$line[$i]}\n";
			next LINE;
		}
	}
	if($new_line_flag == 0){
		$new_line_flag = 1;
	}
	else{
		print OUT "\n";
	}
	print OUT "$line";
}
close IN;
close OUT;


















