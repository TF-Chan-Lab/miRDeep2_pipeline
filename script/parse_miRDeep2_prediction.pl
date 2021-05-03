#!/bin/perl
#Author: ALin
#Purpose: Parse and merge miRDeep2 prediction file and output .fasta file.
#Usage: parse_merge_miRDeep2_prediction.pl <input list> <prefix>

if(@ARGV != 2){
	print "Usage: parse_merge_miRDeep2_prediction.pl <input list> <prefix>\n";
	exit;
}

use Bio::SeqIO;
use Bio::Seq;
use Switch;
use strict;

my $in = shift;
my $prefix = shift;
my $mature = $prefix . "_mature.fa";
my $pre = $prefix . "_pre.fa";
my $pre_gff = $prefix . "_pre.gff";

open(IN, "<$in") or die "Cannot open $in!\n";
open (GFF, ">$pre_gff") or die "Cannot create $pre_gff!\n";
my $mature_out = Bio::SeqIO->new(-file => ">$mature", -format => 'fasta', );
my $pre_out = Bio::SeqIO->new(-file => ">$pre", -format => 'fasta');

my %mature_seq = ();
my %mature_pre = ();
my %mature_miRBase = ();
my %mature_coor = ();
my $mature_count = 1;

my %mature_start = ();
my %mature_end = ();
my %star_start = ();
my %star_end = ();
my %pre_seq = ();
my %pre_miRBase = ();
my $pre_count = 1;



while(<IN>){
	chomp;
	my $file = $_;
	open(FILE, "<$file") or die "Cannot open $file!\n";
	INNER: while(<FILE>){
		chomp;
		my $line = $_;
		my @line = split('\t', $line);
		
		unless($line[8] eq 'yes' || $line[8] eq 'no'){
			next INNER;
		}

		my $temp_seq = uc($line[13]);
		$temp_seq =~ s/U/T/g;
		$line[13] = $temp_seq;
		$temp_seq = uc($line[14]);
		$temp_seq =~ s/U/T/g;
		$line[14] = $temp_seq;
		$temp_seq = uc($line[15]);
		$temp_seq =~ s/U/T/g;
		$line[15] = $temp_seq;
		################################################
		#Mature sequence
		unless(exists $mature_seq{$line[13]}){
			$mature_seq{$line[13]} = 1;
		}
		unless(exists $mature_pre{$line[13]}){
			%{$mature_pre{$line[13]}} = ();
		}
		$mature_pre{$line[13]}{$line[15]} = 1;
		unless(exists $mature_miRBase{$line[13]}){
			%{$mature_miRBase{$line[13]}} = ();
		}
		$mature_miRBase{$line[13]}{$line[9]} = 1;
		unless(exists $mature_coor{$line[13]}){
			%{$mature_coor{$line[13]}} = ();
		}
		$mature_coor{$line[13]}{$line[16]} = 1;
		###############################################
		#Precursor sequence
		unless(exists $pre_seq{$line[15]}){
			$pre_seq{$line[15]} = 1;
		}
		my $pre_length = length($line[15]);
		my $mature_length = length($line[13]);
		my $star_length = length($line[14]);
		my $mature_start;
		my $mature_end;
		my $star_start;
		my $star_end;
		#Mature position
		for(my $i = 0; $i <= ($pre_length - $mature_length); $i++){
			my $temp = substr($line[15], $i, $mature_length);
			if($temp eq $line[13]){
				$mature_start = $i + 1;
				$mature_end = $i + $mature_length;
				$i = $pre_length - $mature_length
			}
		}
		if(exists $mature_start{$line[15]}){
			if($mature_start{$line[15]} > $mature_start){
				$mature_start{$line[15]} = $mature_start;
			}
		}
		else{
			$mature_start{$line[15]} = $mature_start;
		}
		if(exists $mature_end{$line[15]}){
			if($mature_end{$line[15]} < $mature_end){
				$mature_end{$line[15]} = $mature_end;
			}
		}
		else{
			$mature_end{$line[15]} = $mature_end;
		}
		#Star position
		for(my $i = 0; $i <= ($pre_length - $star_length); $i++){
			my $temp = substr($line[15], $i, $star_length);
			if($temp eq $line[14]){
				$star_start = $i + 1;
				$star_end = $i + $star_length;
				$i = $pre_length - $star_length
			}
		}
		if(exists $star_start{$line[15]}){
			if($star_start{$line[15]} > $star_start){
				$star_start{$line[15]} = $star_start;
			}
		}
		else{
			$star_start{$line[15]} = $star_start;
		}
		if(exists $star_end{$line[15]}){
			if($star_end{$line[15]} < $star_end){
				$star_end{$line[15]} = $star_end;
			}
		}
		else{
			$star_end{$line[15]} = $star_end;
		}
		unless(exists $pre_miRBase{$line[15]}){
			%{$pre_miRBase{$line[15]}} = 1;
		}
		else{
			$pre_miRBase{$line[15]}{$line[9]} = ();
		}
	}
	close FILE;
}

close IN;

foreach my $mseq (keys %mature_seq){
	my $seq_obj = Bio::Seq->new();
	my $id = $prefix . "_mature_" . $mature_count;
	my $desc = "miRBase=";
	foreach my $miRID (keys %{$mature_miRBase{$mseq}}){
		$desc .= "$miRID,";
	}
	$desc .= "; COOR=";
	foreach my $coor (keys %{$mature_coor{$mseq}}){
		$desc .= "$coor,";
	}
	$desc .= "; Precursor=";
	foreach my $pre (keys %{$mature_pre{$mseq}}){
		$desc .= "$pre,";
	}
	$seq_obj->id($id);
	$seq_obj->seq($mseq);
	$seq_obj->desc($desc);
	$mature_out->write_seq($seq_obj);
	$mature_count++;
}

my $new_line_flag = 0;

foreach my $pseq (keys %pre_seq){
	my $seq_obj = Bio::Seq->new();
	my $id = $prefix . "_precursor_" . $pre_count;
	my $length = length($pseq);
	my $desc = "miRBase=";
	foreach my $miRID (keys %{$pre_miRBase{$pseq}}){
		$desc .= "$miRID,";
	}
	$desc .= "; MATURE_START=$mature_start{$pseq}; MATURE_END=$mature_end{$pseq}; STAR_START=$star_start{$pseq}; STAR_END=$star_end{$pseq}";
	if($new_line_flag == 0){
		$new_line_flag = 1;
	}
	else{
		print GFF "\n";
	}
	print GFF "$id\tmiRDeep2\tpre-miRNA\t1\t$length\t.\t+\t.\tID=$id; $desc\n";
	if($mature_start{$pseq} < $star_start{$pseq}){
		print GFF "$id\tmiRDeep2\tmiRNA\t$mature_start{$pseq}\t$mature_end{$pseq}\t.\t+\t.\tID=$id-5p_mature; $desc\n";
		print GFF "$id\tmiRDeep2\tmiRNA\t$star_start{$pseq}\t$star_end{$pseq}\t.\t+\t.\tID=$id-3p_star; $desc";
	}
	else{
		print GFF "$id\tmiRDeep2\tmiRNA\t$star_start{$pseq}\t$star_end{$pseq}\t.\t+\t.\tID=$id-5p_star; $desc\n";
		print GFF "$id\tmiRDeep2\tmiRNA\t$mature_start{$pseq}\t$mature_end{$pseq}\t.\t+\t.\tID=$id-3p_mature; $desc";
	}
	$seq_obj->id($id);
	$seq_obj->seq($pseq);
	$seq_obj->desc($desc);
	$pre_out->write_seq($seq_obj);
	$pre_count++;
}

close GFF;

