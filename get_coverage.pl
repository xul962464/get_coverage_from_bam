#!/usr/bin/env perl
use strict;										
use warnings;			
use Getopt::Long;		

#######################################################################################

my $BEGIN_TIME=time();	   
my $version="0.0.1";

my $sam_file;
my $ref_seq_file;
my $data_type;
my $outdir;
my $prefix;
my $MIN_Q;
# ------------------------------------------------------------------
#GetOptions
# ------------------------------------------------------------------

GetOptions(		
				"help|?" =>\&USAGE,			
				"s|sam:s"=>\$sam_file,
				"r|ref:s"=>\$ref_seq_file,
				"p|prefix:s"=>\$prefix,
				"t|type:s"=>\$data_type,
				"q|qual:s"=>\$MIN_Q,
				"o:s"=>\$outdir,
) or &USAGE;

# ------------------------------------------------------------------

if(!$sam_file and !$ref_seq_file and !$prefix and !$data_type){
	&USAGE;
}elsif(!$sam_file){
	info("-s/--sam sam/bam file is empty!");
}elsif(!$ref_seq_file){
	info("-r/--ref reference seqence file is empty!");
}elsif(!$prefix){
	info("-p/--prefix prefix is empty!");
}elsif(!$data_type){
	info("-t/--type type is empty!");
}

$data_type = "\U$data_type";

if($data_type ne "RNA" and $data_type ne "DNA"){
	info("type is wrong, must DNA or RNA");
}

# ------------------------------------------------------------------

$outdir //= "01_run";
$MIN_Q //= 40;

#######################################################################################
# ------------------------------------------------------------------
# Main Body                     
# ------------------------------------------------------------------
msg("Begin");

mkdir $outdir unless(-d "$outdir");

#=====================
msg("get reference");

open SEQ,"$ref_seq_file" or die "cannot open reference seqence file!";

$/ = ">";<SEQ>;
my %for_ref_seq;

while(<SEQ>){
	chomp;

	my ($name,$seq) = split(/\n/,$_,2);
	$name = (split/\s+/,$name)[0];
	$seq =~ s/\n|\s//g;
	$for_ref_seq{$name} = $seq;

	unless($for_ref_seq{$name}){
		error("reference $name is not exists seqence.exit.\n");
		exit;
	}
}
$/ = "\n";

#=====================
#sam file is sorted?
open HEAD,"samtools view -H $sam_file |" or die "sam file cannot open!";

my $has_header = 0;

while(<HEAD>){
	if(/\@HD.*SO:(\S+)/){
		$has_header = 1;
		
		if($1 eq "unsorted"){
			error("your sam/bam file is not sorted!");
			exit;
		}
	}
}

close HEAD;

if($has_header){
	1;
}else{
	error("your sam/bam file has no header");
	exit;
}

#=====================
msg("read sam file");

open IN, "samtools view -F 4 $sam_file | " or die "sam file cannot open!";

my %map_depth;

while(<IN>){

	my ($gene_name,$map_pos,$ciga,$reads,$reads_quality,$map_q) = (split/\t/,$_)[2,3,5,9,10,4];

#	my $reads_len = length $reads;
#
#	if($reads_len < $min_reads_length){
#		next;
#	}
	
	if(!exists $for_ref_seq{$gene_name}){
		error("$gene_name is not in reference seqence file!");
		exit;
	}

	next if($ciga =~ /P|=|X/);
	
	if($data_type eq "RNA"){
		next if($map_q < $MIN_Q);	#min Q
		next if($ciga =~ /^(\d+)S/ and $1 > 4 and $map_pos > 3);
		next if($ciga =~ /\D(\d+?)S$/ and $1 > 4 and $map_pos < (length($for_ref_seq{$gene_name}) - 3));
		next if($ciga =~ /^(\d+)S/ and $1 > 4 and $map_pos > 3);
		next if($ciga =~ /^(\d+)M/ and $1 < 17);
		next if($ciga =~ /(\d+)M$/ and $1 < 17);
	}

	my @arry_reads = split//,$reads;
	my @arry_ciga;

	for my $c($ciga =~ /(\d+[MDISH])/g){
		
		my ($nu,$type) = $c =~ /(\d+)(.)/;
#		print "$nu $type\n";
		next if($type eq "H");

		push @arry_ciga,$type for(1..$nu);
	}

	my $i = 0;
	my $now_map_pos = $map_pos - 1;

	$i++ while($map_depth{$gene_name}->[$map_pos][$i]);

	my $now_read_pos = 0;

	for my $c(0..@arry_reads-1){

		if($arry_ciga[$c] =~ /S|I/){
			$now_map_pos--;
		}elsif($arry_ciga[$c] eq "D"){
			$now_read_pos--;
			$map_depth{$gene_name}->[$now_map_pos][$i] = "-";
		}elsif($arry_ciga[$c] eq "M"){
			$map_depth{$gene_name}->[$now_map_pos][$i] = ($now_read_pos == 0 or $now_map_pos == $map_pos - 1) ? "\l$arry_reads[$now_read_pos]" : $arry_reads[$now_read_pos];
		}else{
			die;
		}

		$now_map_pos++;
		$now_read_pos++;
	}
}
close IN;

#=====================
msg("print mapped result");

open MAP,"> $outdir/$prefix\_$data_type\_mapped_reads.txt" or die"$!";

for my $gene_name(sort keys %map_depth){

	#output reference
	print MAP ">$gene_name:\n$for_ref_seq{$gene_name}\n";

	my @ref_seq = split(//,$for_ref_seq{$gene_name});

	#get max depth
	my $max_depth = 0;

	for my $i(0..@{$map_depth{$gene_name}}-1){
		next if(! defined($map_depth{$gene_name}->[$i]));
		$max_depth = @{$map_depth{$gene_name}->[$i]} > $max_depth ? @{$map_depth{$gene_name}->[$i]} : $max_depth;
	}

	#print
	for my $y(0..$max_depth-1){
		for my $x(0..@{$map_depth{$gene_name}}-1){
			if($map_depth{$gene_name}->[$x][$y]){
				$ref_seq[$x] eq "\U$map_depth{$gene_name}->[$x][$y]" ? print MAP "$map_depth{$gene_name}->[$x][$y]" : print MAP "$map_depth{$gene_name}->[$x][$y]";
			}else{
				print MAP " ";
			}
		}
		print MAP "\n";
	}
	print MAP "-" x length($for_ref_seq{$gene_name}),"\n";
}
close MAP;

#=====================
msg("print consensus result");

open CON,"> $outdir/$prefix\_$data_type\_mapped_consensus.txt";

for my $gene_name(sort keys %map_depth){
	
	print CON ">$gene_name:\n$for_ref_seq{$gene_name}\n";

	my @ref_seq = split(//,$for_ref_seq{$gene_name});
	my $max_depth = 0;

	for my $i(0..@{$map_depth{$gene_name}}-1){
		next if(! defined($map_depth{$gene_name}->[$i]));
		$max_depth = @{$map_depth{$gene_name}->[$i]} > $max_depth ? @{$map_depth{$gene_name}->[$i]} : $max_depth;
	}

	for my $y(0..$max_depth-1){
		for my $x(0..@{$map_depth{$gene_name}}-1){
			if($map_depth{$gene_name}->[$x][$y]){
				$ref_seq[$x] eq "\U$map_depth{$gene_name}->[$x][$y]" ? print CON "*" : print CON "$map_depth{$gene_name}->[$x][$y]";
			}else{
				print CON " ";
			}
		}
		print CON "\n";
	}
}
close CON;

#=====================
msg("count the base coverage of each site");

my %for_depth;

for my $gene_name(sort keys %map_depth){
	for my $i(0..@{$map_depth{$gene_name}}-1){

		next if(! defined($map_depth{$gene_name}->[$i]));

		for my $j(0..@{$map_depth{$gene_name}->[$i]}){
			if($map_depth{$gene_name}->[$i][$j]){
				$for_depth{$gene_name}{$i+1}->{"\U$map_depth{$gene_name}->[$i][$j]"}++;
			}
		}
	}
}

#=====================
msg("print the base coverage of each site");

open OUT,"> $outdir/$prefix\_$data_type\_coverage.txt";
print OUT "Gene\tData\tPos\tRef\tA\tC\tG\tT\n";

for my $gene_name(sort keys %map_depth){

	my @ref_seq = split(//,$for_ref_seq{$gene_name});

	for my $pos(sort {$a <=> $b} keys %{$for_depth{$gene_name}}){
		
		my $ref_index = $pos - 1;

		print OUT "$gene_name\t$data_type\t$pos\t$ref_seq[$ref_index]";
		
		for my $base(qw/A C G T/){
			$for_depth{$gene_name}{$pos}->{$base} //= 0;
		
			print OUT "\t$for_depth{$gene_name}{$pos}->{$base}";
		}
		print OUT "\n";
	}
}

#=====================
msg("done");
#=====================
#sub

sub error{
	print "\n";
	my $info = shift;
	print "#" x 50 , "\n";
	print "#wrong: $info \n";
	print "#" x 50 , "\n";
}

sub info{
	my $info = shift;
	print "\n";
	print "=" x 50 , "\n";
	print "=warn: $info \n";
	print "=" x 50 , "\n";
	&USAGE;
}

sub msg{
	my $info = shift;
	print GetTime()." $info...\n";
}

sub GetTime {		
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("[%4d-%02d-%02d %02d:%02d:%02d]", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {         
	my $usage=<<"USAGE";
	
Program: $0
Version: $version
Contact: xul<xul\@genepioneer.com> <1275875706\@qq.com>
Description:
	
	get reference seqence depth for every site and every base

Usage:
		perl $0 \\
				-s test.bam \\
				-r reference.fa \\
				-p test \\
				-t DNA \\
				-o 01_run
  Options:

	-s --sam	<string>	
		sam/bam file [must be sorted]

	-r --ref	<string>
		reference seqence file [fasta format]

	-p --prefix	<string>
		prefix
	
	-t --type	<string>
		input data type: DNA or RNA
	
	-q --qual	<int>
		Minimum mapping quality
		default: 40;

	-o --outdir	[string]
		default 01_run:
			prefix_datatype_coverage.txt
			prefix_datatype_mapped_consensus.txt
			prefix_datatype_mapped_reads.txt

	-h	help

USAGE
	print $usage;
	exit;
}
