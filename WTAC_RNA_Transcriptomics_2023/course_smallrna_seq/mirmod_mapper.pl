#!/usr/bin/env perl
$|=1;

use strict; 

my $upload_dir = "."; 
my $species = shift(@ARGV);
my $outfile = shift(@ARGV);
my @filenames=();
my $edits=0;
my $depth_threshold=0;

foreach my $a (@ARGV){
	if ($a =~ /--edits/){
		$edits=1;
	} else {
		if ($a =~ /--depth=(\d+)/){
			$depth_threshold=$1;
		} else {
			push(@filenames,$a);
		}
	}
}


if (!@ARGV){
	print "\nMirMod Mapper v1.0 (Enright Lab - EMBL EBI (c) 2012\n";
	print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
	print "Usage: mirmod_mapper.pl hairpins.fasta outputfile.counts (files to process)\n\n";
	exit(1);
} else {
	if (!$species){
		die "Error: please specify an output file for the analysis";
	}
	if (!-e $species){
		die "Error: hairpin fasta file $ARGV[0] not found or not readable";
	}
	if (!-e "$species.folds"){
		die "Error: $species.folds file containing RNAfold folds for each hairpin not present or readable";
	}
}

if (!@filenames){
	die "Error: must have specified at least one file to process";
} else {
	print "Will process the following files: @filenames\n";
}


open (FILEOUT4,">$outfile");
open (FILEOUT5,">mappings.txt");

print "\n";
print "Will store final tabular counts in: $outfile\n";
print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
my $i=0;



my %mirs=();
my %reads=();
my %depths=();
my %rawdepth=();
my %maxdepth=();
my %totalreads=();
my %allreads=();
my %alldepths=();
my %threepreads=();
my %fivepreads=();
my %mmreads=();
my %mirseqs=();
my %totaldepths=();

print "\nProcessing files:\n\n";

print "Reading Hairpin lengths\n";

my $id;
my %hairpinseq=();
my %foldseq=();
my %hairpinlength=();
open(FILE,$species);
while(<FILE>){
	chomp;
	if (/^>(\S+)/){
		$id=$1;
	} else {
		$hairpinseq{$id}.=$_;
		$hairpinlength{$id}=length($hairpinseq{$id});
	}
}

open(FILE,"$species.folds");
while(<FILE>){
        chomp;
        if (/^>(\S+)/){
                $id=$1;
        } else {
                $foldseq{$id}.=$_;
        }
}


foreach my $filename(@filenames){
print "Species: $species Filename:$filename\n";

my $aln_file="$outfile\_$filename" . ".alignments.html";
my $mat_file="$outfile\_$filename" . ".mat";

print "Storing This Samples Alignments in: $aln_file\n";
print "Storing This Samples Matrices in: $mat_file\n";

if ($edits){
	print "Will consider 5', 3' modifications and edits (mode enabled)\n";
} else {
	print "Not considering 5', 3' modifications and edits (mode disabled, use --edits to enable)\n";
}

if ($depth_threshold > 0){
	print "Will filter reads with depth less than $depth_threshold\n";
}

open (FILEOUT1,">$aln_file");
open (FILEOUT2,">$mat_file");
print FILEOUT2 "\tA\tC\tG\tU\n";

my %alignments=();
my %avg_len;
my %avg_len_count;

my %four_way_len;
my %four_way_len_count;

undef(%rawdepth);

my $depth=0;

%mirseqs=();
my $mid="";

if ($edits){
print "Storing Sequences\n";
open (FILE,"gunzip -c $upload_dir/$filename |");
while(<FILE>){
chomp;
if(/^>(\S+)/){
        $mid=$1;
} else {
        $mirseqs{$mid}.=$_;
	$mirseqs{$mid}=~ s/T/U/g;
}

}
print "done\n";
}

my $cat_function="cat";
print "Processing $upload_dir/$filename\n";
print "gunzip -dcf $upload_dir/$filename | ./clean.pl |blastall -p blastn -e 0.001 -d $species -m 8 -a 8 -v 1 -b 1 2>\n";
open(PROC,"gunzip -dcf $upload_dir/$filename | ./clean.pl |blastall -p blastn -e 0.001 -d $species -m 8 -a 8 -v 1 -b 1 2> $filename.blast.stderr | ");

while(<PROC>){
		chomp();
		my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev)=split("\t",$_);
		my $fail=0;
		my $sample_reads++;

                print FILEOUT5 "$filename\t$id\t$match\n";


                my $nice_id = $id;
                $nice_id =~ s/_x(\d+)/ Depth:$1 Modification:/g;
		my $this_depth=$1;

		$allreads{$filename}++;
                $alldepths{$filename}+=$this_depth;

                my $mirlength=length($mirseqs{$id});

		# Remove reverse complement hits (sometimes hit the passenger)
		if ($se < $ss){
			$fail=1;
		}

		# Remove 'edge effect' hits where the miR maps off the boundary of a precursor
		if (($ss-$qs) < 0){
			$fail=1;
		}
		if ($edits){
		if (($ss+($mirlength - $qs)) > $hairpinlength{$match}) {
			$fail=1;
		}

		#Exclude low depth hits
		if ($this_depth < $depth_threshold){
			$fail=1;
		}
		}


		if (!$fail){

		my $arm="";

		if ((($ss+$se)/2) <= ($hairpinlength{$match}/2)){
			$arm="5p";
		} else {
			$arm="3p";
		}

		my $aligned="";

		if ($edits){
			for ($i=0;$i<($ss-$qs);$i++){
                       	 $aligned .= " ";
                	}	

		$aligned .= $mirseqs{$id} . "\t$nice_id";
		}

		my $initial_id=$match;
		$match=$match . "-$arm";

		#print ">>$_\n";
		#print ">>$mirseqs{$id} $mirlength\n";
		my $five_nt="";
		my $three_nt="";	
		if ($edits){
		if ($mm < $mirlength){
			if ($qs != '1'){
				$five_nt=substr($mirseqs{$id},0,$qs-1);
			}
			if ($qe < $mirlength){
				$three_nt=substr($mirseqs{$id},$qe,$mirlength-$qe);
			}
		}
		}


		#print "$five_nt $three_nt\n";
		if ($edits ==0){
			$five_nt="";
			$three_nt="";
		}
		my $modification = "";
		if ($id=~/(\d+)$/){
			$depth=$1;
		} else {
			$depth=1;
		}
		$rawdepth{$match}+=$depth;
		if ($five_nt){
			$depths{$filename}{$match."_nont_5p_$five_nt"}+=$depth;
			$mirs{$match."_nont_5p_$five_nt"}=1;
			$modification = "_nont_5p_$five_nt";
			if ($maxdepth{$match."_nont_5p_$five_nt"} < $depth){
				$maxdepth{$match."_nont_5p_$five_nt"} = $depth;
			}
		}
		elsif ($three_nt){
			$depths{$filename}{$match."_nont_3p_$three_nt"}+=$depth;
			$mirs{$match."_nont_3p_$three_nt"}=1;
			$modification = "_nont_3p_$three_nt";
			if ($maxdepth{$match."_nont_3p_$three_nt"} < $depth){
			                                $maxdepth{$match."_nont_3p_$three_nt"} = $depth;
							                        }
		} else {
			$mirs{$match}=1;
		        $reads{$filename}{$match}++;
		        $depths{$filename}{$match}+=$depth;
			$modification = "no_modification";
			if ($depths{$filename}{$match} > $maxdepth{$match}){
			        $maxdepth{$match}=$depths{$filename}{$match};
				}
		}

		if ($edits){
		$alignments{$match}.=$aligned . " " . $modification . "\n";
                $hairpinseq{$match}=$hairpinseq{$initial_id};
		$foldseq{$match}=$foldseq{$initial_id};
		$avg_len{$match}+=$mirlength;
		$avg_len_count{$match}++;	

	
		if(($modification =~ /no_modification/) || ($modification =~ /_nont_3p_[U]+$/)){
			$four_way_len{$match}+=$mirlength;
                	$four_way_len_count{$match}++;
		}


		}

		$totalreads{$filename}++;
		$totaldepths{$filename}+=$depth;
	} else {
		#FAIL
	}
}


print "mirBase\tUnique\tDepth\n";
print "TOTAL:\t$totalreads{$filename}\t$totaldepths{$filename}\n";
print "INITIAL:\t$allreads{$filename}\t$alldepths{$filename}\n";
print "---------------------------------\n";

print FILEOUT1 "<HTML>\n";
print FILEOUT1 "<BODY BGCOLOR=\"#ECE0F8\">\n";
print FILEOUT1 "<FONT FACE=\"Arial\">\n";
print FILEOUT1 "<H2>MicroRNA Modification Report for $filename</H2>\n";
print FILEOUT1 "<BR>Reads with less than depth $depth_threshold excluded from this analysis\n";
print FILEOUT1 "<BR><BR><BR><BR>\n";


foreach my $th (sort special2(keys(%alignments))){
	print FILEOUT1 "<HR>\n";
	print FILEOUT1 "<H3>$th</H3>\n";
	print FILEOUT1 "Depth: <b>$rawdepth{$th}</b> reads\n";
	print FILEOUT1 "<H4>Aligned Reads on Precursor</H4>\n";
	print FILEOUT1 "<pre style=\"width: 1024px; background-color: FFFFFF; border: 1px solid #bebab0; -webkit-border-radius: 10px; -khtml-border-radius: 10px; -moz-border-radius: 10px; border-radius: 10px;\">\n";
	print FILEOUT1 formatseq("$hairpinseq{$th}\n");
	print FILEOUT1 "$foldseq{$th}\n";
	print FILEOUT1 formatseq("$alignments{$th}\n");
	print FILEOUT1 "</pre>\n";

	my $avg;
	my $four_way_avg;

	if ($avg_len_count{$th} != 0){
		$avg=$avg_len{$th}/$avg_len_count{$th};
	}
	if ($four_way_len_count{$th} != 0){
		$four_way_avg=$four_way_len{$th}/$four_way_len_count{$th};
	}

	print FILEOUT1 "<b>$th: Average Length: $avg</b>\n";
	print FILEOUT1 "<b>$th: FourWay Length: $four_way_avg</b>\n";

	my @matrix=();
	my @histo=();
	my @modif=();

	my $hlength=length($hairpinseq{$th});
	my @array=split("\n",$alignments{$th});

	for (my $j=0;$j<=$hlength;$j++){
		$histo[$j]=0;
		$modif[$j]=0;
	}

	for (my $i=0;$i<=$#array;$i++){
		my $string=$array[$i];
		$string=~ s/(\s+\S+)\s+.*/$1/g;
		for (my $j=0;$j<=$hlength;$j++){
				if (substr($string,$j,1) eq 'A'){
					$matrix[$j][1] ++;
					if (substr($hairpinseq{$th},$j,1) ne 'A'){
						$modif[$j]++;
					} else {
						$histo[$j]++;
					}
				}
				 if (substr($string,$j,1) eq 'C'){
                                        $matrix[$j][2] ++;
					if (substr($hairpinseq{$th},$j,1) ne 'C'){
                                                $modif[$j]++;
                                        } else {
                                                $histo[$j]++;
                                        }
                                }
				 if (substr($string,$j,1) eq 'G'){
                                        $matrix[$j][3] ++;
					if (substr($hairpinseq{$th},$j,1) ne 'G'){
                                                $modif[$j]++;
                                        } else {
                                                $histo[$j]++;
                                        }
                                }
				 if (substr($string,$j,1) eq 'U'){
                                        $matrix[$j][4] ++;
					if (substr($hairpinseq{$th},$j,1) ne 'U'){
                                                $modif[$j]++;
                                        } else {
                                                $histo[$j]++;
                                        }
                                }
		}
	}

	my $max_height=0;
	for (my $i=0;$i<=$#histo;$i++){
		if (($histo[$i]+$modif[$i])>$max_height){
			$max_height=($histo[$i]+$modif[$i]);
		}
	}

	print FILEOUT1 "<H4>Distribution across Precursor</H4>\n";
	print FILEOUT1 "<TABLE style=\'padding:0px; background-color:white; border:1px solid black; border-spacing:0; border-collapse:collapse;\'><TR>\n";
	for (my $i=0;$i<=$#histo;$i++){
		my $height1=($histo[$i]/$max_height)*100;
		my $height2=($modif[$i]/$max_height)*100;
		print FILEOUT1 "<TD>\n";
		print FILEOUT1 "<div style=\'position:relative; height:100px; width:6px\'>\n";
    		print FILEOUT1 "<div style=\'background-color:blue; height:" . $height1 . "px; position:absolute; bottom:0px; width:6px \' />\n";
		print FILEOUT1 "<div style=\'background-color:red; height:" . $height2 . "px;  position:absolute; bottom:". $height1 ."px; width:6px \' />\n";
		print FILEOUT1 "</div>\n";
		print FILEOUT1 "</TD>\n";
	}
	print FILEOUT1 "</TR></TABLE>\n";


	for (my $i=0;$i<=$hlength;$i++){
		print FILEOUT2 "$th" ."_"."$i:\t";
		for (my $j=1;$j<=4;$j++){
			if (!$matrix[$i][$j]){
				print FILEOUT2 "0\t";
			}
			print FILEOUT2 "$matrix[$i][$j]\t";
		}
		print FILEOUT2 "\n";
	}


}




print FILEOUT1 "</HTML>\n";
}

foreach my $fz(@filenames){
	print FILEOUT4 "\t$fz";
}
print FILEOUT4 "\n";

foreach my $mir(sort special(keys(%mirs))){
	print FILEOUT4 "$mir";
	foreach my $fz(@filenames){
		if (!$depths{$fz}{$mir}){
			$depths{$fz}{$mir}=0;
		}
		print FILEOUT4 "\t$depths{$fz}{$mir}";
	}
	print FILEOUT4 "\n";
}

print FILEOUT4 "totaldepth";
foreach my $fz(@filenames){
	print FILEOUT4 "\t$totaldepths{$fz}";
}
print FILEOUT4 "\n";

sub special {
return($maxdepth{$b} <=> $maxdepth{$a});
}

sub special2 {
return($rawdepth{$b} <=> $rawdepth{$a});
}

sub formatseq{

my $seq=$_[0];

$seq=~ s/A/<font bgcolor=\"black\" color=\"red\">A<\/font>/g;
$seq=~ s/U/<font bgcolor=\"black\" color=\"orange\">U<\/font>/g;
$seq=~ s/C/<font bgcolor=\"black\" color=\"blue\">C<\/font>/g;
$seq=~ s/G/<font bgcolor=\"black\" color=\"purple\">G<\/font>/g;
$seq=~ s/N/<font bgcolor=\"black\" color=\"grey\">N<\/font>/g;

return($seq);

}
