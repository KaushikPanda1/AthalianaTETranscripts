#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $depthFile = shift@ARGV; #the file containing samtools depth for each position of the genome
my $gffFile = shift@ARGV; #the file containing transcript gff file either of genes, or TEs or anything in the format mRNA, then exons listed.
my $factor = shift@ARGV; #How many nucleotides before the splice site do you want to include for calculation? farther you go more strict you are in calculation.

#Define variables here
my %depth; #Hash that stores depth for all position in an array whose key is the chromosome name
my %spliceFactor; #store the final printable data splicefactor
my $exonCounter = 0; #to constantly average splicing factor
my $lastAnn = "mRNA"; #memory of what was the last annotation mRNA or exon;
my ($splice5,$splice5Tmp,$splice3,$splice3Tmp,$splice,$transcript,$transcriptTmp,$denomTmp,$denom,$strand)=(0,0,0,0,0,"","",0,0,"");
my ($correctionFactor,$numerator)=(0,0);
#my $factor = 3; #How many nucleotides before the splice site do you want to include for calculation? farther you go more strict you are in calculation.
#First read the depth file and store all that information in a hash
open(my $fh, '<', $depthFile) or die "Cannot open the depth file:$!\n";
while (my $line=<$fh>){

	chomp($line);
	$line=~/(\S+)\s(\d+)\s(\d+)/;
	#print "$1\t$2\t$3\n";
	push @{$depth{$1}}, $3;
}
close $fh;
#print "$exonCounter\t$splice5\t$splice5Tmp\t$splice3\t$splice3Tmp\n";
#Next read the transcript gff file and calculate 5' and 3' splice factor
open($fh, '<', $gffFile) or die "Cannot open the gffFile:$!\n";
while (my $line=<$fh>){

	chomp($line);
	$line=~/(\S+)\s\S+\s(\S+)\s(\d+)\s(\d+)\s\S+\s(\S+).*(transcript_id-"\S+").*/;
	#print "$1\t$2\t$3\t$4\t$5\t$6\n";
	#print "${$depth{$1}}[$3-1]\t${$depth{$1}}[$3]\t${$depth{$1}}[$4]\t${$depth{$1}}[$4+1]\n";
	#print "$transcript\t$splice5\t$splice3\n";
	if($2 eq "mRNA"){

		#print "$exonCounter\t$splice5\t$splice5Tmp\t$splice3\t$splice3Tmp\n";
		if($exonCounter > 1){

			${$spliceFactor{$transcript}}{'splice5'} = $splice5;
			${$spliceFactor{$transcript}}{'splice3'} = $splice3;
			${$spliceFactor{$transcript}}{'splice'} = ($splice5 + $splice3)/2;
			${$spliceFactor{$transcript}}{'strand'} = $strand;
		}
		elsif($exonCounter == 1){

			delete $spliceFactor{$transcript};
		}
		#reset for a new transcript
		$lastAnn = "mRNA";
		$splice5=0;
		$splice5Tmp=0;
		$denomTmp=0;$correctionFactor=0;$numerator=0;
		$denom=0;
		$splice3=0;
		$splice=0;
		$exonCounter=0;
		$transcript = $6;
		$strand = $5;
		if(!exists $spliceFactor{$transcript}){
			${$spliceFactor{$transcript}}{'splice5'} = 0;
			${$spliceFactor{$transcript}}{'splice3'} = 0;
			${$spliceFactor{$transcript}}{'splice'} = 0;
			${$spliceFactor{$transcript}}{'strand'} = $strand;
		}
	}
	elsif($2 eq "exon"){

		if($lastAnn eq "mRNA"){
			
			$denomTmp = ${$depth{$1}}[$4-1-$factor];#ratio of depth at end minus depth of that which exceeded end (end+1) over total depth at end-1 position
			if($denomTmp<${$depth{$1}}[$4-1]){$denomTmp=${$depth{$1}}[$4-1];}
			$correctionFactor = ${$depth{$1}}[$4] - ${$depth{$1}}[$4+$factor];if($correctionFactor<0){$correctionFactor=0;}
			$numerator = ${$depth{$1}}[$4-1] - $correctionFactor;if($numerator<0){$numerator=0;}
			if($denomTmp<=0){
				if($5 eq "+"){$splice5Tmp=0;}
				elsif($5 eq "-"){$splice3Tmp=0;}
			}
			else{
				if($5 eq "+"){$splice5Tmp = $numerator/$denomTmp;}
				elsif($5 eq "-"){$splice3Tmp = $numerator/$denomTmp;}
			}
			$exonCounter++;
			$lastAnn = "exon";
		}
		elsif($lastAnn eq "exon"){

			$splice5 = $splice5*($exonCounter - 1);
			$splice3 = $splice3*($exonCounter - 1);
			$exonCounter++;
			$denomTmp = ${$depth{$1}}[$4-1-$factor];if($denomTmp<${$depth{$1}}[$4-1]){$denomTmp=${$depth{$1}}[$4-1];}
			$denom=${$depth{$1}}[$3-1+$factor];if($denom<${$depth{$1}}[$3-1]){$denom=${$depth{$1}}[$3-1];}
			#$denom = ($denomTmp, ${$depth{$1}}[$3])[$denomTmp > ${$depth{$1}}[$3]];
			#Check what should be the denominator; $set denom
			if ($5 eq "+"){
				$splice5 = $splice5Tmp + $splice5;
				if($denomTmp<=0){$splice5Tmp=0;}
				else{
					$correctionFactor = ${$depth{$1}}[$4] - ${$depth{$1}}[$4+$factor];if($correctionFactor<0){$correctionFactor=0;}
					$numerator = ${$depth{$1}}[$4-1] - $correctionFactor;if($numerator<0){$numerator=0;}
					$splice5Tmp = $numerator/$denomTmp;
				}
				if($denom<=0){$splice3=0;}
				else{
					$correctionFactor = ${$depth{$1}}[$3-2] - ${$depth{$1}}[$3-2-$factor];if($correctionFactor<0){$correctionFactor=0;}
					$numerator = ${$depth{$1}}[$3-1] - $correctionFactor;if($numerator<0){$numerator=0;}
					$splice3 = ($numerator/$denom) + $splice3;
				}
				$splice5 = $splice5/($exonCounter - 1);
				$splice3 = $splice3/($exonCounter - 1);
			}
			elsif($5 eq "-"){
				$splice3 = $splice3Tmp + $splice3;
				if($denomTmp<=0){$splice3Tmp=0;}
				else{
					$correctionFactor = ${$depth{$1}}[$4] - ${$depth{$1}}[$4+$factor];if($correctionFactor<0){$correctionFactor=0;}
					$numerator = ${$depth{$1}}[$4-1] - $correctionFactor;if($numerator<0){$numerator=0;}
					$splice3Tmp = $numerator/$denomTmp;
				}
				if($denom<=0){$splice5=0;}
				else{
					$correctionFactor = ${$depth{$1}}[$3-2] - ${$depth{$1}}[$3-2-$factor];if($correctionFactor<0){$correctionFactor=0;}
					$numerator = ${$depth{$1}}[$3-1] - $correctionFactor;if($numerator<0){$numerator=0;}
					$splice5 = ($numerator/$denom) + $splice5;
				}
				$splice5 = $splice5/($exonCounter - 1);
				$splice3 = $splice3/($exonCounter - 1);
			}
		}
	}
}
if($exonCounter > 1){

	${$spliceFactor{$transcript}}{'splice5'} = $splice5;
	${$spliceFactor{$transcript}}{'splice3'} = $splice3;
	${$spliceFactor{$transcript}}{'splice'} = ($splice5 + $splice3)/2;
	${$spliceFactor{$transcript}}{'strand'} = $strand;
}
elsif($exonCounter == 1){

	delete $spliceFactor{$transcript};
}
close $fh;

#print Dumper(%spliceFactor);
foreach my $key (sort keys %spliceFactor){

	print "$key\t${$spliceFactor{$key}}{'splice5'}\t${$spliceFactor{$key}}{'splice3'}\t${$spliceFactor{$key}}{'splice'}\n";
}

#print ${$depth{"Chr1"}}[3630],"\n";
#print "${$depth{"Chr1"}}[3628]\t${$depth{"Chr1"}}[3629]\t${$depth{"Chr1"}}[3630]\n";
#print Dumper(\%depth);
