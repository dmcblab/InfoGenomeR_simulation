#!/usr/bin/perl -w
use Bio::DB::Fasta;
#use strict;
#use warnings;
#######/DAS_Storage5/leeyh/MELTv2.1.3/me_refs/1KGP_Hg19/########
my $alu_seq = uc Bio::DB::Fasta->new('ALU.fa')->seq("ALU");
my $L1_seq = uc Bio::DB::Fasta->new('LINE1.fa')->seq("LINE1");
my $sva_seq = uc Bio::DB::Fasta->new('SVA.fa')->seq("SVA");
my $numt_seq = uc Bio::DB::Fasta->new('hs37d5_MT.fa')->seq("MT");
###############
my @args;
my $search_length=10000;
my $db = Bio::DB::Fasta->new('/DATA1/AITL/hg19/ucsc.hg19.fasta');


        $hap=$ARGV[0];
	$j=$ARGV[1];
	$IND=$ARGV[2];



	if($j != 23){
		$chr=$j;
	}else{
		$chr="X";
	}

	my $filename = "./${IND}.g${hap}.fa.${chr}";
	open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
	my $seq = $db->seq("chr$chr");
	open(my $index1, '>', "./${IND}.g${hap}.fa.${chr}.delta") or die  "Could not open file";
        open(my $true_SVs, '>', "./${IND}.g${hap}.fa.${chr}.SVs") or die  "Could not open file";

	print $index1 "new\tnew_delta\n";
	
	$i=0;
	open($vcf, '<', "${IND}.g${hap}.vcf.${chr}");
	while(my $line = <$vcf>){
		chomp $line;
		next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
		my @line = split(/\t+/, $line);
		if($line[7] eq "m"){
		        print $index1 "${line[1]+$i}\t";
			print $index1 -$i;
			print $index1 "\n";
			substr($seq, $line[1]+$i-1, length($line[2]))=$line[5];
			$i=$i+length($line[5])-length($line[2]);
		}else{
			if($line[7] eq "alu" || $line[7] eq "sva" || $line[7] eq "L1"){
				my @me_info = split(/\,/, $line[5]);
				my $me_seq = "";
				if($line[7] eq "alu"){$me_seq = $alu_seq;}elsif($line[7] eq "sva"){$me_seq=$sva_seq;}else{$me_seq=$L1_seq;}
 	                        print $index1 "${line[1]+$i}\t";
        	                print $index1 -$i;
	                        print $index1 "\n";
				my $me_seq_sub="";
				if($me_info[3] eq "+"){
					$me_seq_sub=substr($me_seq,$me_info[1]-1,$me_info[2]-$me_info[1]+1);
				}else{
                                        $me_seq_sub=substr($me_seq,$me_info[1]-1,$me_info[2]-$me_info[1]+1);
 		                        $me_seq_sub =~ tr/ACGTacgt/TGCAtgca/;
                		        $me_seq_sub = scalar reverse $me_seq_sub;
				}
	                        substr($seq, $line[1]+$i-1,1)=$line[2].$me_seq_sub;
	                        $i=$i+length($me_seq_sub);
			}elsif($line[7] eq "numt"){
                                my @me_info = split(/\,/, $line[5]);
                                my $me_seq = $numt_seq;
                                print $index1 "${line[1]+$i}\t";
                                print $index1 -$i;
                                print $index1 "\n";
                                my $me_seq_sub="";
                                if($me_info[3] eq "+"){
                                        $me_seq_sub=substr($me_seq,$me_info[1]-1,$me_info[2]-$me_info[1]+1);
                                }else{
                                        $me_seq_sub=substr($me_seq,$me_info[1]-1,$me_info[2]-$me_info[1]+1);
                                        $me_seq_sub =~ tr/ACGTacgt/TGCAtgca/;
                                        $me_seq_sub = scalar reverse $me_seq_sub;
                                }
                                substr($seq, $line[1]+$i-1, $line[3]-$line[1])=$line[2].$me_seq_sub;
                                $i=$i+1+length($me_seq_sub)-$line[3]+$line[1];
			}elsif($line[7] eq "dup" || $line[7] eq "del"){
				if($line[7] eq "del"){
					print $true_SVs "3to5\t${hap}\t$line[0]\t$line[1]\t${hap}\t$line[0]\t$line[3]\t$line[6]\t", "${line[1]+$i}\n";
				}elsif($line[7] eq "dup"){
                                        print $true_SVs "5to3\t${hap}\t$line[0]\t$line[1]\t${hap}\t$line[0]\t$line[3]\t$line[6]\t", "${line[1]+$i}\n";
				}
                                print $index1 "${line[1]+$i}\t";
                                print $index1 -$i;
                                print $index1 "\n";
				my $cn_seq="";
				my $dup_seq="";
                                if($line[7] eq "dup"){
					$dup_seq=substr($seq, $line[1]+$i-1, $line[3]-$line[1]+1);
					for(my $dup_i = 0; $dup_i <= $line[6];$dup_i++){
						$cn_seq = $cn_seq.$dup_seq;
					}
				}else{
					$cn_seq = substr($seq, $line[1]+$i-1,1).substr($seq, $line[3]+$i-1,1);
				}
				substr($seq, $line[1]+$i-1, $line[3]-$line[1]+1)=$cn_seq;
				$i=$i+length($cn_seq)-($line[3]-$line[1]+1);
			}elsif($line[7] eq "inv"){
				my $inv_seq_sub="";
				print $true_SVs "3to3\t${hap}\t$line[0]\t$line[1]\t${hap}\t$line[0]\t$line[3]\t$line[6]\t", "${line[1]+$i}\n";
                                print $true_SVs "5to5\t${hap}\t$line[0]\t$line[1]\t${hap}\t$line[0]\t$line[3]\t$line[6]\t", "${line[1]+$i}\n";
                                print $index1 "${line[1]+$i}\t";
                                print $index1 -$i;
                                print $index1 "\n";
                                $inv_seq_sub=substr($seq,$line[1]+$i, $line[3]-$line[1]);
                                $inv_seq_sub =~ tr/ACGTacgt/TGCAtgca/;
                                $inv_seq_sub = scalar reverse $inv_seq_sub;
				substr($seq, $line[1]+$i, $line[3]-$line[1])=$inv_seq_sub;
			}
		}
	}
	close($vcf);
	##}
	
	
	print $fh ">${hap}\_${chr}\n";
	while( my $chunk = substr($seq, 0, 80, "")){
	          print $fh "$chunk\n";
	}



