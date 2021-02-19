use POSIX;
use Bio::DB::Fasta;
use Switch;
open(my $map_log, '>', "map_log") or die;
my $G1000_dir="/DASstorage6/leeyh/1000G_SNP_indel_SV_integrated";
my $low_map_line=`cat /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.wig.filtered | wc -l`;

my $PLOIDY=2;
my $PLOIDY_delta_1=3;
my $PLOIDY_delta_2=2;
my $intra_trans_n_min=4;
my $intra_trans_n_max=40; ## should be even
my $intra_trans_size_MIN=1e6;
#my $intra_trans_min = 1e6;
#my $intra_trans_max_delta = 2e6;
my $offset_to_3_end=100;
my $number_of_WGD=0;
my @chr_start;
my @chr_end;
my @chr_cumul;
my @chr_hap;
my @chr_number;
my $chr_index;
my $deletion_bridge=1000;
my $IND=$ARGV[0];
my $mappability_bias=0.2;
my $max_map_score=0.5;
open(my $log, '>', 'logs') or die;
open(my $germ_results, '>', 'germ_results') or die;
my $number_of_germ=2000;
my $number_of_somatic=100;
my @chr_fa;
my $germline_chrs;
for(my $chrs=0; $chrs<=22; $chrs++){
        my $chr_name;
        if($chrs != 22){
                $chr_name=$chrs+1;
        }else{
                $chr_name="X";
        }
        $chr_fa[0][$chrs]= Bio::DB::Fasta->new("${G1000_dir}/${IND}.g1.fa.${chr_name}");
        $chr_fa[1][$chrs]= Bio::DB::Fasta->new("${G1000_dir}/${IND}.g2.fa.${chr_name}");
}


for(my $chr_index = 0; $chr_index <=22; $chr_index++){
	my $chr_name;
	if($chr_index != 22){
		$chr_name=$chr_index+1;
	}else{
		$chr_name="X";
	}	
	my $db = Bio::DB::Fasta->new("${G1000_dir}/${IND}.g1.fa.${chr_name}");
	my $seq = $db->seq("1_${chr_name}");
	my $chr_length=length($seq);
	$chr_start[$chr_index][0]=1;
	$chr_end[$chr_index][0]=$chr_length;
	$chr_hap[$chr_index][0]=1;
	$chr_number[$chr_index][0]=$chr_name;
        my $db = Bio::DB::Fasta->new("${G1000_dir}/${IND}.g2.fa.${chr_name}");
        my $seq = $db->seq("2_${chr_name}");
        my $chr_length=length($seq);
        $chr_start[23+$chr_index][0]=1;
        $chr_end[23+$chr_index][0]=$chr_length;
        $chr_hap[23+$chr_index][0]=2;
        $chr_number[23+$chr_index][0]=$chr_name;


#	$chr_start[$chr_index][1]=10001;
#	$chr_end[$chr_index][1]=20000;
#	$chr_start[$chr_index][2]=20001;
#	$chr_end[$chr_index][2]=30000;
#	$chr_start[$chr_index][3]=30001;
#	$chr_end[$chr_index][3]=40000;
}

for(my $chr_index = 0; $chr_index <= $#chr_start; $chr_index++){

	$chr_cumul[$chr_index][0]=$chr_end[$chr_index][0]-$chr_start[$chr_index][0]+1;
	for($i = 1; $i <=$#{@chr_end[$chr_index]}; $i++){
	        @{@chr_cumul[$chr_index]}[$i]=@{@chr_cumul[$chr_index]}[$i-1]+@{@chr_end[$chr_index]}[$i]-@{@chr_start[$chr_index]}[$i]+1;
	}
}


open(my $fh, '>', 'SVs') or die "Could not open file 'SVs' $!";

my %SV_hash;
%SV_hash = (0=>"dup", 1=>"del", 2=>"simple_inv");

my @line;
for(my $SVs = 0; $SVs < $number_of_germ+$number_of_somatic; $SVs++){
	@line =();
	
	if($SVs < $number_of_germ ){ ####### germline$$$
		@line[0]="initiation";
		@line[3]=0;
		@line[1]=0;
		@line[2]=inf;
		my $N_test=1;
		my $mappability=1;
		while($line[1]<=0 ||  $N_test==1 || @line[1]<1 || @line[2] > $chr_cumul[@line[3]][$#{@chr_cumul[@line[3]]}] || ($SVs <= floor($number_of_germ*$mappability_bias) && $mappability>$max_map_score) ||($SVs > floor($number_of_germ*$mappability_bias) && $mappability<$max_map_score)){
			@line = ();
                        @line[0]=SV_roulette_germ(rand(1));
			my $c=`Rscript rbeta_germline.R $line[0]`;my @c = ();@c=split(/\s+/, $c);
			my $SV_size = @c[1];
			if($SVs <= floor($number_of_germ*$mappability_bias)){
				my @low_map_random_coor;
				($low_map_random_coor[0], $low_map_random_coor[1], $low_map_random_coor[2])=low_map_select();
				my( $new_chr_index , $new_chr_coor)=old_fa_coor_to_new_fa_coor($low_map_random_coor[0], $low_map_random_coor[1], $low_map_random_coor[2], $SV_size);
				@line[3]=$new_chr_index;
				@line[1]=$new_chr_coor;
			}else{
					@line[3]=int(rand($#chr_start+1));
					@line[1]=1+int(rand($chr_cumul[@line[3]][$#{@chr_cumul[@line[3]]}]-$offset_to_3_end-$SV_size));
			}
                        @line[2]=@line[1]+$SV_size-1;
                        $N_test = N_sequence_test_1_chr();
			if($line[1]>0){
	                        $mappability=Mappability_cal();
			}
		}
		print $fh "@line[0]\t@line[1]\t@line[2]\t@line[3]\n";
	}else{

                @line[0]="initiation";
                @line[5]=0;
                @line[1]=0;
                @line[2]=inf;
                @line[6]=0;
                @line[3]=0;
                @line[4]=inf;
		my $N_test=1;
		my $mappability=1;
		my $while_iter=1;
		my $uniform_test=0;
                while($uniform_test ==0 && $line[1]<= 0 || (($line[0] eq "balanced_translocation" ||$line[0] eq "unbalanced_translocation" || $line[0] eq "intra_trans_2" || $line[0] eq "insertion") && $line[3] <= 0) || ($SVs <= $number_of_germ+floor($number_of_somatic*$mappability_bias) && $mappability>$max_map_score) || ($SVs > $number_of_germ+floor($number_of_somatic*$mappability_bias) && $mappability<$max_map_score) || $N_test==1 || 
			(( $line[0] eq "insertion" || @line[0] eq "balanced_translocation" || @line[0] eq "unbalanced_translocation" || @line[0] eq "intra_trans_2") && (@line[5] == @line[6] || @line[1]<1 || @line[2] > $chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}] || @line[3]<1 || @line[4] > $chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}])) 
			||
		       (($line[0] ne "insertion" && @line[0] ne "balanced_translocation" && @line[0] ne "unbalanced_translocation" && @line[0] ne "intra_trans_2") && (@line[1]<1 || @line[2] > $chr_cumul[@line[3]][$#{@chr_cumul[@line[3]]}]))){

			if($SVs <= $number_of_germ+floor($number_of_somatic*$mappability_bias)){$mappability=0;}else{$mappability=1}## just of dup,del,inv###

	                @line = ();
			if($while_iter==1){
	                	@line[0]=SV_roulette(rand(1));
			}else{
				@line[0]=$prev_SV_type;
			}
			$while_iter=$while_iter+1;
			$prev_SV_type=$line[0];
			switch(@line[0]){
				case "foldback_3" {
	 	                       @line[3]=int(rand($#chr_start+1));
	        	               my $SV_size = 2000+int(rand(9500));
	              		       @line[1]=1+int(rand($chr_cumul[@line[3]][$#{@chr_cumul[@line[3]]}]-$offset_to_3_end-$SV_size));
		                       @line[2]=@line[1]+$SV_size-1;
				}
				case "insertion"{
                                       @line[5]=int(rand($#chr_start+1));
                                       @line[6]=int(rand($#chr_start+1));
                                       my $c=`Rscript rbeta_somatic.R insertion`;my @c = ();@c=split(/\s+/, $c);
                                       my $SV_size = @c[1];
                                       @line[1]=1+int(rand($chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}]-$offset_to_3_end-$SV_size));
                                       @line[2]=@line[1]+$SV_size-1;
                                       $c=`Rscript rbeta_somatic.R insertion`;my @c = ();@c=split(/\s+/, $c);
                                       $SV_size = @c[1];
                                       @line[3]=1+int(rand($chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}]-$offset_to_3_end-$SV_size));
                                       @line[4]=@line[3]+$SV_size-1;
	
				}
				case "balanced_translocation" {
	                               @line[5]=int(rand($#chr_start+1));
	                               @line[6]=int(rand($#chr_start+1));
	                               my $SV_size = 1000+int(rand(10000));
	                               @line[1]=1+int(rand($chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}]-$offset_to_3_end-$SV_size));
	                               @line[2]=@line[1]+$SV_size-1;
	                               my $SV_size = 1000+int(rand(10000));
	                               @line[3]=1+int(rand($chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}]-$offset_to_3_end-$SV_size));
	                               @line[4]=@line[3]+$SV_size-1;
				}
                                case "unbalanced_translocation" {
                                       @line[5]=int(rand($#chr_start));
                                       @line[6]=int(rand($#chr_start));
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[1]=1+int(rand($chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}]-$offset_to_3_end-$SV_size));
                                       @line[2]=@line[1]+$SV_size-1;
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[3]=1+int(rand($chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}]-$offset_to_3_end-$SV_size));
                                       @line[4]=@line[3]+$SV_size-1;
                                }
				case "chromosome_deletion"{
					@line[3]=int(rand($#chr_start));
					@line[1]=$chr_start[@line[3]][0];
					@line[2]=$chr_end[@line[3]][$#{@chr_end[@line[3]]}];
				}
				case "chromosome_amplification"{
                                        @line[3]=int(rand($#chr_start+1));
                                        @line[1]=$chr_start[@line[3]][0];
                                        @line[2]=$chr_end[@line[3]][$#{@chr_end[@line[3]]}];
				
				}
				case "WGD"{
					@line[3]=0;
					@line[1]=1;
					@line[2]=1;
				}
				case "intra_trans"{
                                        @line[3]=int(rand($#chr_start+1));
                                        @line[1]=1;
                                        @line[2]=1;
				}
				case "intra_trans_2"{
                                       @line[5]=int(rand($#chr_start));
                                       @line[6]=int(rand($#chr_start));
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[1]=1+int(rand($chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}]-$offset_to_3_end-$SV_size));
                                       @line[2]=@line[1]+$SV_size-1;
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[3]=1+int(rand($chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}]-$offset_to_3_end-$SV_size));
                                       @line[4]=@line[3]+$SV_size-1;
				}
				else{
					my $c=`Rscript rbeta_somatic.R $line[0]`;my @c = ();@c=split(/\s+/, $c);
	                               my $SV_size = @c[1];
					if($SVs <= $number_of_germ+floor($number_of_somatic*$mappability_bias)){
						my @low_map_random_coor;
						($low_map_random_coor[0], $low_map_random_coor[1], $low_map_random_coor[2])=low_map_select();
						my( $new_chr_index , $new_chr_coor)=old_fa_coor_to_new_fa_coor($low_map_random_coor[0], $low_map_random_coor[1], $low_map_random_coor[2], $SV_size);
						@line[3]=$new_chr_index;
						@line[1]=$new_chr_coor;
					}else{
						@line[3]=int(rand($#chr_start+1));
						@line[1]=1+int(rand($chr_cumul[@line[3]][$#{@chr_cumul[@line[3]]}]-$offset_to_3_end-$SV_size));
					}
					@line[2]=@line[1]+$SV_size-1;
					if($line[1]>0){
	 		                        $mappability=Mappability_cal();
					}
				}
			}
			if($line[0] eq "foldback_3" || $line[0] eq "dup" || $line[0] eq "del" || $line[0] eq "simple_inv"){
				$N_test=N_sequence_test_1_chr();
			}elsif($line[0] eq "balanced_translocation" || $line[0] eq "unbalanced_translocation" || $line[0] eq "intra_trans_2" || $line[0] eq "insertion"){
				$N_test=N_sequence_test_2_chr();
			}else{
				$N_test=0;
			}

			if($line[0] eq "chromosome_deletion" || $line[0] eq "chromosome_amplification"){
				my $selected_chr=representative_chromosome(\@{@chr_start[$line[3]]}, \@{@chr_end[$line[3]]},\@{@chr_number[$line[3]]});
 	                  	my @where_chrs;
		                for(my $wgd_i=0; $wgd_i <=$#chr_start; $wgd_i++){
		                	my $representative_chr = representative_chromosome(\@{@chr_start[$wgd_i]}, \@{@chr_end[$wgd_i]},\@{@chr_number[$wgd_i]});
		                        if($representative_chr == $selected_chr){
		                        	$where_chrs[$#where_chrs+1]=$wgd_i;
		                        }
		                }
				if(($line[0] eq "chromosome_deletion" && $#where_chrs+1 < $PLOIDY_delta_1) || $line[0] eq "chromosome_amplification" && $#where_chrs+1 > $PLOIDY_delta_1){
					$uniform_test=0;
				}else{
					$uniform_test=1;
				}

			}elsif(@line[0] eq "balanced_translocation" || @line[0] eq "unbalanced_translocation" || $line[0] eq "intra_trans_2"){
                               my $selected_chr=representative_chromosome(\@{@chr_start[$line[5]]}, \@{@chr_end[$line[5]]},\@{@chr_number[$line[5]]});
                                my @where_chrs1;
                                for(my $wgd_i=0; $wgd_i <=$#chr_start; $wgd_i++){
                                        my $representative_chr = representative_chromosome(\@{@chr_start[$wgd_i]}, \@{@chr_end[$wgd_i]},\@{@chr_number[$wgd_i]});
                                        if($representative_chr == $selected_chr){
                                                $where_chrs1[$#where_chrs1+1]=$wgd_i;
                                        }
                                }
                                my $selected_chr=representative_chromosome(\@{@chr_start[$line[6]]}, \@{@chr_end[$line[6]]},\@{@chr_number[$line[6]]});
                                my @where_chrs2;
                                for(my $wgd_i=0; $wgd_i <=$#chr_start; $wgd_i++){
                                        my $representative_chr = representative_chromosome(\@{@chr_start[$wgd_i]}, \@{@chr_end[$wgd_i]},\@{@chr_number[$wgd_i]});
                                        if($representative_chr == $selected_chr){
                                                $where_chrs2[$#where_chrs2+1]=$wgd_i;
                                        }
                                }
				if($#where_chrs1 +1 < $PLOIDY_delta_2 || $#where_chrs1 +1 < $PLOIDY_delta_2){
					$uniform_test=0;
				}else{
					$uniform_test=1;
				}
			}else{
				$uniform_test=1;
			}



		}
                if(@line[0] eq "balanced_translocation" || @line[0] eq "unbalanced_translocation" || $line[0] eq "insertion"){
			print $fh "@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\n";
		}elsif($line[0] eq "intra_trans_2" || $line[0] eq "intra_trans"){
		        print $fh "@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t";
			######################### field 8 is printed at the sub function########################
		}
		else{
			print $fh "@line[0]\t@line[1]\t@line[2]\t@line[3]\n";
		}
	}


	if (@line[0] eq "dup"){
		duplication();
	}
        if (@line[0] eq "del"){
		deletion();
	}
        if (@line[0] eq "simple_inv"){
		simple_inversion();
	}
        if (@line[0] eq "foldback_3"){
		foldback_inversion_3();
	}
        if (@line[0] eq "insertion"){
                insertion();
        }

        if (@line[0] eq "balanced_translocation"){
		balanced_translocation();
	}
	if (@line[0] eq "unbalanced_translocation"){
		balanced_translocation();
		if(rand(1)<0.5){
			chromosome_deletion(@line[5]);
		}else{
                        chromosome_deletion(@line[6]);
		}
	}
        if (@line[0] eq "chromosome_deletion"){
		chromosome_deletion(@line[3]);
	}
        if (@line[0] eq "chromosome_amplification"){
                chromosome_amplification(@line[3]);
        }
        if (@line[0] eq "WGD"){
		my @wgd_chrs_index; ### chromosomes for duplication 
	
		for(my $chr_i_in_wgd=0; $chr_i_in_wgd<=22; $chr_i_in_wgd++){
			my @where_chrs;	
			for(my $wgd_i=0; $wgd_i <=$#chr_start; $wgd_i++){
				my $representative_chr = representative_chromosome(\@{@chr_start[$wgd_i]}, \@{@chr_end[$wgd_i]},\@{@chr_number[$wgd_i]});
				if($representative_chr == $chr_i_in_wgd){
					$where_chrs[$#where_chrs+1]=$wgd_i;
				}
			}
			if($#where_chrs != -1){
				@wgd_chrs_index[$#wgd_chrs_index+1]=$where_chrs[int(rand($#where_chrs+1))];
			}
		}
		
		for( my $wgd_chrs_index_i = 0; $wgd_chrs_index_i <=$#wgd_chrs_index; $wgd_chrs_index_i++){
 	                chromosome_amplification($wgd_chrs_index[$wgd_chrs_index_i]);
		}

        }
	if(@line[0] eq "intra_trans"){
		if(rand(1) < 0.5){
			chromosome_amplification(@line[3]);
		}
		intra_trans();
	}
        if(@line[0] eq "intra_trans_2"){
		if(rand(1) < 0.5){
			chromosome_amplification(@line[5]);
			chromosome_amplification(@line[6]);
		}
               balanced_translocation();
                if(rand(1)<0.5){
                        chromosome_deletion(@line[5]);
			$line[3]=$line[6];
                }else{
                        chromosome_deletion(@line[6]);
			$line[3]=$line[5];
                }
                intra_trans();
        }

	if($SVs ==  $number_of_germ - 1){
		$germline_chrs=$#chr_start;
		for(my $chr_index = 0; $chr_index <=$#chr_start; $chr_index++){
		        print $germ_results "\n\n$chr_index\n";
		        open(my $germ_SVs_files, '>', "germline_SVs.results.${chr_index}") or die "Could not open file 'SVs' $!";
 		       for($i = 0; $i <= $#{@chr_start[$chr_index]}; $i++){
		              print $germ_SVs_files "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
		              print $germ_results "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
        		}
		}
	}
}

for(my $chr_index = 0; $chr_index <=$#chr_start; $chr_index++){
	print "\n\n$chr_index\n";
        open(my $SVs_files, '>', "SVs.results.${chr_index}") or die "Could not open file 'SVs' $!";

	for($i = 0; $i <= $#{@chr_start[$chr_index]}; $i++){
              print $SVs_files "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
	      print "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
	}
}

`Rscript chromosomes_graph.R $#chr_start`;

`perl results_to_fa.pl $#chr_start $IND $G1000_dir`;
`perl germ_results_to_fa.pl $germline_chrs $IND $G1000_dir`;


sub SV_roulette_germ {
 my $SV_type;
  if($_[0]<17/(2788+37+17)){
        $SV_type="dup";
  }
  elsif($_[0]<(17+2788)/(2788+37+17)){
        $SV_type="del";
  }
  elsif($_[0]<1){
        $SV_type="simple_inv";
  }
}

#dup 0.17 del 0.17 inv 0.02 balanced translocation 0.035 unbalanced translocation 0.035 chromosome amplication 0.02 chromosome deletion 0.02 intra_trans 0.015 intra_trans_2 0.01 ins 0.07
#my $SV_prior=0.4+0.4+0.05+0.05+0.05+0.05+0.05;
sub SV_roulette {
  my $SV_type;
  if($_[0]<0.17/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="dup";
  }
  elsif($_[0]<(0.17+0.17)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="del";
  }
  elsif($_[0]<(0.17+0.17+0.02)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="simple_inv";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="balanced_translocation";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035+0.035)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="unbalanced_translocation";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035+0.035+0.02)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="chromosome_amplification";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035+0.035+0.02+0.02)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="chromosome_deletion";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
	$SV_type="intra_trans";
  }
  elsif($_[0]<(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
        $SV_type="intra_trans_2";
  }elsif($_[0]<(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)/(0.17+0.17+0.02+0.035+0.035+0.02+0.02+0.015+0.01+0.07)*0.99){
	$SV_type="insertion";
  }

  elsif($_[0]<1){
        if($number_of_WGD <$PLOIDY-2){
                $SV_type="WGD";
                $number_of_WGD=$number_of_WGD+1;
        }else{
                $SV_type="chromosome_amplification";
        }
  }
  return $SV_type;
}


sub duplication {
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                my @hap=();
                my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};
         
	        my $last, $split_5, $split_3, $SV_5, $SV_3;

                my $SV_5_searched=0;
                my $SV_3_searched=0;
                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @line[1]=@start[$SV_5]+@line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @line[1]=@start[$SV_5]+@line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @line[2]=@start[$SV_3]+@line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @line[2]=@start[$SV_3]+@line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }
                if(@line[1]==@start[$SV_5]){
                        $split_5=0;
                }else{
                        $split_5=1;
                }
                if(@line[2]==@end[$SV_3]){
                        $split_3=0;
                }else{
                        $split_3=1;
                }
                $last=$#start;
                $dup_size = $SV_3-$SV_5+1;

                @dup_start[0 .. ($SV_3-$SV_5)]= @start[$SV_5 .. $SV_3];
                @dup_end[0 .. ($SV_3-$SV_5)]= @end[$SV_5 .. $SV_3];
                @dup_start[0]=@line[1];
                @dup_end[($SV_3-$SV_5)]=@line[2];
                @dup_hap[0 .. ($SV_3-$SV_5)]= @hap[$SV_5 .. $SV_3];
                @dup_number[0 .. ($SV_3-$SV_5)]= @number[$SV_5 .. $SV_3];

                if($split_5==1){
                        @segment_start[0]=@start[$SV_5];
                        @segment_end[0]=@line[1]-1;
                        @segment_hap[0]=@hap[$SV_5];
                        @segment_number[0]=@number[$SV_5];
                }
                if($split_3==1){
                        @segment_start[0+2*($SV_3-$SV_5)+$split_5+2]=@line[2]+1;
                        @segment_end[0+2*($SV_3-$SV_5)+$split_5+2]=@end[$SV_3];
                        @segment_hap[0+2*($SV_3-$SV_5)+$split_5+2]=@hap[$SV_3];
                        @segment_number[0+2*($SV_3-$SV_5)+$split_5+2]=@number[$SV_3];
                }

                @segment_start[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_start[0 .. ($SV_3-$SV_5)];
                @segment_start[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_start[0 .. ($SV_3-$SV_5)];
                @segment_end[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_end[0 .. ($SV_3-$SV_5)];
                @segment_end[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_end[0 .. ($SV_3-$SV_5)];
                @segment_hap[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_hap[0 .. ($SV_3-$SV_5)];
                @segment_hap[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_hap[0 .. ($SV_3-$SV_5)];
                @segment_number[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_number[0 .. ($SV_3-$SV_5)];
                @segment_number[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_number[0 .. ($SV_3-$SV_5)];

                $segment_length=$#segment_start+1;

                if($last>$SV_3){
                        @start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
                        @end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
                        @hap[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @hap[$SV_3+1 .. $last];
                        @number[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @number[$SV_3+1 .. $last];
                }

                @start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
                @end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];
                @hap[$SV_5 .. $SV_5+$segment_length-1]=@segment_hap[0 .. $segment_length-1];
                @number[$SV_5 .. $SV_5+$segment_length-1]=@segment_number[0 .. $segment_length-1];

               @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;
                @{@chr_hap[@line[3]]}=@hap;
                @{@chr_number[@line[3]]}=@number;

}

sub intra_trans {
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
		my @dup_hap=();my @dup_number=();my @inv_hap=();my @inv_number=();
                my @hap=();
                my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};
		my $intra_trans_n; my $intra_trans_size_min;
		
		my $trans_while_iter=0;
		do{
                	$intra_trans_n = $intra_trans_n_min+int(rand($intra_trans_n_max));
                	$intra_trans_size_min = int($cumul[$#cumul]/$intra_trans_n);
			$trans_while_iter=$trans_while_iter+1;
		}while($intra_trans_size_min < $intra_trans_size_MIN || $trans_while_iter < 100);
		
		if($intra_trans_size_min < $intra_trans_size_MIN){
			return;
		}

		print $fh "$intra_trans_n\n";

		my @new_start=();my @new_end=(); my @new_hap=(); my @new_number=(); my @new_cumul=();
		for(my $trans_i = 0; $trans_i < $intra_trans_n;$trans_i++){
			@dup_start=();@dup_end=();@dup_hap=();@dup_number=(); @inv_start=();@inv_end=();@inv_hap=();@inv_number=();
			do{
				$intra_trans_size= $intra_trans_size_min+int(rand(2*$intra_trans_size_min));
				$line[1] = int(rand($cumul[$#cumul]-$intra_trans_size-100))+1;
				$line[2] = $line[1]+$intra_trans_size;
			}while(N_sequence_test_1_chr() || $line[1]>$line[2]);
			my $last, $split_5, $split_3, $SV_5, $SV_3;
			my $SV_5_searched=0;
			my $SV_3_searched=0;
			for( $j = 0; $j <= $#cumul; $j++){
				if($SV_5_searched == 0 && @line[1] <= @cumul[$j]){
				       $SV_5 = $j;
				       if($SV_5 != 0){
					       @line[1]=@start[$SV_5]+@line[1]-@cumul[$SV_5-1]-1;
					}else{
					       @line[1]=@start[$SV_5]+@line[1]-1;
					}
					$SV_5_searched=1;
				}
				if($SV_3_searched == 0 && @line[2] <= @cumul[$j]){
				       $SV_3 = $j;
				       if($SV_3 != 0){
					       @line[2]=@start[$SV_3]+@line[2]-@cumul[$SV_3-1]-1;
					}else{
					       @line[2]=@start[$SV_3]+@line[2]-1;
					}
					$SV_3_searched=1;
				}

			}
			if(@line[1]==@start[$SV_5]){
				$split_5=0;
			}else{
				$split_5=1;
			}
			if(@line[2]==@end[$SV_3]){
				$split_3=0;
			}else{
				$split_3=1;
			}
			$last=$#start;
			$dup_size = $SV_3-$SV_5+1;

			@dup_start[0 .. ($SV_3-$SV_5)]= @start[$SV_5 .. $SV_3];
			@dup_end[0 .. ($SV_3-$SV_5)]= @end[$SV_5 .. $SV_3];
 	               @dup_start[0]=@line[1];
	                @dup_end[($SV_3-$SV_5)]=@line[2];
	                @dup_hap[0 .. ($SV_3-$SV_5)]= @hap[$SV_5 .. $SV_3];
	                @dup_number[0 .. ($SV_3-$SV_5)]= @number[$SV_5 .. $SV_3];

			if( rand(1) < 0.5){
		                for(my $k=0; $k <= $SV_3-$SV_5; $k++){
		                        @inv_start[$k]=-@dup_end[$SV_3-$SV_5-$k];
		                        @inv_end[$k]=-@dup_start[$SV_3-$SV_5-$k];
		                        @inv_hap[$k]=@dup_hap[$SV_3-$SV_5-$k];
		                        @inv_number[$k]=@dup_number[$SV_3-$SV_5-$k];
				}
				@dup_start=@inv_start;
				@dup_end=@inv_end;
				@dup_hap=@inv_hap;
				@dup_number=@inv_number;
			}
			@new_start[($#new_start+1) .. ($SV_3-$SV_5+$#new_start+1)]=@dup_start[0 .. ($SV_3-$SV_5)];
                        @new_end[($#new_end+1) .. ($SV_3-$SV_5+$#new_end+1)]=@dup_end[0 .. ($SV_3-$SV_5)];
                        @new_hap[($#new_hap+1) .. ($SV_3-$SV_5+$#new_hap+1)]=@dup_hap[0 .. ($SV_3-$SV_5)];
                        @new_number[($#new_number+1) .. ($SV_3-$SV_5+$#new_number+1)]=@dup_number[0 .. ($SV_3-$SV_5)];

		}
                @new_cumul[0]=@new_end[0]-@new_start[0]+1;
                for(my $new_cumul_i = 1; $new_cumul_i <=$#new_end; $new_cumul_i++){
                        @new_cumul[$new_cumul_i]=@new_cumul[$new_cumul_i-1]+@new_end[$new_cumul_i]-@new_start[$new_cumul_i]+1;
                }

                @{@chr_start[@line[3]]}=@new_start;
                @{@chr_end[@line[3]]}=@new_end;
                @{@chr_cumul[@line[3]]}=@new_cumul;
                @{@chr_hap[@line[3]]}=@new_hap;
                @{@chr_number[@line[3]]}=@new_number;

}

sub insertion {
                my @start1=(); my @end1=(); my @cumul1=(); my @start2=(); my @end2=(); my @cumul2=(); my @hap1=(); my @hap2=(); my @number1=(); my @number2=();
                my @segment1_start=(); my @segment1_end=(); my @segment1_hap=(); my @segment1_number=(); my @segment2_start=(); my @segment2_end=(); my @segment2_hap=(); my @segment2_number=(); my @segment1_cumul=();
                @start1=@{@chr_start[@line[5]]};
                @end1=@{@chr_end[@line[5]]};
                @cumul1=@{@chr_cumul[@line[5]]};
                @hap1=@{@chr_hap[@line[5]]};
                @number1=@{@chr_number[@line[5]]};

                @start2=@{@chr_start[@line[6]]};
                @end2=@{@chr_end[@line[6]]};
                @cumul2=@{@chr_cumul[@line[6]]};
                @hap2=@{@chr_hap[@line[6]]};
                @number2=@{@chr_number[@line[6]]};

                my $SV_5_searched1=0;
                my $SV_5_searched2=0;
                my $SV_3_searched1=0;
                my $SV_3_searched2=0;
                for( $j = 0; $j <= $#cumul1; $j++){
                        if($SV_5_searched1 == 0 && @line[1] <= @cumul1[$j]){
                               $SV_5_1 = $j;
                               if($SV_5_1 != 0){
                                       @line[1]=@start1[$SV_5_1]+@line[1]-@cumul1[$SV_5_1-1]-1;
                                }else{
                                       @line[1]=@start1[$SV_5_1]+@line[1]-1;
                                }
                                $SV_5_searched1=1;
                        }
                        if($SV_3_searched1 == 0 && @line[2] <= @cumul1[$j]){
                               $SV_3_1 = $j;
                               if($SV_3_1 != 0){
                                       @line[2]=@start1[$SV_3_1]+@line[2]-@cumul1[$SV_3_1-1]-1;
                                }else{
                                       @line[2]=@start1[$SV_3_1]+@line[2]-1;
                                }
                                $SV_3_searched1=1;
                        }
                }
                for( $j = 0; $j <= $#cumul2; $j++){
                        if($SV_5_searched2 == 0 && @line[3] <= @cumul2[$j]){
                               $SV_5_2 = $j;
                               if($SV_5_2 != 0){
                                       @line[3]=@start2[$SV_5_2]+@line[3]-@cumul2[$SV_5_2-1]-1;
                                }else{
                                       @line[3]=@start2[$SV_5_2]+@line[3]-1;
                                }
                                $SV_5_searched2=1;
                        }
                        if($SV_3_searched2 == 0 && @line[4] <= @cumul2[$j]){
                               $SV_3_2 = $j;
                               if($SV_3_2 != 0){
                                       @line[4]=@start2[$SV_3_2]+@line[4]-@cumul2[$SV_3_2-1]-1;
                                }else{
                                       @line[4]=@start2[$SV_3_2]+@line[4]-1;
                                }
                                $SV_3_searched2=1;
                        }
                }
                @segment1_start[0 .. $SV_5_1] = @start1[0 .. $SV_5_1];
                @segment1_hap[0 .. $SV_5_1] = @hap1[0 .. $SV_5_1];
                @segment1_number[0 .. $SV_5_1] = @number1[0 .. $SV_5_1];
                @segment1_end[0 .. $SV_5_1] = @end1[0 .. $SV_5_1];
		@segment1_end[$SV_5_1]=$line[1];

                @segment1_start[$SV_5_1+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2] = @start2[$SV_5_2 .. $SV_3_2];
                @segment1_hap[$SV_5_1+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2] = @hap2[$SV_5_2 .. $SV_3_2];
                @segment1_number[$SV_5_1+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2] = @number2[$SV_5_2 .. $SV_3_2];
                @segment1_end[$SV_5_1+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2] = @end2[$SV_5_2 .. $SV_3_2];
		@segment1_start[$SV_5_1+1]=$line[3];
		@segment1_end[$SV_5_1+1+$SV_3_2-$SV_5_2]=$line[4];

                @segment1_start[$SV_5_1+1+$SV_3_2-$SV_5_2+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2+1+$#start1-$SV_3_1] = @start1[$SV_3_1 .. $#start1];
                @segment1_hap[$SV_5_1+1+$SV_3_2-$SV_5_2+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2+1+$#start1-$SV_3_1] = @hap1[$SV_3_1 .. $#start1];
                @segment1_number[$SV_5_1+1+$SV_3_2-$SV_5_2+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2+1+$#start1-$SV_3_1] = @number1[$SV_3_1 .. $#start1];
                @segment1_end[$SV_5_1+1+$SV_3_2-$SV_5_2+1 .. $SV_5_1+1+$SV_3_2-$SV_5_2+1+$#start1-$SV_3_1] = @end1[$SV_3_1 .. $#start1];
		@segment1_start[$SV_5_1+1+$SV_3_2-$SV_5_2+1]=$line[2];

                @segment1_cumul[0]=@segment1_end[0]-@segment1_start[0]+1;
                for($i = 1; $i <=$#segment1_end; $i++){
                        @segment1_cumul[$i]=@segment1_cumul[$i-1]+@segment1_end[$i]-@segment1_start[$i]+1;
                }

               @{@chr_start[@line[5]]}=@segment1_start;
                @{@chr_end[@line[5]]}=@segment1_end;
                @{@chr_cumul[@line[5]]}=@segment1_cumul;
                @{@chr_hap[@line[5]]}=@segment1_hap;
                @{@chr_number[@line[5]]}=@segment1_number;

}
sub deletion {
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                my @hap=();
                my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};

                my $last, $split_5, $split_3, $SV_5, $SV_3;
                my $SV_5_searched=0;
                my $SV_3_searched=0;
                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @line[1]=@start[$SV_5]+@line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @line[1]=@start[$SV_5]+@line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @line[2]=@start[$SV_3]+@line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @line[2]=@start[$SV_3]+@line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }

               if(@line[1]==@start[$SV_5]){
                        $split_5=0;
                }else{
                        $split_5=1;
                }
                if(@line[2]==@end[$SV_3]){
                        $split_3=0;
                }else{
                        $split_3=1;
                }
                $last=$#start;

                if($split_5==1){
                        @segment_start[0]=$start[$SV_5];
                        @segment_end[0]=@line[1]-1;
                        @segment_hap[0]=@hap[$SV_5];
                        @segment_number[0]=@number[$SV_5];

                }
                if($split_3==1){
                        @segment_start[$split_5+$split_3-1]=@line[2]+1;
                        @segment_end[$split_5+$split_3-1]=@end[$SV_3];
                        @segment_hap[$split_5+$split_3-1]=@hap[$SV_3];
                        @segment_number[$split_5+$split_3-1]=@number[$SV_3];

                }
                $segment_length=$#segment_start+1;
                if($last>$SV_3){
                        @start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
                        @end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
                        @hap[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @hap[$SV_3+1 .. $last];
                        @number[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @number[$SV_3+1 .. $last];

                 }
                if($SV_5+$segment_length + $last - $SV_3 -1 < $last){
                        delete @start[$SV_5+$segment_length + $last - $SV_3 .. $last];
                        delete @end[$SV_5+$segment_length + $last - $SV_3 .. $last];
                        delete @cumul[$SV_5+$segment_length + $last - $SV_3 .. $last];
                        delete @hap[$SV_5+$segment_length + $last - $SV_3 .. $last];
                        delete @number[$SV_5+$segment_length + $last - $SV_3 .. $last];

                }
                @start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
                @end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];
                @hap[$SV_5 .. $SV_5+$segment_length-1]=@segment_hap[0 .. $segment_length-1];
                @number[$SV_5 .. $SV_5+$segment_length-1]=@segment_number[0 .. $segment_length-1];

               @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;
                @{@chr_hap[@line[3]]}=@hap;
                @{@chr_number[@line[3]]}=@number;

}

sub simple_inversion{
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                my @hap=();
                my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};

                my $last, $split_5, $split_3, $SV_5, $SV_3;

                my $SV_5_searched=0;
                my $SV_3_searched=0;
                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @line[1]=@start[$SV_5]+@line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @line[1]=@start[$SV_5]+@line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @line[2]=@start[$SV_3]+@line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @line[2]=@start[$SV_3]+@line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }
                if(@line[1]==@start[$SV_5]){
                        $split_5=0;
                }else{
                        $split_5=1;
                }
                if(@line[2]==@end[$SV_3]){
                        $split_3=0;
                }else{
                        $split_3=1;
                }
                $last=$#start;
                $dup_size = $SV_3-$SV_5+1;

                @dup_start[0 .. ($SV_3-$SV_5)]= @start[$SV_5 .. $SV_3];
                @dup_end[0 .. ($SV_3-$SV_5)]= @end[$SV_5 .. $SV_3];
                @dup_hap[0 .. ($SV_3-$SV_5)]= @hap[$SV_5 .. $SV_3];
                @dup_number[0 .. ($SV_3-$SV_5)]= @number[$SV_5 .. $SV_3];

                @dup_start[0]=@line[1];
                @dup_end[($SV_3-$SV_5)]=@line[2];

                if($split_5==1){
                        @segment_start[0]=@start[$SV_5];
                        @segment_end[0]=@line[1]-1;
                        @segment_hap[0]=@hap[$SV_5];
                        @segment_number[0]=@number[$SV_5];

                }
                if($split_3==1){
                        @segment_start[0+$SV_3-$SV_5+$split_5+$split_3]=@line[2]+1;
                        @segment_end[0+$SV_3-$SV_5+$split_5+$split_3]=@end[$SV_3];
                        @segment_hap[0+$SV_3-$SV_5+$split_5+$split_3]=@hap[$SV_3];
                        @segment_number[0+$SV_3-$SV_5+$split_5+$split_3]=@number[$SV_3];

                }

                for(my $k=0; $k <= $SV_3-$SV_5; $k++){
                        @inv_start[$k]=-@dup_end[$SV_3-$SV_5-$k];
                        @inv_end[$k]=-@dup_start[$SV_3-$SV_5-$k];
                        @inv_hap[$k]=@dup_hap[$SV_3-$SV_5-$k];
                        @inv_number[$k]=@dup_number[$SV_3-$SV_5-$k];

                }
                @segment_start[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_start[0 .. ($SV_3-$SV_5)];
                @segment_end[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_end[0 .. ($SV_3-$SV_5)];
                @segment_hap[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_hap[0 .. ($SV_3-$SV_5)];
                @segment_number[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_number[0 .. ($SV_3-$SV_5)];


                $segment_length=$#segment_start+1;

                if($last>$SV_3){
                        @start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
                        @end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
                        @hap[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @hap[$SV_3+1 .. $last];
                        @number[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @number[$SV_3+1 .. $last];

                }

                @start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
                @end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];
                @hap[$SV_5 .. $SV_5+$segment_length-1]=@segment_hap[0 .. $segment_length-1];
                @number[$SV_5 .. $SV_5+$segment_length-1]=@segment_number[0 .. $segment_length-1];

                @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;
                @{@chr_hap[@line[3]]}=@hap;
                @{@chr_number[@line[3]]}=@number;

}

sub foldback_inversion_3{
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                my @hap=();
                my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};

                my $last, $split_5, $split_3, $SV_5, $SV_3;

                my $SV_5_searched=0;
                my $SV_3_searched=0;
                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @line[1]=@start[$SV_5]+@line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @line[1]=@start[$SV_5]+@line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @line[2]=@start[$SV_3]+@line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @line[2]=@start[$SV_3]+@line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }
                if(@line[1]==@start[$SV_5]){
                        $split_5=0;
                }else{
                        $split_5=1;
                }
                if(@line[2]==@end[$SV_3]){
                        $split_3=0;
                }else{
                        $split_3=1;
                }

               $last=$#start;

               if($SV_3 < $last){
                        delete @start[$SV_3+1 .. $last];
                        delete @end[$SV_3+1 .. $last];
                        delete @cumul[$SV_3+1 .. $last];
                        delete @hap[$SV_3+1 .. $last];
                        delete @number[$SV_3+1 .. $last];
               }
                @end[$SV_3]=$line[2];

                @segment_start[0 .. $SV_5] = @start[0 .. $SV_5];
                @segment_end[0 .. $SV_5] = @end[0 .. $SV_5];
                @segment_end[$SV_5] = @line[1];
                @segment_hap[0 .. $SV_5] = @hap[0 .. $SV_5];
                @segment_number[0 .. $SV_5] = @number[0 .. $SV_5];


               for(my $k=0; $k <= $SV_5; $k++){
                        @inv_start[$k]=-@segment_end[$SV_5-$k];
                        @inv_end[$k]=-@segment_start[$SV_5-$k];
                        @inv_hap[$k]=@segment_hap[$SV_5-$k];
                        @inv_number[$k]=@segment_number[$SV_5-$k];

                }

                @start[$SV_3+1 .. $SV_3+$#inv_start+1] = @inv_start[0 .. $#inv_start];
                @end[$SV_3+1 .. $SV_3+$#inv_end+1] = @inv_end[0 .. $#inv_end];
                @hap[$SV_3+1 .. $SV_3+$#inv_end+1] = @inv_hap[0 .. $#inv_end];
                @number[$SV_3+1 .. $SV_3+$#inv_end+1] = @inv_number[0 .. $#inv_end];

                @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;
                @{@chr_hap[@line[3]]}=@hap;
                @{@chr_number[@line[3]]}=@number;

}

sub balanced_translocation {
                my @start1=(); my @end1=(); my @cumul1=(); my @start2=(); my @end2=(); my @cumul2=(); my @hap1=(); my @hap2=(); my @number1=(); my @number2=();
                my @segment1_start=(); my @segment1_end=(); my @segment1_hap=(); my @segment1_number=(); my @segment2_start=(); my @segment2_end=(); my @segment2_hap=(); my @segment2_number=();
                @start1=@{@chr_start[@line[5]]};
                @end1=@{@chr_end[@line[5]]};
                @cumul1=@{@chr_cumul[@line[5]]};
                @hap1=@{@chr_hap[@line[5]]};
                @number1=@{@chr_number[@line[5]]};

                @start2=@{@chr_start[@line[6]]};
                @end2=@{@chr_end[@line[6]]};
                @cumul2=@{@chr_cumul[@line[6]]};
                @hap2=@{@chr_hap[@line[6]]};
                @number2=@{@chr_number[@line[6]]};

                $last1=$#start1;
                $last2=$#start2;

                my $SV_5_searched1=0;
                my $SV_5_searched2=0;
                my $SV_3_searched1=0;
                my $SV_3_searched2=0;
                for( $j = 0; $j <= $#cumul1; $j++){
                        if($SV_5_searched1 == 0 && @line[1] <= @cumul1[$j]){
                               $SV_5_1 = $j;
                               if($SV_5_1 != 0){
                                       @line[1]=@start1[$SV_5_1]+@line[1]-@cumul1[$SV_5_1-1]-1;
                                }else{
                                       @line[1]=@start1[$SV_5_1]+@line[1]-1;
                                }
                                $SV_5_searched1=1;
                        }
                        if($SV_3_searched1 == 0 && @line[2] <= @cumul1[$j]){
                               $SV_3_1 = $j;
                               if($SV_3_1 != 0){
                                       @line[2]=@start1[$SV_3_1]+@line[2]-@cumul1[$SV_3_1-1]-1;
                                }else{
                                       @line[2]=@start1[$SV_3_1]+@line[2]-1;
                                }
                                $SV_3_searched1=1;
                        }
                }
                for( $j = 0; $j <= $#cumul2; $j++){
                        if($SV_5_searched2 == 0 && @line[3] <= @cumul2[$j]){
                               $SV_5_2 = $j;
                               if($SV_5_2 != 0){
                                       @line[3]=@start2[$SV_5_2]+@line[3]-@cumul2[$SV_5_2-1]-1;
                                }else{
                                       @line[3]=@start2[$SV_5_2]+@line[3]-1;
                                }
                                $SV_5_searched2=1;
                        }
                        if($SV_3_searched2 == 0 && @line[4] <= @cumul2[$j]){
                               $SV_3_2 = $j;
                               if($SV_3_2 != 0){
                                       @line[4]=@start2[$SV_3_2]+@line[4]-@cumul2[$SV_3_2-1]-1;
                                }else{
                                       @line[4]=@start2[$SV_3_2]+@line[4]-1;
                                }
                                $SV_3_searched2=1;
                        }
                }
                @segment1_start[0 .. $#start2-$SV_3_2 ] = @start2[$SV_3_2 .. $#start2];
                @segment1_end[0 .. $#start2-$SV_3_2 ] = @end2[$SV_3_2 .. $#start2];
                @segment1_hap[0 .. $#start2-$SV_3_2 ] = @hap2[$SV_3_2 .. $#start2];
                @segment1_number[0 .. $#start2-$SV_3_2 ] = @number2[$SV_3_2 .. $#start2];
                @segment1_start[0] = @line[4];

                @segment2_start[0 .. $#start1-$SV_3_1 ] = @start1[$SV_3_1 .. $#start1];
                @segment2_end[0 .. $#start1-$SV_3_1 ] = @end1[$SV_3_1 .. $#start1];
                @segment2_hap[0 .. $#start1-$SV_3_1 ] = @hap1[$SV_3_1 .. $#start1];
                @segment2_number[0 .. $#start1-$SV_3_1 ] = @number1[$SV_3_1 .. $#start1];
                @segment2_start[0] = @line[2];

                @end1[$SV_5_1]=@line[1];

                @start1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_start[0 .. $#segment1_start];
                @end1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_end[0 .. $#segment1_end];
                @hap1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_hap[0 .. $#segment1_end];
                @number1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_number[0 .. $#segment1_end];

                if($last1>$#segment1_start+$SV_5_1+1){
                        delete @start1[$#segment1_start+$SV_5_1+2..$last1];
                        delete @end1[$#segment1_end+$SV_5_1+2..$last1];
                        delete @cumul1[$#segment1_end+$SV_5_1+2..$last1];
                        delete @hap1[$#segment1_end+$SV_5_1+2..$last1];
                        delete @number1[$#segment1_end+$SV_5_1+2..$last1];

                }
                @end2[$SV_5_2]=@line[3];
                @start2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_start[0 .. $#segment2_start];
                @end2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_end[0 .. $#segment2_end];
                @hap2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_hap[0 .. $#segment2_end];
                @number2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_number[0 .. $#segment2_end];
                if($last2>$#segment2_start+$SV_5_2+1){
                        delete @start2[$#segment2_start+$SV_5_2+2..$last2];
                        delete @end2[$#segment2_end+$SV_5_2+2..$last2];
                        delete @cumul2[$#segment2_end+$SV_5_2+2..$last2];
                        delete @hap2[$#segment2_end+$SV_5_2+2..$last2];
                        delete @number2[$#segment2_end+$SV_5_2+2..$last2];

                }

                @cumul1[0]=@end1[0]-@start1[0]+1;
                for($i = 1; $i <=$#end1; $i++){
                        @cumul1[$i]=@cumul1[$i-1]+@end1[$i]-@start1[$i]+1;
                }
                @{@chr_start[@line[5]]}=@start1;
                @{@chr_end[@line[5]]}=@end1;
                @{@chr_cumul[@line[5]]}=@cumul1;
                @{@chr_hap[@line[5]]}=@hap1;
                @{@chr_number[@line[5]]}=@number1;

                @cumul2[0]=@end2[0]-@start2[0]+1;
                for($i = 1; $i <=$#end2; $i++){
                        @cumul2[$i]=@cumul2[$i-1]+@end2[$i]-@start2[$i]+1;
                }

                @{@chr_start[@line[6]]}=@start2;
                @{@chr_end[@line[6]]}=@end2;
                @{@chr_cumul[@line[6]]}=@cumul2;
                @{@chr_hap[@line[6]]}=@hap2;
                @{@chr_number[@line[6]]}=@number2;

}

sub chromosome_deletion{

@{@chr_start[$_[0]]}=@{@chr_start[$#chr_start]};
@{@chr_end[$_[0]]}=@{@chr_end[$#chr_end]};
@{@chr_cumul[$_[0]]}=@{@chr_cumul[$#chr_cumul]};
@{@chr_hap[$_[0]]}=@{@chr_hap[$#chr_hap]};
@{@chr_number[$_[0]]}=@{@chr_number[$#chr_number]};

delete @chr_start[$#chr_start];
delete @chr_end[$#chr_end];
delete @chr_cumul[$#chr_cumul];
delete @chr_hap[$#chr_hap];
delete @chr_number[$#chr_number];

}

sub chromosome_amplification{

@{$chr_start[$#chr_start+1]}=@{@chr_start[$_[0]]};
@{$chr_end[$#chr_end+1]}=@{@chr_end[$_[0]]};
@{$chr_cumul[$#chr_cumul+1]}=@{@chr_cumul[$_[0]]};
@{$chr_hap[$#chr_hap+1]}=@{@chr_hap[$_[0]]};
@{$chr_number[$#chr_number+1]}=@{@chr_number[$_[0]]};

}


sub representative_chromosome{
	my ($re_chr_start,$re_chr_end,$re_chr_name) = @_;
	my @re_chr_start=@{$re_chr_start};
        my @re_chr_end=@{$re_chr_end};
        my @re_chr_name=@{$re_chr_name};
	my @re;
	for( my $re_i=0; $re_i <=22; $re_i++){
		$re[$re_i]=0;
	}

	for (my $re_i=0; $re_i <=$#re_chr_start; $re_i++){
		if($re_chr_name[$re_i] ne "X"){
			$re[$re_chr_name[$re_i]-1]=$re[$re_chr_name[$re_i]-1]+abs($re_chr_start[$re_i]-$re_chr_end[$re_i]);
		}else{
                        $re[22]=$re[22]+abs($re_chr_start[$re_i]-$re_chr_end[$re_i]);
		}
	}
	my $max_chr;
	my $max=0;
	for( my $re_i=0; $re_i <=$#re; $re_i++){
		if($re[$re_i] > $max){
			$max_chr=$re_i;
			$max=$re[$re_i];
		}
	}
	
	return $max_chr;
}

sub N_sequence_test_1_chr{
                my @start=();my @end=();my @cumul=();my @hap=();my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};
		my $j;
                my $SV_5, $SV_3;
                my $SV_5_searched=0;
                my $SV_3_searched=0;
		my @N_sequence_line=@line;

                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @N_sequence_line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @N_sequence_line[1]=@start[$SV_5]+@N_sequence_line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @N_sequence_line[1]=@start[$SV_5]+@N_sequence_line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @N_sequence_line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @N_sequence_line[2]=@start[$SV_3]+@N_sequence_line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @N_sequence_line[2]=@start[$SV_3]+@N_sequence_line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }
		$N_sequence_line[1]=abs $N_sequence_line[1];
                $N_sequence_line[2]=abs $N_sequence_line[2];

		my $N_chr1=$number[$SV_5];
		if($number[$SV_5] eq "X"){
			$N_chr1=22;
		}else{
			$N_chr1=$N_chr1-1;
		}
		my $N_seq1=$chr_fa[$hap[$SV_5]-1][$N_chr1]->seq("$hap[$SV_5]_$number[$SV_5]", $N_sequence_line[1]=>$N_sequence_line[1]);


                my $N_chr2=$number[$SV_3];
                if($number[$SV_3] eq "X"){
                        $N_chr2=22;
                }else{
			$N_chr2=$N_chr2-1;
		}
                my $N_seq2=$chr_fa[$hap[$SV_3]-1][$N_chr2]->seq("$hap[$SV_3]_$number[$SV_3]", $N_sequence_line[2]=>$N_sequence_line[2]);

		if($N_seq1 eq "N" || $N_seq2 eq "N"){
			print $log "@line[0]\t@line[1]\t@line[2]\t@line[3]\t$N_seq1\t$N_seq2\n";
			return 1;
		}else{
                        print $log "@line[0]\t@line[1]\t@line[2]\t@line[3]\t$N_seq1\t$N_seq2\n";
			return 0;
		}
}
sub N_sequence_test_2_chr{
                my @start1=(); my @end1=(); my @cumul1=(); my @start2=(); my @end2=(); my @cumul2=(); my @hap1=(); my @hap2=(); my @number1=(); my @number2=();
                @start1=@{@chr_start[@line[5]]};
                @end1=@{@chr_end[@line[5]]};
                @cumul1=@{@chr_cumul[@line[5]]};
                @hap1=@{@chr_hap[@line[5]]};
                @number1=@{@chr_number[@line[5]]};

                @start2=@{@chr_start[@line[6]]};
                @end2=@{@chr_end[@line[6]]};
                @cumul2=@{@chr_cumul[@line[6]]};
                @hap2=@{@chr_hap[@line[6]]};
                @number2=@{@chr_number[@line[6]]};

                my @N_sequence_line=@line;
                my $SV_5_searched1=0;
                my $SV_5_searched2=0;
                my $SV_3_searched1=0;
                my $SV_3_searched2=0;

                for( $j = 0; $j <= $#cumul1; $j++){
                        if($SV_5_searched1 == 0 && @line[1] <= @cumul1[$j]){
                               $SV_5_1 = $j;
                               if($SV_5_1 != 0){
                                       @N_sequence_line[1]=@start1[$SV_5_1]+@N_sequence_line[1]-@cumul1[$SV_5_1-1]-1;
                                }else{
                                       @N_sequence_line[1]=@start1[$SV_5_1]+@N_sequence_line[1]-1;
                                }
                                $SV_5_searched1=1;
                        }
                        if($SV_3_searched1 == 0 && @line[2] <= @cumul1[$j]){
                               $SV_3_1 = $j;
                               if($SV_3_1 != 0){
                                       @N_sequence_line[2]=@start1[$SV_3_1]+@N_sequence_line[2]-@cumul1[$SV_3_1-1]-1;
                                }else{
                                       @N_sequence_line[2]=@start1[$SV_3_1]+@N_sequence_line[2]-1;
                                }
                                $SV_3_searched1=1;
                        }
                }
                for( $j = 0; $j <= $#cumul2; $j++){
                        if($SV_5_searched2 == 0 && @line[3] <= @cumul2[$j]){
                               $SV_5_2 = $j;
                               if($SV_5_2 != 0){
                                       @N_sequence_line[3]=@start2[$SV_5_2]+@N_sequence_line[3]-@cumul2[$SV_5_2-1]-1;
                                }else{
                                       @N_sequence_line[3]=@start2[$SV_5_2]+@N_sequence_line[3]-1;
                                }
                                $SV_5_searched2=1;
                        }
                        if($SV_3_searched2 == 0 && @line[4] <= @cumul2[$j]){
                               $SV_3_2 = $j;
                               if($SV_3_2 != 0){
                                       @N_sequence_line[4]=@start2[$SV_3_2]+@N_sequence_line[4]-@cumul2[$SV_3_2-1]-1;
                                }else{
                                       @N_sequence_line[4]=@start2[$SV_3_2]+@N_sequence_line[4]-1;
                                }
                                $SV_3_searched2=1;
                        }
                }
                $N_sequence_line[1]=abs $N_sequence_line[1];
                $N_sequence_line[2]=abs $N_sequence_line[2];
                $N_sequence_line[3]=abs $N_sequence_line[3];
                $N_sequence_line[4]=abs $N_sequence_line[4];

                my $N_chr1=$number1[$SV_5_1];
                my $N_chr2=$number1[$SV_3_1];

                if($number1[$SV_5_1] eq "X"){
                        $N_chr1=22;
                }else{
			$N_chr1=$N_chr1-1;
		}
                if($number1[$SV_3_1] eq "X"){
                        $N_chr2=22;
                }else{
                        $N_chr2=$N_chr2-1;
		}
                my $N_seq1=$chr_fa[$hap1[$SV_5_1]-1][$N_chr1]->seq("$hap1[$SV_5_1]_$number1[$SV_5_1]", $N_sequence_line[1]=>$N_sequence_line[1]);
                my $N_seq2=$chr_fa[$hap1[$SV_3_1]-1][$N_chr2]->seq("$hap1[$SV_3_1]_$number1[$SV_3_1]", $N_sequence_line[2]=>$N_sequence_line[2]);

                my $N_chr3=$number2[$SV_5_2];
                my $N_chr4=$number2[$SV_3_2];
                if($number2[$SV_5_2] eq "X"){
                        $N_chr3=22;
                }else{
                        $N_chr3=$N_chr3-1;
		}
                if($number2[$SV_3_2] eq "X"){
                        $N_chr4=22;
                }else{
                        $N_chr4=$N_chr4-1;
		}
                my $N_seq3=$chr_fa[$hap2[$SV_5_2]-1][$N_chr3]->seq("$hap2[$SV_5_2]_$number2[$SV_5_2]", $N_sequence_line[3]=>$N_sequence_line[3]);
                my $N_seq4=$chr_fa[$hap2[$SV_3_2]-1][$N_chr4]->seq("$hap2[$SV_3_2]_$number2[$SV_3_2]", $N_sequence_line[4]=>$N_sequence_line[4]);



                if($N_seq1 eq "N" || $N_seq2 eq "N" || $N_seq3 eq "N" || $N_seq4 eq "N"){
                        return 1;
                }else{
                        return 0;
                }

}

sub Mappability_cal{
                my @start=();my @end=();my @cumul=();my @hap=();my @number=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};
                @hap=@{@chr_hap[@line[3]]};
                @number=@{@chr_number[@line[3]]};
                my $j;
                my $SV_5, $SV_3;
                my $SV_5_searched=0;
                my $SV_3_searched=0;
                my @N_sequence_line=@line;
		print $map_log "$N_sequence_line[1]\t$N_sequence_line[2]\n";
                for( $j = 0; $j <= $#cumul; $j++){
                        if($SV_5_searched == 0 && @N_sequence_line[1] <= @cumul[$j]){
                               $SV_5 = $j;
                               if($SV_5 != 0){
                                       @N_sequence_line[1]=@start[$SV_5]+@N_sequence_line[1]-@cumul[$SV_5-1]-1;
                                }else{
                                       @N_sequence_line[1]=@start[$SV_5]+@N_sequence_line[1]-1;
                                }
                                $SV_5_searched=1;
                        }
                        if($SV_3_searched == 0 && @N_sequence_line[2] <= @cumul[$j]){
                               $SV_3 = $j;
                               if($SV_3 != 0){
                                       @N_sequence_line[2]=@start[$SV_3]+@N_sequence_line[2]-@cumul[$SV_3-1]-1;
                                }else{
                                       @N_sequence_line[2]=@start[$SV_3]+@N_sequence_line[2]-1;
                                }
                                $SV_3_searched=1;
                        }

                }
                $N_sequence_line[1]=abs $N_sequence_line[1];
                $N_sequence_line[2]=abs $N_sequence_line[2];
		print $map_log "$SV_5_searched\t$SV_3_searched\t$cumul[$#cumul]\n";
                print $map_log "$hap[$SV_5]\t$number[$SV_5]\t$hap[$SV_3]\t$number[$SV_3]\t", "$N_sequence_line[1]\t$N_sequence_line[2]\n";

#                print $map_log "$N_sequence_line[1]\t$N_sequence_line[2]\n";
#		if($N_sequence_line[1] < 0){
#			open(my $negative_log, '>', 'negative_log') or die;
#			for(my $print_i=0;$print_i<=$#cumul;$print_i++){
#				print $negative_log "$start[$print_i]\t$end[$print_i]\t$cumul[$print_i]\t$hap[$print_i]\t$number[$print_i]\n";
#			}
#		}

		my $deltafile;
		open($deltafile, '<', "${G1000_dir}/${IND}.g$hap[$SV_5].fa.$number[$SV_5].delta") or die;
		my @delta_line;
		my $delta_line = <$deltafile>;
		my $delta=0;
		my $last_found=0;
		while($delta_line = <$deltafile>){
		        @delta_line = ();
		        chomp $delta_line;
		        @delta_line = split(/\t+/, $delta_line);
			if($N_sequence_line[1]<=$delta_line[0]){
				$N_sequence_line[1]=$N_sequence_line[1]+$delta;
				$last_found=1;
				last;
			}
			$delta = $delta_line[1];
		}close($deltafile);
		if($last_found==0){
                        $N_sequence_line[1]=$N_sequence_line[1]+$delta;
		}
	#	print "$hap[$SV_3]\t$number[$SV_3]\n";
                open($deltafile, '<', "${G1000_dir}/${IND}.g$hap[$SV_3].fa.$number[$SV_3].delta") or die;
                $delta_line = <$deltafile>;
                $delta=0;
                $last_found=0;
                while($delta_line = <$deltafile>){
                        @delta_line = ();
                        chomp $delta_line;
                        @delta_line = split(/\t+/, $delta_line);
                        if($N_sequence_line[2]<=$delta_line[0]){
                                $N_sequence_line[2]=$N_sequence_line[2]+$delta;
                                $last_found=1;
                                last;
                        }
                        $delta = $delta_line[1];
                }close($deltafile);
                if($last_found==0){
                        $N_sequence_line[2]=$N_sequence_line[2]+$delta;
                }
		my $bed1;
		my $bed2;
		open($bed1, '>', "map.bed1") or die;
		print $bed1 "chr$number[$SV_5]\t", $N_sequence_line[1]-500, "\t", $N_sequence_line[1]+500, "\n";
                #print "chr$number[$SV_5]\t", $N_sequence_line[1]-500, "\t", $N_sequence_line[1]+500, "\n";
		$map_score1=`bwtool extract bed map.bed1 /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.bigWig /dev/stdout | awk 'BEGIN{sum=0;e=0;}{n=split(\$5,f,","); for(i=1;i<=n;i++){if(f[i]=="NA"){sum=sum+1;e=e+1;}else{sum=sum+f[i];e=e+1;}}}END{print sum/e}'`;
		close($bed1);
                open($bed2, '>', "map.bed2") or die;
                print $bed2 "chr$number[$SV_3]\t", $N_sequence_line[2]-500, "\t", $N_sequence_line[2]+500, "\n";
                #print "chr$number[$SV_3]\t", $N_sequence_line[2]-500, "\t", $N_sequence_line[2]+500, "\n";
                $map_score2=`bwtool extract bed map.bed2 /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.bigWig /dev/stdout | awk 'BEGIN{sum=0;e=0;}{n=split(\$5,f,","); for(i=1;i<=n;i++){if(f[i]=="NA"){sum=sum+1;e=e+1;}else{sum=sum+f[i];e=e+1;}}}END{print sum/e}'`;
		close($bed2);
                print $map_log "$hap[$SV_5]\t$number[$SV_5]\t$hap[$SV_3]\t$number[$SV_3]\t", "$N_sequence_line[1]\t$N_sequence_line[2]\n";
		print $map_log "${map_score1}";
		print $map_log "${map_score2}";
		print $map_log ($map_score1+$map_score2)/2, "\n";
		return ($map_score1+$map_score2)/2;
}
sub low_map_select{
        my $random_line=1+int(rand($low_map_line));
        #print "$random_line\n";
        $random_line=`sed '$random_line!d' /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.wig.filtered`;
        my @random_line = split(/\t+/, $random_line);
        my $random_hap=1+int(rand(2));
        my $deltafile;
        my $delta_line;
        my @delta_line;
        open($deltafile, '<', "${G1000_dir}/${IND}.g$random_hap.fa.$random_line[0].delta") or die;
        $delta_line = <$deltafile>;
        my $delta=0;
        my $last_found=0;
        while($delta_line = <$deltafile>){
                @delta_line = ();
                chomp $delta_line;
                @delta_line = split(/\t+/, $delta_line);
                if($random_line[1]<=$delta_line[0]+$delta_line[1]){
                        $random_line[1]=$random_line[1]-$delta;
                        $last_found=1;
                        last;
                }
                $delta = $delta_line[1];
        }close($deltafile);
        if($last_found==0){
                $random_line[1]=$random_line[1]-$delta;
        }
        return ($random_hap, $random_line[0], $random_line[1]);
}
sub old_fa_coor_to_new_fa_coor{
	my ($old_fa_hap, $old_fa_chr, $old_fa_coor, $old_SV_size)=(@_);
#	print "$old_fa_hap\t$old_fa_chr\t$old_fa_coor\n";
	my @where_is;
	for(my $new_fa_i=0; $new_fa_i<=$#chr_start; $new_fa_i++){
		for(my $new_fa_j=0; $new_fa_j<=$#{@chr_cumul[$new_fa_i]}; $new_fa_j++){
			if($old_fa_hap == $chr_hap[$new_fa_i][$new_fa_j] && $old_fa_chr == $chr_number[$new_fa_i][$new_fa_j] 
				&& $old_fa_coor >= (abs $chr_start[$new_fa_i][$new_fa_j] >= abs $chr_end[$new_fa_i][$new_fa_j] ? abs $chr_end[$new_fa_i][$new_fa_j] : abs $chr_start[$new_fa_i][$new_fa_j]) && $old_fa_coor <= (abs $chr_start[$new_fa_i][$new_fa_j] >= abs $chr_end[$new_fa_i][$new_fa_j] ? abs $chr_start[$new_fa_i][$new_fa_j] : abs $chr_end[$new_fa_i][$new_fa_j])){
				$where_is[$#where_is+1][0]=$new_fa_i;
				if($chr_end[$new_fa_i][$new_fa_j] > 0 ){
                                	$where_is[$#where_is][1]=$chr_cumul[$new_fa_i][$new_fa_j]-($chr_end[$new_fa_i][$new_fa_j]-$old_fa_coor);
				}else{
                                        $where_is[$#where_is][1]=$chr_cumul[$new_fa_i][$new_fa_j]-($old_fa_coor - (abs $chr_end[$new_fa_i][$new_fa_j]));
				}
			}
		}
	}
        my $where_is_random_i=int(rand($#where_is+1));

	if($#where_is<0 || $chr_cumul[$where_is[$where_is_random_i][0]][$#{@chr_cumul[$where_is[$where_is_random_i][0]]}]<=$old_SV_size+$where_is[$where_is_random_i][1]){
                ####################### not exits ==> simply seleect first element"#
                my $low_map_select_random_coor;
		my $low_map_select_random_chr_index;
                do{
			#print "here?";
			$low_map_select_random_chr_index=int(rand($#chr_start+1));
			$low_map_select_random_coor=1+int(rand($chr_cumul[$low_map_select_random_chr_index][$#{@chr_cumul[$low_map_select_random_chr_index]}]-$offset_to_3_end-$old_SV_size));
		}while($low_map_select_random_coor<=0);
		return ($low_map_select_random_chr_index, $low_map_select_random_coor);
	}else{
		return ($where_is[$where_is_random_i][0], $where_is[$where_is_random_i][1]);
	}
}
