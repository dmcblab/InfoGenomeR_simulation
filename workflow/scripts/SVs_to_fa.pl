use POSIX;
use Bio::DB::Fasta;
use feature qw( switch );
no warnings qw( experimental::smartmatch );
open(my $map_log, '>', "map_log") or die;
my $lib="InfoGenomeR_simulation/simulation";
my $G1000_dir="indv_select_high_coverage";
my $low_map_line=0;
#my $low_map_line=`cat /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.wig.filtered | wc -l`;

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
my @delta_cumul;
my @chr_hap;
my @chr_number;
my $chr_index;
my $deletion_bridge=1000;
my $IND=$ARGV[0];
my $mappability_bias=0.2;
my $max_map_score=0.5;
open(my $log, '>', 'logs') or die;
open(my $germ_results, '>', 'germ_results') or die;
my $number_of_germ=0;
my $number_of_somatic=1;
my @chr_fa;
my $germline_chrs=45;

for(my $chrs=0; $chrs<=45; $chrs++){
	$delta_cumul[$chrs]=0;
}

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



my %SV_hash;
%SV_hash = (0=>"dup", 1=>"del", 2=>"simple_inv");

for(my $chr_index = 0; $chr_index <=$#chr_start; $chr_index++){
	print $germ_results "\n\n$chr_index\n";
	open(my $germ_SVs_files, '>', "germline_SVs.results.${chr_index}") or die "Could not open file 'SVs' $!";
       for($i = 0; $i <= $#{@chr_start[$chr_index]}; $i++){
	      print $germ_SVs_files "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
	      print $germ_results "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
	}
}


open($SVs, '<', 'SVs');


my @line_sampleCoord;
my @line;
while(my $line = <$SVs>){
	chomp $line;
	next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
	@line = split(/\t+/, $line);
	@line_sampleCoord=@line; ## origin
	diploid_liftover();

	if (@line[0] eq "dup"){
		adjust_delta_from_cumul();
		duplication();
	}
        if (@line[0] eq "del"){
		adjust_delta_from_cumul();
		deletion("del");
	}
	if (@line[0] eq "aoh"){
		adjust_delta_from_cumul();
		deletion("aoh");
	}
	### for not messing up total chromosome buckets
	## loss first from last chr
        if (@line[0] eq "chr_loss"){
                chromosome_deletion(@line[3]);
        }
        if (@line[0] eq "chr_gain"){
               chromosome_amplification(@line[3]);
        }

	if(@line[0] eq "chr_aoh"){
		chromosome_aoh(@line[3]);
	}
	if(@line[0] eq "arm_gain"){
		arm_amplification(@line[3]);
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

`perl $lib\/results_to_fa.pl $#chr_start $IND $G1000_dir`;
`perl $lib\/germ_results_to_fa.pl $germline_chrs $IND $G1000_dir`;

sub diploid_liftover {
	my $delta_f;
	my $hap_idx;
	my $chr_idx;
        if(@line[3]<23){
                $hap_idx=1;
                $chr_idx=@line[3]+1;
        }else{
                $hap_idx=2;
                $chr_idx=@line[3]-22;
        }
        my $chr_name;
        if($chr_idx==23){
                $chr_name="X";
        }else{
                $chr_name=$chr_idx;
        }


	for(my $coord_idx=1; $coord_idx<=2; $coord_idx++){ # 1 for start 2 for end
		open($delta_f, '<', "$G1000_dir/${IND}.g${hap_idx}.fa.${chr_name}.delta") or die "Could not open file 'delta' $!";
		my $d_line;
		my $found=0;

		while($d_line = <$delta_f>){
			chomp $d_line;
			next if $d_line =~ /new/;
			@d_line = split(/\t+/, $d_line);
			if(@d_line[0]>=@line[$coord_idx]){
				my $liftover_coord=@d_line[1];
				@line[$coord_idx]=@line[$coord_idx]-$liftover_coord;
				$found=1;
				last;
			}
		}
		if($found!=1){
			my $liftover_coord=@d_line[1];
			@line[$coord_idx]=@line[$coord_idx]-$liftover_coord;
		}
		close($delta_f);
	}
}

sub adjust_delta_from_cumul {
	my $sv_size=@line[2]-@line[1]+1;

	@line[1]+=$delta_cumul[@line[3]];
	@line[2]+=$delta_cumul[@line[3]];
        if (@line[0] eq "dup"){
		$delta_cumul[@line[3]]+=$sv_size;
        }
        if (@line[0] eq "del"){
		$delta_cumul[@line[3]]-=$sv_size;
        }
	if(@line[0] eq "aoh"){
		@line_swp=@line;
	        @line=@line_sampleCoord;
		if(@line[3]<23){
			@line[3]=@line[3]+23;
		}else{
			@line[3]=@line[3]-23;
		}
		diploid_liftover();
		my $new_sv_size=@line[2]-@line[1]+1;
		@line=@line_swp;
		$delta_cumul[@line[3]]+=$new_sv_size-$sv_size;
	}
}


sub duplication {
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @start=();my @end=();my @cumul=();
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

sub deletion {
		my $isAOH=0;
		if($_[0] eq "aoh"){
			$isAOH=1;
		}
                my @segment_start=();my @segment_end=();my @start=();my @end=();my @cumul=();
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

		if($isAOH){
			# a segment from the other haplotype
                        # hap swp
			my @line_swp=@line;
			@line=@line_sampleCoord;

			my $hap_swp=1;
                        if(@line[3]<23){
				@line[3]=@line[3]+23;
				$hap_swp=2;
                        }else{
                        	@line[3]=@line[3]-23;
                        }
                        diploid_liftover();	

			$seg_idx=$split_5+1;
			@segment_start[$seg_idx .. $segment_length]=@segment_start[$split_5 .. $segment_length+1];
			@segment_end[$seg_idx .. $segment_length]=@segment_end[$split_5 .. $segment_length+1];
			@segment_hap[$seg_idx .. $segment_length]=@segment_hap[$split_5 .. $segment_length+1];
			@segment_number[$seg_idx .. $segment_length]=@segment_number[$split_5 .. $segment_length+1];

			@segment_start[$split_5]=@line[1];
			@segment_end[$split_5]=@line[2];
			@segment_hap[$split_5]=$hap_swp;
			@segment_number[$split_5]=$chr_number[@line_swp[3]][0];
			
			@line=@line_swp;
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

sub chromosome_aoh{
	my $hap_idx;
	if($_[0] > 22){
		$hap_idx=$_[0]-23;
	}else{
                $hap_idx=$_[0]+23;

	}
@{@chr_start[$_[0]]}=@{@chr_start[$hap_idx]};
@{@chr_end[$_[0]]}=@{@chr_end[$hap_idx]};
@{@chr_cumul[$_[0]]}=@{@chr_cumul[$hap_idx]};
@{@chr_hap[$_[0]]}=@{@chr_hap[$hap_idx]};
@{@chr_number[$_[0]]}=@{@chr_number[$hap_idx]};
}

sub arm_amplification{

if(@line_sampleCoord[1] == 1){
	@{$chr_start[$#chr_start+1]}=1;
}else{
	@{$chr_start[$#chr_start+1]}=@line[1];
}

@{$chr_end[$#chr_end+1]}=@line[2];
@{$chr_cumul[$#chr_cumul+1]}=@line[2]-@line[1]+1;
@{$chr_hap[$#chr_hap+1]}=@{@chr_hap[@line[3]]};
@{$chr_number[$#chr_number+1]}=@{@chr_number[@line[3]]};

}
