my $intra_trans_n_min=4;
my $intra_trans_n_max=50; ## should be even
my $intra_trans_size_MIN=1e6;
#my $intra_trans_min = 1e6;
#
use POSIX;
use Bio::DB::Fasta;
use Switch;
open(my $map_log, '>', "map_log") or die;
my $G1000_dir="/DASstorage6/leeyh/1000G_SNP_indel_SV_integrated";
my $low_map_line=`cat /DASstorage6/leeyh/BG_job/2017.10.9.simulation_scripts_NGS/wgEncodeCrgMapabilityAlign100mer.wig.filtered | wc -l`;
my $PLOIDY=3;
my $intra_trans_n=20; ## should be even
my $offset_to_3_end=100;
my $number_of_WGD=0;
my @chr_start;
my @chr_end;
my @chr_cumul;
my @chr_hap;
my @chr_number;
my $chr_index;
my $deletion_bridge=1000;
my $IND=1762;
my $mappability_bias=0.2;
my $max_map_score=0.5;
open(my $log, '>', 'logs') or die;
open(my $germ_results, '>', 'germ_results') or die;
my $number_of_germ=3000;
my $number_of_somatic=200;
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


#       $chr_start[$chr_index][1]=10001;
#       #       $chr_end[$chr_index][1]=20000;
#       #       $chr_start[$chr_index][2]=20001;
#       #       $chr_end[$chr_index][2]=30000;
#       #       $chr_start[$chr_index][3]=30001;
#       #       $chr_end[$chr_index][3]=40000;
#       }
}

for(my $chr_index = 0; $chr_index <= $#chr_start; $chr_index++){

        $chr_cumul[$chr_index][0]=$chr_end[$chr_index][0]-$chr_start[$chr_index][0]+1;
        for($i = 1; $i <=$#{@chr_end[$chr_index]}; $i++){
                @{@chr_cumul[$chr_index]}[$i]=@{@chr_cumul[$chr_index]}[$i-1]+@{@chr_end[$chr_index]}[$i]-@{@chr_start[$chr_index]}[$i]+1;
        }
}

my @line=();
                                       @line[5]=int(rand($#chr_start+1));
                                       @line[6]=int(rand($#chr_start+1));
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[1]=1+int(rand($chr_cumul[@line[5]][$#{@chr_cumul[@line[5]]}]-$offset_to_3_end-$SV_size));
                                       @line[2]=@line[1]+$SV_size-1;
                                       my $SV_size = 1000+int(rand(10000));
                                       @line[3]=1+int(rand($chr_cumul[@line[6]][$#{@chr_cumul[@line[6]]}]-$offset_to_3_end-$SV_size));
                                       @line[4]=@line[1]+$SV_size-1;
               balanced_translocation();
                if(rand(1)<0.5){
                        chromosome_deletion(@line[5]);
                        $line[3]=$line[6];
			print "deleted $line[5]";
                }else{
                        chromosome_deletion(@line[6]);
                        $line[3]=$line[5];
			print "deleted $line[6]";
                }
		print "intra $line[3]\n";
                intra_trans();




for(my $chr_index = 0; $chr_index <=$#chr_start; $chr_index++){
        print "\n\n$chr_index\n";

        for($i = 0; $i <= $#{@chr_start[$chr_index]}; $i++){
              print "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\t$chr_hap[$chr_index][$i]\t$chr_number[$chr_index][$i]\n";
        }
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
			print "$intra_trans_size_min\n";
                }while($intra_trans_size_min < $intra_trans_size_MIN || $trans_while_iter < 100);

                if($intra_trans_size_min < $intra_trans_size_MIN){
                        return;
                }


                my @new_start=();my @new_end=(); my @new_hap=(); my @new_number=(); my @new_cumul=();
                for(my $trans_i = 0; $trans_i < $intra_trans_n;$trans_i++){
                        @dup_start=();@dup_end=();@dup_hap=();@dup_number=(); @inv_start=();@inv_end=();@inv_hap=();@inv_number=();
                        do{
				print"why loop";
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


