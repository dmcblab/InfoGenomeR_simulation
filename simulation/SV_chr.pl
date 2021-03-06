my @chr_start;
my @chr_end;
my @chr_cumul;
my $chr_index;
my $deletion_bridge=1000;
for(my $chr_index = 0; $chr_index <=3; $chr_index++){
	$chr_start[$chr_index][0]=1;
	$chr_end[$chr_index][0]=100000000;
#	$chr_start[$chr_index][1]=10001;
#	$chr_end[$chr_index][1]=20000;
#	$chr_start[$chr_index][2]=20001;
#	$chr_end[$chr_index][2]=30000;
#	$chr_start[$chr_index][3]=30001;
#	$chr_end[$chr_index][3]=40000;
}

for(my $chr_index = 0; $chr_index <=3; $chr_index++){

	$chr_cumul[$chr_index][0]=$chr_end[$chr_index][0]-$chr_start[$chr_index][0]+1;
	for($i = 1; $i <=$#{@chr_end[$chr_index]}; $i++){
	        @{@chr_cumul[$chr_index]}[$i]=@{@chr_cumul[$chr_index]}[$i-1]+@{@chr_end[$chr_index]}[$i]-@{@chr_start[$chr_index]}[$i]+1;
	}
}

open($SVs, '<', 'SVs');

while(my $line = <$SVs>){
        chomp $line;
        next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        my @line = split(/\t+/, $line);


	if (@line[0] eq "dup"){
 	        my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
        	@start=@{@chr_start[@line[3]]};
	        @end=@{@chr_end[@line[3]]};
	        @cumul=@{@chr_cumul[@line[3]]};

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

		if($split_5==1){
			@segment_start[0]=@start[$SV_5];
			@segment_end[0]=@line[1]-1;
		}
		if($split_3==1){
			@segment_start[0+2*($SV_3-$SV_5)+$split_5+2]=@line[2]+1;
			@segment_end[0+2*($SV_3-$SV_5)+$split_5+2]=@end[$SV_3];
		}
		@segment_start[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_start[0 .. ($SV_3-$SV_5)];
		@segment_start[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_start[0 .. ($SV_3-$SV_5)];
		@segment_end[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@dup_end[0 .. ($SV_3-$SV_5)];
		@segment_end[0+$split_5+$SV_3-$SV_5+1 .. 0+2*($SV_3-$SV_5)+$split_5+1]=@dup_end[0 .. ($SV_3-$SV_5)];

		$segment_length=$#segment_start+1;

		if($last>$SV_3){
			@start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
			@end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
		}

		@start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
		@end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];


		@cumul[0]=@end[0]-@start[0]+1;
		for($i = 1; $i <=$#end; $i++){
			@cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
		}

		@{@chr_start[@line[3]]}=@start;
		@{@chr_end[@line[3]]}=@end;
 	        @{@chr_cumul[@line[3]]}=@cumul;



	}

        if (@line[0] eq "del"){
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};

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
		}
                if($split_3==1){
                        @segment_start[$split_5+$split_3-1]=@line[2]+1;
                        @segment_end[$split_5+$split_3-1]=@end[$SV_3];
                }
		$segment_length=$#segment_start+1;
                if($last>$SV_3){
                        @start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
                        @end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
                }
		if($SV_5+$segment_length + $last - $SV_3 -1 < $last){
			delete @start[$SV_5+$segment_length + $last - $SV_3 .. $last];
			delete @end[$SV_5+$segment_length + $last - $SV_3 .. $last];
		}
                @start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
                @end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];
                @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;


	}
        if (@line[0] eq "simple_inv"){
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};

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

                if($split_5==1){
                        @segment_start[0]=@start[$SV_5];
                        @segment_end[0]=@line[1]-1;
                }
                if($split_3==1){
                        @segment_start[0+$SV_3-$SV_5+$split_5+$split_3]=@line[2]+1;
                        @segment_end[0+$SV_3-$SV_5+$split_5+$split_3]=@end[$SV_3];
                }
		
		for(my $k=0; $k <= $SV_3-$SV_5; $k++){
			@inv_start[$k]=-@dup_end[$SV_3-$SV_5-$k];
			@inv_end[$k]=-@dup_start[$SV_3-$SV_5-$k];
		}
                @segment_start[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_start[0 .. ($SV_3-$SV_5)];
                @segment_end[0+$split_5 .. 0+$split_5+$SV_3-$SV_5]=@inv_end[0 .. ($SV_3-$SV_5)];
		
                $segment_length=$#segment_start+1;

                if($last>$SV_3){
                        @start[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @start[$SV_3+1 .. $last];
                        @end[$SV_5+$segment_length .. $SV_5+$segment_length + $last - $SV_3 -1 ] = @end[$SV_3+1 .. $last];
                }

		@start[$SV_5 .. $SV_5+$segment_length-1]=@segment_start[0 .. $segment_length-1];
                @end[$SV_5 .. $SV_5+$segment_length-1]=@segment_end[0 .. $segment_length-1];
                @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;


	}
        if (@line[0] eq "foldback_3"){
                my @dup_start=();my @dup_end=();my @segment_start=();my @segment_end=();my @inv_start=();my @inv_end=();my @start=();my @end=();my @cumul=();
                @start=@{@chr_start[@line[3]]};
                @end=@{@chr_end[@line[3]]};
                @cumul=@{@chr_cumul[@line[3]]};

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
               }
		@end[$SV_3]=$line[2];

		@segment_start[0 .. $SV_5] = @start[0 .. $SV_5];
                @segment_end[0 .. $SV_5] = @end[0 .. $SV_5];
		@segment_end[$SV_5] = @line[1];
 
               for(my $k=0; $k <= $SV_5; $k++){
                        @inv_start[$k]=-@segment_end[$SV_5-$k];
                        @inv_end[$k]=-@segment_start[$SV_5-$k];
                }

		@start[$SV_3+1 .. $SV_3+$#inv_start+1] = @inv_start[0 .. $#inv_start];
                @end[$SV_3+1 .. $SV_3+$#inv_end+1] = @inv_end[0 .. $#inv_end];
                @cumul[0]=@end[0]-@start[0]+1;
                for($i = 1; $i <=$#end; $i++){
                        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
                }

                @{@chr_start[@line[3]]}=@start;
                @{@chr_end[@line[3]]}=@end;
                @{@chr_cumul[@line[3]]}=@cumul;



	}
        if (@line[0] eq "balanced_translocation"){
		my @start1=(); my @end1=(); my @cumul1=(); my @start2=(); my @end2=(); my @cumul2=();
   	        @start1=@{@chr_start[@line[5]]};
	        @end1=@{@chr_end[@line[5]]};
	        @cumul1=@{@chr_cumul[@line[5]]};
	        @start2=@{@chr_start[@line[6]]};
	        @end2=@{@chr_end[@line[6]]};
	        @cumul2=@{@chr_cumul[@line[6]]};
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
		@segment1_start[0] = @line[4];

                @segment2_start[0 .. $#start1-$SV_3_1 ] = @start1[$SV_3_1 .. $#start1];
                @segment2_end[0 .. $#start1-$SV_3_1 ] = @end1[$SV_3_1 .. $#start1];
                @segment2_start[0] = @line[2];


		@end1[$SV_5_1]=@line[1];
		@start1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_start[0 .. $#segment1_start];
		@end1[$SV_5_1+1 .. $#segment1_start+$SV_5_1+1] = @segment1_end[0 .. $#segment1_end];
		if($last1>$#segment1_start+$SV_5_1+1){
                        delete @start1[$#segment1_start+$SV_5_1+2..$last1];
                        delete @end1[$#segment1_end+$SV_5_1+2..$last1];

		}		
                @end2[$SV_5_2]=@line[3];
                @start2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_start[0 .. $#segment2_start];
                @end2[$SV_5_2+1 .. $#segment2_start+$SV_5_2+1] = @segment2_end[0 .. $#segment2_end];

                if($last2>$#segment2_start+$SV_5_2+1){
                        delete @start2[$#segment2_start+$SV_5_2+2..$last2];
                        delete @end2[$#segment2_end+$SV_5_2+2..$last2];
                }

                @cumul1[0]=@end1[0]-@start1[0]+1;
                for($i = 1; $i <=$#end1; $i++){
                        @cumul1[$i]=@cumul1[$i-1]+@end1[$i]-@start1[$i]+1;
                }
                @{@chr_start[@line[5]]}=@start1;
                @{@chr_end[@line[5]]}=@end1;
                @{@chr_cumul[@line[5]]}=@cumul1;

                @cumul2[0]=@end2[0]-@start2[0]+1;
                for($i = 1; $i <=$#end2; $i++){
                        @cumul2[$i]=@cumul2[$i-1]+@end2[$i]-@start2[$i]+1;
                }

                @{@chr_start[@line[6]]}=@start2;
                @{@chr_end[@line[6]]}=@end2;
                @{@chr_cumul[@line[6]]}=@cumul2;

	}

}
close($SVs);

for(my $chr_index = 0; $chr_index <=3; $chr_index++){
	print "$chr_index\n";
	for($i = 0; $i <= $#{@chr_start[$chr_index]}; $i++){
	      print "$chr_start[$chr_index][$i]\t$chr_end[$chr_index][$i]\t$chr_cumul[$chr_index][$i]\n";
	}
}
