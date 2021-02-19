my @start;
my @end;
@start[0]=1;
@end[0]=10000;
@start[1]=10001;
@end[1]=20000;
@start[2]=20001;
@end[2]=30000;
@start[3]=30001;
@end[3]=40000;

@cumul[0]=@end[0]-@start[0]+1;
for($i = 1; $i <=$#end; $i++){
        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
}


open($SVs, '<', 'SVs');

while(my $line = <$SVs>){
        chomp $line;
        next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        my @line = split(/\t+/, $line);

        my @dup_start=();
	my @dup_end=();
	my @segment_start=();
	my @segment_end=();
	my @inv_start=();
	my @inv_end=();
	if (@line[0] eq "dup"){
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


	}

        if (@line[0] eq "del"){
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

	}
        if (@line[0] eq "simple_inv"){
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


	}
        if (@line[0] eq "foldback_3"){
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



	}


@cumul[0]=@end[0]-@start[0]+1;
for($i = 1; $i <=$#end; $i++){
        @cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
}

}


close($SVs);
@cumul=();
@cumul[0]=@end[0]-@start[0]+1;

for($i = 1; $i <=$#end; $i++){
	@cumul[$i]=@cumul[$i-1]+@end[$i]-@start[$i]+1;
}
for($i = 0; $i <= $#start; $i++){
      print "@start[$i]\t@end[$i]\t@cumul[$i]\n";

}

for($i = 0; $i <= $#segment_start; $i++){

#  print "$segment_start[$i]\t$segment_end[$i]\n";
}
