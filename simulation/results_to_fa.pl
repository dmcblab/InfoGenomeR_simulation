use Bio::DB::Fasta;
use warnings;
my $final_index_of_chr=$ARGV[0];
my $IND=$ARGV[1];
my $G1000_dir=$ARGV[2];

open(my $fh, '>', 'simulated.fa') or die "Could not open file 'SVs' $!";

my @chr_fa;
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

for(my $result_i =0; $result_i <=$final_index_of_chr; $result_i++){

	my $new_fa="";
	open(my $SVs, '<', "SVs.results.${result_i}");
	while(my $line = <$SVs>){
       		chomp $line;
       		next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        	my @line = split(/\t+/, $line);
		my $partial_chr_index;
		if($line[4] eq "X"){
			$partial_chr_index=22;
		}else{
			$partial_chr_index=$line[4]-1;
		}
		if($line[0]>0){
			$partial_fa=$chr_fa[$line[3]-1][$partial_chr_index]->seq("$line[3]_$line[4]",$line[0]=>$line[1]);
		}else{
                        $partial_fa=$chr_fa[$line[3]-1][$partial_chr_index]->seq("$line[3]_$line[4]",-$line[1]=>-$line[0]);
			$partial_fa =~ tr/ACGTacgt/TGCAtgca/;
			$partial_fa = scalar reverse $partial_fa;

		}
		$new_fa="$new_fa$partial_fa"

	}
	close($SVs);

        print $fh ">${result_i}\n";
        while( my $chunk = substr($new_fa, 0, 80, "")){
                  print $fh "$chunk\n";
        }
}


