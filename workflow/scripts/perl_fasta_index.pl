#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Fasta;
use File::Basename; 

my $humandb=$ARGV[0];
my $ref_version=$ARGV[1];

sub force_index {
    my ($file) = @_;
    my $dirname = dirname($file);
    my $basename = basename($file);
    my @index_files = glob("$dirname/$basename.index");

    foreach my $index_file (@index_files) {
        if (-e $index_file) {
            unlink $index_file or warn "Failed to delete $index_file: $!";
        }
    }
}

my @fasta_files = (
    "${humandb}/mobileElements/ALU.fa",
    "${humandb}/mobileElements/LINE1.fa",
    "${humandb}/mobileElements/SVA.fa",
    "${humandb}/${ref_version}/ref/${ref_version}.fa",
);

foreach my $file (@fasta_files) {
    print "forced indexing $file\n";
    force_index($file);
}
