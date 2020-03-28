#!/usr/bin/perl -w

##
# This script converts a min code formatted graph to DMTL format.
# Eg:
# 0 1 1 0 1:1 2 1 0 4:2 0 4 0 1:2 3 4 0 6:3 4 6 0 2:3 5 6 0 4:3 6 6 0 6:6 0 6 0 1
# to
# t # 0
# v 0 1
# ....
##

my $in_f=$ARGV[0];  # File containing min dfs codes for a bunch of graphs
open IN, "$in_f";

my $curr_graph=0;

while(<IN>) {

  my %nodes=();
  my @edges=();

  chomp;
  my @prts=split(':', $_);

  # Each edge in the graph.
  foreach my $e (@prts) {

    my @p=split(' ', $e);  # <src_id> <dest_id> <src_lbl> <edge_lbl> <dest_lbl>
    if(!exists $nodes{$p[0]}) {
      $nodes{$p[0]}=$p[2];
    }
    
    if(!exists $nodes{$p[1]}) {
      $nodes{$p[1]}=$p[4];
    }
   
    push(@edges, "e $p[0] $p[1] $p[3]"); 
  }

  # Print
  print "t # $curr_graph\n";
 
  # nodes 
  @sorted_nodes = sort {$a <=>$b} keys %nodes;
  foreach my $v (@sorted_nodes) {
    print "v $v $nodes{$v}\n";
  }

  # Edges 
  foreach my $e (@edges) {
    print "$e\n";
  }

  $curr_graph++;
}

close(IN);
