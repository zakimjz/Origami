#!/usr/bin/perl -w

my $sup=1;

my $cmd="./graph_test -i ../data/GRAPH_large.dat -s $sup > out";

while(1) {
  my @lines = `$cmd`;

  # my $full_file = join('', @lines);
  open IN, "out";
  undef $/;
  my $full_file = <IN>;
  close(IN);

  $full_file =~ m/Support: (\d+)\n\nVAT for the max graph:.*?Tid = (\d+).*?Pattern size = (\d+)$/gs;
  my $sup = $1;
  my $tid = $2;
  my $sz = $3;
  my $num_edges;

  if($sup > $sup) {  # If support > $sup, not good.
    print "Error!!\n";
    print "Support = $sup\n";
    exit;
  } else {
  
    my $cmd_2 ="../data/draw_graph.pl ../data/GRAPH_large.dat $tid";
    my @lines_2 = `$cmd_2`;
    my $full_file_2 = join('', @lines_2);

    $full_file_2 =~ m/num edges = (\d+)$/gs;
    $num_edges = $1;
  }
  
  if($sz != $num_edges) {
   print "Error!!\n";
   print "sz = $sz, num_edges = $num_edges\n";
   exit;
  }
}  
