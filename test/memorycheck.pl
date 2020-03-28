#!/usr/bin/perl

if (!defined($kidpid = fork())) {
  # fork returned undef, so failed
  die "Cannot fork: $!";
}
elsif ($kidpid == 0) {
  # fork returned 0, so this branch is child
  system("top | grep graph_test > _mem_profile");
  # if exec fails, fall through to the next statement
}
else {
  # fork returned non-zero
  # this branch is patent
  system("./graph_test -i ../data/GRAPH_large.dat -s 50 -tm 5");
  # waitpid($kidpid, 0);
  @processes = system ("ps -u alhasan | grep top | awk '{print \$2}'");
  print @processes;
  # kill $kidpid;
}
