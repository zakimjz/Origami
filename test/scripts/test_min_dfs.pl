#!/usr/bin/perl -w

my $dmtl_pats_f=$ARGV[0];
my $rand_pats_f=$ARGV[1];

my %dmtl_hash=();

# Read the DMTL pats.
open D, "$dmtl_pats_f";
while(<D>) {
  chomp;

  my $pat = $1 if $_ =~ /Min DFS code= (.*)$/;

  if(!exists $dmtl_hash{$pat}) {
    $dmtl_hash{$pat} = 1;
  }
}
close(D);

# Read the Randomly generated pats.
open R, "$rand_pats_f";
while(<R>) {
  chomp;

  my $pat = $1 if $_ =~ /Min DFS code= (.*)$/;

  if(!exists $dmtl_hash{$pat}) {
    print "Pattern $pat not in DMTL freq patterns\n";
  }
}
close(R);

