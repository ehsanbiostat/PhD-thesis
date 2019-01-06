#!/usr/local/bin/perl -w

use strict;

open (OUT, "> GO_separated_combination.txt");
print OUT "a=list()\n";
for(my $i=1;$i<148;$i++){
print OUT "load(\"AG_GO_separated$i.RData\")\n";
print OUT "q=data.frame(w[1],w[2],w[3])\n";
print OUT "a[$i]=list(q)\n";
}
print OUT "sep=a\n";
print OUT "save(sep, file=\"AG_GO_separated.RData\")\n";
close OUT;

