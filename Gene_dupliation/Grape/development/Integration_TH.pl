
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 11; $i++){
	
	my $jobfile = "Integ".$i.".r";
     open(OUT,">$jobfile");
     print OUT "


load(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData\")
load(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score99/A.score.RData\")
load(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData\")

i = ".$i."
## Search orthologous modules for a pair of paralogous modules
RR = table(G.score[(G.score[, \"orthologous\"] == Result[i , \"Gene\"] |
              	G.score[, \"orthologous\"] == Result[i , \"Paralogous\"]), \"Gene\"])
RR = as.numeric(names(RR[RR > 1]))


save(RR, file = \"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/".$i.".RData\")

"
	;
	close OUT;
}





