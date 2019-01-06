#!/usr/local/bin/perl -w

use strict;


use strict;
print "module load gridengine\n";
my $i = 0;
while ($i <= 7179){
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | awk '{print \$2}' > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | gawk '{print \$3}' | sed 's/[^0-9]*//g' >> /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`awk '{ sum += \$1 } END { print sum }' /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt`;
open (MYFILE, '</ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt') or die "$!";
my $line = <MYFILE>;

if (not defined $line) {
	 $line = 0
}
if ($line < 150) {
	$i = $i + 1;

  my $rfile="SA".$i."_Permut.SD_V2_Ind.r";	
  my $jobfile="job".$i."_T_Permut_Ind.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Perm$i.txt\n";
  close OUT;
  `qsub -S /bin/bash -l hostname="cyclone1|blizzard2|blizzard1|tanith|crunch1|cyclone2|cyclone3|cyclone4|cyclone5|cyclone6|cyclone7|cyclone8|hurricane1|hurricane2" -l h_vmem=1.5G $jobfile`;

}
# print $line;
print "Job number is $i\n";
}
