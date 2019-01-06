#!/usr/local/bin/perl -w
	
use strict;


use strict;
print "module load gridengine\n";
my $i = 1021;
while ($i <= 1607){
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | awk '{print \$2}' > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | gawk '{print \$3}' | sed 's/[^0-9]*//g' >> /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`awk '{ sum += \$1 } END { print sum }' /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt`;
open (MYFILE, '</ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt') or die "$!";
my $line = <MYFILE>;

if (not defined $line) {
         $line = 0
}
if ($line < 190) {
        $i = $i + 1;

	my $rfile="Random".$i."_99.r";	
	my $jobfile="job".$i."_Randomization99.sh";
	open(OUT,">$jobfile");
	print OUT "module load R\n";
	print OUT "R --vanilla < $rfile > out_Randomization99$i.txt\n";
	close OUT;
	`qsub -S /bin/bash -pe serial 8 -l mem_free=8G,h_vmem=1G $jobfile`;


}
# print $line;
print "Job number is $i\n";
}

