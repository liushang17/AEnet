use warnings;
use strict;

opendir I,$ARGV[0] or die "$!";
foreach(readdir I){
        chomp;
        next if($_ eq "\." or $_ eq "\.\.");
	my $file = "$ARGV[0]/$_";
	my $num = 0;
        open I1,$file or die "$!";
	my $clu1 = $1 if($_ =~ /(\d+)_\d+.xls/);
        my $clu2 = $1 if($_ =~ /\d+_(\d+).xls/);
        my $genename = $` if($_ =~ /\.\d+_\d+.xls/);

        foreach my $line(<I1>){
                chomp($line);
                next if($line =~ /logFC/);
	        my @word = split " ",$line;
                if(abs($word[1]) >=0 and $word[5] <= 0.05){
			print "$word[0]\t$word[1]\t$word[5]\t$clu1\t$clu2\t$genename\n";	
                }
        }
        close I1;
}
closedir I;
