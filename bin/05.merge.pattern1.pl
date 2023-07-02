use warnings;
use strict;

my %cell;
my $cellname = "cell";
open O,">$ARGV[1]" or die "$!";
opendir I,$ARGV[0] or die "$!";
foreach(readdir I){
	chomp;
	next if($_ eq "\." or $_ eq "\.\.");
	if($_ =~ /.pattern.v1.xls/){
	my $file2 = "$ARGV[0]/$_";
	if(-e $file2){
	open I2,$file2 or die "$!";
#	print "$file2\n";
	foreach my $line1(<I2>){
		chomp($line1);
		if($line1 =~ /cell/){
			my @word = split " ",$line1;
			shift (@word);
			my $linetmp = join "\t",@word;
			$cellname .= "\t$linetmp";
		}else{
			my @word = split " ",$line1;
			my $tmp = $word[0];
			shift (@word);
			my $linetmp = join "\t",@word;
			$cell{$tmp} .= "\t$linetmp";
		}
	}
	close I2;
	}
	}
}
closedir I;
close O;

open O1,">$ARGV[1]" or die "$!";
print O1 "$cellname\n";
foreach(keys %cell){
	print O1 "$_$cell{$_}\n";
}
close O1;
