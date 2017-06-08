#!/user/bin/perl

$seq = "HapMap_6_chr22";
$geno = "HGDP_238_chr22";
$newseq = $seq."_subset";
$newgeno = $geno."_subset";

my @idx1=(0,1);
my @idx2=(0,1);

open(SITE0, $seq.".site") || die "can't open the sitefile of the input data: $!";
open(SITE1, ">".$newseq.".site") || die "can't open the output sitefile1: $!";
open(SITE2, ">".$newgeno.".site") || die "can't open the output sitefile2: $!";

$line = <SITE0>;
print SITE1 $line;
print SITE2 $line;
$i = 2;
while($line = <SITE0>){
	if(rand()<0.8){
		print SITE1 $line;
		push(@idx1, $i); 
	}
	if(rand()<0.8){
		print SITE2 $line;
		push(@idx2, $i); 
	}	
	$i = $i+1;
}
close SITE0;
close SITE1;
close SITE2;


open(SEQ1, $seq.".seq") || die "can't open the seqfile of the input data: $!";
open (SEQ2, ">".$newseq.".seq") || die "can't open the output seqfile: $!";
while($line = <SEQ1>){
	chomp($line);
	@tmp = split("\t", $line);
	@seq = @tmp[@idx1];
	print SEQ2 join("\t", @seq)."\n";
}
close SEQ2;
close SEQ1;


open(GENO1, $geno.".geno") || die "can't open the genofile of the input data: $!";
open (GENO2, ">".$newgeno.".geno") || die "can't open the output genofile: $!";
while($line = <GENO1>){
	chomp($line);
	@tmp = split("\t", $line);
	@seq = @tmp[@idx2];
	print GENO2 join("\t", @seq)."\n";
}
close GENO2;
close GENO1;
