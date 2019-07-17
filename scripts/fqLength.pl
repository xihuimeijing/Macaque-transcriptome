#!perl
while ($l1=<>){
chomp $l1;
$l1=~s/^@//;
$l2=<>;
$l3=<>;
$l4=<>;
$b=length($l2);
print "$l1\t$b\n";
}
