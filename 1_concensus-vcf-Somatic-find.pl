my $coord = "Somatic_Variant/concensus-list1";
my $vcf   = "Somatic_Variant/concensus-vcf";
my $out   = "Somatic_Variant/Final-Somatic-ALL-CONCENSUS.vcf";

my %coordhash;
open(coord, $coord);
while (my $line = <coord>) {
  chomp $line; 
  my $chr = (split(/\s+/,$line))[0];
  my $pos = (split(/\s+/,$line))[1];
  $coordhash{$chr}{$pos} = 1;
}
close coord;

open(out, ">$out");
open(vcf, $vcf);
while (my $line = <vcf>) {
  chomp $line;
  if ($line =~ /^chr/) {
    my $chr = (split(/\s+/,$line))[0];
    my $pos = (split(/\s+/,$line))[1];
    if ($coordhash{$chr}{$pos} == 1) {
      print out ("$line\n");
      $coordhash{$chr}{$pos} += 1;
    }
  }
}
close vcf;
close out; 
undef %coordhash;