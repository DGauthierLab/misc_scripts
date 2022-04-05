#!/sw/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
#use IPC::System::Simple;
#use autodie qw(:all);


my ($pair1, $pair2, $shrimp_genomefile,$bedtools_genomefile);

GetOptions(
           '-1:s' => \$pair1,
           '-2:s' => \$pair2,
           'r|reference:s'  => \$shrimp_genomefile,
           'b|bedtools:s' => \$bedtools_genomefile,
          );

print ("Enter bacterium prefix:\n");
my $prefix = <STDIN>;
chomp $prefix;

print ("Reference strain?\n");
my $refstrain = <STDIN>;
chomp $refstrain;

print ("Which chromosome?\n");
my $chromo = <STDIN>;
chomp $chromo;

## alternate entry method
#print ("I expect to find the following files:\n");
#print ("~/Volumes/500GB/marinum_paired_files/$prefix/$prefix.trim.1.fastq\n");
#print ("~/Volumes/500GB/marinum_paired_files/$prefix/$prefix.trim.2.fastq\n");
#
#
#print ("Where is the reference genome?\n");
#my $shrimp_genomefile = <STDIN>;
#chomp $shrimp_genomefile;
#
#print ("Where is the genomefile for Bedtools?\n");
#my $bedtools_genomefile = <STDIN>;
#chomp $bedtools_genomefile;

#my $pair1 = "~/Volumes/500GB/marinum_paired_files/$prefix/$prefix.trim.1.fastq";
#my $pair2 = "~/Volumes/500GB/marinum_paired_files/$prefix/$prefix.trim.2.fastq";

system ("mkdir ${prefix}_ref_${refstrain}");
system ("gmapper -N 20 -p opp-in --qv-offset 33 --progress 0 -1 $pair1 -2 $pair2 $shrimp_genomefile > ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sam"); 
print ("SHRiMP mapping completed!\n");
system("samtools view -b -S -o ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.bam ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sam");
print (".bam file created!\n");
system("samtools sort -o ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bam -T $prefix ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.bam");
print (".bam file sorted!\n");
system ("bedtools genomecov -d -ibam ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bam -g $bedtools_genomefile > ${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bed");
print ("genome coverage estimation complete!\n");

open(my $OUTPUT, ">", "${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bed.cat") or die $!;

open(my $INPUT,"<", "${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bed") or die $!;

while (my $row = <$INPUT>) {
    chomp $row;
    
    $row =~ s/.*\t(.*\t.*)/$chromo\t$1/;
    (my $a, my $b, my $c) = split("\t", $row);
    if ($c < 5) {$row =~ s/(.*\t.*)\t.*/$1\t0/}
           else {if ($c > 5){$row =~ s/(.*\t.*)\t.*/$1\t1/}
                          }
#print("$row\n");
print $OUTPUT "$row\n";
}
close $OUTPUT;
print ("binning completed!\n");
my $chr;
my $pos;
my $cov_type;  

my $chr_prev;
my $pos_start;
my $pos_end;
my $cov_type_prev=-1;  


 open (FILE, "${prefix}_ref_${refstrain}/${prefix}_ref_${refstrain}.map.sorted.bed.cat");
 open(OUTFILE, ">", "heatmaps/${prefix}_ref_${refstrain}.heatmap.circos") or die $!;
 while (<FILE>) {#o_while
  chomp;
  ($chr, $pos, $cov_type) = split("\t");
  if($cov_type == $cov_type_prev)   # 1st round -- evaluates false.  go to else.
                                    # 2nd round -- reads cov type, compares.  if same, bumps pos_end to current position
  {  #O if
      $pos_end = $pos; #pos_start stays the same as this loops. pos_end advances with every line
      $chr_prev = $chr; #chr_prev is current chromosome name
  }
      elsif($cov_type_prev == -1) {                     #1st round --evaluates true
            $pos_start=$pos_end=$pos;                    # $pos_start = 1, $pos_end=1, $pos= 1
            $chr_prev = $chr;                           # $chr_prev = mm1, $chr = mm1
            $cov_type_prev=$cov_type;                   # $cov_type_prev = 1 $cov_type=1
            }         
      elsif(($cov_type_prev != -1) && (($pos_start - $pos_end)) <= -750) {   # O do # 1st round -- FALSE, otw print the line 
            print OUTFILE "$chr_prev\t$pos_start\t$pos_end\t$cov_type_prev\n"; # 
            $pos_start=$pos_end=$pos; # next three lines reset variables before next line read.  only happens when coverage varies
            $chr_prev = $chr;
            $cov_type_prev=$cov_type;
            }# cov_type_prev sets to cov_type
      elsif(($cov_type_prev != -1) && (($pos_start - $pos_end) > -750)){
            $cov_type_prev=$cov_type;
            }# cov_type_prev sets to cov_type
      }
  
 print OUTFILE "$chr_prev\t$pos_start\t$pos\t$cov_type_prev\n";
 close (FILE);
 
 
