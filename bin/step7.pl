#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
my %opts;

unless (getopts("s:f:p:h:m:a", \%opts))
{
        printhelp();
        print "Error: some options are not set properly!\n";
        exit;
}

# haplotypes 2 sample -------------

my $h2f = "";
if ($opts{"f"})
{
    $h2f = $opts{"f"};
}
$h2f =~ s/\[//ig;
$h2f =~ s/\]//ig;
$h2f =~ s/\s+//ig;
my (@h2f) = split(/(\w+?.\w+.txt)/, $h2f);

# haplotypes seq file -------------

my $hf = "";
if ($opts{"h"})
{
    $hf = $opts{"h"};
}
$hf =~ s/\[//ig;
$hf =~ s/\]//ig;
$hf =~ s/\s+//ig;
my (@hf) = split(/(\w+?.\w+.txt)/, $hf);

# samples -------------------------------
my $samples = "";
if ($opts{"s"})
{
    $samples = $opts{"s"};
}
$samples =~ s/\[//ig;
$samples =~ s/\]//ig;
$samples =~ s/\s+//ig;
my (@samples) = split ",", $samples;

# primers -------------------------------
my $primers = "";
if ($opts{"p"})
{
    $primers = $opts{"p"};
}
$primers =~ s/\[//ig;
$primers =~ s/\]//ig;
$primers =~ s/\s+//ig;
my (@genes) = split ",", $primers;

my $haplotype_reads_ratio = 0.2;
if ($opts{"a"}) 
{
	$haplotype_reads_ratio = $opts{"a"};
}

my $runmsa = "";
if (defined $opts{"m"}) 
{
	my @t = split ":", $opts{"m"};
	if ($t[0]=~/clustalw/i) 
	{
		$runmsa = "clustalw -INFILE=xxxxx -TYPE=DNA -OUTORDER=INPUT -OUTPUT=$t[1] -QUIET";
	}
	elsif ($t[0]=~/clustalo/i) 
	{
		my $my_wrap = "";
		if ($t[2])
		{
			$my_wrap = " --wrap=$t[2] ";
		}
		
		$runmsa = "clustalo --threads=8 -i xxxxx --seqtype=DNA --output-order=input-order --outfmt=$t[1] -o xxxxx.alignment $my_wrap --force";
	}
	else
	{
		print "Error: the parameter for -m is not right. It should by clustalw:FORMAT or clustalo:FORMAT:length.\n";
		exit;
	}
}

my %sample2gt; 
my %haplotypefrq;

print "Step 7. Identifying top 2 haplotypes from each sample.\n";
LOOPi7:foreach my $gene (@genes)
{
    my $haplotype2samplefile;
    my $haplotype2seqfile;

    foreach my $hf (@hf)
    {
        my $file1 = basename($hf);
        if ($file1 eq "$gene.haplotypes.txt") {
            $haplotype2seqfile = $hf;
            last;
        }
    }

    foreach my $h2f (@h2f)
    {
        my $file2 = basename($h2f);
        if ($file2 eq "$gene.haplotype2sample.txt") {
            $haplotype2samplefile = $h2f;
            last;
        }
    }

    unless (-e $haplotype2samplefile) 
    {
        next LOOPi7;
    }

    my $filesize = -s $haplotype2samplefile;
    if  ($filesize == 0) 
    {
        next LOOPi7;
    }

    ### read in the haplotype2sample file
    my %tbt = ();
    open (IN, "$haplotype2samplefile") || die "Cannot open $haplotype2samplefile \n";
    my $line = <IN>;
    chomp $line;
    my @thessamples = split "\t", $line;
    shift @thessamples;
    while (<IN>) 
    {
        chomp;
        my ($haplotype, @counts) = split "\t";
        for (my $i=0; $i<=$#counts; $i++) 
        {
            $tbt{$thessamples[$i]}{$haplotype} = $counts[$i];
        }
    }
    close IN;

    LOOPi7_1:foreach  my $sample(@samples) 
    {
        if (exists $tbt{$sample}) 
        {
            my %haplotype2counts = %{$tbt{$sample}};
            my @tophaplotypes = reverse sort {$haplotype2counts{$a} <=> $haplotype2counts{$b}} keys %haplotype2counts;
            my ($a1, $a2);
            if (@tophaplotypes>=2) 
            {
                ($a1, $a2) = @tophaplotypes[0,1];
            }
            elsif (@tophaplotypes==1) 
            {
                $a1 = $tophaplotypes[0];
                $a2 = "";
                $haplotype2counts{""} = 0;
            }
            else
            {
                next LOOPi7_1;
            }

            if ($haplotype2counts{$a1}<=2) 
            {
                $sample2gt{$gene}{$sample} = "./.:0";
            }
            elsif (($a2 eq "") || ($haplotype2counts{$a2}==0) )
            {
                $sample2gt{$gene}{$sample} = "$a1/$a1:$haplotype2counts{$a1}";
                $haplotypefrq{$gene}{$a1} += 2;
            }
            elsif ($haplotype2counts{$a2}/$haplotype2counts{$a1}<$haplotype_reads_ratio) 
            {
                $sample2gt{$gene}{$sample} = "$a1/$a1:$haplotype2counts{$a1}";
                $haplotypefrq{$gene}{$a1} += 2;
            }
            else
            {
                $sample2gt{$gene}{$sample} = "$a1/$a2:$haplotype2counts{$a1},$haplotype2counts{$a2}";
                $haplotypefrq{$gene}{$a1} += 1;
                $haplotypefrq{$gene}{$a2} += 1;
            }
            #print "$gene\t$sample\t$sample2gt{$gene}{$sample}\n";          
        }
    }
    my $fastafile = "$gene.haplotypes.fa";

    open OUT, ">$fastafile";
    open IN, "$haplotype2seqfile";
    my %haplotype2seq = ();
    while (<IN>) 
    {
        chomp;
        my ($a, $c, $s) = split "\t";
        if (exists $haplotypefrq{$gene}{$a}) 
        {
            print OUT ">$a Reads:$c Samples: $haplotypefrq{$gene}{$a}\n$s\n";
            $haplotype2seq{$a} = $s;
        }
    }
    close IN;
    close OUT;
    if ($runmsa=~/\w/) 
    {
        my $t= $runmsa;
        $t=~s/xxxxx/$fastafile/g;
        system ($t);
    }
}

open OUT, ">hap_genotype.txt";
print OUT "Locus\tHaplotypes\t";
print OUT (join "\t", @samples);
print OUT "\n";
my @existing_genes = ();
LOOPi7_2:foreach my $gene (@genes) 
{
    unless ($haplotypefrq{$gene}) 
    {
        next LOOPi7_2;
    }

    my %frq = %{$haplotypefrq{$gene}};
    my @haplotypes = reverse sort {$frq{$a} <=> $frq{$b}} keys %frq;
    my $total=0;
    foreach my $haplotype (@haplotypes) 
    {
        $total += $frq{$haplotype};
    }
    my $haplotypestr = "";
    if ($total==0) 
    {
        next LOOPi7_2;
    }
    push @existing_genes,  $gene;
    foreach my $haplotype (@haplotypes) 
    {
        my $frq =  sprintf("%.2f", $frq{$haplotype}/$total);
        $haplotypestr .= "$haplotype($frq);";
    }

    print OUT "$gene\t$haplotypestr";
    foreach my $sample (@samples) 
    {
        if (exists $sample2gt{$gene}{$sample}) 
        {
            print OUT "\t", $sample2gt{$gene}{$sample}; 
        }
        else
        {
            print OUT "\t", "./.:0"; 
        }
    }
    print OUT "\n";
}
close OUT;


open OUT, ">hap_genotype_matrix.txt";
print OUT "\t";
print OUT (join "\t", @existing_genes);
print OUT "\n";
foreach my $sample (@samples) 
{
    print OUT "$sample";
    foreach my $gene (@existing_genes)
    {
        if (exists $sample2gt{$gene}{$sample}) 
        {
            my $gt =  $sample2gt{$gene}{$sample};
            $gt=~s/:.+//;
            $gt=~s/\//|/;
            print OUT "\t$gt"; 
        }
        else
        {
            print OUT "\tNA"; 
        }

    }
    print OUT "\n";
}
close OUT;

open OUT, ">tag_genotypes.txt";
print OUT "Gene\tAlleles\t";
print OUT (join "\t", @samples);
print OUT "\n";


LOOPi7_3:foreach my $gene (@genes) 
{
    unless ($haplotypefrq{$gene}) 
    {
        next LOOPi7_3;
    }
    my %frq = %{$haplotypefrq{$gene}};
    my @haplotypes = reverse sort {$frq{$a} <=> $frq{$b}} keys %frq;
    my $total=0;
    foreach my $haplotype (@haplotypes) 
    {
        $total += $frq{$haplotype};
    }
    my $haplotypestr = "";
    if ($total==0) 
    {
        next LOOPi7_3;
    }
    foreach my $haplotype (@haplotypes) 
    {
        my $frq =  sprintf("%.2f", $frq{$haplotype}/$total);
        $haplotypestr .= "$haplotype($frq);";
    }

    print OUT "$gene\t$haplotypestr";
    foreach my $sample (@samples) 
    {
        if (exists $sample2gt{$gene}{$sample}) 
        {
            print OUT "\t", $sample2gt{$gene}{$sample}; 
        }
        else
        {
            print OUT "\t", "./.:0"; 
        }
    }
    print OUT "\n";
}
close OUT;
                                
print "Identifying top 2 haplotypes finished.\n\n";
