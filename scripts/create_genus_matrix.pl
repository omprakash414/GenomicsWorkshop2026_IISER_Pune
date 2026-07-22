use strict;

my %sample_hash;
my %genus_hash;
my %sample_genus_hash;

open(LIST,"$ARGV[0]");
while(<LIST>)
{
	my $FF = $_;
	chomp $FF;
	my $sample_file = $FF;
	my $sample_name = $FF;
	$sample_name =~ s/.spingo.out.txt//g;
	if(not exists $sample_hash{$sample_name})
	{
		$sample_hash{$sample_name} = 1;
	}
	open(SAMPLE_FILE,"$sample_file");
	while(<SAMPLE_FILE>)
	{
		my $FF1 = $_;
		chomp $FF1;
		my @a = split(/\t/,$_);
		if($a[5] >= 0.70)
		{
			my @g_array = split("_",$a[4]);
			my $genus_name = $g_array[0];
			if(not exists $genus_hash{$genus_name})
			{
				$genus_hash{$genus_name} = 1;
			}

			if(not exists $sample_genus_hash{$sample_name}{$genus_name})
			{
				$sample_genus_hash{$sample_name}{$genus_name} = 1;
			}
			else
			{
				$sample_genus_hash{$sample_name}{$genus_name} += 1;
			}
		}
	}
}

foreach my $ckeys (keys %genus_hash)
{
	print "\t$ckeys";
}
print "\n";

foreach my $ckeys (keys %sample_hash)
{
	print "$ckeys";
	foreach my $ckeys1 (keys %genus_hash)
	{
		if(not exists $sample_genus_hash{$ckeys}{$ckeys1})
                {
                	$sample_genus_hash{$ckeys}{$ckeys1} = 0;
                }
		print "\t$sample_genus_hash{$ckeys}{$ckeys1}";
	}
	print "\n";
}

