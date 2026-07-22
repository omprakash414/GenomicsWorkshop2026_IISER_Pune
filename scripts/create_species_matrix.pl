use strict;

my %sample_hash;
my %species_hash;
my %sample_species_hash;

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
		if($a[$#a] >= 0.70)
		{
			my $species_name = $a[$#a - 1];
			if(not exists $species_hash{$species_name})
			{
				$species_hash{$species_name} = 1;
			}

			if(not exists $sample_species_hash{$sample_name}{$species_name})
			{
				$sample_species_hash{$sample_name}{$species_name} = 1;
			}
			else
			{
				$sample_species_hash{$sample_name}{$species_name} += 1;
			}
		}
	}
}

foreach my $ckeys (keys %species_hash)
{
	print "\t$ckeys";
}
print "\n";

foreach my $ckeys (keys %sample_hash)
{
	print "$ckeys";
	foreach my $ckeys1 (keys %species_hash)
	{
		if(not exists $sample_species_hash{$ckeys}{$ckeys1})
                {
                	$sample_species_hash{$ckeys}{$ckeys1} = 0;
                }
		print "\t$sample_species_hash{$ckeys}{$ckeys1}";
	}
	print "\n";
}

