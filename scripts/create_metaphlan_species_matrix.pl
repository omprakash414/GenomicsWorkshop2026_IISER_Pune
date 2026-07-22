use strict;
my $file_list = $ARGV[0];
my $out_file = $ARGV[1];

my %sample_hash;
my %species_hash;
my %sample_species_hash;

open(FILE_LIST,"$file_list");
while(<FILE_LIST>)
{
	my $FF = $_;
	chomp $FF;
	my @a = split("\t",$FF);
	my $file_name = $a[0];
	my $sample_id = $a[1];
	if(not exists $sample_hash{$sample_id})
	{
		$sample_hash{$sample_id} = 1;
	}
	#print "$file_name\n";
	open(SAMPLE_FILE,"$file_name");
	while(<SAMPLE_FILE>)
	{
		my $FF1 = $_;
		chomp $FF1;
		unless($FF1 =~ /^#/)
		{
			#print "$FF1\n";
			my @a = split(/\t/,$FF1);
			my $lineage = $a[0];
			my $abundance = $a[2];
			my @lineage_array = split(/\|/,$lineage);
			if($lineage_array[$#lineage_array] =~ /s__/)
			{
				my $species_name = $lineage_array[$#lineage_array];
				$species_name =~ s/s__//g;
				if(not exists $species_hash{$species_name})
				{	
					$species_hash{$species_name} = 1;
				}
				if(not exists $sample_species_hash{$sample_id}{$species_name})
				{
					$sample_species_hash{$sample_id}{$species_name} = $abundance;
				}
			}
		}
	}
}

open(OUT,">$out_file");
foreach my $ckeys (sort keys %species_hash)
{
	print OUT "\t$ckeys";
}
print OUT "\n";

foreach my $ckeys (keys %sample_hash)
{
	print OUT "$ckeys";
	foreach my $ckeys1 (sort keys %species_hash)
	{
		if(not exists $sample_species_hash{$ckeys}{$ckeys1})
		{
			$sample_species_hash{$ckeys}{$ckeys1} = 0;
		}
		print OUT "\t$sample_species_hash{$ckeys}{$ckeys1}";
	}
	print OUT "\n";
}

