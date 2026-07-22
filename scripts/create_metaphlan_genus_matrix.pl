use strict;
my $file_list = $ARGV[0];
my $out_file = $ARGV[1];

my %sample_hash;
my %genus_hash;
my %sample_genus_hash;

open(FILE_LIST, "$file_list");
while (<FILE_LIST>) 
{
    my $FF = $_;
    chomp $FF;
    my @a = split("\t", $FF);
    my $file_name = $a[0];
    my $sample_id = $a[1];
    if(not exists $sample_hash{$sample_id})
    {
		$sample_hash{$sample_id} = 1;
	}

    open(SAMPLE_FILE, "$file_name");
    while (<SAMPLE_FILE>) 
    {
        my $FF1 = $_;
        chomp $FF1;
        unless ($FF1 =~ /^#/) 
        { 
            # Skip header lines
            my @a = split(/\t/, $FF1);
            my $lineage = $a[0];
            my $abundance = $a[2];
            my @lineage_array = split(/\|/, $lineage);
            if ($lineage_array[$#lineage_array] =~ /g__/) # Extract genus-level taxonomy
            { 
                my $genus_name = $lineage_array[$#lineage_array];
                $genus_name =~ s/g__//g; # Remove 'g__'
                if(not exists $genus_hash{$genus_name})
				{	
					$genus_hash{$genus_name} = 1;
				}
				if(not exists $sample_genus_hash{$sample_id}{$genus_name})
				{
					$sample_genus_hash{$sample_id}{$genus_name} = $abundance;
				}                
            }
        }
    }
}

open(OUT, ">$out_file");
foreach my $ckeys (sort keys %genus_hash) 
{
    print OUT "\t$ckeys";
}
print OUT "\n";

foreach my $ckeys (keys %sample_hash) 
{
    print OUT "$ckeys";
    foreach my $ckeys1 (sort keys %genus_hash)
    {
		if(not exists $sample_genus_hash{$ckeys}{$ckeys1})
		{
			$sample_genus_hash{$ckeys}{$ckeys1} = 0;
		}
		print OUT "\t$sample_genus_hash{$ckeys}{$ckeys1}";
	}
    print OUT "\n";
}

