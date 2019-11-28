#usage perl count.hm.reads.v3.pl $samfile $asfile $outfile 
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use Bio::DB::Sam;

use strict;
use warnings;

my $flanking = 150;
my $swindow = 15;

sub count_hm_bin{
	my $infile = shift;
	my $asfile = shift;
	my $outfile = shift;
	my %aligned_pos;

	##processing HM bam file ####
	my $bam = Bio::DB::Bam -> open($infile);
	my $header = $bam->header;
	my $target_names = $header->target_name;
	while(my $align = $bam->read1){
		my $hit = reads_hit($align);
		my $chr = $target_names->[$align->tid];
		my $start = $align->start;
		print "$chr\t$start\n";
		exit;
	}
	close FILE;

	###processing AS file###
	open(FILE, $asfile) or die;
	open(OUT, "> $outfile") or die;
	while(<FILE>){
		chomp;
		s/"//g;
		if(/^ID/){
			print OUT "$_\tchip_left\tchip_right\n";
			next;
		}
		my @line = split(/\t/,$_);
		my $sexon_left = $line[5] - $flanking;
		my $sexon_right = $line[5] + $flanking;
		my $eexon_left = $line[6] - $flanking;
		my $eexon_right = $line[6] + $flanking;
		my $left_count = "";
		my $right_count = "";
		while($sexon_left < $sexon_right){
			my $count = 0;
			for(my $i=$sexon_left; $i<$sexon_left+$swindow; $i++){
				if(exists $aligned_pos{"$line[3]\t$i"}){
					$count = $count + $aligned_pos{"$line[3]\t$i"};	
				}
			}
			$left_count = "$left_count\,$count";
			$sexon_left = $sexon_left + $swindow;
		}
		while($eexon_left < $eexon_right){
			my $count = 0;
			for(my $i=$eexon_left; $i<$eexon_left+$swindow; $i++){
                        	if(exists $aligned_pos{"$line[3]\t$i"}){
                                	$count = $count + 1;
                        	}
                	}
                	$right_count = "$right_count\,$count";
                	$eexon_left = $eexon_left + $swindow;
		}
		$left_count =~ s/^,//; $right_count =~ s/^,//;
		print OUT "$_\t$left_count\t$right_count\n";
		#exit;
	}
	close FILE;
	close OUT;
}


# please replace the inidividual file with your own file names
$samfile = "sam file fore each sample";
$asfile = "rmatsClassified for each sample"
$outfile = "output file names"

&count_hm_bin("$samfile", "$asfile", "$outfile");
#my $as_dir = "/project/eheller_itmat_lab/HM_splicing/data/encode/mouse/embryo/rMAST/tissue_different_time_points/";
#my $sam_dir = "/project/eheller_itmat_lab/HM_splicing/data/encode/mouse/embryo/chip_seq/all_tissue/sam/";
#my $out_dir = "/project/eheller_itmat_lab/HM_splicing/data/encode/mouse/embryo/HM_signal/chip_signal/timepoints/forebrain/";

#my $base = "11.5";
#open(IN, "$ARGV[0].difftimepoints.list.txt") or die;
#while(<IN>){
#	chomp;
#	my $samfile = $_;
#	my $asfile = $samfile;
#	$asfile =~ s/mixed\.//;
#	$asfile =~ s/day.*//;
#	$asfile = "$asfile\.vs.$base\.rmast.out.SE";
#	my $outfile = "$samfile.hm.signal";
#	&count_hm_bin("$sam_dir$samfile", "$as_dir$asfile", "$out_dir$outfile");
#	exit;
#}
#close IN;
