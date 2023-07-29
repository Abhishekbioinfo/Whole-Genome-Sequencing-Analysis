my %fileData;

sub process_files {

  my $file = $_[0];
  my $name = $_[1];
  
  if ($name ne "dbnsfp42a") {
    open(in, $file);
    while (my $line = <in>) {
      chomp $line; 
      my $id  = (split(/\s+/,$line))[1];
      my $chr = (split(/\s+/,$line))[2];
      my $pos = (split(/\s+/,$line))[3];
      $fileData{$name}->{$chr}{$pos} = $id;
    }
    close in;  
  } else {
    my @search = qw(6 8 12 14 15 17 18 20 21 23 24 26 27 29 30 32 33 75 77 78 80 82 35 37 38 40 92 94 98 100 102 104 106 108 44 46 110 111 112 123);
    open(in, $file);
    while (my $line = <in>) {
      chomp $line;
      my $info;
      my $chr = (split(/\t/,$line))[0];
      my $pos = (split(/\t/,$line))[1]; 
      for(my $i = 0; $i <= $#search; $i++) {
        my $j  = $search[$i]-1;
        $info .= (split(/\t/,$line))[$j]."\t";
      }
      $fileData{$name}->{$chr}{$pos} = $info;
    }
    close in; 
  }
}

process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_avsnp150_dropped", "snp");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_cosmic96_coding_dropped", "cosmic");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_esp6500siv2_all_dropped", "esp");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_tmcsnpdb_dropped", "pgt2");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_clinvar_20220320_dropped", "clinvar");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_ALL.sites.2015_08_dropped", "g1000");
process_files("Somatic_Annovar/Variant-annovar.annotate.avinput.hg38_dbnsfp42a.hg38_multianno.txt", "dbnsfp42a");


my $header = "Line\tMutation_Types\tGene\tTranscriptID\tExon_Pos\tcDNA_Change\tAA_Change\tCHR\tLOCATION_ST\tLOCATION_END\tREF-Al\tALT-AL\tNA\tGenotype\tTotal_Depth\tCHR\tLOCATION_ST\tNA\tREF-Al\tALT-AL\tGenotype\tNA\tGenotype\tTotal_Depth\tRef_Depth\tVariant_Depth\tVariant_Freq(%)\tsnp\tcosmic\tesp\tTmcSNP\tclinvar\tMAF\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\tMutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score\tFATHMM_pred\tPROVEAN_score\tPROVEAN_pred\tVEST4_score\tCADD_raw\tCADD_phred\tDANN_score\tfathmm-MKL_coding_score\tfathmm-MKL_coding_pred\tMetaSVM_score\tMetaSVM_pred\tMetaLR_score\tMetaLR_pred\tintegrated_fitCons_score\tintegrated_confidence_value\tGERP_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphastCons100way_vertebrate\tphastCons30way_mammalian\tSiPhy_29way_logOdds\tM-CAP_score\tM-CAP_pred\tInterpro_domain\tGTEx_V8_gene\tGTEx_V8_tissue\tOtherinfo\n";


print ("$header");
open(file, "Somatic_Annovar/Variant-annovar.annotate.avinput.exonic_variant_function_non-synon");
while (my $line = <file>) {
  chomp $line;
  my $chr = (split(/\t/,$line))[7];
  my $pos = (split(/\t/,$line))[8];
  my $snp       = $fileData{"snp"}->{$chr}{$pos} eq "" ? 0 : $fileData{"snp"}->{$chr}{$pos};
  my $cosmic    = $fileData{"cosmic"}->{$chr}{$pos} eq "" ? 0 : $fileData{"cosmic"}->{$chr}{$pos}; 
  my $esp       = $fileData{"esp"}->{$chr}{$pos} eq "" ? 0 : $fileData{"esp"}->{$chr}{$pos};    
  my $pgt2      = $fileData{"pgt2"}->{$chr}{$pos} eq "" ? 0 : $fileData{"pgt2"}->{$chr}{$pos};  
  my $clinvar   = $fileData{"clinvar"}->{$chr}{$pos} eq "" ? 0 : $fileData{"clinvar"}->{$chr}{$pos};
  my $g1000     = $fileData{"g1000"}->{$chr}{$pos} eq "" ? 0 : $fileData{"g1000"}->{$chr}{$pos};
  my $dbnsfp42a = $fileData{"dbnsfp42a"}->{$chr}{$pos} eq "" ? "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" : $fileData{"dbnsfp42a"}->{$chr}{$pos};

  #print ("$chr\t$pos\t$snp\t$cosmic\t$esp\t$pgt2\t$clinvar\t$g1000\t$dbnsfp42a\n")
	print ("$line\t$snp\t$cosmic\t$esp\t$pgt2\t$clinvar\t$g1000\t$dbnsfp42a\n")
}
close file;
