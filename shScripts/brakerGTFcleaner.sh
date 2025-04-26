awk 'BEGIN{OFS="\t"} 
$3 == "transcript" {
  if ($9 !~ /transcript_id/) {
    match($9, /g[0-9]+\.t[0-9]+/, tid);
    if (tid[0] != "") {
      $9 = $9 "; transcript_id \"" tid[0] "\"; gene_id \"" tid[0] "\";";
    }
  }
}
{ print }' braker.gtf > braker_fixed.gtf

# braker creates some incomplete annotations sometimes, if you want to include these
# use thisscript to add the necessary parameters for a proper GTF file

# !!!!IMPORTANT CAVEAT!!!!
# If you do include these, know that they may be incomplete annotations and that
# inclusion may increase the noise of annotations and lower overall quality
