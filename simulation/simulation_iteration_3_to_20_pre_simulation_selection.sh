#!/bin/bash
IND=1762;
Rscript true_SV_sets_dup_check.R
Rscript true_SV_sets_germline_dup_check.R
Rscript true_SV_sets_somatic.R
perl true_SV_sets_to_which_ref_coor.pl $IND true_SV_sets_somatic_dup_checked true_SV_sets_somatic_dup_checked.refcoor
perl true_SV_sets_to_which_ref_coor.pl $IND true_SV_sets_dup_checked true_SV_sets_dup_checked.refcoor
perl germline_initial_SVs_changed.pl $IND  | awk '{if($4-$7>1000 || $4-$7 < -1000) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >   germline_initial_SVs_changed_size1000
cat true_SV_sets_dup_checked.refcoor  germline_initial_SVs_changed_size1000 > true_SV_sets_dup_checked.refcoor.germline_initial_SVs_changed_size1000 


perl SVs.results_to_ref_coor.pl $IND germline
perl SVs.results_to_ref_coor.pl $IND total

total=`ls SVs.results.* | awk 'BEGIN{max=0}{split(\$0,f,"."); if(max<f[3]){max=f[3]}}END{print max}'`;
germ_total=`ls germline_SVs.results.* | awk 'BEGIN{max=0}{split(\$0,f,"."); if(max<f[3]){max=f[3]}}END{print max}'`;

Rscript results_to_cnv.R $total total
Rscript results_to_cnv.R $germ_total germline
Rscript results.cnv.germ_initial_edit.R total
Rscript results.cnv.germ_initial_edit.R germline
Rscript results.cnv.somatic.R

