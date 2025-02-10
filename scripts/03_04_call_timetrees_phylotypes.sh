# 03 timetree

## ZM
# gag
Rscript 03_timetrees.R "popart" "global_aln_naive" "gag" "results/rds/popart_md_subype_art_status.rds" > 03_timetrees_zm_gag.log 2>&1 &
# pol
Rscript 03_timetrees.R "popart" "global_aln_naive" "pol" "results/rds/popart_md_subype_art_status.rds" > 03_timetrees_zm_pol.log 2>&1 &
# env
Rscript 03_timetrees.R "popart" "global_aln_naive" "env" "results/rds/popart_md_subype_art_status.rds" > 03_timetrees_zm_env.log 2>&1 &
# wgs
Rscript 03_timetrees.R "popart" "global_aln_naive" "wgs" "results/rds/popart_md_subype_art_status.rds" > 03_timetrees_zm_wgs.log 2>&1 &
# V3 loop
Rscript 03_timetrees.R "popart" "global_aln_naive" "v3" "results/rds/popart_md_subype_art_status.rds" > 03_timetrees_zm_v3.log 2>&1 &

# 04 timetree2

## ZM
# gag
Rscript 03_timetrees2.R "popart" "global_aln_naive" "gag" "0.000001" > 03_timetrees2_zm_gag.log 2>&1 &
# pol
Rscript 03_timetrees2.R "popart" "global_aln_naive" "pol" "0.00001" > 03_timetrees2_zm_pol.log 2>&1 &
# env
Rscript 03_timetrees2.R "popart" "global_aln_naive" "env" "0.000001" > 03_timetrees2_zm_env.log 2>&1 &
# wgs
Rscript 03_timetrees2.R "popart" "global_aln_naive" "wgs" "0.2" > 03_timetrees2_zm_wgs.log 2>&1 &
# v3
Rscript 03_timetrees2.R "popart" "global_aln_naive" "v3" "0.2" > 03_timetrees2_zm_v3.log 2>&1 &