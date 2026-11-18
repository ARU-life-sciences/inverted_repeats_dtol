# 0. merge the files together 

csvtk join -t -f species taxon_working latest_assemblies_hdr.txt > TWLA

csvtk replace -t -f species -p "_" -r " " ../out/summary_stats.tsv > summary_stats_spaces.tsv

csvtk join -t -f species TWLA summary_stats_spaces.tsv > TWLASS

csvtk replace -t -f species -p "_" -r " " assembly_quality.tsv > assembly_quality_spaces.tsv

csvtk join -t -f species TWLASS assembly_quality_spaces.tsv > TWLASSAQ

csvtk replace -t -f species -p "_" -r " " otol_tip_to_species.tsv > otol_tip_to_species_spaces.tsv

csvtk join -t -f species TWLASSAQ.tsv ../phylo/otol_tip_to_species_spaces.tsv > TWLASSAQO


