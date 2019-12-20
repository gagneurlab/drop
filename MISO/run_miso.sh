#!/bin/bash

# first recalculate gff index
rm -rf gff_data/
index_gff --index ./settings/new_splice_site_events.gff ./gff_data/

# create sashimi plots
sashimi_plot --plot-event "TIMMDC1" gff_data/ settings/plot_timmdc1_new_exon.conf   --output-dir plots/
