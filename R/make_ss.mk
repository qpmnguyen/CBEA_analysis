
samp_eval_fdr.qs: ss_typeIerror_pwr.R fdr_sim/simulation_*.qs
	Rscript ss_typeIerror_pwr.R --ncores 18 --

