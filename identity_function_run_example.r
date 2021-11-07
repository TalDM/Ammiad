#########################
### load the function ###

source("identity_function.r")


### run it with specific files: ###  

ammiad_summarize_vcf(vcf_path = "845.vcf",
					map_csv_path = "map.csv",
					output_csv_path = "output.csv",
					identity_threshold = 0.981)


### assign colors to DGGs ### 

ammiad_colorize_by_groups(summarized_table_path="output.csv",
    colors_table_path="DGG_colors.csv",
		output_csv_path="output_with_colors.csv")


#########################
