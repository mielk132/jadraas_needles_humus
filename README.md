### jadraas_experiment ###

#Mielke et al., 2025 https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.70316

#organice bags (pine needle litter and humus) incubated in two rounds (Set A and B) from 2017-2018 and 2018-2019.

---------------
#datafiles: with info on column names

---------------
#community.txt all samples (background, seq replicates (e.g. dna extracts 60,61), missing samples, neg controls, even mock)

---------------
#metadata.txt all 342 samples (background = initial substrates, seq replicates (e.g. dna extracts 60,61), missing samples, neg controls, even mock)

#Sample	
#Unique_Sample	
#order_for_unique_id	
#treatment = control, disturbance control, shrub removal, combined shrub removal and pine root exclusion (trenching), pine root exclusion (trenching)
#shrubs =  whether shrubs are present (1) or absent (0)	
#pine =  whether pine roots are present (1) absent (0)	
#block	substrate	
#start_date	= start of incuvation in the field
#end_date	= end of incubation in the experiment 
#label = color of bag to tell whether it was incubation for 5 months (blue) or 17 months (green)
#set = whether the incubations began in 2017 (set 1 or A) or 2018 (set 2 or B)
#incubation	= length of incubatiuon
#original_substrate_mass_g	
#deployed_bag_mass_g	
#bag_mass_g	
#final_bag_mass_g	
#final_substrate_g	
#substrate_mass_difference_g	= difference between 17 and 5 months
#percent_loss	
#percent_mass_remaining
#DNA_Extraction_ID	= labels on tubes for DNA extraction
#copies_DNA_per_g_substrate	= number of ITS2 copies estimated by qPCR per g of substrate (whether humus or pine needles)
#DNA_ug_per_g_substrate	= amount of DNA extracted as measured by the Nanodrop
#ITS_copies_per_ng_DNA	= number of ITS2 copies estimated by qPCR per ng of DNA
#tagNum	= primer tags for reverse and forward primers
#pool	- pacBio sequencing pool (1,2,3, or 4)
#notes	
#scata_ID = unique sample ID with matching tag numbers

---------------
#taxa.txt ordered by original global abundances after sequencing and scata

#scata	= OTU cluster
#Taxon = name for matching
#Sci_Name = unique name
#Rank = Original abundance order
#mock = true or false if in mock community (10 sequences with different bp lengths)
#Taxonomy: Kingdom Phylum	Subphylum	Class	Order	Family	Genus	Species abv.
#guild = ecto white rot, ecto other, sap white rot, sap other, ericoid, ericoid-ecto, moulds and yeasts, unknown, plants
#emf = yes (1) or no (0)
#other categories not essential for analyses: primary_lifestyle_phylum, Decay_type, Decay_Lifestyle	
#guild_notes = literature or sequence references used to determine lifestyle, white rot capacity and guild
#reference = UNITE Species Hypothesis (1.5%)

#raw data PRJNA834027 in the Sequence Read Archive at NCBI
