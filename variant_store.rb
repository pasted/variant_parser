class VariantStore
		attr_accessor :variants
		
		def initialize(variants)
			self.variants = variants
  	end
  	    
    def select_variants(special_gene_symbols, special_gene_ids)
    	selected_variants = Array.new
  		not_selected_variants = Array.new
  		
  		#Loop through variants 
  		self.variants.each do |this_variant|
  				
  				selected = false
  				#Check if the transcript is correct
  				if ( special_gene_symbols.include?(this_variant.gene) || special_gene_ids.include?("#{this_variant.gene_id}") )
  					
  						if ['DM', 'DM?', 'FTV', 'R'].include?(this_variant.hgmd_sub_category)
								# DP, DFP, FP, FTV, DM?, DM, R
								#1.	Select all HGMD variants marked as DM
								this_variant.reason_for_selection = "HGMD sub-category"
								selected_variants.push(this_variant)
								selected = true
							elsif this_variant.var_type == 'substitution' && this_variant.filter_vcf == 'PASS'
								#2.	Select substitutions that contain ‘PASS’ in Filter(VCF) field; select all indels
								
								if ['missense', 'nonsense', 'start loss', 'stop loss'].include?(this_variant.coding_effect)
									#3.	Select all coding non-synonymous variants
									this_variant.reason_for_selection = "Coding effect"
									selected_variants.push(this_variant)
									selected = true
								else
										if [-3..3].include?(this_variant.distance_nearest_splice_site)
											#This should also catch canonical splice site variants
											#4.	Select all synonymous variants within 3bp from splice site
											this_variant.reason_for_selection = "Distance to splice site"
											selected_variants.push(this_variant)
											selected = true
										elsif ( (-50..10).include?(this_variant.distance_nearest_splice_site.to_i) ) && this_variant.nearest_splice_site_change && ( ((this_variant.nearest_splice_site_change != nil) || (this_variant.nearest_splice_site_change.empty? == false)) ) && (this_variant.nearest_splice_site_change != 0)
											#5.	Select all variants that have values in the "nearestSSChange" or "localSpliceEffect" fields
											this_variant.reason_for_selection = "Nearest splice site change"
											selected_variants.push(this_variant)
											selected = true
										elsif ( (-50..10).include?(this_variant.distance_nearest_splice_site.to_i) ) && this_variant.local_splice_effect && ((this_variant.local_splice_effect != nil) || (this_variant.local_splice_effect.empty? == false))
											#Possible values: Cryptic Donor Strongly Activated, Cryptic Donor Weakly Activated, Cryptic Acceptor Strongly Activated, Cryptic Acceptor Weakly Activated, New Donor Site, New Acceptor Site
											this_variant.reason_for_selection = "Local splice effect"
											selected_variants.push(this_variant)
											selected = true
										elsif ( ['upstream', '5\'UTR', '3\'UTR', 'downstream'].include?(this_variant.var_location) ) && ( (-50..10).include?(this_variant.distance_nearest_splice_site.to_i) )
										#upstream, 5'UTR, exon, intron, 3'UTR, downstream AND -50 to +10 of a known splice site
										#6.	Select all variants with 'varLocation' of '3_UTR', '5_UTR', 'Upstream' and 'Downstream'
											this_variant.reason_for_selection = "Variant location"
											selected_variants.push(this_variant)
											selected = true
										end # Rest of selectors
								
							 end # coding_effect
						
						elsif ( this_variant.var_type != 'substitution' ) && ( ['exon'].include?(this_variant.var_location) )
							#2. Select all Indels : var_type = ['duplication', 'insertion', 'deletion', 'delins'] AND with a variant_location of 'exon'	
  
							this_variant.reason_for_selection = "Exonic Indel"
							selected_variants.push(this_variant)
							selected = true
						elsif ( this_variant.var_type != 'substitution' ) && ( (-50..10).include?(this_variant.distance_nearest_splice_site.to_i) )
							#2. Select all Indels : var_type = ['duplication', 'insertion', 'deletion', 'delins'] AND within -50 to +10 of a known splice site
							
							this_variant.reason_for_selection = "Indel"
							selected_variants.push(this_variant)
							selected = true

						end # HGMD sub-cat / var_type

  					
  				end
  				if selected == false
  					not_selected_variants.push(this_variant)
  				end
  		
  		end#variants
  		
  			tmp_store = [selected_variants, not_selected_variants]

  		return tmp_store
  		
  	end

end
