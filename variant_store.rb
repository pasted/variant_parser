class VariantStore
		require 'ruby-progressbar'
		require 'digest/bubblebabble'
		
		attr_accessor :variants, :hashed_variants
		attr_accessor :gene_list_variant_count
		
		def initialize(variants)
			self.variants = variants
  	end
  	
  	def convert_variants_to_hash
  		hashed_variants = Array.new
  		self.variants.each do |this_variant|
  			variant_hash = {}
  			variant_hash = this_variant.to_hash
  			hashed_variants.push(variant_hash)
  		end
  		return hashed_variants
  	end
  	
  	def collapse_variants(variants)

  		outer_tmp_variants = variants
  		inner_tmp_variants = variants
  		result_tmp_variants = Array.new
  		result_hash = Hash.new
  		processed_uids = Array.new
  		
  		outer_tmp_variants.each_with_index do |this_variant, outer_index|
  				
  				
  				inner_tmp_variants.each_with_index do |variant_to_check, inner_index|
  					if (this_variant.chromosome == variant_to_check.chromosome) && (this_variant.position == variant_to_check.position) && (this_variant.nearest_splice_site_change == variant_to_check.nearest_splice_site_change) && (this_variant.local_splice_effect == variant_to_check.local_splice_effect) 
  							
  							hash_key = Digest::SHA256.bubblebabble "#{variant_to_check.chromosome}:#{variant_to_check.position}:#{variant_to_check.nearest_splice_site_change}:#{variant_to_check.local_splice_effect}"
  							if result_hash.has_key?(hash_key)
  								variants_ordered_by_length = result_hash.fetch(hash_key)									
  							else
  								variants_ordered_by_length = Array.new
  							end
  							#check if variants share genomic location and variant data
  							#variant_to_check.collapsed = true 
  							if !processed_uids.include?(variant_to_check.uid) 
  								variants_ordered_by_length.push(variant_to_check)
  								result_hash[hash_key] = variants_ordered_by_length
  								processed_uids.push(variant_to_check.uid)
  							end
  							

  					end#if clause for equivalent genomic location
  					
  				end
  				
  		end#outer loop

  		result_hash.each_pair do |this_key, variants_ordered_by_length|
  				variants_ordered_by_length.sort! { |a,b| a.transcript_length.to_i <=> b.transcript_length.to_i}
  				transcripts = variants_ordered_by_length.map { |v| "#{v.full_transcript}(#{v.assembly})" }
  				transcripts_string = transcripts.join(",")
  				collapsed_variant = variants_ordered_by_length[0]
  				collapsed_variant.transcripts = transcripts_string
  				result_tmp_variants.push(collapsed_variant)  			
  		end

  		return result_tmp_variants
  	end
  	
    def select_variants(special_gene_symbols, special_gene_ids, sample_id, parse_all, select_clinvar, maf_cutoff, research)
    	progressbar = ProgressBar.create(:title => "Processing variants",:format => "%t %c : %C %w", :total => self.variants.length)
    	selected_variants = Array.new
  		not_selected_variants = Array.new
  		self.gene_list_variant_count = 0
  		
  		#Loop through variants 
  		self.variants.each do |this_variant|
  				
  				passed_maf_cutoff = false
  				selected = false 
  				
  				if sample_id
  					genotype_check = this_variant.check_genotype(sample_id)				
  				else
  					genotype_check = true
  				end
  				 				
  				passed_maf_cutoff = this_variant.check_maf_cutoff(maf_cutoff)
  				
  				this_variant.find_highest_maf

  				#Check if the transcript is correct
  				if ( special_gene_symbols.include?(this_variant.gene) || special_gene_ids.include?("#{this_variant.gene_id}") || parse_all ) && genotype_check  && passed_maf_cutoff
  					self.gene_list_variant_count = self.gene_list_variant_count + 1
  					#if this_variant.check_maf_cutoff(maf_cutoff)
  						if ['DM', 'DM?', 'FTV', 'R'].include?(this_variant.hgmd_sub_category)
								# DP, DFP, FP, FTV, DM?, DM, R
								#1.	Select all HGMD variants marked as DM
								this_variant.reason_for_selection = "HGMD sub-category"
								selected_variants.push(this_variant)
								selected = true
								#:clin_var_ids, :clin_var_origins, :clin_var_methods, :clin_var_clin_signifs, :clin_var_review_status, :clin_var_phenotypes
							elsif select_clinvar && (this_variant.clin_var_clin_signifs.contains?("Pathogenic") && this_variant.clin_var_origins.contains?("germline") && (this_variant.clin_var_review_status.contains?("3") || this_variant.clin_var_review_status.contains?("4")) )
								#optional.	Select all ClinVar variants that meet criteria, if flag in command line
								this_variant.reason_for_selection = "ClinVar category"
								selected_variants.push(this_variant)
								selected = true
							elsif this_variant.var_type == 'substitution' && this_variant.filter_vcf == 'PASS'
								#2.	Select substitutions that contain ‘PASS’ in Filter(VCF) field; select all indels
								
								if ['missense', 'stop gain', 'nonsense', 'start loss', 'stop loss'].include?(this_variant.coding_effect)
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
										elsif ( ['upstream', '5\'UTR', '3\'UTR', 'downstream'].include?(this_variant.var_location) ) && ( (-50..10).include?(this_variant.distance_nearest_splice_site.to_i) && research)
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
					#	end#check MAF cutoff
  					
  				end#check gene list
  				if selected == false
  					not_selected_variants.push(this_variant)
  				end
  				progressbar.increment
  		end#variants
  		
  			tmp_store = [selected_variants, not_selected_variants]

  		return tmp_store
  		
  	end

end
