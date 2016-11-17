class CnvStore
		attr_accessor :cnvs
		
		def initialize(cnvs)
			self.cnvs = cnvs
  	end
  	    

    def select_cnvs(special_gene_symbols, special_gene_ids, sample_id, parse_all)

    	selected_cnvs = Array.new
  		not_selected_cnvs = Array.new
  		
  		#Loop through cnvs 
  		self.cnvs.each do |this_cnv|

  				if sample_id != ''
  					genotype_check = this_cnv.check_genotype(sample_id)
  				else
  					genotype_check = true
  				end

  				selected = false
  				#Check if the gene symbol is present
  				if ( special_gene_symbols.include?(this_variant.gene) || special_gene_ids.include?("#{this_variant.gene_id}") || parse_all ) && (genotype_check == true)

								this_cnv.reason_for_selection = "Gene list"
								selected_cnvs.push(this_cnv)
								selected = true
  				end
  				if selected == false
  					not_selected_cnvs.push(this_cnv)
  				end
  		
  		end#cnvs
  		
  			tmp_store = [selected_cnvs, not_selected_cnvs]

  		return tmp_store
  		
  	end

end
