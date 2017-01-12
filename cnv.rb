class Cnv
	 attr_accessor :chrom, :pos, :gene, :transcript, :type, :var_type, :coding_effect
	 attr_accessor :g_dna_start, :g_dna_end, :g_nomen, :exon, :allele_freq
	 attr_accessor :filter, :proband_genotype, :alleles
	 
	def initialize
		self.alleles = Hash.new
	end
	
	def check_genotype(query_sample_id)
			genotype_check = false

			self.alleles.each_pair do |this_sample_id, this_allele|

				if query_sample_id.downcase == this_sample_id.downcase
					
					if (this_allele.gt.eql? './.') || (this_allele.gt.eql? '0/0')
						genotype_check = false
					else
						self.proband_genotype = this_allele.gt
						genotype_check = true
					end
				end
			end
			return genotype_check
	end
	
end
