class VariantParser
  require 'yaml'
  require 'smarter_csv'
  require 'spreadsheet'
  require 'trollop'
  require_relative 'variant'
  require_relative 'allele'
  require_relative 'variant_store'
  require_relative 'symbol_checker'
  
  
  def parse_alamut_file(file_name, sample_ids)
  	options = { :col_sep => "\t" }
  	variant_array = Array.new 
  	
  	if File.exists?(file_name) && ( File.stat(file_name).size > 0 )
  		
  		SmarterCSV.process( file_name, options ) do |csv|
  			this_variant = Variant.new
  			
  			this_variant.chromosome											= csv.first[:chrom]
  			this_variant.position												= csv.first[:pos]
  			this_variant.genomic_dna_start							= csv.first[:gdnastart]
  			this_variant.genomic_dna_end								= csv.first[:gdnaend]
  			this_variant.gene														= csv.first[:gene]
  			this_variant.gene_id												=	csv.first[:geneid]                             					
  			this_variant.full_transcript								= csv.first[:transcript]
  			this_variant.var_type												= csv.first[:vartype]
  			this_variant.coding_effect									= csv.first[:codingeffect]
  			this_variant.var_location										= csv.first[:varlocation]
  			this_variant.genomic_nomen									= csv.first[:gnomen]
  			this_variant.cdna_nomen											= csv.first[:cnomen]
  			this_variant.protein_nomen									= csv.first[:pnomen]
  			this_variant.exon														= csv.first[:exon]
  			this_variant.intron													= csv.first[:intron]
  			this_variant.distance_nearest_splice_site		= csv.first[:distnearestss]
  			this_variant.nearest_splice_site_type				= csv.first[:nearestsstype]
  			this_variant.nearest_splice_site_change			= csv.first[:nearestsschange]
  			this_variant.local_splice_effect						= csv.first[:localspliceeffect]
  			this_variant.rs_id													= csv.first[:rsid]
  			this_variant.rs_validated										= csv.first[:rsvalidated]
  			this_variant.rs_maf													= csv.first[:rsmaf]
  			this_variant.genomes_1000_freq							= csv.first[:"1000g_af"]
  			this_variant.genomes_1000_afr_freq					= csv.first[:"1000g_afr_AF"]
  			this_variant.genomes_1000_sas_freq					= csv.first[:"1000g_sas_AF"]
  			this_variant.genomes_1000_eas_freq					= csv.first[:"1000g_eas_AF"]
  			this_variant.genomes_1000_eur_freq					= csv.first[:"1000g_eur_AF"]
  			this_variant.genomes_1000_amr_freq					= csv.first[:"1000g_amr_AF"]
  			this_variant.exac_all_freq									= csv.first[:exacallfreq]
  			this_variant.exac_afr_freq									= csv.first[:exacafrfreq]
  			this_variant.exac_amr_freq									= csv.first[:exacamrfreq]
  			this_variant.exac_eas_freq									= csv.first[:exaceasfreq]
  			this_variant.exac_sas_freq									= csv.first[:exacsasfreq]
  			this_variant.exac_nfe_freq									= csv.first[:exacnfefreq]
  			this_variant.exac_fin_freq									= csv.first[:exacfinfreq]
  			this_variant.exac_oth_freq									= csv.first[:exacothfreq]
  			this_variant.exac_afr_hmz										= csv.first[:exacafrhmz]
  			this_variant.exac_amr_hmz										= csv.first[:exacamrhmz]
  			this_variant.exac_eas_hmz										= csv.first[:exaceashmz]
  			this_variant.exac_sas_hmz										= csv.first[:exacsashmz]
  			this_variant.exac_nfe_hmz										= csv.first[:exacnfehmz]
  			this_variant.exac_fin_hmz										= csv.first[:exacfinhmz]
  			this_variant.exac_oth_hmz										= csv.first[:exacothhmz]
  			this_variant.exac_filter										= csv.first[:exacfilter]
  			this_variant.exac_read_depth								= csv.first[:exacreaddepth]
  			this_variant.clin_var_ids										= csv.first[:clinvarids]
  			this_variant.clin_var_origins								= csv.first[:clinvarorigins]
  			this_variant.clin_var_methods								= csv.first[:clinvarmethods]
  			this_variant.clin_var_clin_signifs					= csv.first[:clinvarclinsignifs]
  			this_variant.clin_var_review_status					= csv.first[:clinvarreviewstatus]
  			this_variant.clin_var_phenotypes						= csv.first[:clinvarphenotypes]
  			this_variant.cosmic_ids											= csv.first[:cosmicids]
  			this_variant.cosmic_tissues									= csv.first[:cosmictissues]
  			this_variant.cosmic_freqs										= csv.first[:cosmicfreqs]
  			this_variant.cosmic_sample_counts						= csv.first[:cosmicsamplecounts]
  			                                    				
  			this_variant.esp_all_maf										= csv.first[:espallmaf]
  			this_variant.hgmd_phenotype									= csv.first[:hgmdphenotype]
  			this_variant.hgmd_pub_med_id								= csv.first[:hgmdpubmedid]
  			this_variant.hgmd_sub_category							= csv.first[:hgmdsubcategory]
  			this_variant.n_orthos												= csv.first[:northos]
  			this_variant.conserved_orthos								= csv.first[:conservedorthos]
  			this_variant.conserved_dist_species					= csv.first[:conserveddistspecies]
  			this_variant.grantham_dist									= csv.first[:granthamdist]
  			this_variant.agv_gd_class										= csv.first[:agvgdclass]
  			this_variant.sift_prediction								= csv.first[:siftprediction]
  			this_variant.wt_nuc													=	csv.first[:wtnuc]
  			this_variant.var_nuc												=	csv.first[:varnuc]
  			this_variant.ins_nucs												=	csv.first[:insnucs]
  			this_variant.del_nucs												= csv.first[:delnucs]
  			this_variant.filter_vcf											= csv.first[:"filter_(vcf)"]
  			
  			sample_ids.each do |sample_id|
  				sample_id_lc = sample_id.downcase
  				this_allele = Allele.new
  				this_allele.ad											= csv.first[:"ad_(#{sample_id_lc})"]
  				this_allele.dp											= csv.first[:"dp_(#{sample_id_lc})"]
  				this_allele.gq											= csv.first[:"gq_(#{sample_id_lc})"]
  				this_allele.gt											= csv.first[:"gt_(#{sample_id_lc})"]
  				this_allele.pl											= csv.first[:"pl_(#{sample_id_lc})"]
  				this_variant.alleles.store("#{sample_id}", this_allele)
  			end
  			
  			#this_variant.parse_transcript()
  			variant_array.push(this_variant)
  		end
  	else
  		puts "ERROR :: variants file has no content :: SAMPLE ID : #{sample_id_lc}"
  		if File.exists?(file_name)
  			puts "File exists? #{File.exists?(file_name)}"
  			puts "File size? #{File.stat(file_name).size}"
  		else
  			puts "File exists? #{File.exists?(file_name)}"
  		end
  	end
  	return variant_array
  end
  
  def parse_sample_ids(variants_filepath)
  	
  	header = File.foreach(variants_filepath).first
  	header_array = header.split(/\t/)
  	genotype_fields = header_array.select{ |e| e=~/GT/ }
  	sample_ids = genotype_fields.collect{ |e| e.gsub(/GT |\(|\)/, '') }
  	sample_ids.map{|id| id.downcase!}
  	return sample_ids
	end
	
	def parse_gene_list(genes_filepath)
		
		gene_symbol_list = Array.new
		File.foreach(genes_filepath) do |line|
  	  	gene_symbol_list.push(line.strip)
		end
		
		return gene_symbol_list
	end
	
	def parse_hpo_gene_list(hpo_genes_filepath)
		
		gene_id_list = Array.new
		File.foreach(hpo_genes_filepath) do |line|
			if !( line.match('Export for') || line.match('Associated diseases') || line.match('Total') )
				gene_symbol = line.split('","').first.delete('"').split(' ').first
				hgnc_id = line.split('","').first.delete('"').split(' ')[1].delete('(').delete(')')
				gene_id_list.push(hgnc_id)
				
			end
			
		end
		gene_id_list.sort!
		return gene_id_list
	end
	
	def variants_to_excel(results, workbook)
		this_sheet = workbook.create_worksheet :name => "Selected variants"
  	row_number = 0
		header_array = Array.new

		results[0].each do |this_selected_variant|
			
			if header_array.empty?
				header_array = this_selected_variant.variable_order
				header_array.map!{ |element| element.to_s }
				vcf_headers = this_selected_variant.print_alleles.map {|var| var[0]}
				header_array = header_array + vcf_headers
				header_array.each do |this_header|
					this_sheet.row(row_number).push this_header  
 				end
 			end
 			row_number = row_number + 1
 		
 			#output variables to spreadsheet in given order
 			this_selected_variant.variable_order.map {|var|  this_sheet.row(row_number).push  "#{this_selected_variant.send(var)}" }
 			this_selected_variant.print_alleles.map {|var| this_sheet.row(row_number).push "#{var[1]}" }
 		end
 		return workbook
	end
	
  opts = Trollop::options do
  	opt :variants, "Filepath to Alamut file to parse.", :type => String
  	opt :genes, "Filepath to text file with list of valid HGVS gene symbols - one symbol per line.", :type => String
  	opt :ontology_genes, "Filepath to CSV export from HPO web browser.", :type => String
  	opt :check_genes, "Check gene symbols against HGNC - only process those which are the current symbol"
  	opt :excel_output, "Export results to an Excel spreadsheet."
  	opt :proband_sample, "Proband sample ID to include in output", :type => String
  end
  
  #Trollop::die :alamut_file, "Alamut file must exist." unless File.exist?(opts[:alamut_file]) if opts[:alamut_file]

  variants_filepath = opts[:variants]
  genes_filepath = opts[:genes]
  hpo_genes_filepath = opts[:hpo_genes]
  excel_output = opts[:excel_output]
  proband_sample_id = opts[:proband_sample]
  
  
  if variants_filepath && (genes_filepath || hpo_genes_filepath)

  	gene_symbol_list = []
  	gene_id_list = []
  	parser = VariantParser.new()
  	
  	if hpo_genes_filepath
  		gene_id_list = parser.parse_hpo_gene_list(hpo_genes_filepath)
  	elsif genes_filepath
  		gene_symbol_list = parser.parse_gene_list(genes_filepath)
  	end
  	
  	if opts[:check_genes]
  		checked_symbols = Array.new
  		this_symbol_checker = SymbolChecker.new	
  		checked_gene_symbol_list = this_symbol_checker.check_gene_symbols(gene_symbol_list)
  		gene_symbol_list = checked_gene_symbol[0]
  	end
  	
  	sample_ids = parser.parse_sample_ids(variants_filepath)
  	variants = parser.parse_alamut_file(variants_filepath, sample_ids)
  	
  	variant_store = VariantStore.new(variants)
  	results = variant_store.select_variants(gene_symbol_list, gene_id_list, proband_sample_id)
  	puts results[0].length
  	
  	if excel_output
  		this_book = Spreadsheet::Workbook.new
  		this_book = parser.variants_to_excel(results, this_book)
  		this_book.write "results/#{Time.now.strftime("%Y-%m-%d-%H%M%S")}_#{ENV['USER']}.xls"
  	else
  		puts "#{results.length} variants selected"
  	end

  else
  	 puts "Check variants file and gene symbols file, try --help for further help."
  end
end
