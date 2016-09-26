class SymbolChecker
	require 'net/http'
	require 'uri'
	require 'rubygems'
	require 'json'
	require 'pp'
	
	def check_gene_symbols(gene_symbols_array)
		server="http://rest.genenames.org"
		confirmed_gene_symbols = Array.new
		missing_gene_symbols = Array.new

		gene_symbols_array.each do |this_symbol|
			get_path = "/fetch/symbol/#{this_symbol}"
			url = URI.parse(server)
			http = Net::HTTP.new(url.host, url.port)
			 
			request = Net::HTTP::Get.new(get_path, {'Accept' => 'application/json'})
			response = http.request(request)
			 
			if response.code != "200"
			 puts "Invalid response: #{response.code}"
			 puts response.body
			 exit
			end
			
			result = JSON.parse(response.body)
			
			if result.fetch("response").fetch("numFound") > 0
				confirmed_gene_symbols.push("#{this_symbol}")
			else
				missing_gene_symbols.push("#{this_symbol}")
			end
			#pp result
		end
		
		results[confirmed_gene_symbols, missing_gene_symbols]
		return results
	end

end
