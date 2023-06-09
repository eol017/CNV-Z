#Julia script to analyse ICR96 dataset Fx-$id.csv file by filtering
#Usage julia FilterICR96.jl 

using CSV, DataFrames, Statistics

for id in ["17296", "17297", "17298", "17299", "17300", "17301", "17302", "17303", "17304", "17305", "17306", "17307", "17308", "17309", "17311", "17312", "17313", "17314", "17315", "17316", "17317", "17318", "17319", "17320", "17321", "17322", "17323", "17324", "17326", "17327", "17328", "17329", "17330", "17331", "17332", "17333", "17334", "17335", "17336", "17337", "17338", "17339", "17340", "17341", "17342", "17343", "17356", "17358", "17359", "17360", "17361", "17362", "17363", "17364", "17365", "17366", "17367", "17368", "17369", "17370", "17371", "17372", "17373", "17374", "17375", "17376", "17377", "17378", "17379", "17380", "17381", "17382", "17383", "17384", "17385", "17386", "17387", "17388", "17389", "17390", "17391", "17392", "17393", "17394", "17395", "17396", "17397", "17398", "17399", "17400", "17401"]
        
	df = DataFrame(CSV.File("Fx-$id.csv"))
        CSV.write("Fx-$id.filter2.3.csv",filter([:copynumber,:zscore] => (y,z) -> (y < 1.2 || y > 2.8) && abs(z) >= 2.3,df))

	dft = DataFrame(CSV.File("hglft_ICR96target.bed",header=false, delim='\t'))
	len = dft[:,3]-dft[:,2]
	dft.len = len
	idx = collect(1:nrow(dft))
	dft.idx = idx
	hits = zeros(Int64,nrow(dft))
	dft.hits = hits
	zscore = zeros(Float64,nrow(dft))
	dft.zscore = zscore 
	
	f1 = open("Fx-$id.filter2.3.csv")
	lines = readlines(f1)
	println("\nSampleID :",id,"\nlength of lines: ",length(lines),"\nnrow(dft): ",nrow(dft))

	for i in 2:length(lines)
	    line = split(lines[i],",")
	    chrom = line[1]
	    pos = parse(Int64,line[2])
	    for j in 1:nrow(dft)
	        if dft[j,1]==chrom && dft[j,2]<pos && dft[j,3]>=pos
           		dft[j,8] = dft[j,8]+1
        	end
    	    end
	end

	for i in 1:nrow(dft)
	    Z = Float64[]
	    for j in  2:length(lines)
		line = split(lines[j],",")
	    	chrom = line[1]
	    	pos = parse(Int64,line[2])
		score = parse(Float64,line[8])
		if dft[i,1]==chrom && dft[i,2]<pos && dft[i,3]>=pos
		    push!(Z, score)
		end
	    end
	    dft[i,:zscore] = mean(Z)
	end     
	
	dff = filter(:hits => x -> x>0,dft)
	println(dff)

	CSV.write("Fx-$id.res.2.3.csv",dff)
	close(f1)
end
