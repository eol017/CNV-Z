#Usage: julia -t 8 CNV-Z.jl <sampleID> [<sampleID> <sampleID> ... ]
#Adjust the number of threads, and paths to targets and bamfiles below if necessary

using CSV, DataFrames, Statistics, Base.Threads

println(ARGS)

Threads.@threads for i in ARGS
	open("mdepth-$i.txt","w") do io
		write(io,read(`samtools depth -a -b hglft_ICR96target.bed $i.sorted.bam`))
	end
end

for i in ARGS
	open("mdepth-$i.h.txt","w") do io
		write(io,read(`cat header.txt mdepth-$i.txt`))
	end
end

Files = []

for id in ARGS
       	   push!(Files, string("mdepth-",id,".h.txt"))
      end

Fx = [DataFrame(CSV.File(f)) for f in Files]

function normgene(chr,pos)::Bool
	    ichr = chr != "chrX"
	    ipos = pos > 0	
           #ichr = chr == "chr1"			#MTOR
           #ipos = pos >= 11107485 && pos <=11259409	#MTOR
           ichr && ipos
       end

for i in 1:length(Files)
           Fx[i].prop = Fx[i][!,:Depth]/sum(filter([:Chr, :Position] => normgene,Fx[i])[!,:Depth])
	end

M = Matrix{Float64}(undef, size(Fx[1],1), length(Fx))

for i in 1:length(Fx)
           M[:,i] = Fx[i][:,:prop]
       end

for i in 1:size(M,2)
	A = Float64[]
	for j in 1:size(M,1)
           push!(A, mean(M[j,:]))
       	end	
	Fx[i].mean = A
    end

for i in 1:size(M,2)
	B = Float64[]
	for j in 1:size(M,1)
           push!(B, std(M[j,:]))
       	end	
	Fx[i].std = B
    end

for i in 1:size(M,2)
	C = Float64[]
        for j in 1:size(M,1)
           push!(C, Fx[i][j,:Depth]/Fx[i][j,:prop]*Fx[i][j,:mean])
        end
        Fx[i].expDepth = C
    end

for i in 1:size(M,2)
	D = Float64[]
        for j in 1:size(M,1)
           push!(D, (Fx[i][j,:prop]-Fx[i][j,:mean])/Fx[i][j,:std])
        end
        Fx[i].zscore = D
    end

for i in 1:size(M,2)
	E = Float64[]
        for j in 1:size(M,1)
           push!(E, 2*Fx[i][j,:prop]/Fx[i][j,:mean])
        end
        Fx[i].copynumber = E
    end

#println(Fx)

for i in 1:length(ARGS)
           name = ARGS[i]
           CSV.write("Fx-$name.csv", bufsize = 2^24, Fx[i])
       end


