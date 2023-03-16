# CNV-Z

CNV-Z is a collection of Julia- and R-based scripts that we developed to identify and visualize CNVs in NGS data. An algorithm for copy number estimation based on samtools depth was incorporated. The output generated is a .csv file, which is then filtered to call CNVs. A simple filtering algorithm is also provided. The user can choose stringency level to find harmony between sensitivity and specific-ity. To visualize the CNVs, an Rscript is also provided, where sample and gene is selected and a plot generated displaying copy number and z-score, which are plotted for easy visual confirmation of filtering results.

The CNV-Z.jl script takes mapped bam files as input and uses samtools depth to determine read-depth in each position of the target bed file provided by the user. 

The copy number (CN) is estimated by the ratio of normalized read depth of the sample and mean read depth at each position as given below. By default, CN~2 indicates a diploid state, CN< 1.2 indicates deletion, and CN>2.8 indicates amplification.

CN=2*(Normalized read depth (X))/(Mean read depth (μ))

The significance of a CN change is determined using the Z-score. The Z-score calculates the standard deviation of the sample coverage from the mean coverage of the background data set. Here X denotes normalized read depth, µ denotes mean read depth and σ denotes the standard deviation.

Z-score =   (X-µ)/σ 

In the filtering julia script Filter.jl filtering is set as fil-ter([:copynumber,:zscore] => (y,z) -> (y < 1.2 || y > 2.8) && abs(z) >= 2.3,df), i.e. filtering for a Z-score corresponding to a probability of less than approximately 0.01 and a copy number value of either < 1.2 (indi-cating deletion) or > 2.8 (indicating gain). Depending on the number of samples and quality of mapping other filtering strategies can easily be implemented by changing these variables. The filtering script also gener-ates a table in .csv format, which includes the number of hits within each deviating exon and their mean Z-score

For visualization, the included Rscript generates a scatterplot with esti-mated copy number and z-score for every single position included in the target bed file. 

