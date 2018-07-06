library(kohonen)
library(fBasics)
rm(list=ls(all=TRUE))
graphics.off()
set.seed(42)
default_par <- par()
graphics.off()

### PARAMETERS START

### Columns of interest! Comment out the columns you don't wish to include in
### the analysis
columns_to_use <- c(
    #"locus_ref",
    #"contig",
    #"function.",
    "flowers",                   #Flowers
    "suspension_culture",        #SC YE 12h
    "suspension_culture.1",      #SC YE 24h
    "suspension_culture.2",      #SC Ctrol 24h
    "sterile_seedings",          #SS Ctrol
    "sterile_seedings.1",        #SS MJ 12d
    "sterile_seedings.2",        #SS MJ 5d
    "suspension_culture.3",      #SC YE 6h
    "suspension_culture.4",      #SC MJ 6h
    "suspension_culture.5",      #SC MJ 12h
    "suspension_culture.6",      #SC MJ 24h
    "suspension_culture.7",      #SC MJ 0h
    "mature_leaf",               #Mature leaf
    "immature_leaf",             #Immature leaf
    "stem",                      #Stem
    "root",                      #Root
    "hairy_root",                #TDCi HR
    "hairy_root.1",              #RebH F HR
    "wild_type_hairy_root",      #WT HR
    "widl_type_hairy_root",      #WT HR MJ 0h
    "wild_type_hairy_root.1",    #WT HR MJ 24h
    "hairy_roots",               #TDCi HR MJ 0h
    "hairy_roots.1"              #TDCi HR MJ 24h
)

### Number of SOM clusters
tot <- 500

### Whether to generate a PDF or not
generate_pdf <- TRUE

### The PDF file name
pdf_filename <- 'diagnostic_all_tissues.pdf'

### The file name of the total output table
all_som_table_filename <- 'all_som_results.csv'

### The file name of the total output table
subset_som_table_filename <- 'subset_som_results.csv'

### PARAMETERS END

### Import the expression data, which is in a TAB separated .csv file
### stringsAsFactors=FALSE is required as otherwise R tries to optimise its
### memory usage by using awkward to use factors
raw_exp_data <- read.table('c.roseus.transcriptome.csv', sep='\t', header=TRUE,
    stringsAsFactors=FALSE)

### Uses regular expression matching to get the locus ID from the starting
### column
locus_ref <- as.numeric(regmatches(raw_exp_data[,1],
    regexpr("[0-9]+", raw_exp_data[,1], perl=TRUE)))

### Attaches the locus IDs to the table
raw_exp_data <- cbind(locus_ref, raw_exp_data)
rownames(raw_exp_data) <- locus_ref
rm(locus_ref)

### Converts the column selection into a numeric vector
column_subset_indices <- which(colnames(raw_exp_data) %in% columns_to_use)

### Remove expression profiles that have zero expression in more than half of
### samples
zero_count_vec <- apply(raw_exp_data[,column_subset_indices], 1,
    function(x) table(as.numeric(x))['0'])
raw_exp_data <- raw_exp_data[
    which(zero_count_vec < (length(column_subset_indices)/2) |
    is.na(zero_count_vec)),]
rm(zero_count_vec)

### Mean centre and unit variance normalisation for each gene
### Mean centre and unit variance normalisation is not required for each
### experimental tissue as that normalisation has already been carried out
MeanCentreUnitVariance <- function(num_vector){
    if (sd(num_vector)==0) { return(num_vector-mean(num_vector)) }
    else { return( (num_vector-mean(num_vector))/sd(num_vector) ) }
}

normed_data <- t(apply(raw_exp_data[,column_subset_indices], 1,
    MeanCentreUnitVariance))
norm_exp_data <- raw_exp_data
norm_exp_data[,column_subset_indices] <- normed_data
rm(normed_data)

### SOM stuff starts
print('Start SOM')

### Make the ratio of sides the same as the ratio between the first and second
### eigenvalues for the data. This helps the SOM to capture the variation
### present in the data
eigenvalues <- sort(svd(as.matrix(norm_exp_data[,column_subset_indices]))$d,
    decreasing=TRUE)
rat <- eigenvalues[1] / eigenvalues[2]
grid_y <- sqrt(tot / rat)
grid_y <- c(floor(grid_y), ceiling(grid_y))
grid_y <- grid_y[which(grid_y%%2==0)]
grid_x <- floor(rat * grid_y)

### Actually perform the SOM analysis
som.norm_exp_data <- som(as.matrix(norm_exp_data[,column_subset_indices]),
    grid=somgrid(grid_x, grid_y, 'hexagonal'), toroidal=TRUE)

### HIGHLIGHTING GENES OF INTEREST START

### Make all gene traces transparent initially
col_codes <- rep('transparent', length(norm_exp_data$locus_ref))

### Make known loci visible
known_loci <- c(1389, 17400, 10318, 4962, 2176, 8420, 2708, 2522, 2046)
col_codes[which(norm_exp_data$locus_ref %in% known_loci)] <- 'black'

### Using the function column, try and work out if transporters tend to be
### clustered together
transporter_rows <- unique(c(
    grep('efflux', norm_exp_data$function., ignore.case=TRUE),
    grep('pump', norm_exp_data$function., ignore.case=TRUE),
    grep('transport', norm_exp_data$function., ignore.case=TRUE),
    grep('drug', norm_exp_data$function., ignore.case=TRUE),
    grep('pdr', norm_exp_data$function., ignore.case=TRUE)
))
# col_codes[transporter_rows] <- 'black'

### Reads in a table of genes implicated in the monoterpene indole pathway
imp_contigs <- read.table('implicated_genes.csv', sep='\t',
    stringsAsFactors=FALSE, header=TRUE)

# col_codes[which(norm_exp_data$locus_ref %in% imp_contigs$locus_ref)] <- 'black'
#
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$localization=='IPAP cells',]$locus_ref)] <- 'red'
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$localization=='Epidermis',]$locus_ref)] <- 'white'
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$localization=='MLI cells',]$locus_ref)] <- 'blue'
#
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$abbreviation=='T16H',]$locus_ref)] <- rainbow(5)[1]
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$abbreviation=='OMT',]$locus_ref)] <- rainbow(5)[2]
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$abbreviation=='NMT',]$locus_ref)] <- rainbow(5)[3]
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$abbreviation=='D4H',]$locus_ref)] <- rainbow(5)[4]
# col_codes[which(norm_exp_data$locus_ref %in%
#     imp_contigs[imp_contigs$abbreviation=='DAT',]$locus_ref)] <- rainbow(5)[5]

### HIGHLIGHTING GENES OF INTEREST END

### PLOT DIAGNOSTIC PLOTS START

if (generate_pdf) {
    pdf(pdf_filename, onefile=TRUE, paper='a4r', width=11, height=8)
}

### The quality plot is a measure of how similar the gene expression traces are
### to the cluster trace to which they are assigned
plot(som.norm_exp_data, type='quality')
qual.dists <- sapply(split(som.norm_exp_data$distances,
    som.norm_exp_data$unit.classif), mean)
qual.dist.threshold <- quantile(qual.dists)['25%']
add.cluster.boundaries(som.norm_exp_data,
    as.numeric(qual.dists < qual.dist.threshold) + 1)

if (!generate_pdf){x11()}

### The neighbour plot is a measure of how similar the cluster traces are to the
### neighbouring six clusters
nhbrdist <- unit.distances(som.norm_exp_data$grid, som.norm_exp_data$toroidal)
nhbrdist[nhbrdist > 1.05] <- NA
for (i in 2:nrow(nhbrdist)) {
    for (j in 1:(i - 1)) {
        if (!is.na(nhbrdist[i, j]))
        nhbrdist[i, j] <- nhbrdist[j, i] <- dist(
            som.norm_exp_data$codes[c(i, j), ])
    }
}
neigh.dists <- colSums(nhbrdist, na.rm = TRUE)
rm(nhbrdist)
neigh.dist.threshold <- quantile(neigh.dists)['25%']

plot(som.norm_exp_data, type='dist.neighbours', palette.name=greyPalette)
add.cluster.boundaries(som.norm_exp_data,
    as.numeric(neigh.dists < neigh.dist.threshold) + 1)

if (!generate_pdf){x11()}

### Clusters in the bottom 25% of the two distance measures above are regarded
### as 'high quality', in that they represent the genes mapped to them well and
### are similar to their neighbouring clusters, meaning that you be confident
### looking in neighbouring clusters for genes of interest.
### Here, these high quality clusters are coloured green
bgcols <- rep('white',nrow(som.norm_exp_data$codes))
cluster_overlap <- intersect(which(neigh.dists<neigh.dist.threshold),
    which(qual.dists<qual.dist.threshold))
bgcols[cluster_overlap] <- 'green'

### Plot the cluster traces
plot(som.norm_exp_data,type='codes', bgcol=bgcols, codeRendering='lines')

if (!generate_pdf){x11()}

### Plot where the genes of interest (as defined in the HIGHLIGHTING GENES OF
### INTEREST section) are located in the map
plot(som.norm_exp_data,type='mapping', pchs=1, cex=0.8, col=col_codes,
    bgcol=bgcols, keepMargins=TRUE)

if (!generate_pdf){x11()}

### This section of code can be used to find out cluster numbers
# while (TRUE) {
# plot(som.norm_exp_data,type='mapping', pchs=1, cex=0.8, col=col_codes,
#     bgcol=bgcols, keepMargins=TRUE)
# cell <- identify(som.norm_exp_data, n=1)
# print(cell)
# }

### Count the number of transporters (as defined in the HIGHLIGHTING GENES OF
#### INTEREST section) and see where they are located in the map
transporter_sums <- sapply(split(as.numeric(1:nrow(norm_exp_data) %in%
    transporter_rows), som.norm_exp_data$unit.classif), sum)
plot(som.norm_exp_data,type='property',
     property=transporter_sums,
     main='Transporter Count')
add.cluster.boundaries(som.norm_exp_data, as.numeric(
    1:max(som.norm_exp_data$unit.classif) %in% cluster_overlap) + 1)

if (generate_pdf){dev.off()}

### A function to write tables of the gene cluster assignments
write_som_table <- function(output_filename, rows_to_write)
{
    table_subset <- norm_exp_data[rows_to_write, 1:3]

    transporter_y_n <- rep('n', nrow(norm_exp_data))
    transporter_y_n[transporter_rows] <- 'y'

    table_subset <- cbind(
        table_subset,
        neuron_number=som.norm_exp_data$unit.classif[rows_to_write],
        imp_contig=c('n','y')[as.numeric(
            table_subset$locus_ref %in% imp_contigs$locus_ref) + 1],
        transporter=transporter_y_n[rows_to_write]
    )

    table_subset <- cbind(table_subset, raw_exp_data[
        match(table_subset$locus_ref, raw_exp_data$locus_ref),
        which(colnames(raw_exp_data) %in% c('mature_leaf', 'immature_leaf'))])

    write.table(table_subset, file=output_filename, sep='\t', row.names=FALSE)
}
write_som_table(all_som_table_filename, 1:nrow(norm_exp_data))

### This section of code allows you to select clusters of interest to allow for
### further investigation. It asks on the command line for the number of
### clusters, then will open an interactive plot allowing you to select the
### clusters you would like more information for
print('How many interesting clusters?')
num_int_clust <- as.numeric(readline())
if (is.na(num_int_clust) | num_int_clust==0) {graphics.off()} else {
    if (!generate_pdf){x11()}
    plot(som.norm_exp_data,type='codes', bgcol=bgcols, codeRendering='lines')
    interesting_clusters <- identify(som.norm_exp_data,n=num_int_clust)
    graphics.off()

    bgcols <- rep('white',nrow(som.norm_exp_data$codes))
    bgcols[interesting_clusters] <- rep('green', length(interesting_clusters))

    ### Plot a codes plot with the selected clusters highlighted
    plot(som.norm_exp_data,type='codes', bgcol=bgcols)

    x11()

    ### Plot a line graph of the cluster traces for the selected clusters
    maximum_expr = max(som.norm_exp_data$codes[interesting_clusters,])
    minimum_expr = min(som.norm_exp_data$codes[interesting_clusters,])
    par(mar=c(7.1,4.1,4.1,2.1))
    plot(0, type='n', xlim=c(1,ncol(som.norm_exp_data$codes)),
        ylim=c(minimum_expr, maximum_expr), main='Interesting Clusters',
        xaxt='n', xlab='', ylab='Expression Level')
    axis(1, at=1:ncol(som.norm_exp_data$codes),
        colnames(som.norm_exp_data$codes), las=2, cex.axis=0.7)

    for (cluster_id in interesting_clusters) {
        lines(1:ncol(som.norm_exp_data$codes),
            som.norm_exp_data$codes[cluster_id,])
    }
    par(default_par$mar)

    ### Write a subset of the table corresponding to the clusters of interest
    write_som_table(subset_som_table_filename,
        which(som.norm_exp_data$unit.classif %in% interesting_clusters))
}
