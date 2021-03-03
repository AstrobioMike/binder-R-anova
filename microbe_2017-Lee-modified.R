### ENVIRONMENT ###

library("phyloseq")
library("tidyr")
library("ggplot2")
library("viridis")
library("DESeq2")
library("dendextend")
library("ggrepel")
library("vegan")

### READING IN DATA ###

otu_counts_tab <- read.table("data/combined_otu_counts_tab.txt", header=TRUE, row.names=1, sep="\t")
head(otu_counts_tab)
dim(otu_counts_tab)
colnames(otu_counts_tab)
head(rownames(otu_counts_tab))

taxonomy_tab <- read.table("data/combined_taxonomy_tab.txt", header=TRUE, row.names=1, sep="\t", na.strings="<NA>")
head(taxonomy_tab)
dim(taxonomy_tab)
colnames(taxonomy_tab)

sample_info_tab <- read.table("data/sample_info_tab.txt", header=T, check.names=F, row.names=1, sep="\t")
head(sample_info_tab)
dim(sample_info_tab)
colnames(sample_info_tab)

### TAXONOMIC SUMMARIES ###

counts_phy <- otu_table(otu_counts_tab, taxa_are_rows=TRUE)
tax_phy <- tax_table(as.matrix(taxonomy_tab))
sample_phy <- sample_data(sample_info_tab)
aqua_phy <- phyloseq(counts_phy, tax_phy, sample_phy)


















phyla_counts_tab <- otu_table(tax_glom(aqua_phy, taxrank="Phylum")) # this is summarizing our counts table at the phylum level, it will add together all counts of the same phylum in the same sample
rownames(phyla_counts_tab) <- as.vector(tax_table(tax_glom(aqua_phy, taxrank = "Phylum"))[,2])


  # we also have to account for sequences that weren't assigned any taxonomy even at the phylum level
    # these came into R as 'NAs' in the taxonomy table, but their counts are still in the count table
      # so we can get that value for each sample by substracting the column sums of this new table (that has everything that had a phylum assigned to it) from the column sums of the starting count table (that has all representative sequences)
unclassified_counts <- colSums(counts_phy) - colSums(phyla_counts_tab)
  # and we'll add this row to our phylum count table:
phyla_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_counts)

colSums(otu_counts_tab) - colSums(phyla_counts_tab) # making sure we have everything (if we missed any the difference in that sample wouldn't be 0)

  # now we'll remove the Proteobacteria, so we can next add them back in broken down by class
phyla_counts_no_proteo_tab <- phyla_counts_tab[!row.names(phyla_counts_tab) %in% "Proteobacteria", , drop=F]

  # making count table broken down by class (all classes at this point)
class_counts_tab <- data.frame(otu_table(tax_glom(aqua_phy, taxrank="Class")))

  # making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(aqua_phy, taxrank="Class"))

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)

  # making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"])

  # changing the row names like above so that they correspond to the taxonomy,
    # rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$Class)

  # making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ]

  # there are also possibly some some sequences that were resolved to the level
    # of Proteobacteria, but not any further, and therefore would be missing from
    # our class table
    # we can find the sum of them by subtracting the proteo class count table
    # from just the Proteobacteria row from the original phylum-level count table
unclassified_proteos <- phyla_counts_tab[row.names(phyla_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

#   # and adding that to our proteo count table:
# proteo_counts_tab <- rbind(proteo_counts_tab, "Unclassified_Proteos"=unclassified_proteos)

  # now combining the tables:
major_taxa_counts_tab <- rbind(phyla_counts_no_proteo_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria" = unclassified_proteos)

  # and making sure again that we're not missing anything:
colSums(otu_counts_tab) - colSums(major_taxa_counts_tab) # great, all 0s

    ## right now these are all counts and not "normalized" across samples
      # what we mean by that here, is some samples may have a lot of sequences, and some a lot less, so comparing the raw counts like this could be misleading
        # let's take a peek at the total sequences in each:
colSums(major_taxa_counts_tab) # some have a lot, like T2_Tank_2 with 213,421 sequences, compared to say T1_Tank_1 with only 55,306 sequences
  # so here since we are just summarizing taxonomy proportions, we can turn these into percents of the total for each sample, and then compare those percents

  # first let's peek at the starting table again real quick:
head(major_taxa_counts_tab)

  # making proportions (again, don't get bogged down in the code being confusing, it has to be that way at first! and again, i don't have most of this memorized, i google and figure out what i need to do at the time, or i look at old scripts I've written)
major_taxa_proportions_tab <- round(apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100), 2)

  # now let's peek at the proportions tab
head(major_taxa_proportions_tab) # cool
    # but note, they don't all add exactly to 100 because i used the round() function above to keep the table clean
colSums(major_taxa_proportions_tab) # this is negligible, and could be avoided just by not rounding...
major_taxa_proportions_tab_not_rounded <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
    # just now they are not as easy to visually assess
head(major_taxa_proportions_tab_not_rounded)

  # if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
  # we see there are currently 36 rows, which might be a little busy for a summary figure (you'll see in a second colors are not friendly when you need many different shades)
    # many of these taxa make up a very small percentage though, so we're going to filter some of the lower abundance bugs out
      # this is a completely arbitrary decision solely to ease visualization and intepretation, entirely up to your data and you
        # here, we'll only keep rows (taxa) that make up greater than 3% in any individual sample
temp_filtered_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 3, ])
  # checking how many we have that were above this threshold
dim(temp_filtered_major_taxa_proportions_tab) # now we have 10, much more manageable for an overview figure

  # though each of the filtered taxa made up less than 3% alone in each individual sample, together they can be more than 3% in an individual sample
    # so we're going to add a row called "Other" that keeps track of how much we filtered out (which will also keep our totals at ~100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filtered_major_taxa_proportions_tab)
    # let's look at how much this "Other" category makes up in our samples:
sort(filtered_proportions, decreasing=TRUE) # the most is about 5% in sample T2_Tank_10, not too big of a deal, this is up to the researcher too though

  # here we're combining our filtered table with this "Other" row to keep track of those we filtered out
filtered_major_taxa_proportions_tab <- rbind(temp_filtered_major_taxa_proportions_tab, "Other"=filtered_proportions)

  # let's look at our table
dim(filtered_major_taxa_proportions_tab) # now we have 11 rows because we added the "Other"
filtered_major_taxa_proportions_tab # and we can see the whole table like this

  # for convencience, we can write this table out of R, to a file, so we have it saved for the future and in case we wanted to send it to someone or do stuff with it in excel or whatever
write.table(filtered_major_taxa_proportions_tab, "major_taxa_proportions_tab.tsv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

  # now if we check our files again:
list.files() # we see that the "major_taxa_proportions_tab.txt" is there too

    ## now that we have our summary table, we can finally make some visualizations with it
# making copy to manipulate for plottingm
filtered_major_taxa_proportions_tab_for_plotting <- filtered_major_taxa_proportions_tab

  # and add a column of the taxa names so that it is within the table, rather than just as row names
filtered_major_taxa_proportions_tab_for_plotting$Major_Taxa <- row.names(filtered_major_taxa_proportions_tab_for_plotting)

  # now we'll transform the table into "long" format, also known as "narrow", format (you can read a bit about this here: https://en.wikipedia.org/wiki/Wide_and_narrow_data#Narrow), but we'll also look at the change ourselves
filtered_major_taxa_proportions_tab_for_plotting.g <- gather(filtered_major_taxa_proportions_tab_for_plotting, Sample, Proportion, -Major_Taxa)

  # take a look at the old table ("wide" format) and compare it with the old one
head(filtered_major_taxa_proportions_tab_for_plotting) # this is a two dimensional table like we are used to seeing with rows, columns, and cells connected to the row and column they fall within
head(filtered_major_taxa_proportions_tab_for_plotting.g) # this is "long" format, where it's not as concise or intuitive, but it can hold more different types of information
  # manipulating tables like this is something you may need to do frequently in R

  # one thing we're going to mention here quickly are that different "classes" or types of objects in R have different properties (e.g. a number is different than text)
    # we can check the classes of data in our table here like this
str(filtered_major_taxa_proportions_tab_for_plotting.g)
      # which tells us the whole thing is a data frame, and in it, there are 3 columns "Major_Taxa", "Sample", and "Proportion"
        # right now Major_Taxa is what we're going to use to print, and by default it would just order them alphabetically
          # but I want "Other" and "Unclassified" to be on the bottom, so I'm going to change it from "character" to "factor", and give it a special order
            # "Factors" in R are like character data in the sense that they just hold text, but you can also give them an order
unique(sort(filtered_major_taxa_proportions_tab_for_plotting.g$Major_Taxa)) # here's what the unique values in this column look like before
filtered_major_taxa_proportions_tab_for_plotting.g$Major_Taxa <- factor(filtered_major_taxa_proportions_tab_for_plotting.g$Major_Taxa, levels=c("Alphaproteobacteria", "Bacteroidetes", "Betaproteobacteria", "Cyanobacteria", "Deinococcus-Thermus", "Firmicutes", "Gammaproteobacteria", "Planctomycetes", "Verrucomicrobia", "Other", "Unclassified"))

unique(sort(filtered_major_taxa_proportions_tab_for_plotting.g$Major_Taxa)) # and here's what it looks like after changing it a factor
str(filtered_major_taxa_proportions_tab_for_plotting.g) # now it is called a "factor" with 9 levels

  # we need to do the same thing for the samples so that they appear in the appropriate order, as right now the order goes Tank_1 to Tank_10, then Tank_2 (computer numeric order is not how we're used to thinking at first)
unique(sort(filtered_major_taxa_proportions_tab_for_plotting.g$Sample))
        # instead of writing out all 32 samples, i can provide them in a variable too, or part of a table
          # remember our sample_info_table
sample_info_tab
            # that has them in the order I want, so i can use the row names of that like this:
rownames(sample_info_tab) # just running that prints them to the screen just like the c() funciton did above, so i'm going to put that in the code here
filtered_major_taxa_proportions_tab_for_plotting.g$Sample <- factor(filtered_major_taxa_proportions_tab_for_plotting.g$Sample, levels=rownames(sample_info_tab))
unique(sort(filtered_major_taxa_proportions_tab_for_plotting.g$Sample)) # and now they're in the order I want! (where Tank_1 is followed by Tank_2, rather than Tank_10)

  # and now we're ready to make some summary figures with our wonderfully constructed table
    ## a good color scheme can be hard to find, i included the viridis package here because sometimes it's been really helpful for me, though this is not demonstrated in all of the following :/
      ### and IMPORTANTLY, again, don't get bogged down in the code! same deal, i don't have this memorized either and i look it up each time i write something :)
# stacked bar charts for each taxon per sample:
ggplot(filtered_major_taxa_proportions_tab_for_plotting.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_viridis(discrete=TRUE) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

  # here is with default colors:
ggplot(filtered_major_taxa_proportions_tab_for_plotting.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

  ## see? also ugly, (and not color-blind friendly, which the viridis package is)

    # looking at that, nothing too huge jumps out at me. The EM1 and the Lambda look very diffferent from all the other samples, which is good
      # but it can be hard to make out from this image which color those bars belong to
        # we can look at the original table like this:
filtered_major_taxa_proportions_tab # and see that EM1 has 79% Firmicutes and Lambda has 82% Gammaproteobacteria
    # we could also look at just a part of the table we are interested in - see here for more on subsetting tables like this: https://astrobiomike.github.io/R/basics#subsetting-tables
filtered_major_taxa_proportions_tab[, c("T3_EM1", "T3_Lambda")]

      # we can also specify the colors manually to try to make it a little better, but this is hard and usually doesn't help unless you only need a few:
ggplot(filtered_major_taxa_proportions_tab_for_plotting.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values=c("blue", "blueviolet", "darkgreen", "dodgerblue", "deeppink", "lightsteelblue", "seagreen2", "brown1", "cadetblue4", "lightyellow3", "black")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples") #### note here, if you get a weird Error message about "non-finite location..." or something, run "dev.off()", just like that, then try the plot again (i don't know why that happens but it's a bug somewhere)
          # this is still ugly, but you can at least work on your own color scheme by using the color names listed here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
            # if it's on that sheet then entering the name of the color in quotes like i did above will work
    # now we can at least tell a couple things from the figure:
      # T3_Tank_6 has a much greater proportion of "Unclassified" sequences than the other samples
      # T3_EM1 is dominated by Firmicutes
      # T3_Lambda is dominated by Gammaproteobacteria

  #### you can also make circle plots of all if you'd like:
ggplot(filtered_major_taxa_proportions_tab_for_plotting.g, aes(x="", y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=1, stat="identity") +
  scale_fill_manual(values=c("blue", "blueviolet", "darkgreen", "dodgerblue", "deeppink", "lightsteelblue", "seagreen2", "brown1", "cadetblue4", "lightyellow3", "black")) +
  coord_polar(theta = "y", direction = -1) +
  facet_wrap(~Sample) +
  theme_void() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position="right") +
  guides(fill=guide_legend(title="Major Taxa"))

  # or individually with something like this (this one I'm specifying sample "T1_Tank_1"):
ggplot(filtered_major_taxa_proportions_tab_for_plotting.g[filtered_major_taxa_proportions_tab_for_plotting.g$Sample == "T1_Tank_1", ], aes(x="T1_Tank_1", y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=1, stat="identity") +
  scale_fill_manual(values=c("blue", "blueviolet", "darkgreen", "dodgerblue", "deeppink", "lightsteelblue", "seagreen2", "brown1", "cadetblue4", "lightyellow3", "black")) +
  coord_polar(theta="y", direction = -1) +
  theme_void() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position="right") +
  guides(fill=guide_legend(title="Major Taxa")) +
  ggtitle("T1 Tank 1") + theme(plot.title = element_text(hjust=0.5))
        ## to modify plots like this, i would search google for whatever i want to do
          ## for instance, imagine we want to add percents to the pie wedges, by searching google for: "r ggplot add percents to pie chart", the first page that comes up for me is this one: https://stackoverflow.com/questions/41338757/adding-percentage-labels-on-pie-chart-in-r
            ## and on that page there is an example of how to add these, so i would then look at that and modify my code until it worked. If i got stuck, I'd find another page (this is really all we bioinformaticians do all day - please don't tell anyone that though)


  # or a boxplot type like so, thought it's maybe a little busy for this many samples:

ggplot(filtered_major_taxa_proportions_tab_for_plotting.g, aes(Sample, Proportion)) +
    geom_boxplot(fill=NA, outlier.color=NA) +
    geom_jitter(aes(color=Major_Taxa), size=2, width=0.15, height=0) +
    scale_color_manual(values=c("blue", "blueviolet", "darkgreen", "dodgerblue", "deeppink", "lightsteelblue", "seagreen2", "brown1", "cadetblue4", "lightyellow3", "black")) +
    theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
    labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  # those are just a few examples of how you can visualize taxonomic summaries of samples, but again, google is your friend, look for things and try whatever you want, and manipulate the above code to try different things


#### BETA-DIVERSITY ####
  # "beta-diversity" involves comparing samples to other samples, an individual sample cannot have a beta-diversity metric by itself, these metrics require at least 2 samples because they only exist for one sample when it is compared to something else
    # there is some more information about this here: https://astrobiomike.github.io/amplicon/workflow_ex#beta-diversity
      # and we will be following the code there so refer to that for further explanation of the following

  # normalizing for sample depth (like we did for the taxonomic summaries above by calculating proportions (percents), but now that we're going to be using these metrics for beta-diversity, there are better ways to normalize)
    # we're using something called variance stabilizing transformation - this is pretty far in the weeds, and you certainly don't have to dive all the way into it, I did once, then forgot most and trust the statisticians that tell me to use this
      # but there is some more info and links to some literature here in the "Normalizing" section: https://astrobiomike.github.io/amplicon/workflow_ex#beta-diversity
  # first we need to make a DESeq2 object, and we are going back to using our full count table from the start (not summarized by taxonomy, and not proportions)
deseq_counts <- DESeqDataSetFromMatrix(otu_counts_tab, colData = sample_info_tab, design = ~Tank) # we have to include the colData (which is our sample info table) and a "design" argument because they are required by deseq2, but we aren't doing anything that uses that information here so it doesn't matter what we put as the "design"
  # and this is how we run the normalization
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

  # and this is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
head(vst_trans_count_tab) # note these values aren't interpretable the way the counts or proportions were, because they can be negative now, they are how most useful this way for comparing samples like we're going to do

  # and calculating our Euclidean distance matrix (a distance matrix in this context is something that relates every individual sample to every other individual sample)
euc_dist <- dist(t(vst_trans_count_tab))
    # it's not easy to see the whole thing here, but if you print this object out and look at the top you'll see something like this:
      # where T1_Tank_1 has a distance of 213.0185 to T1_Tank_2, but a distance of 193.0769 to T1_Tank_3
        # this really means in this sort of "comparison-space" of these samples (that is built based on the sequences we recovered from their microbial communities), sample T1_Tank_1 is more similar (closer to) T1_Tank_3 than it is to T1_Tank_2
#             T1_Tank_1 T1_Tank_2 ...
# T1_Tank_2   213.0185
# T1_Tank_3   193.0769  206.2715
# ...
        # this has converted the counts of sequence abundances we had for each sample into a, well, "distance" between each pair of two samples
          # this concept should absolutely be confusing at first, but will make some more sense when we look at things visualized

  ## hiearchical clustering
    # one way to visualize these distances is through what's known as hierarchical clustering, it's too far off the rails here to go into that
      # but in general, this just means grouping the closest two samples together into one group, then grouping that group with the next closest samples and so on
        # here's the code to do it
euc_clust <- hclust(euc_dist, method="ward.D2")
    # and here's how we can plot it:
plot(euc_clust)
        # so quick notes on how to interpret this:
          # - the groupings of samples matter and the point at which samples are connected matter (the horizontal lines connecting the samples beneath them)
            # - starting from the very top, you can break things up however you'd like, imagine moving a horizonal line downwards from the 450 at the top
              # - at first there is nothing, then there are two vertical lines going down (2 groups)
                # so we can say, for instance, there are 2 main groups:
                  # 1) the group with T3_EM1 and T3_Lambda
                  # 2) everything else
                    # this tells us that T3_EM1 and the T3_Lambda samples are more different from all the other "real" samples, than any of the "real" samples are to other "real" samples
                      # continuing to look at T3_EM1 and T3_Lambda, the point at which the horizontal line that joins them sits is their distance from each other: which is a little bit under 300 (this is a "euclidean distance", but it doesn't really have units)
                        # let's look at this in the distance matrix to confirm (side note - I didn't know I needed to use data.frame(as.matrix()) like this in here, i got there through trial and error just trying to get the info i wanted, that's normal too; try it other ways and see if it works or not; when things like this happen, sometimes i think about why and sometimes i just move on when i know it's working how i want, just depends on your internal, ever-shifting balance between "learning things" and "getting things done"):
data.frame(as.matrix(euc_dist))["T3_EM1", ] # this shows T3_EM1's distance to all other samples (including itself which is 0)
data.frame(as.matrix(euc_dist))["T3_EM1", "T3_Lambda", drop=FALSE] # this shows just T3_EM1's distance to T3_Lambda
                      # the two lowest samples are the most similar based on our euclidean distance here
                        # looking at our plot again, the two lowest hanging samples are T2_Tank_8 and T2_Tank_10, which appear to be connected at a little over 100
data.frame(as.matrix(euc_dist))["T2_Tank_8", ] # T2_Tank_8's distances to all others
data.frame(as.matrix(euc_dist))["T2_Tank_8", "T2_Tank_10", drop=FALSE] # T2_Tank_8's distance to T2_Tank_10: 105.5607, this is the smallest distance between any combination of samples (meaning they are the most similar in the entire dataset)

                # hopefully that helps a little on how to interpret this, but we can chat about it more :)
    # now, onto making this more useful, we can color the sample names (known as "leaves" in a dendrogram, the places where two samples are joined are called "nodes")
  # we could use plot() on the euc_clust object like we did above, but i like to change them to dendrograms (dendrogram "objects" in R) for two reasons:
    # 1) it's easier to color the dendrogram plot by groups
    # 2) if you want you can rotate clusters with the rotate() function of the dendextend package

  # so here let's make one that colors the "leaves" (samples) by which timepoint they come from
euc_dend_color_by_timepoint <- as.dendrogram(euc_clust, hang=0.1)
  # making a sample info table we can use to create our needed color vectors
sample_info_tab_with_color_groups <- sample_info_tab
  # looking at this table again first, there are no columns specifying any colors yet
sample_info_tab_with_color_groups
colnames(sample_info_tab_with_color_groups)
  # giving a color depending upon which timepoint a sample is from (this will be confusing at first, we can chat about it, but it just takes some time playing with it to eventually realize this is one of the things that makes R awesome; it's covered some here: https://astrobiomike.github.io/R/basics#subsetting-by-conditional-statements and here: https://astrobiomike.github.io/R/more_indexing)
sample_info_tab_with_color_groups$color_by_timepoint[sample_info_tab_with_color_groups$Time_point == "T1"] <- "darkgreen"
sample_info_tab_with_color_groups$color_by_timepoint[sample_info_tab_with_color_groups$Time_point == "T2"] <- "blue"
sample_info_tab_with_color_groups$color_by_timepoint[sample_info_tab_with_color_groups$Time_point == "T3"] <- "red"

  # now if we look at this table again, we see there is a new column called "color_by_timepoint":
colnames(sample_info_tab_with_color_groups)
sample_info_tab_with_color_groups # note that you could also add these columns in excel or something before reading into R, but that leaves a little more room for human error

  # and we can now use this to color our dendrogram leaves based on timepoint of the samples
labels_colors(euc_dend_color_by_timepoint) <- sample_info_tab_with_color_groups$color_by_timepoint[order.dendrogram(euc_dend_color_by_timepoint)]

plot(euc_dend_color_by_timepoint)
  # and we can apply a bit of auto-sorting which may help visually with the dendsort function
library(dendsort)

euc_dend_color_by_timepoint_sorted <- dendsort(euc_dend_color_by_timepoint, isReverse=TRUE)
plot(euc_dend_color_by_timepoint_sorted)

  # and we can add a label to our y-axis like so:
plot(euc_dend_color_by_timepoint_sorted, ylab="VST Euc. dist.")
      # ok, cool, things are colored as they should be, but doesn't seem to make things all that clear as the time points seem to mix together pretty well regarding their distances to each other
  # let's try coloring by Tank
    # but for this we need 10 colors (there are better/more efficient ways to do this, but since it's a small enough amount to do by hand i'm just going to do it the same way rather than spend time figuring one of those out)
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_1"] <- "darkgreen"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_2"] <- "blue"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_3"] <- "red"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_4"] <- "blueviolet"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_5"] <- "dodgerblue"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_6"] <- "deeppink"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_7"] <- "steelblue"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_8"] <- "seagreen"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_9"] <- "chocolate"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Tank_10"] <- "brown2"
  # and making the Lambda and EM1 samples black
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "EM1"] <- "black"
sample_info_tab_with_color_groups$color_by_tank[sample_info_tab_with_color_groups$Tank == "Lambda"] <- "black"

  # now if we look at this table yet again, we see there is another new column, now called "color_by_tank":
colnames(sample_info_tab_with_color_groups)
sample_info_tab_with_color_groups

  # and we can now use this to color our dendrogram leaves based on the tank each sample comes from
euc_dend_color_by_tank <- as.dendrogram(euc_clust, hang=0.1)
labels_colors(euc_dend_color_by_tank) <- sample_info_tab_with_color_groups$color_by_tank[order.dendrogram(euc_dend_color_by_tank)]

  # auto-sorting
euc_dend_color_by_tank_sorted <- dendsort(euc_dend_color_by_tank, isReverse = TRUE)

  # now let's plot it:
plot(euc_dend_color_by_tank_sorted, ylab="VST Euc. dist.")
        # let's look around a bit:
          # - all of Tank_1's and Tank_9's are together in one cluster (the six of them form a group that no other samples mix in with)
          # - all of Tank_3's are together in their own cluster
          # - all Tank_5's are pretty close together, thrown off a little by a T3_Tank_7, and it seems all Tank_7 timepoints are a little all over the place

        # this is just for visualization and exploration of the data, but down further we'll look at one way to test if any of these we see visually are actually statistically significant

    ## you can add colors here however you want of course, so adding a column to the table that colors things by say Plants or Media would be fine too, maybe give that a shot

  ## ordination
    # another way to visualize this is with an ordiation, see https://astrobiomike.github.io/amplicon/workflow_ex#beta-diversity and the "ordination" section for a little discussion on this
      # overall, ordination is just another way of organizing our samples in space, and then visualizing their locations with respect to each other (it's done differently than hierarchical clustering like we did above, they each provide ways to try to gain insight into your data, they are complementary, not competitors)
        # we will use phyloseq to do this, but we need to make a new phyloseq object using the transformed table we made above:
  # making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <- phyloseq(vst_count_phy, sample_data(sample_info_tab_with_color_groups))

  # generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples (don't worry too much about this right now if this is confusing)

  # here's looking with coloring by timepoints
plot_ordination(vst_physeq, vst_pcoa, color="color_by_timepoint") +
  labs(col="Time_point") + geom_point(size=1) +
  geom_text(aes(label=rownames(sample_info_tab_with_color_groups))) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups$color_by_timepoint[order(sample_info_tab_with_color_groups$Time_point)])) +
  theme(legend.position="none")
      # okay, this is a mess, but we can see something. Like the hiearchical clustering we saw above, the Lambda and EM1 samples are much more different than everything else (which is good and what we'd expect)
        # in this type of figure, have these two relatively "far away samples" included is lowering the overall resolution we have to separate out the other samples, so for this part we're going to drop those samples
  # to do that, we're going to make a table with out Lambda and EM1, and process it through to this point:
colnames(otu_counts_tab) # all column names, we can see the T3_EM1 and T3_Lambda are the last two
  # here is one way to subset our table to only include the real samples, in the subsetting brackets we provide the rows we want first (all here so left blank), then the columns we want (1 through 30): again, see the links above for more on subsetting in R
otu_counts_samples_only_tab <- otu_counts_tab[ , 1:30]
  # and subsampling the sample info table to only include "real" samples
sample_info_tab_with_color_groups_samples_only <- sample_info_tab_with_color_groups[1:30, ]
deseq_counts_samples_only <- DESeqDataSetFromMatrix(otu_counts_samples_only_tab, colData = sample_info_tab_with_color_groups_samples_only, design = ~Tank) # we have to include the colData (which is our sample info table) and a "design" argument because they are required by deseq2, but we aren't doing anything that uses that information here so it doesn't matter what we put as the "design"
  # running the normalization on subset table
deseq_counts_samples_only_vst <- varianceStabilizingTransformation(deseq_counts_samples_only)

# pulling out transformed subset table
vst_trans_count_sample_only_tab <- assay(deseq_counts_samples_only_vst)
head(vst_trans_count_sample_only_tab)

  # making phyloseq object of transformed, subst table
vst_count_samples_only_phy <- otu_table(vst_trans_count_sample_only_tab, taxa_are_rows=T)
  # needed to fix one thing about the Tank column
sample_info_tab_with_color_groups_samples_only$Tank <- factor(as.character(sample_info_tab_with_color_groups_samples_only$Tank), levels=unique(as.character(sample_info_tab_with_color_groups_samples_only$Tank)))

vst_samples_only_physeq <- phyloseq(vst_count_samples_only_phy, sample_data(sample_info_tab_with_color_groups_samples_only))

  # generating PCoA:
vst_samples_only_pcoa <- ordinate(vst_samples_only_physeq, method="MDS", distance="euclidean")
eigen_vals_samples_only <- vst_samples_only_pcoa$values$Eigenvalues

  # now let's see how it looks with EM1 and Lambda removed
plot_ordination(vst_samples_only_physeq, vst_samples_only_pcoa, color="color_by_timepoint") +
  labs(col="Time_point") + geom_point(size=1) +
  geom_text_repel(aes(label=rownames(sample_info_tab_with_color_groups_samples_only))) +
  coord_fixed(sqrt(eigen_vals_samples_only[2]/eigen_vals_samples_only[1])) + ggtitle("PCoA") +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_timepoint[order(sample_info_tab_with_color_groups_samples_only$Time_point)])) +
  theme(legend.position="none")

    # ok things are spread out more, but we saw above that coloring by Tank seemed to be more informative in the hierarchical clustering, so let's try that with ordination and see if that trend holds up

  # coloring by Tank
plot_ordination(vst_samples_only_physeq, vst_samples_only_pcoa, color="color_by_tank") +
  labs(col="Tank") + geom_point(size=1) +
  geom_text_repel(aes(label=rownames(sample_info_tab_with_color_groups_samples_only))) +
  coord_fixed(sqrt(eigen_vals_samples_only[2]/eigen_vals_samples_only[1])) + ggtitle("PCoA") +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_tank[order(sample_info_tab_with_color_groups_samples_only$Tank)])) +
  theme(legend.position="none")
    # eh, still doesn't look all that clean based on this approach, though there are some groupings like all 3 Tank_1's, and all 3 Tank_9's, (as we saw with the hclust above)
      # but looking at the axis labels, we have the amount of variation explained by this separation, and it's pretty small anyway (only 11.3% on the x axis)



#### ALPHA-DIVERSITY ####
  # note - Alpha diversity is iffy to begin with when working with DNA sequencing data, but it is extra specious due to the nature of the sequencing technology used here, so I wouldn't necessarily trust these for much, but here is code similar to what is here just for practice's sake and inclusion: https://astrobiomike.github.io/amplicon/workflow_ex#alpha-diversity

rarecurve(t(otu_counts_samples_only_tab), step=1000, col=sample_info_tab_with_color_groups_samples_only$color_by_timepoint, lwd=2, ylab="OTUs")
abline(v=(min(rowSums(t(otu_counts_samples_only_tab)))))

rarecurve(t(otu_counts_samples_only_tab), step=1000, col=sample_info_tab_with_color_groups_samples_only$color_by_tank, lwd=2, ylab="OTUs")
abline(v=(min(rowSums(t(otu_counts_samples_only_tab)))))

  # making new phyloseq object:
aqua_phy2 <- phyloseq(otu_table(otu_counts_samples_only_tab, taxa_are_rows=T), tax_phy, sample_data(sample_info_tab_with_color_groups_samples_only))

# and now we can call the plot_richness() function on our phyloseq object
plot_richness(aqua_phy2, color="Time_point", measures=c("Chao1", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_timepoint[order(sample_info_tab_with_color_groups_samples_only$Time_point)])) +
  theme(legend.title = element_blank())

plot_richness(aqua_phy2, color="Tank", measures=c("Chao1", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_tank[order(sample_info_tab_with_color_groups_samples_only$Tank)])) +
  theme(legend.title = element_blank())


plot_richness(aqua_phy2, x="Time_point", color="Time_point", measures=c("Chao1", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_timepoint[order(sample_info_tab_with_color_groups_samples_only$Time_point)])) +
  theme(legend.title = element_blank())


plot_richness(aqua_phy2, x="Tank", color="Tank", measures=c("Chao1", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_tank[order(sample_info_tab_with_color_groups_samples_only$Tank)])) +
  theme(legend.title = element_blank())


## i think it'd be best to not report chao1, and use just "observed" richness like this:
plot_richness(aqua_phy2, color="Time_point", measures=c("Observed", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_timepoint[order(sample_info_tab_with_color_groups_samples_only$Time_point)])) +
  theme(legend.title = element_blank())

plot_richness(aqua_phy2, color="Tank", measures=c("Observed", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_tank[order(sample_info_tab_with_color_groups_samples_only$Tank)])) +
  theme(legend.title = element_blank())


plot_richness(aqua_phy2, x="Time_point", color="Time_point", measures=c("Observed", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_timepoint[order(sample_info_tab_with_color_groups_samples_only$Time_point)])) +
  theme(legend.title = element_blank())


plot_richness(aqua_phy2, x="Tank", color="Tank", measures=c("Observed", "Shannon")) +
  scale_color_manual(values=unique(sample_info_tab_with_color_groups_samples_only$color_by_tank[order(sample_info_tab_with_color_groups_samples_only$Tank)])) +
  theme(legend.title = element_blank())



#### ONE WAY TO TEST FOR CORRELATIONS ####
  # this follows from here: https://astrobiomike.github.io/amplicon/workflow_ex#betadisper-and-permutational-anova
euc_dist_samples_only <- dist(t(vst_trans_count_sample_only_tab))

  # looking by tank
anova(betadisper(euc_dist_samples_only, sample_info_tab_with_color_groups_samples_only$Tank)) # 0.9922, assumption for permutational anova test met (see above link)
adonis(euc_dist_samples_only~sample_info_tab_with_color_groups_samples_only$Tank, permutations=9999) # sig. diff, 1e-04, but semi-low R^2 (predictive capability), 0.44

  # let's look at it by timepoint, which from what we saw above, shouldn't be significant
anova(betadisper(euc_dist_samples_only, sample_info_tab_with_color_groups_samples_only$Time_point)) # 0.1021, assumption ok
adonis(euc_dist_samples_only~sample_info_tab_with_color_groups_samples_only$Time_point, permutations=9999) # hm, sig diff still, 2e-04, but even lower R^2, 0.11

  # we can run this test on anything in our sample info table that is categorical
colnames(sample_info_tab_with_color_groups_samples_only)

    # let's try by "Plants"
anova(betadisper(euc_dist_samples_only, sample_info_tab_with_color_groups_samples_only$Plants)) # 0.056, eh, a little close, assumption violated (see link above)

