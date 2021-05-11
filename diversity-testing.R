# couldn't figure out how to install these permanently, so just doing here, only take a couple minutes
devtools::install_github("adw96/breakaway")
devtools::install_github("adw96/DivNet")

library(phyloseq)
library(DivNet)
library(tidyverse)
library(magrittr)

# making helper functions to run divnet and testing following along with tutorial here: https://github.com/adw96/DivNet/blob/master/vignettes/getting-started.Rmd)
# helper functions for mass running things
make_div_obj <- function(phy_obj, wanted_formula) {

    # running divnet
    div_obj <- phy_obj %>% divnet(formula = wanted_formula, ncores = 4)

    return(div_obj)

}

run_div_test <- function(phy_obj, div_obj, target_covariate, target_ref = NULL) {

    # running testing
    estimates <- div_obj$shannon %>% summary %$% estimate
    ses <- sqrt(div_obj$`shannon-variance`)
    sample_tab <- phy_obj %>% sample_data %>% data.frame()
    predictors_ <- sample_tab[, names(sample_tab) %in% target_covariate]

    # there is no built-in way to change the baseline in the contrasts, so modifying the factor levels as suggested on divnet issue here: https://github.com/adw96/DivNet/issues/37
    if ( is.null(target_ref) ) {
        X <- model.matrix(~predictors_)
        baseline <- predictors_ %>% levels %>% head(1)
    } else {
        predictors_ <- predictors_ %>% relevel(ref = target_ref)
        baseline <- predictors_ %>% levels %>% head(1)
        X <- model.matrix(~predictors_)
    }


    res <- betta(estimates, ses, X)

    # global p-value (are any of the groups different)
    global_p <- res$global[2]

    cat("Global p-value: ", global_p, "\n\n", sep = "")

    cat("With '", baseline, "' as the baseline:\n\n", sep = "")

    print(res$table)

    return(res)

}

# doing a test with the data stored with the divnet package and used in the tutorial linked above - it's called Lee in the package because it comes from me :)
data(Lee)

# collapsing at phylum level for testing
lee_phylum <- tax_glom(Lee, taxrank = "Phylum")
lee_phylum %>% sample_data
# making divnet object based on "char" column in sample data
test_div_obj <- make_div_obj(lee_phylum, ~char)

# running test based on "characteristic" in the test data
test_div_res <- run_div_test(lee_phylum, test_div_obj, "char")
# Global p-value: 0
#
# With 'altered' as the baseline:
#
#                        Estimates Standard Errors p-values
# (Intercept)           1.34484669      0.01007056    0.000
# predictors_biofilm   -0.91111064      0.11083815    0.000
# predictors_carbonate -0.05379851      0.13020724    0.679
# predictors_glassy    -0.20930449      0.02420628    0.000
# predictors_water      0.11472333      0.09311783    0.218

    # good, works like it's supposed to, global diff p-value is 0 saying something is diff
    # 'altered' is the baseline and has greater shannon diversity than glassy.
    # run as above, we are only comparing everything against 'altered'
    # here's how we can run it to compare everything against something else, e.g. 'water'

test_div_res2 <- run_div_test(lee_phylum, test_div_obj, "char", target_ref = "water")
# Global p-value: 0
#
# With 'water' as the baseline:
#
#                       Estimates Standard Errors p-values
# (Intercept)           1.4595700      0.01007056    0.000
# predictors_altered   -0.1147233      0.01125231    0.000
# predictors_biofilm   -1.0258340      0.11083815    0.000
# predictors_carbonate -0.1685218      0.13020724    0.196
# predictors_glassy    -0.3240278      0.02420628    0.000

    # and that shows us everything is diff from the water (and lower in diversity), except the carbonate sample
    # okay, helper functions check out based on what the test data set looks like, now moving onto our stuff


# reading in our data (excluding EM1 and Lambda)
otu_counts_tab <- read.table("data/combined_otu_counts_tab.txt", header=TRUE, row.names=1, sep="\t")[,c(1:30)]
taxonomy_tab <- read.table("data/combined_taxonomy_tab.txt", header=TRUE, row.names=1, sep="\t", na.strings="<NA>")
sample_info_tab <- read.table("data/sample_info_tab.txt", header=T, check.names=F, row.names=1, sep="\t")[c(1:30),]

# making phyloseq object
counts_phy <- otu_table(otu_counts_tab, taxa_are_rows=TRUE)
tax_phy <- tax_table(as.matrix(taxonomy_tab))
sample_phy <- sample_data(sample_info_tab)
aqua_phy <- phyloseq(counts_phy, tax_phy, sample_phy)

# grouping at diff ranks
aqua_class_phy <- tax_glom(aqua_phy, taxrank = "Class")
aqua_order_phy <- tax_glom(aqua_phy, taxrank = "Order")
aqua_family_phy <- tax_glom(aqua_phy, taxrank = "Family")

    # factors we'll look at
names(sample_info_tab)

## starting at CLASS level
    # Time_point
class_time_point_div_obj <- make_div_obj(aqua_class_phy, ~Time_point)
class_time_point_T1_base_div_res <- run_div_test(aqua_class_phy, class_time_point_div_obj, "Time_point")
# Global p-value: 0
#
# With 'T1' as the baseline:
#
#                 Estimates Standard Errors p-values
# (Intercept)    1.48577272      0.01873717    0.000
# predictors_T2  0.12386158      0.02641508    0.000
# predictors_T3 -0.06554238      0.04737824    0.167

    # Intercept is T1, as that's the only one not listed
    # this says T1 was sig. diff from T2 (p-value), and T2 had a higher diversity (estimate is positive 0.1238...)
    # let's look with T2 as the base, this will tell us if T2 was sig diff from T3
class_time_point_T2_base_div_res <- run_div_test(aqua_class_phy, class_time_point_div_obj, "Time_point", target_ref = "T2")
# Global p-value: 0
#
# With 'T2' as the baseline:
#
#                Estimates Standard Errors p-values
# (Intercept)    1.6096343      0.01873717        0
# predictors_T1 -0.1238616      0.03211331        0
# predictors_T3 -0.1894040      0.04737824        0

    # this says T2 was sig. diff from T3 (says T1 too, but we already knew that one from above)
    # T3 is 0.1894... lower in estimated shannon diversity than T2, with a p-value < 0.000

    ## so based on time point, T2 has the highest diversity, T2 is sig higher than T1, and T3 is sig lower than T2, T3 is not sig higher than T1
    # let's see if we can plot these with their error bars

class_time_point_div_obj$shannon %>%
    plot(aqua_class_phy, color = "Time_point") +
    xlab("Time point") +
    ylab("Shannon diversity estimate\n(Class level)")

   # since the estimate is based on all of them, they are all the same, so only going to plot one per Time point, rather than all
    # i'm sure there is a better way to do this, but now isn't the time for figuring out better ways, it's the time for doing things, ha
sub_time_point_div_obj_shannon <- class_time_point_div_obj$shannon[c(1,11,21)]
names(sub_time_point_div_obj_shannon) <- c("T1", "T2", "T3")

div_est_vec <- c(unlist(sub_time_point_div_obj_shannon[1])[1] %>% as.numeric(), unlist(sub_time_point_div_obj_shannon[2])[1] %>% as.numeric(), unlist(sub_time_point_div_obj_shannon[3])[1] %>% as.numeric())
ci_vec <- c(unlist(sub_time_point_div_obj_shannon[1])[2] %>% as.numeric() * 2, unlist(sub_time_point_div_obj_shannon[2])[2] %>% as.numeric() * 2, unlist(sub_time_point_div_obj_shannon[3])[2] %>% as.numeric() * 2)
sub_time_plot_tab <- data.frame(Time_point = names(sub_time_point_div_obj_shannon), shannon = div_est_vec, ci = ci_vec)

ggplot(sub_time_plot_tab, aes(x = Time_point, y = shannon, color = Time_point)) +
    geom_errorbar(aes(ymin = shannon - ci, ymax = shannon + ci), width = 0) +
    geom_point() + theme_bw() + coord_cartesian(ylim = c(0,2)) +
    xlab("Time point") + ylab("Shannon diversity estimate\n(Class level)") +
    theme(legend.position = "none")


## will look at others too, but haven't yet, e.g. others we have in our sample table
names(sample_info_tab)

  # and can look at other rank summarizations, meaning other than class, probably shouldn't go lower than family though, and we staying higher is generally better because of the high error rate in nanopore sequencing
  
