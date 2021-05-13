# couldn't figure out how to install these permanently, so just doing here, only take a couple minutes
devtools::install_github("adw96/breakaway", upgrade = "never")
devtools::install_github("adw96/DivNet", upgrade = "never")

library(phyloseq)
library(DivNet)
library(tidyverse)
library(magrittr)

# making helper functions to run divnet and testing following along with tutorial here: https://github.com/adw96/DivNet/blob/master/vignettes/getting-started.Rmd)

make_div_obj <- function(phy_obj, wanted_formula) {

    # running divnet
    div_obj <- phy_obj %>% divnet(formula = wanted_formula)

    return(div_obj)

}

run_div_test <- function(phy_obj, div_obj, target_ref = NULL) {

    # getting covariate
    target_covariate <- attributes(div_obj$X)$contrasts %>% names

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

    cat("With '", baseline, "' as the baseline:\n\n", sep = "")

    tab <- res$table %>% data.frame

    # p-values are calculated as: 2*(1-pnorm(abs(Estimate/SE))) (see betta function)
    # I think getting 0s is just due to how many digits R is working with. It seems
    # pnorm(5.32672)  = 0.9999999 and pnorm(5.32673) = 1. So if our abs(Estimate/SE)
    # is greater than 5.32672, we’re going to get a p-value of 0. So setting 0's to
    # < 1e-8 for what is printed out

    new_p_vec <- vector()

    for ( value in tab$p.values ) {
        if ( value == 0 ) {
            new_p_vec <- c(new_p_vec, "< 1e-8")
        } else {
            new_p_vec <- c(new_p_vec, value)
        }
    }

    tab$p.values <- new_p_vec

    print(tab)

    return(res)

}

plot_DivNet_diversity <- function(div_obj, sample_tab = sample_info_tab) {

    # getting covariate
    target_covariate <- attributes(div_obj$X)$contrasts %>% names

    # getting unique treatment types
    unique_treatment_types <- sample_tab %>% pull(target_covariate) %>% unique %>% as.character

    unique_samples_needed <- c()

    for ( type in unique_treatment_types ) {

        sample_to_add <- sample_tab %>% filter(!!as.symbol(target_covariate) == type) %>% head(1) %>% row.names()

        unique_samples_needed <- c(unique_samples_needed, sample_to_add)
    }

    sub_div_obj_shannon <- div_obj$shannon[c(unique_samples_needed)]

    div_est_vec <- c()
    ci_vec <- c()

    for ( i in 1:length(unique_samples_needed) ) {
        div_est_vec <- c(div_est_vec, unlist(sub_div_obj_shannon[i])[1] %>% as.numeric())
        ci_vec <- c(ci_vec, unlist(sub_div_obj_shannon[i])[2] %>% as.numeric() * 2)
    }

    sub_plot_tab <- data.frame(target_covariate = unique_treatment_types, shannon = div_est_vec, ci = ci_vec)

    # getting max value for plot y-axis
    max_y <- sub_plot_tab$shannon %>% max %>% ceiling

    plot <- ggplot(sub_plot_tab, aes(x = target_covariate, y = shannon, color = target_covariate)) +
        geom_errorbar(aes(ymin = shannon - ci, ymax = shannon + ci), width = 0) +
        geom_point() + theme_bw() +
        xlab(target_covariate) + ylab("Shannon diversity estimate\n(Class level)") +
        theme(legend.position = "none") + coord_cartesian(ylim = c(0,max_y))

    return(plot)

}


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
# aqua_order_phy <- tax_glom(aqua_phy, taxrank = "Order")
# aqua_family_phy <- tax_glom(aqua_phy, taxrank = "Family")

    # factors we'll look at (we'll need these specific names exactly as used below to tell the programs what to look at)
names(sample_info_tab)

    ## Time_point

        # first function creates the DivNet object holding the estimated diversities
        # we need to give it the phyloseq object, that's the first argument here
        # and the second starts with the tilde '~', which means "as a function of", and we then give it the column name from our sample info table, here doing 'Time_point'
class_time_point_div_obj <- make_div_obj(aqua_class_phy, ~Time_point)

        # this second function does the actual hypothesis test of whether there is a difference or not
        # we need to give it the same phyloseq object (first argument)
        # then the DivNetobject we just made (second argument)
        # then last we need to pass a string this time (meaning in quotes), that specifies what we're testing, the same "Time_point" column we specified above (just needs to be done with quotes this time and not before)
class_time_point_T1_base_div_res <- run_div_test(aqua_class_phy, class_time_point_div_obj)

# With 'T1' as the baseline:

#                 Estimates Standard.Errors p.values
# (Intercept)    1.48609385      0.01642516   < 1e-8
# predictors_T2  0.12242525      0.04933613    0.013
# predictors_T3 -0.06621529      0.02088762    0.002

    # Intercept is T1, as that's the only one not listed
    # this says T1 was sig. diff from T2 (p-value of 0.013), and T2 had a higher diversity (estimate is positive 0.122...)
    # and says T3 was sig. diff than T1, lower diversity (This changed a little since last time as I modified how we are running the code)

    # let's look with T2 as the base, this will tell us if T2 was sig diff from T3
class_time_point_T2_base_div_res <- run_div_test(aqua_class_phy, class_time_point_div_obj, "Time_point", target_ref = "T2")

# With 'T2' as the baseline:
# 
#                Estimates Standard.Errors p.values
# (Intercept)    1.6085191      0.01642516   < 1e-8
# predictors_T1 -0.1224252      0.03156335   < 1e-8
# predictors_T3 -0.1886405      0.02088762   < 1e-8


    # this says T2 was sig. diff from T3 (says T1 too, but we already knew that one from above)
    # T3 is 0.188... lower in estimated shannon diversity than T2, with a p-value < 1e-8

    ## so based on time point, T2 has the highest diversity, T2 is sig higher than T1, and T3 is sig lower than T2, T3 is not sig higher than T1
    # we can plot these with confidence intervals with this helper function now :)
plot_DivNet_diversity(class_time_point_div_obj)


## now looking at lighting
names(sample_info_tab)

## lighting

class_lighting_div_obj <- make_div_obj(aqua_class_phy, ~Lighting)
class_lighting_div_res <- run_div_test(aqua_class_phy, class_lighting_div_obj)

#                    Estimates Standard.Errors p.values
# (Intercept)         1.420162      0.01193468   < 1e-8
# predictors_Natural  0.136010      0.02141959   < 1e-8

    # says natural has greater diversity 

# plotting
plot_DivNet_diversity(class_lighting_div_obj)
  
    # looking at counts of plants in each lighting type
sample_info_tab %>% group_by(Lighting, Plants) %>% summarize(n = n())
#   Lighting   Plants       n
#   <fct>      <fct>    <int>
# 1 Artificial Basil        6
# 2 Artificial Peas         6
# 3 Artificial Tomatoes     3
# 4 Natural    Basil        6
# 5 Natural    Peas         3
# 6 Natural    Tomatoes     6


## plants
class_plants_div_obj <- make_div_obj(aqua_class_phy, ~Plants)
class_plants_div_res <- run_div_test(aqua_class_phy, class_plants_div_obj)

# With 'Basil' as the baseline:

#                     Estimates Standard.Errors p.values
# (Intercept)         1.3583161      0.01947881   < 1e-8
# predictors_Peas     0.1378274      0.04539046    0.002
# predictors_Tomatoes 0.2520203      0.03179659   < 1e-8

# running one with Peas as baseline so we can see if peas and tomatoes are diff
class_plants_div_res_peas_baseline <- run_div_test(aqua_class_phy, class_plants_div_obj, target_ref = "Peas")
# With 'Peas' as the baseline:
# 
#                      Estimates Standard.Errors p.values
# (Intercept)          1.4961435      0.01947881   < 1e-8
# predictors_Basil    -0.1378274      0.02934699   < 1e-8
# predictors_Tomatoes  0.1141929      0.03179659   < 1e-8

plot_DivNet_diversity(class_plants_div_obj)

#### i think i need to double check with Amy (developer of DivNet) that there isn't a problem with the tutorial i'm following, because this seems a little off to me that they are almost all significant

    ## for now, let's finish figuring out how we can split the phyloseq object based on plant type

#### splitting into plant-specific phyloseq objects
aqua_class_phy
    # ah, wonderfully there is already a function built for this :)
aqua_class_peas_phy <- subset_samples(aqua_class_phy, Plants == "Peas")
aqua_class_basil_phy <- subset_samples(aqua_class_phy, Plants == "Basil")
aqua_class_tomatoes_phy <- subset_samples(aqua_class_phy, Plants == "Tomatoes")

    # looking at timepoint of just the peas
class_timepoints_peas_div_obj <- make_div_obj(aqua_class_peas_phy, ~Time_point)
class_timepoints_peas_div_res <- run_div_test(aqua_class_peas_phy, class_timepoints_peas_div_obj)
    # we need to specify a subset sample table too for this plotting function to work properly (that's the second argument here now)
plot_DivNet_diversity(class_timepoints_peas_div_obj, sample_tab = data.frame(sample_data(aqua_class_peas_phy)))

    # just the basil
class_timepoints_basil_div_obj <- make_div_obj(aqua_class_basil_phy, ~Time_point)
class_timepoints_basil_div_res <- run_div_test(aqua_class_basil_phy, class_timepoints_basil_div_obj)
plot_DivNet_diversity(class_timepoints_basil_div_obj, sample_tab = data.frame(sample_data(aqua_class_basil_phy)))


