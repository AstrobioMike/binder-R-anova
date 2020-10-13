library(tidyverse)
library(ggpubr)
library(car)
library(rstatix)

# reading in data
tab <- read.table("Plant-Growth-2017.tsv", sep="\t", header=TRUE)

  # checking what "type" the columns are stored as
str(tab)
    # looks okay as read in here, all are factors except the AvgGrowth column

    ### anova in R on non-diversity stuff ### generally following along with my go-to starting place for anova in R: https://www.datanovia.com/en/lessons/anova-in-r/

## starting with just looking at each column by itself ##

##### Tank #####
levels(tab$Tank)
  # it doesn't really matter for us, but we'll order this factor anyway
tab$Tank <- factor(tab$Tank, levels = c("TANK1", "TANK2", "TANK3", "TANK4", "TANK5", "TANK6", "TANK7", "TANK8", "TANK9", "TANK10"))
levels(tab$Tank)

  # summary stats
tab %>% group_by(Tank) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
   # looks like some pretty big differences in some by eye based on just the averages, like tank 1 and tank 7 at ~0.5

  # vis.
ggboxplot(tab, x = "Tank", y = "AvgGrowth")
    # looks like some pretty high variance in those that looked high given their boxplot spreads relative to the others

  # outlier check
tab %>% group_by(Tank) %>% identify_outliers(AvgGrowth)
    # seems 1 marked as extreme from tank 9, that would be something that is > 3 * the interquartile range (range between 25th percentile and the 75th percentile)
    # as noted on the page, we'll run it both with and without that one and see if it changes the outcome of the anova before we worry about it more

  # normality check
tank_model <- lm(AvgGrowth ~ Tank, data = tab)
ggqqplot(residuals(tank_model))
shapiro_test(residuals(tank_model))
    # doesn't look to satisfy the assumption of a normal distribution, looking by each tank

ggqqplot(tab, "AvgGrowth", facet.by = "Tank")
    # doesn't look that bad for any individually
tab %>% group_by(Tank) %>% shapiro_test(AvgGrowth)
    # only Tank 9 is lower than 0.05 by a bit, i think we're close enough to meet the assumption of normality, but will also run a kruskal-wallis on this one like suggested if not meeting normality assumption

  # homogeneity of variance
plot(tank_model, 1)
    # looks like much larger spread on the far right, let's see what the quantitative test says
tab %>% levene_test(AvgGrowth ~ Tank)
    # that comes up pretty low p-value, saying variances are not equal, so we don't meet that assumption
    # page notes we can use a Welch one-way anova with welch_anova_test(), so we'll look at that and the kruskal-wallis anova

  # welch anova (doesn't assume equal variances)
welch_tank_anova_result <- tab %>% welch_anova_test(AvgGrowth ~ Tank)
    # this has p value of 0.0002
  # let's look at kruskal-wallis ()
kruskal_wallis_tank_anova_result <- tab %>% kruskal_test(AvgGrowth ~ Tank)
    # that has a p-value of 0.0002 also

  ## so i'd say the tanks alone show a difference in AvgGrowth
  # post hoc (which ones are different?)
tank_pwc <- tab %>% tukey_hsd(AvgGrowth ~ Tank)
tank_pwc %>% filter(p.adj <= 0.10) # based on adj. p <= 0.10
    # looks like Tank 1 vs all of them except Tank 7
    # tank 7 vs all except tank 1
    # let's look at the overall plot again with this info in mind and see if it looks like it makes sense:
ggboxplot(tab, x = "Tank", y = "AvgGrowth")
    # 1 and 7 were the two highest, with similar medians (black horizontal bars in the middle)
        # looks sensible to me

    ## so differences in growth do seem to correlate with which tank they were in
    # looking at what else varies the same between tanks 1 and 7 just by eye on the table just out of curiousity,
        # tanks 1 and 7 both used artificial lighting and were both peas
        # tanks 2 and 3 both were inoculant, tomatoes, and used wicks

  # lastly just taking a look at the results if we removed that one value that was an outlier
    # there may be a better way, but here's one way we can remove it
tab_no_tank_outlier <- tab %>% filter(!(Tank == "TANK9" & AvgGrowth == 0.36))
welch_tank_no_outlier_anova_result <- tab_no_tank_outlier %>% welch_anova_test(AvgGrowth ~ Tank)
    # this has p value of 0.0002
  # let's look at kruskal-wallis ()
kruskal_wallis_tank_no_outlier_anova_result <- tab_no_tank_outlier %>% kruskal_test(AvgGrowth ~ Tank)
        # and both are still basically the same, so i wouldn't worry about that one outlier

    # overall, tanks 1 and 7 were different from all, and both held peas

##### InoGrav #####
  # now doing all the same, but for this column
levels(tab$InoGrav)

  # summary stats
tab %>% group_by(InoGrav) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
   # look pretty similar by eye based on averages and standard deviations

  # vis.
ggboxplot(tab, x = "InoGrav", y = "AvgGrowth")
    # look pretty similar

  # outlier check
tab %>% group_by(InoGrav) %>% identify_outliers(AvgGrowth)
    # seems to mark quite a few as extreme
    # again, we'll run it both with and without them

  # normality check
InoGrav_model <- lm(AvgGrowth ~ InoGrav, data = tab)
ggqqplot(residuals(InoGrav_model))
shapiro_test(residuals(InoGrav_model))
    # also doesn't look to satisfy the assumption of a normal distribution, looking by each tank

ggqqplot(tab, "AvgGrowth", facet.by = "InoGrav")
    # same deal, doesn't look that great by group
tab %>% group_by(InoGrav) %>% shapiro_test(AvgGrowth)
    # both have pretty low p-values, i'd say we don't meet the assumption of normality
    # going to move to the kruskal wallis test now

kruskal_wallis_InoGrav_anova_result <- tab %>% kruskal_test(AvgGrowth ~ InoGrav)
    # that has a p-value of 0.832, so InoGrav doesn't seem to be significant by itself

  # making a table without extreme outliers (ok forced to write something to do it this time as there were a bunch, ha)
    # the code may look a bit confusing, i had to try a few things to get it to work, this is not intuitive for me either, but it is using the rule the program is using mentioned above, that something is an "extreme" outlier if it is greater than 3 times the interquartile range above the 3rd quartile, or 3 times the interquartile range below the 1st quartile, see ?identify_outliers()
tab_no_InoGrav_outliers <- tab %>% group_by(InoGrav) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()

kruskal_wallis_InoGrav_no_outliers_anova_result <- tab_no_InoGrav_outliers %>% kruskal_test(AvgGrowth ~ InoGrav)
    # p-value of 0.629, so still not significant

  ## seeing if it passes the normality assumption when outliers removed
no_outliers_InoGrav_model <- lm(AvgGrowth ~ InoGrav, data = tab_no_InoGrav_outliers)
ggqqplot(residuals(no_outliers_InoGrav_model))
shapiro_test(residuals(no_outliers_InoGrav_model))
  # it does, moving forward with other tests
  # homogeneity of variance
plot(no_outliers_InoGrav_model, 1)
    # looks okay
no_outliers_InoGrav_model %>% levene_test(AvgGrowth ~ InoGrav)
    # good on homoskedasticity
  # seeing what a regular anova looks like with outliers removed
no_outliers_InoGrav_anova_results <- tab_no_InoGrav_outliers %>% anova_test(AvgGrowth ~ InoGrav)
    # comes back with p-value of 0.836

 # moving on to next done alone

##### Lighting #####

  # summary stats
tab %>% group_by(Lighting) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
   # looks like there might be a difference, though Artificial has a high standard deviation that might cause it to not show up clearly

  # vis.
ggboxplot(tab, x = "Lighting", y = "AvgGrowth")
    # yeah, look close based on medians, but large variance for the artificial category

  # outlier check
tab %>% group_by(Lighting) %>% identify_outliers(AvgGrowth)
    # there are 4 extreme, will run both with and without them as done above

  # normality check
Lighting_model <- lm(AvgGrowth ~ Lighting, data = tab)
ggqqplot(residuals(Lighting_model))
shapiro_test(residuals(Lighting_model))
    # doesn't satisfy normality assumption, checking by each group
tab %>% group_by(Lighting) %>% shapiro_test(AvgGrowth)
    # no good, moving on to kruskal wallis

kruskal_wallis_Lighting_anova_result <- tab %>% kruskal_test(AvgGrowth ~ Lighting)
    # that has a p-value of 0.072, let's take a look if we remove the outliers
tab_no_Lighting_outliers <- tab %>% group_by(Lighting) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()

kruskal_wallis_Lighting_no_outliers_anova_result <- tab_no_Lighting_outliers %>% kruskal_test(AvgGrowth ~ Lighting)
    # p-value of 0.138, so still not significant, and actually increased

  # trying normality test with outliers removed to see if we can furthur pursue standard anova that way
no_outliers_Lighting_model <- lm(AvgGrowth ~ Lighting, data = tab_no_Lighting_outliers)
ggqqplot(residuals(no_outliers_Lighting_model))
shapiro_test(residuals(no_outliers_Lighting_model))
  # still doesn't pass normality test, so can't try standard anova on it

  # seems Lighting alone doesn't show a significant difference
    # moving on to Plant

##### Plant #####

  # summary stats
tab %>% group_by(Plant) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
   # looks like a pretty big difference for Peas, but much higher variance too

  # vis.
ggboxplot(tab, x = "Plant", y = "AvgGrowth")

  # outlier check
tab %>% group_by(Plant) %>% identify_outliers(AvgGrowth)
    # none marked extreme

  # normality check
Plant_model <- lm(AvgGrowth ~ Plant, data = tab)
ggqqplot(residuals(Plant_model))
shapiro_test(residuals(Plant_model))
    # doesn't look to satisfy the assumption of a normal distribution, checking by group

ggqqplot(tab, "AvgGrowth", facet.by = "Plant")
    # doesn't look that bad for any individually
tab %>% group_by(Plant) %>% shapiro_test(AvgGrowth)
    # no good for tomatoes or peas, moving on to kruskal wallis

kruskal_wallis_Plant_anova_result <- tab %>% kruskal_test(AvgGrowth ~ Plant)
    # 2e-6 p-value, running posthoc test

  # post hoc (which ones are different?)
Plant_pwc <- tab %>% tukey_hsd(AvgGrowth ~ Plant)
    # looks like Peas are the standouts, different from basil and different from tomatos, which is what the initial plot looked like:
ggboxplot(tab, x = "Plant", y = "AvgGrowth")
    # i don't know if that is just a consequence of the peas happening to grow more, that would depend on if the average growth values were somehow normalized across the 3 plant types (like growth relative to final size or something)

  # moving on to next column

##### Media #####
  # summary stats
tab %>% group_by(Media) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
    # coconut is lower, others have higher variance too though
  # vis.
ggboxplot(tab, x = "Media", y = "AvgGrowth")

  # outlier check
tab %>% group_by(Media) %>% identify_outliers(AvgGrowth)
    # several are marked as extreme, will again look at with and without them removed as done above

  # normality check
Media_model <- lm(AvgGrowth ~ Media, data = tab)
ggqqplot(residuals(Media_model))
shapiro_test(residuals(Media_model))
    # doesn't look to satisfy the assumption of a normal distribution, looking by each tank

ggqqplot(tab, "AvgGrowth", facet.by = "Media")
tab %>% group_by(Media) %>% shapiro_test(AvgGrowth)
    # doesn't meet assumption of normality, moving to kruskal wallis

kruskal_wallis_Media_anova_result <- tab %>% kruskal_test(AvgGrowth ~ Media)
    # that has a p-value of 0.443, no sig differences

    # taking a look with outliers removed
tab_no_Media_outliers <- tab %>% group_by(Media) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()

kruskal_wallis_Media_no_outliers_anova_result <- tab_no_Media_outliers %>% kruskal_test(AvgGrowth ~ Media)
    # p-value of 0.570, so still not significant

  # trying normality test with outliers removed to see if we can furthur pursue standard anova that way
no_outliers_Media_model <- lm(AvgGrowth ~ Media, data = tab_no_Media_outliers)
ggqqplot(residuals(no_outliers_Media_model))
shapiro_test(residuals(no_outliers_Media_model))
  # still doesn't pass normality test, so can't try standard anova on it

  # seems Media alone didn't make a difference
    # on to next

##### Water #####
  # summary stats
tab %>% group_by(Water) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
  # vis.
ggboxplot(tab, x = "Water", y = "AvgGrowth")
      # look pretty similar
  # outlier check
tab %>% group_by(Water) %>% identify_outliers(AvgGrowth)
    # several are marked as extreme, will again look at with and without them removed as done above

  # normality check
Water_model <- lm(AvgGrowth ~ Water, data = tab)
ggqqplot(residuals(Water_model))
shapiro_test(residuals(Water_model))
    # doesn't look to satisfy the assumption of a normal distribution, looking by each tank

ggqqplot(tab, "AvgGrowth", facet.by = "Water")
tab %>% group_by(Water) %>% shapiro_test(AvgGrowth)
    # doesn't meet assumption of normality, moving to kruskal wallis

kruskal_wallis_Water_anova_result <- tab %>% kruskal_test(AvgGrowth ~ Water)
    # that has a p-value of 0.134, no sig differences

    # taking a look with outliers removed
tab_no_Water_outliers <- tab %>% group_by(Water) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()

kruskal_wallis_Water_no_outliers_anova_result <- tab_no_Water_outliers %>% kruskal_test(AvgGrowth ~ Water)
    # p-value of 0.133

  # trying normality test with outliers removed to see if we can furthur pursue standard anova that way
no_outliers_Water_model <- lm(AvgGrowth ~ Water, data = tab_no_Water_outliers)
ggqqplot(residuals(no_outliers_Water_model))
shapiro_test(residuals(no_outliers_Water_model))
  # still doesn't pass normality test, so can't try standard anova on it


## Testing for interactions (following along with the Two-way anova section: https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova)
# based on how this is noted there, "If an interaction effect does not exist, main effects could be reported.", so i'm thinking of this as for us to, for example, look at the Tank and Plant factors together, since they both correlate (Tanks 1 and 7 were highest growth and both had peas)

  ## i'm not sure we can do any two-way or more because i don't think any individual groupings encompass all other factor types, and i think that's needed

