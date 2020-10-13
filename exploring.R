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
  # i don't think we can include Tank, because each tank doesn't have combinations of all factors, which i also think means we can't "Block" the Tank factor out of the comparison. I read about "Blocking" from this page: https://stat.ethz.ch/~meier/teaching/anova/block-designs.html
  # i also read about "blocking" on this page: https://online.stat.psu.edu/stat503/book/export/html/646 that says we can try analysis of covariance when it's not possible for us to block, so we could maybe try that (https://www.datanovia.com/en/lessons/ancova-in-r/), but i didn't have any luck yet and didn't get to spend too much time with it. i'm not sure what's up with that yet
    # we also need to keep in mind that although the Tanks show a large difference, it tracks with the Plant type (peas), so maybe some tests can be done on just basil and tomatoes (which would eliminate the large Tank/Peas effect), or the AvgGrowth can maybe somehow be normalized based on total size or total change or something


  # anyway, to start, i think we can look at InoGrav and Lighting
tab %>% group_by(InoGrav, Lighting) %>% summarize(n=n())
#   InoGrav   Lighting       n
#   <fct>     <fct>      <int>
# 1 Gravel    Artificial    20
# 2 Gravel    Natural       23
# 3 Inoculant Artificial    35
# 4 Inoculant Natural       33

  # InoGrav and Plant
tab %>% group_by(InoGrav, Plant) %>% summarize(n=n())
#   InoGrav   Plant        n
#   <fct>     <fct>    <int>
# 1 Gravel    Basil       23
# 2 Gravel    Peas         9
# 3 Gravel    Tomatoes    11
# 4 Inoculant Basil       12
# 5 Inoculant Peas        22
# 6 Inoculant Tomatoes    34

  # InoGrav and Media
tab %>% group_by(InoGrav, Media) %>% summarize(n=n())

#   InoGrav   Media       n
#   <fct>     <fct>   <int>
# 1 Gravel    Coconut    23
# 2 Gravel    Float      20
# 3 Inoculant Clay       45
# 4 Inoculant Float      23

  # InoGrav and Water
tab %>% group_by(InoGrav, Water) %>% summarize(n=n())
#   InoGrav   Water     n
#   <fct>     <fct> <int>
# 1 Gravel    Pumps    11
# 2 Gravel    Wicks    32
# 3 Inoculant Pumps    45
# 4 Inoculant Wicks    23

  # Lighting and Plant
tab %>% group_by(Lighting, Plant) %>% summarize(n=n())
#   Lighting   Plant        n
#   <fct>      <fct>    <int>
# 1 Artificial Basil       12
# 2 Artificial Peas        20
# 3 Artificial Tomatoes    23
# 4 Natural    Basil       23
# 5 Natural    Peas        11
# 6 Natural    Tomatoes    22

  # Lighting and Media
tab %>% group_by(Lighting, Media) %>% summarize(n=n())
#   Lighting   Media       n
#   <fct>      <fct>   <int>
# 1 Artificial Clay       23
# 2 Artificial Coconut    11
# 3 Artificial Float      21
# 4 Natural    Clay       22
# 5 Natural    Coconut    12
# 6 Natural    Float      22

  # Lighting and Water
tab %>% group_by(Lighting, Water) %>% summarize(n=n())
#   Lighting   Water     n
#   <fct>      <fct> <int>
# 1 Artificial Pumps    23
# 2 Artificial Wicks    32
# 3 Natural    Pumps    33
# 4 Natural    Wicks    23

  # Plant and Media is missing a Coconut for Peas (those could probably be split apart if wanted, meaning getting read of Peas to look at the 3 media in basil and tomatoes, and/or getting rid of coconut to look at clay and float across all 3 plants)
tab %>% group_by(Plant, Media) %>% summarize(n=n())
#   Plant    Media       n
#   <fct>    <fct>   <int>
# 1 Basil    Clay       12
# 2 Basil    Coconut    12
# 3 Basil    Float      11
# 4 Peas     Clay       22
# 5 Peas     Float       9
# 6 Tomatoes Clay       11
# 7 Tomatoes Coconut    11
# 8 Tomatoes Float      23

  # Plant and Water
tab %>% group_by(Plant, Water) %>% summarize(n=n())
#   Plant    Water     n
#   <fct>    <fct> <int>
# 1 Basil    Pumps    23
# 2 Basil    Wicks    12
# 3 Peas     Pumps    22
# 4 Peas     Wicks     9
# 5 Tomatoes Pumps    11
# 6 Tomatoes Wicks    34

  # Media and Water, i don't think can be done together because Pumps is missing a coconut
tab %>% group_by(Media, Water) %>% summarize(n=n())


### starting with InoGrav and Lighting ###
tab %>% group_by(InoGrav, Lighting) %>% summarize(n=n())
  # taking a look
tab %>% group_by(InoGrav, Lighting) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
ggboxplot(tab, x = "InoGrav", y = "AvgGrowth", color = "Lighting")
  # checking for outliers
tab %>% group_by(InoGrav, Lighting) %>% identify_outliers(AvgGrowth)
    # quite a few, but all Tank related, so i don't think we should necessarily remove them...

  # normality test
InoGrav_Lighting_model <- lm(AvgGrowth ~ InoGrav * Lighting, data = tab)
ggqqplot(residuals(InoGrav_Lighting_model))
shapiro_test(residuals(InoGrav_Lighting_model))
    # yea yea, looking by groups
tab %>% group_by(InoGrav, Lighting) %>% shapiro_test(AvgGrowth)
    # still no good, but when we include Tanks as a group it mostly is ok
tab %>% group_by(InoGrav, Lighting, Tank) %>% shapiro_test(AvgGrowth)
    # i'm just not sure how to incorporate that in the anova, anytime its included, all other factors are NA (because they are represented across different tanks)


  # after some more poking, i found this "robust anova" program (https://rdrr.io/cran/walrus/man/ranova.html), that i think we can use without worrying about the assumptions as much
 # can be installed in the web-based R environment like so:

install.packages("walrus")
library(walrus)
ranova(tab, "AvgGrowth", c("InoGrav"), ph=TRUE) # 0.73
ranova(tab, "AvgGrowth", c("Lighting"), ph=TRUE) # 0.16
ranova(tab, "AvgGrowth", c("Plant"), ph=TRUE) # 0.0003 (each different)
ranova(tab, "AvgGrowth", c("Media"), ph=TRUE) # 0.44
ranova(tab, "AvgGrowth", c("Water"), ph=TRUE) # 0.33

ranova(tab, "AvgGrowth", c("InoGrav", "Lighting"), ph=TRUE) # interaction 0.06
ranova(tab, "AvgGrowth", c("InoGrav", "Plant"), ph=TRUE) # interaction 0.25
ranova(tab, "AvgGrowth", c("InoGrav", "Media"), ph=TRUE) # not possible due to missing combinations
ranova(tab, "AvgGrowth", c("InoGrav", "Water"), ph=TRUE) # interaction 0.13

ranova(tab, "AvgGrowth", c("Lighting", "Plant"), ph=TRUE) # interaction 0.10
ranova(tab, "AvgGrowth", c("Lighting", "Media"), ph=TRUE) # interaction 0.22
ranova(tab, "AvgGrowth", c("Lighting", "Water"), ph=TRUE) # interaction 0.93

ranova(tab, "AvgGrowth", c("Plant", "Media"), ph=TRUE) # not possible
ranova(tab, "AvgGrowth", c("Plant", "Water"), ph=TRUE) # 0.34

ranova(tab, "AvgGrowth", c("Media", "Water"), ph=TRUE) # not possible

  # i would just rely on those results (as i don't know better if or how we can "correct"/take-account of the Tank thing that makes it hard to look across the other treatments as you noted)
  ## and/or i would try running some things on each plant individually, that will help get rid of the overlapping problem between Tanks and Plants
  # e.g. making tables of each plant
basil_tab <- tab %>% filter(Plant == "Basil")
peas_tab <- tab %>% filter(Plant == "Peas")
tomatoes_tab <- tab %>% filter(Plant == "Tomatoes")

  # and then i'd try some of what we did at the top, but working on one plant at a time and seeing if things are cleaner that way



