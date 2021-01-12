### Mike's added-in comments start with 3 #'s like this

library(tidyverse)
library(ggpubr)
library(car)
library(rstatix)
tab <- read.csv(file = 'Plant-Growth-2018-3.csv')
str(tab)
tab %>% group_by(Tank.) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
ggboxplot(tab, x = "Tank.", y = "AvgGrowth")

  # trying it by soil and water first since that's an important factor in this one
tab %>% group_by(SoilWater) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
tab %>% group_by(SoilWater) %>% identify_outliers(AvgGrowth)
  # there are extreme outliers, but I'll run it without them first
soilwater_model <- lm(AvgGrowth ~ SoilWater, data = tab)
ggqqplot(residuals(soilwater_model))
shapiro_test(residuals(soilwater_model))
ggqqplot(tab, "AvgGrowth", facet.by = "SoilWater")
tab %>% group_by(SoilWater) %>% shapiro_test(AvgGrowth)
plot(soilwater_model, 1)
tab %>% levene_test(AvgGrowth ~ SoilWater)
  # since it's not normally distributed, I'm gonna do the kruskal wallis
kruskal_wallis_soilwater_anova_result <- tab %>% kruskal_test(AvgGrowth ~ SoilWater)
  # ey! p-value of 1.66 e-07!
    ### nice! i'm adding in a boxplot figure of these two groups
ggboxplot(tab, x = "SoilWater", y = "AvgGrowth")

  # trying it by inoculant and gravel
tab %>% group_by(InoGrav) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
### these look much more similar than the Soil/Water did

tab %>% group_by(InoGrav) %>% identify_outliers(AvgGrowth)
inograv_model <- lm(AvgGrowth ~ InoGrav, data = tab)
ggqqplot(residuals(inograv_model))
shapiro_test(residuals(inograv_model))
ggqqplot(tab, "AvgGrowth", facet.by = "InoGrav")
tab %>% group_by(InoGrav) %>% shapiro_test(AvgGrowth)
plot(inograv_model, 1)
tab %>% levene_test(AvgGrowth ~ InoGrav)
  # since it's not normally distributed, I'm gonna do the kruskal wallis
kruskal_wallis_InoGrav_anova_result <- tab %>% kruskal_test(AvgGrowth ~ InoGrav)
  # p-value is .292, not significant.

  # trying by tank to check
tab %>% group_by(Tank.) %>% get_summary_stats(AvgGrowth, type = "mean_sd")  
tab %>% group_by(Tank.) %>% identify_outliers(AvgGrowth)
  # no extreme outliers
tank._model <- lm(AvgGrowth ~ Tank., data = tab)
ggqqplot(residuals(tank._model))
shapiro_test(residuals(tank._model))
  # not normal :(
### this isn't anything to be sad about :) some data just aren't, haha

ggqqplot(tab, "AvgGrowth", facet.by = "Tank.")
tab %>% group_by(Tank.) %>% shapiro_test(AvgGrowth)
  # only tank 6 is normal, so basically shouldn't assume this for any of them, right?
### yeah that's what i think too, especially for the purposes of doing a test across all of them
plot(tank._model, 1)
  # still don't know exactly what this means though haha
### haha, me either. I think we just want any points plotted to be randomly distributed, so like the fact that they are tighter on the left side then they are on the right side is no bueno
tab %>% levene_test(AvgGrowth ~ Tank.)
kruskal_wallis_Tank._anova_result <- tab %>% kruskal_test(AvgGrowth ~ Tank.)
  # p-value is .000744, looking at individual comparisons:
tank._pwc <- tab %>% tukey_hsd(AvgGrowth ~ Tank.)
tank._pwc %>% filter(p.adj <= 0.10)
  # It doesn't seem to be consistently one tank. Pairs : 1/2, 1/5, 1/6, 2/4, 3/4, 4/5, 4/6.
### i'm interested to see what the tank tests look like with outliers removed

  # now without outliers for soil and water
tab_no_SoilWater_outliers <- tab %>% group_by(SoilWater) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
kruskal_wallis_SoilWater_no_outliers_anova_result <- tab_no_SoilWater_outliers %>% kruskal_test(AvgGrowth ~ SoilWater)
  # it's even MORE significant! lovely. I'm just going to do the same tests you did to check if we can run the ANOVA too.
no_outliers_SoilWater_model <- lm(AvgGrowth ~ SoilWater, data = tab_no_SoilWater_outliers)
ggqqplot(residuals(no_outliers_SoilWater_model))
shapiro_test(residuals(no_outliers_SoilWater_model))
  # nah.. still not normal :( 
plot(no_outliers_SoilWater_model, 1)
no_outliers_SoilWater_model %>% levene_test(AvgGrowth ~ SoilWater)
  # wait and now it doesn't meet the variance parameter either?
### hmm, seems the 3 outliers dropped were all soil:
tab %>% group_by(SoilWater) %>% identify_outliers(AvgGrowth) %>% filter(is.extreme == TRUE)
### and looking at the previous residuals vs fitted plot with the soil outliers included:
plot(soilwater_model, 1)
### vs this one:
plot(no_outliers_SoilWater_model, 1)
### we can see the high points dropped from the left column (soil aparently), make the two distributions look more different now


welch_SoilWater_no_outlier_anova_result <- tab_no_SoilWater_outliers %>% welch_anova_test(AvgGrowth ~ SoilWater)
  # also significant, so I guess we're good either way then.
### yea, looks to me like there is a pretty clear difference in soil vs water :) 
### i'll also say we can use a regular welch "t-test" here, since there are only these 2 groups, but i think the results will be the same, let's see
   ### the t_test function by default uses the Welch version (suitable for non-equal variances)
welch_SoilWater_no_outlier_ttest_result <- tab_no_SoilWater_outliers %>% t_test(AvgGrowth ~ SoilWater)

### an outlier-free boxplot with our:
ggboxplot(tab_no_SoilWater_outliers, x = "SoilWater", y = "AvgGrowth", title = "Welch p-value = 2.82e-9")

### i'll also note something i just came across, we may want to lean less on the shapiro test for normality when we have this many data points
### this page (https://www.datanovia.com/en/lessons/t-test-in-r/#assumptions-and-preliminary-tests) notes this in the middle:
   # "Note that, if your sample size is greater than 50, the normal QQ plot is preferred because at larger sample sizes the Shapiro-Wilk test becomes very sensitive even to a minor deviation from normality."
### now we are really in the realm of just not having enough experience looking at these to say whether something is "normal" or not, but based on the figure on that page and what this one looks like:
ggqqplot(residuals(no_outliers_SoilWater_model))
### i might be inclined to say it looks kinda normal...

### last thing i'm going to check is a wilcoxon/mann-whitney test (a non-parametric t-test that is more robust to non-normality, just for the heck of it) i looked at it quickly here: https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/
wilcox_res_no_outliers_soilwater <- wilcox.test(AvgGrowth ~ SoilWater, tab_no_SoilWater_outliers) # still tiny p-value


  # trying no outliers for InoGrav too
tab_no_InoGrav_outliers <- tab %>% group_by(InoGrav) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_InoGrav_model <- lm(AvgGrowth ~ InoGrav, data = tab_no_InoGrav_outliers)
ggqqplot(residuals(no_outliers_InoGrav_model))
### this one looks a bit uglier than the last one (regarding us having no experience at looking at many of these:
ggqqplot(residuals(no_outliers_SoilWater_model))

shapiro_test(residuals(no_outliers_InoGrav_model))
  # not normal
plot(no_outliers_InoGrav_model, 1)
no_outliers_InoGrav_model %>% levene_test(AvgGrowth ~ InoGrav)
  # meets the variance criteria, so kruskal wallis
kruskal_wallis_InoGrav_no_outliers_anova_result <- tab_no_InoGrav_outliers %>% kruskal_test(AvgGrowth ~ InoGrav)
  # still not significant
### trying the wilcoxon test, to see if it's the same 
wilcox_res_no_outliers_InoGrav <- wilcox.test(AvgGrowth ~ InoGrav, tab_no_InoGrav_outliers) # yep, came up the same


  # now I'm going to try the interactions thing (though I'm going into it a little blind)
### you and me both!

tab %>% group_by(InoGrav, SoilWater) %>% summarize(n=n())
tab %>% group_by(InoGrav, Tank.) %>% summarize(n=n())
tab %>% group_by(Tank., SoilWater) %>% summarize(n=n())
  # InoGrav and SoilWater
tab %>% group_by(InoGrav, SoilWater) %>% summarize(n=n())
tab %>% group_by(InoGrav, SoilWater) %>% get_summary_stats(AvgGrowth, type = "mean_sd")
ggboxplot(tab, x = "InoGrav", y = "AvgGrowth", color = "SoilWater")
tab %>% group_by(InoGrav, SoilWater) %>% identify_outliers(AvgGrowth)
  # all of the extreme outliers are in tank 4, the rest in tank 1, but I'm not gonna remove them for now
InoGrav_SoilWater_model <- lm(AvgGrowth ~ InoGrav * SoilWater, data = tab)
ggqqplot(residuals(InoGrav_SoilWater_model))
shapiro_test(residuals(InoGrav_SoilWater_model))
tab %>% group_by(InoGrav, SoilWater) %>% shapiro_test(AvgGrowth)
  # not normal overall, only one is Inoc/Water is normal
tab %>% group_by(InoGrav, SoilWater, Tank.) %>% shapiro_test(AvgGrowth)
  # when factoring in tank, only tank4/gravel/water and tank1/inoc/soil have p below .05
 
  # okay now doing the robust anova
install.packages("walrus")
library(walrus)
ranova(tab, "AvgGrowth", c("InoGrav"), ph=TRUE)
ranova(tab, "AvgGrowth", c("SoilWater"), ph=TRUE) 
ranova(tab, "AvgGrowth", c("Tank."), ph=TRUE) 
 # again, the soil and water are very significantly different. there are a few tank combos that are, but no tanks are different from all others.
ranova(tab, "AvgGrowth", c("InoGrav", "SoilWater"), ph=TRUE)  # there's something about this that the program didn't like, not sure what it is though.
### hmm, i think maybe it needs one of the factors to have more than 2 groups, which isn't the case when we are looking at this one
tab %>% group_by(InoGrav, SoilWater) %>% summarize(n=n())

ranova(tab, "AvgGrowth", c("InoGrav", "Tank."), ph=TRUE)  # it says incomplete design
### ah, this makes sense because we don't have each tank in each 
tab %>% group_by(InoGrav, Tank.) %>% summarize(n=n())

ranova(tab, "AvgGrowth", c("Tank.", "SoilWater"), ph=TRUE) # okay idk what's happening, I think it's doing them separately. 
tab %>% group_by(SoilWater, Tank.) %>% summarize(n=n())
### this one has a soil and a water in each tank, so i think that's why it works
### it's still checking each individually, but it's checking them with additional information going into the model
### meaning, if we look at just soilwater alone again
ranova(tab, "AvgGrowth", c("SoilWater"), ph=TRUE) 
### the p-value is < .0000001
### if we look at the soil/water p-value from the combined one:
ranova(tab, "AvgGrowth", c("Tank.", "SoilWater"), ph=TRUE)
### we get a p-value of 0.0001171
### i don't know which would be the better way to go

### i say let's forget the ranova and worrying about interactions business, since I haven't been able to dig in and figure it out either
### this 2018 dataset is simpler in design fortunately, so no biggie here, and i think we're good sticking with the single contrasts using the z-standardization in the 2017 dataset


  # okay I think I'm going to stop here for now. 

tab %>% group_by(InoGrav, Tank.) %>% summarize(n=n())


