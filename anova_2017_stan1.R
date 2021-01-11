library(tidyverse)
library(ggpubr)
library(car)
library(rstatix)
tab <- read.csv(file = 'Plant-Growth-2017-Standardized.csv')
str(tab)
tab %>% group_by(Tank.) %>% get_summary_stats(StanGrowth, type = "mean_sd")
ggboxplot(tab, x = "Tank.", y = "StanGrowth")
tab %>% group_by(Tank.) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(Tank.) %>% identify_outliers(StanGrowth)
    #there is one extreme outlier in tank 9, i'll ignore it for now

    ## Tanks ##
tank._model <- lm(StanGrowth ~ Tank., data = tab)
ggqqplot(residuals(tank._model))
shapiro_test(residuals(tank._model))  
ggqqplot(tab, "StanGrowth", facet.by = "Tank.")
tab %>% group_by(Tank.) %>% shapiro_test(StanGrowth) #all normal except for tank 9. I'll do both anova and kruskal wallis
plot(tank._model, 1)
tab %>% levene_test(StanGrowth ~ Tank.) #okay assumption not met, so welch

welch_tank._anova_result <- tab %>% welch_anova_test(StanGrowth ~ Tank.) 
kruskal_wallis_tank._anova_result <- tab %>% kruskal_test(StanGrowth ~ Tank.)
    #basically said that the tanks themselves have significant variance, probably because of tanks 5 and 7
tank._pwc <- tab %>% tukey_hsd(StanGrowth ~ Tank.)
tank._pwc %>% filter(p.adj <= 0.10)
    # yeah so basically tanks 5 and 7 are going a bit crazy

    #Tanks no Outliers
tab_no_tank._outliers <- tab %>% group_by(Tank.) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_tank._model <- lm(StanGrowth ~ Tank., data = tab_no_tank._outliers)
ggqqplot(residuals(no_outliers_tank._model))
shapiro_test(residuals(no_outliers_tank._model)) #still not normal..
plot(no_outliers_tank._model, 1)
no_outliers_tank._model %>% levene_test(AvgGrowth ~ Tank.) #nope
kruskal_wallis_tank._no_outliers_anova_result <- tab_no_tank._outliers %>% kruskal_test(StanGrowth ~ Tank.) #significant
welch_tank._no_outlier_anova_result <- tab_no_tank._outliers %>% welch_anova_test(AvgGrowth ~ Tank.) #significant


    ## InoGrav ##
tab %>% group_by(InoGrav) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(InoGrav) %>% identify_outliers(StanGrowth) # one tank 5 xtreme, two tank 7 xtreme
InoGrav_model <- lm(StanGrowth ~ InoGrav, data = tab)
ggqqplot(residuals(InoGrav_model))
shapiro_test(residuals(InoGrav_model))  
ggqqplot(tab, "StanGrowth", facet.by = "InoGrav")
tab %>% group_by(InoGrav) %>% shapiro_test(StanGrowth) #inoculant not normal, but probably tank 5 and 7 striking again.
plot(InoGrav_model, 1)
tab %>% levene_test(StanGrowth ~ InoGrav) #assumption not met, so welch

welch_InoGrav_anova_result <- tab %>% welch_anova_test(StanGrowth ~ InoGrav) 
kruskal_wallis_InoGrav_anova_result <- tab %>% kruskal_test(StanGrowth ~ InoGrav)
        #nada

    #InoGrav no Outliers
tab_no_InoGrav_outliers <- tab %>% group_by(InoGrav) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_InoGrav_model <- lm(StanGrowth ~ InoGrav, data = tab_no_InoGrav_outliers)
ggqqplot(residuals(no_outliers_InoGrav_model))
shapiro_test(residuals(no_outliers_InoGrav_model)) #still not normal..
plot(no_outliers_InoGrav_model, 1)
no_outliers_InoGrav_model %>% levene_test(AvgGrowth ~ InoGrav) #homegeneity secured
kruskal_wallis_InoGrav_no_outliers_anova_result <- tab_no_InoGrav_outliers %>% kruskal_test(StanGrowth ~ InoGrav) # not significant
welch_InoGrav_no_outlier_anova_result <- tab_no_InoGrav_outliers %>% welch_anova_test(AvgGrowth ~ InoGrav) # not significant (just running the same code for all of them haha)



        ## Plant ##
tab %>% group_by(Plant) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(Plant) %>% identify_outliers(StanGrowth) # one tank 7 xtreme
Plant_model <- lm(StanGrowth ~ Plant, data = tab)
ggqqplot(residuals(Plant_model))
shapiro_test(residuals(Plant_model))  
ggqqplot(tab, "StanGrowth", facet.by = "Plant")
tab %>% group_by(Plant) %>% shapiro_test(StanGrowth) #basil and peas not normal, but probably tank 5 and 7 striking again.
plot(Plant_model, 1)
tab %>% levene_test(StanGrowth ~ Plant) #assumption not met, so welch

welch_Plant_anova_result <- tab %>% welch_anova_test(StanGrowth ~ Plant) 
kruskal_wallis_Plant_anova_result <- tab %>% kruskal_test(StanGrowth ~ Plant)
        # both have p < .05!
Plant_pwc <- tab %>% tukey_hsd(StanGrowth ~ Plant)
Plant_pwc %>% filter(p.adj <= 0.10)
        #not sure if I needed to do this, but the significant difference is between basil and tomatoes

    #Plant no outliers

tab_no_Plant_outliers <- tab %>% group_by(Plant) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_Plant_model <- lm(StanGrowth ~ Plant, data = tab_no_Plant_outliers)
ggqqplot(residuals(no_outliers_Plant_model))
shapiro_test(residuals(no_outliers_Plant_model)) #still not normal..
plot(no_outliers_Plant_model, 1)
no_outliers_Plant_model %>% levene_test(AvgGrowth ~ Plant) #nope
kruskal_wallis_Plant_no_outliers_anova_result <- tab_no_Plant_outliers %>% kruskal_test(StanGrowth ~ Plant) #significant
welch_Plant_no_outlier_anova_result <- tab_no_Plant_outliers %>% welch_anova_test(AvgGrowth ~ Plant) #very significant


        ## Lighting ##
tab %>% group_by(Lighting) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(Lighting) %>% identify_outliers(StanGrowth) # nothing extreme
Lighting_model <- lm(StanGrowth ~ Lighting, data = tab)
ggqqplot(residuals(Lighting_model))
shapiro_test(residuals(Lighting_model))  
ggqqplot(tab, "StanGrowth", facet.by = "Lighting")
tab %>% group_by(Lighting) %>% shapiro_test(StanGrowth) #artificial not normal, but probably tank 5 and 7 striking again.
plot(Lighting_model, 1)
tab %>% levene_test(StanGrowth ~ Lighting) #assumption not met, so welch

welch_Lighting_anova_result <- tab %>% welch_anova_test(StanGrowth ~ Lighting) 
kruskal_wallis_Lighting_anova_result <- tab %>% kruskal_test(StanGrowth ~ Lighting)
        #welch shows significance, but kruskal is not.

    #Lighting no outliers

tab_no_Lighting_outliers <- tab %>% group_by(Lighting) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_Lighting_model <- lm(StanGrowth ~ Lighting, data = tab_no_Lighting_outliers)
ggqqplot(residuals(no_outliers_Lighting_model))
shapiro_test(residuals(no_outliers_Lighting_model)) #still not normal..
plot(no_outliers_Lighting_model, 1)
no_outliers_Lighting_model %>% levene_test(AvgGrowth ~ Lighting) #nope
kruskal_wallis_Lighting_no_outliers_anova_result <- tab_no_Lighting_outliers %>% kruskal_test(StanGrowth ~ Lighting) #not significant :(
welch_Lighting_no_outlier_anova_result <- tab_no_Lighting_outliers %>% welch_anova_test(AvgGrowth ~ Lighting) #significant

        ## Media ##
tab %>% group_by(Media) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(Media) %>% identify_outliers(StanGrowth) # no extreme outliers
Media_model <- lm(StanGrowth ~ Media, data = tab)
ggqqplot(residuals(Media_model))
shapiro_test(residuals(Media_model))  
ggqqplot(tab, "StanGrowth", facet.by = "Media") #wow clay looks wack, and when I look at the variables, it's tank 5 and 7 again!
tab %>% group_by(Media) %>% shapiro_test(StanGrowth) #clay and float not normal, so it's not just tanks 5 and 7 this time
plot(Media_model, 1)
tab %>% levene_test(StanGrowth ~ Media) #assumption not met, so welch

welch_Media_anova_result <- tab %>% welch_anova_test(StanGrowth ~ Media) 
kruskal_wallis_Media_anova_result <- tab %>% kruskal_test(StanGrowth ~ Media)
        #welch signficant, kruskal not

    #Media no outliers

tab_no_Media_outliers <- tab %>% group_by(Media) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_Media_model <- lm(StanGrowth ~ Media, data = tab_no_Media_outliers)
ggqqplot(residuals(no_outliers_Media_model))
shapiro_test(residuals(no_outliers_Media_model)) #still not normal..
plot(no_outliers_Media_model, 1)
no_outliers_Media_model %>% levene_test(AvgGrowth ~ Media) #nope
kruskal_wallis_Media_no_outliers_anova_result <- tab_no_Media_outliers %>% kruskal_test(StanGrowth ~ Media) # not significant
welch_Media_no_outlier_anova_result <- tab_no_Media_outliers %>% welch_anova_test(AvgGrowth ~ Media) #not significant (now both of them aren't)


        ## Water ##
tab %>% group_by(Water) %>% get_summary_stats(StanGrowth, type = "mean_sd")
tab %>% group_by(Water) %>% identify_outliers(StanGrowth) # no extreme outliers
Water_model <- lm(StanGrowth ~ Water, data = tab)
ggqqplot(residuals(Water_model))
shapiro_test(residuals(Water_model))  
ggqqplot(tab, "StanGrowth", facet.by = "Water")
tab %>% group_by(Water) %>% shapiro_test(StanGrowth) #pumps not normal, but probably tank 5 and 7 striking again.
plot(Water_model, 1)
tab %>% levene_test(StanGrowth ~ Water) #assumption not met, so welch

welch_Water_anova_result <- tab %>% welch_anova_test(StanGrowth ~ Water) 
kruskal_wallis_Water_anova_result <- tab %>% kruskal_test(StanGrowth ~ Water)
        #nada, but welch is close.

    # Water no outliers

tab_no_Water_outliers <- tab %>% group_by(Water) %>% filter(!(AvgGrowth > ((3 * IQR(AvgGrowth)) + quantile(AvgGrowth, probs = 0.75))) | (AvgGrowth < ((3 * IQR(AvgGrowth)) - quantile(AvgGrowth, probs = 0.25)))) %>% data.frame()
no_outliers_Water_model <- lm(StanGrowth ~ Water, data = tab_no_Water_outliers)
ggqqplot(residuals(no_outliers_Water_model))
shapiro_test(residuals(no_outliers_Water_model)) #still not normal..
plot(no_outliers_Water_model, 1)
no_outliers_Water_model %>% levene_test(AvgGrowth ~ Water) #nope
kruskal_wallis_Water_no_outliers_anova_result <- tab_no_Water_outliers %>% kruskal_test(StanGrowth ~ Water) # not significant
welch_Water_no_outlier_anova_result <- tab_no_Water_outliers %>% welch_anova_test(AvgGrowth ~ Water) #not significant (both p>0.9)


        ## just lookin ##
ggboxplot(tab, x = "Tank.", y = "StanGrowth")
ggboxplot(tab, x = "InoGrav", y = "StanGrowth")
ggboxplot(tab, x = "Media", y = "StanGrowth")
ggboxplot(tab, x = "Plant", y = "StanGrowth")
ggboxplot(tab, x = "Lighting", y = "StanGrowth")
ggboxplot(tab, x = "Water", y = "StanGrowth")





