ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("VCA","ggpubr","multcomp","ggridges","emmeans","broom","ggplot2", "dr4pl", "tidyverse","ggplot2", "drc", "lme4", "lsmeans", "plyr", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "ggsci","scales")
ipak(packages)


setwd(dir="C:/Users/mbreu/OneDrive - Michigan State University/desthio exp")

##Comparison Between Prothioconazole and Prothioconazole -Desthio in Poison-Plate Mycelial Growth Assays of Fusarium graminearum
##Plant health Progress 2022 doi.org/10.1094/PHP-06-21-0087-RS


rawrep1<-read_delim(file = "desthioRep1_091820.csv",delim=",",na=".")
rawrep2<-read_delim(file = "desthioRep2_091920.csv",delim=",",na=".")
rawrep3<-read_delim(file = "desthioRep3_092020.csv",delim=",",na=".")

View(screen1)
screen1<-rawrep1 %>% 
  filter(ppm==1) %>% 
  gather(key="Day", value="length",day5, day7) %>% 
  group_by(isolate,chemistry,ppm,Day) %>% 
  spread(measurement_rep,length) %>% 
  transmute(avg_length=(one+two)/2) %>% 
  transmute(avg_length_correct=(avg_length - 5)) %>% 
  mutate(run="1")
screen2<-rawrep2 %>% 
  filter(ppm==1) %>% 
  gather(key="Day", value="length",day5, day7) %>% 
  group_by(isolate,chemistry,ppm,Day) %>% 
  spread(measurement_rep,length) %>% 
  transmute(avg_length=(one+two)/2) %>% 
  transmute(avg_length_correct=(avg_length - 5)) %>% 
  mutate(run="2")

screen3<-rawrep3 %>% 
  filter(ppm==1) %>% 
  gather(key="Day", value="length",day5, day7) %>% 
  group_by(isolate,chemistry,ppm,Day) %>% 
  spread(measurement_rep,length) %>% 
  transmute(avg_length=(one+two)/2) %>% 
  transmute(avg_length_correct=(avg_length - 5)) %>% 
  mutate(run="3")
View(screen3)


totalscreen<-bind_rows(screen3,screen2,screen1)
View(totalscreen)
Day5<-totalscreen %>% 
  filter(Day=="day5") %>% 
  spread(chemistry,avg_length_correct)

Avg5<-Day5 %>% 
  group_by(isolate) %>% 
  summarise(desthioavg=mean(desthio, na.rm=TRUE),prothavg=mean(prothioconazole,na.rm=TRUE),prolineavg=mean(proline,na.rm=TRUE))
View(Avg5)

Totalsum<-Avg5 %>% 
  summarise(desthio=mean(desthioavg, na.rm=TRUE),proth=mean(prothavg,na.rm=TRUE),proline=mean(prolineavg,na.rm=TRUE))
View(Totalsum)


Day7<-totalscreen %>% 
  filter(Day=="day7") %>% 
  spread(chemistry,avg_length_correct)

cor.test(Avg5$desthioavg, Avg5$prothavg,method="pearson")
cor.test(Avg5$prolineavg, Avg5$prothavg,method="pearson")
cor.test(Avg5$prolineavg, Avg5$desthioavg,method="pearson")


###Scatterplots with correlation

a<-ggscatter(Avg5, x="prolineavg",y="prothavg",
             #color = "black", shape = 21, size = 3, # Points color, shape and size
             add = "reg.line",  # Add regressin line
             cor.method = "pearson",
             #add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE, # Add confidence interval
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             #ylim=c(3,16),
             #xlim=c(3,16),
             cor.coeff.args = list(method = "pearson", label.x = 3),#, label.sep = "\n")
             ggtheme=theme_bw()
)
b<-ggscatter(Avg5, x="desthioavg",y="prothavg",
             #color = "black", shape = 21, size = 3, # Points color, shape and size
             add = "reg.line",  # Add regressin line
             cor.method = "pearson",
             #add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE, # Add confidence interval
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             #ylim=c(3,16),
             #xlim=c(0,2),
             cor.coeff.args = list(method = "pearson"),#, label.sep = "\n")
             ggtheme=theme_bw()
)
c<-ggscatter(Avg5, x="desthioavg",y="prolineavg",
             #color = "black", shape = 21, size = 3, # Points color, shape and size
             add = "reg.line",  # Add regressin line
             cor.method = "pearson",
             #add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE, # Add confidence interval
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             #ylim=c(3,16),
             #xlim=c(0,2),
             cor.coeff.args = list(method = "pearson"),#, label.sep = "\n")
             ggtheme=theme_bw()
)

CorrFig <- ggarrange(a, b, c,
                     labels = c("A", "B", "C"),
                     ncol =3, nrow = 1)

ggsave(CorrFig, filename = "CrossSensitivity.JPEG",  bg = "transparent", dpi=600,units="in", height=5.5, width=12.5)

####subsetting data for EC50 determination##
#####

fullset<-c("24A","57G", "73S","89B","27A","21A","KMISO2-2-17")

Run3<-rawrep3 %>% 
  filter( isolate %in% fullset) %>% 
  filter (chemistry!="pda") %>% 
  filter (chemistry!="proline") %>% 
  select (-day7) %>% 
  mutate (run="3")
View(rawrep3)
Run1<-rawrep1 %>% 
  filter( isolate %in% set1) %>% 
  filter (chemistry!="pda") %>% 
  filter (chemistry!="proline") %>% 
  select (-day7) %>% 
  mutate (run="1")

Run2<-rawrep2 %>% 
  filter( isolate %in% fullset) %>% 
  filter (chemistry!="pda") %>% 
  filter (chemistry!="proline") %>% 
  select (-day7) %>% 
  mutate (run="2")


DataEC50<-bind_rows(Run1,Run2,Run3)

prothedited<-DataEC50 %>% 
  filter(chemistry=="prothioconazole" | ppm==0) %>% 
  spread(key="measurement_rep", value=day5) %>% 
  mutate(average_length=(((prothedited$one+prothedited$two)/2)-5))  %>% 
  group_by(isolate) %>% 
  mutate(relative=average_length/average_length[ppm==0])
View(proth)

desthioedited<-DataEC50 %>% 
  filter(chemistry=="desthio" | ppm==0) %>% 
  spread(key="measurement_rep", value=day5) %>% 
  mutate(average_length=(((desthioedited$one+desthioedited$two)/2)-5))  %>% 
  group_by(isolate) %>% 
  mutate(relative=average_length/average_length[ppm==0])


growthcurve.prothLL3<-drm(100 * relative ~ ppm, data = prothedited, curveid = isolate, fct = LL.3(), na.action = na.omit)
LL3tableproth50abs<-ED(growthcurve.prothLL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95)

growthcurve.desthioLL3<-drm(100 * relative ~ ppm, data = desthioedited, curveid = isolate, fct = LL.3(), na.action = na.omit)
LL3tabledesthio50abs<-ED(growthcurve.desthioLL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95)

isolateName<-c("21A","24A","27A","57G","73S","89B","Kmiso2217")
save<-as_tibble(LL3tableproth50abs)
savedes<-as_tibble(LL3tabledesthio50abs)

save$isolate=isolateName
savedes$isolate=isolateName

save$chem="prothioconazole"
savedes$chem="desthio"

newEc50<-bind_rows(save,savedes) %>% 
  select(-Upper,-Lower,-"Std. Error") %>% 
  spread(key="chem",value=Estimate)
View(newEc50)
load(EnvStats)


t.test(newEc50$desthio , newEc50$prothioconazole , paired = TRUE, alternative = "two.sided")

write.csv(newEc50,file="EC50desthio.csv")

cor.test(newEc50$desthio,newEc50$prothioconazole,method="spearman")

