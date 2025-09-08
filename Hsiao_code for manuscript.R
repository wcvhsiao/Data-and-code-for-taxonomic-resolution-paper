rm(list=ls())

#################################
###Draw sst map around Taiwan###
################################

if (!require("rerddap")) install.packages("rerddap"); library("rerddap")
if (!require("akima")) install.packages("akima"); library("akima")
if (!require("dplyr")) install.packages("dplyr"); library("dplyr")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("mapdata")) install.packages("mapdata"); library("mapdata")
if (!require("ncdf4")) install.packages("ncdf4"); library("ncdf4")
if (!require("plot3D")) install.packages("plot3D"); library("plot3D")

browse('NOAA_DHW_monthly')

longitude = c(119, 124)
latitude = c(21, 27)
time = c("1994-8-16","2023-7-16")
sstInfo <- info('NOAA_DHW_monthly')
DHWSST <- griddap(sstInfo, latitude = c(latitude),longitude = c(longitude), time = time, fields = 'sea_surface_temperature')
mycolor <- colors$temperature
data<-DHWSST$data
df_aggregated <- aggregate(sea_surface_temperature ~ latitude + longitude, data = data, mean)


if (!require('rnaturalearth')) install.packages('rnaturalearth'); library("rnaturalearth")
# Load country shape data
world <- ne_countries(scale = "large", returnclass = "sf")


windows()
ggplot() +
  geom_tile(data = df_aggregated, aes(x = longitude, y = latitude, fill = sea_surface_temperature)) +
  scale_fill_gradientn(colours = mycolor, na.value = NA) +
  geom_sf(data = world,color = NA) +
  coord_sf(xlim = longitude, ylim = latitude, expand = FALSE)


################################
###Measure demographic rates###
###############################

#Load in data for all three genera

area_poci <- read.csv("set working directory/poci_area.csv")
area_pori<-read.csv('set working directory/pori_area.csv')
area_acro<-read.csv('set working directory/acro_area.csv')

#Remove data in Pocillopora for unidentified individuals 
#(including recruits and failed sequenced data)
area_poci <- area_poci[-which(area$Species == ""), ]

#Calculate percentage growth rate and Arithmetric mean radius
growth_2023<-(area$Jul2023-area$Jul2022)/area$Jul2022
area$amr<-(area$Jul2023/area$Jul2023_perim)-(area$Jul2022/area$Jul2022_perim)
area<-cbind(area,growth_2023,amr)


amr_acro<-(area_acro$Jul2023/area_acro$Jul2023_perim)-(area_acro$Jul2022/area_acro$Jul2022_perim)
amr_pori<-(area_pori$Jul2023/area_pori$Jul2023_perim)-(area_pori$Jul2022/area_pori$Jul2022_perim)
genus<-rep(c("poci",'acro','pori'),c(length(amr_poci),length(amr_acro),length(amr_pori)))
Genus_gro<-rbind(amr_poci,amr_acro,amr_pori)
Genus_gro<-cbind(genus,Genus_gro)

ggplot(area,aes(x=Haplotype,y=amr),fill=area$species)+
  +   geom_boxplot()

ggplot(Genus_gro,aes(x=genus,y=amr),fill=Genus_gro$genus)+
  +   geom_boxplot()

#Calculate survival rate and present as binary (1=alive,0=dead)
area$survival<-area$Jul2023
require(dplyr)
area <- area %>%
  mutate(survival = ifelse(survival == 'NA',0,1))
area$survival[is.na(area$survival)] <- 0  

area$survival<-jitter(area$survival, amount=0.05)


# Survival curve of different Pocillopora species in accordance to their size 
p<-ggplot(area,aes(log10(Jul2022),survival,color=area$Species))+ 
  geom_smooth(method="glm",method.args=list(family='binomial'),formula=y~x, alpha=0.2,linewidth=2,aes(fill=Species))+ 
  geom_point(position=position_jitter(height=0.03,width=0),alpha=0.8)+ 
  scale_color_manual(values=c('darkblue',"deepskyblue",'beige', "tan",'black','lightgreen', "coral"))+
  scale_fill_manual(values=c('darkblue',"deepskyblue",'beige', "tan",'black','lightgreen', "coral"))+
  xlab("Log(area)")+ylab("Pr(survived)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

#ANOVA testing for vital rates difference among taxonomic resolution/region 
#percentage growth rate
gro_anova<-aov(area$growth_2023~area$Region)
summary(gro_anova)  
TukeyHSD(gro_anova)

#AMR
amr_anova_hp<-aov(area$amr~area$Haplotype)
summary(amr_anova_hp)
TukeyHSD(amr_anova_hp)
amr_anova_sp<-aov(area$amr~area$Species)
summary(amr_anova_sp)
TukeyHSD(amr_anova_sp)

#survival
sur_anova_sp<-aov(area$survival~area$Species)
summary(sur_anova_sp)
TukeyHSD(sur_anova_sp)

sur_anova_reg<-aov(area$survival~area$Region)
summary(sur_anova_reg)
TukeyHSD(sur_anova_reg)

sur_anova_re<-aov(area$survival~area$reeftype)
summary(sur_anova_re)
TukeyHSD(sur_anova_re)

#Summary descriptive states for growth on haplotypes/species/genus 

area_summary_hp <- area %>% 
  filter(!is.na(growth_2023))%>%
  group_by(Haplotype) %>%
  summarise(mean(growth_2023), min(growth_2023), max(growth_2023),sd(growth_2023))

area_summary_sp <- area %>% 
  filter(!is.na(growth_2023))%>%
  group_by(Species) %>%
  summarise(mean(growth_2023), min(growth_2023), max(growth_2023),sd(growth_2023))

area$survival<-jitter(area$survival)
area_survival_hp <- area%>% 
  filter(!is.na(survival))%>%
  group_by(Haplotype) %>%
  summarise(mean(survival), sd(survival))

area_survival_sp <- area%>% 
  filter(!is.na(survival))%>%# Get min by group
  group_by(Species) %>%
  summarise(mean(survival), sd(survival))




##HCPC and FAMD##
FAMD<-read.csv('genus_traits.csv')
FAMD$Reproductive.strategy<-as.factor(FAMD$Reproductive.strategy)
FAMD$Species<-as.factor(FAMD$Species)
FAMD$genus<-as.factor(FAMD$genus)
FAMD$Haplotype<-as.factor(FAMD$Haplotype)
FAMD$morphology<-as.factor(FAMD$morphology)
FAMD$region<-as.factor(FAMD$region)
res.famd <- FAMD(FAMD,graph = FALSE)
get_eigenvalue(res.famd)
res.hcpc <- HCPC(res.famd, graph = FALSE)

res.hcpc1<-HCPC(res.famd,nb.clust=-1,min=3,max=NULL,graph=TRUE)

plot(res.hcpc, choice='3D.map',xlim=c(-2,5),ylim=c(-3,4))

# Plot of variables
fviz_famd_var(res.famd, repel = TRUE)
# Contribution to the first dimension
fviz_contrib(res.famd, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2)

fviz_famd_ind(res.famd, habillage = 'genus',label='none',
              addEllipses = TRUE, ellipse.type = "norm",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)+
  xlim(-4,5)+
  ylim(-3,5)
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_mfa_ind(res.famd, 
             habillage = "genus", # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "norm", 
             repel = TRUE # Avoid text overlapping
) 
