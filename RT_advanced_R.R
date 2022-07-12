library(tidyverse)
tidyverse_logo()
library(datasets)
install.packages("gapminder")
library(gapminder)
attach(iris)
view(iris)
dim(iris)
iris %>% filter(Species=="setosa")
sepal <- select(iris,Sepal.Length,Sepal.Width)
sepal <- iris%>%select(Sepal.Length,Sepal.Width)
##to remove a column
less_sepal <- select(sepal,-Sepal.Length)
less_sepal
##Renaming a column
tidy_sepal <- iris%>%rename(Sepal.Width=sepal.width)
tidy_sepal1 <- iris%>%rename(sepal.width=Sepal.Width)
tidy_sepal1
##filter
iris_virginica <- iris %>%filter(Species=="virginica")
iris_virginica_sepal <- iris %>%filter(Species=="virginica",Sepal.Length>6) 
iris_virginica_sepal
iris_setosa_sepal <- iris%>%filter(Species=="setosa")%>%select(Sepal.Length,Sepal.Width)
iris_setosa_sepal

#Exercise
df_species <- iris%>%filter(Sepal.Length>5)%>%select(Species,Sepal.Width)
df_species

##groupby
##check properties of the dataset
str(iris)
str(iris%>%group_by(Species))
##summarize
###group per species and get the mean of Sepal.Width
gdp_byspecies <- iris%>%group_by(Species)%>%summarise(mean_species=mean(Sepal.Width))
gdp_byspecies

###group per two variable
gdp_byspecies1 <- iris%>%group_by(Species,Petal.Length)%>%summarise(mean_sepal_width=mean(Sepal.Width))
gdp_byspecies1

###
gdp_byspecies2 <- iris %>% group_by(Species, Petal.Length) %>% summarize(mean_sepal_width = mean(Sepal.Width), mean_petal_length = mean(Petal.Length), sd_sepal_length = mean(Sepal.Length))
gdp_byspecies2

###Use group_by(), summarize(), mean(), sd(), min(), max() to calculate the mean, standard deviation, get maximum value, minimum value of each Species’ Sepal.Width
gdp_byspecies3 <- iris %>% group_by(Species) %>% summarise(mean_sepal_width = mean(Sepal.Width),sd_sepal_width=sd(Sepal.Width),min_sepal_width=min(Sepal.Width),max_sepal_width=max(Sepal.Width))
gdp_byspecies3

##mutate--updates or creates new variables of a dataframe
iris_SLMm<-iris %>% mutate(Sepal.Length=Sepal.Length*10) 
iris_sl <- iris %>% mutate(SLMm=Sepal.Length*10) 
###
mutate_iris <- iris%>%group_by(Species)%>% mutate(SPlength = Sepal.Length/Petal.Length)%>%summarise(mean_SPlength=mean(SPlength),sd_SPlength=sd(SPlength),min_SPlength=min(SPlength),max_SPlength=max(SPlength))
mutate_iris

##Create a subset of iris using Sepal.Length>5 as a condition
iris_small <- iris %>% filter(Sepal.Length > 5)

##Create a scatter plot that compares petal width and length 
ggplot(iris_small, aes(x=Petal.Length, y=Petal.Width)) + geom_point()

##Adding color
ggplot(iris_small, aes(x=Petal.Length, y=Petal.Width, color=Species)) + geom_point() 

##Using different size for different variable
ggplot(iris_small, aes(x=Petal.Length, y=Petal.Width, color=Species, size=Sepal.Length)) + geom_point()
  
ggplot(iris_small, aes(x=Petal.Length, y=Petal.Width, color=Species, size=Sepal.Length)) + geom_point()+facet_wrap(~Species)

###Line plot
attach(gapminder)
By_year <- gapminder%>% group_by(year) %>% summarise(medianGdpperCap=median(gdpPercap))
view(By_year)
ggplot(By_year,aes(x=year,y=medianGdpperCap))+geom_line()

ggplot(By_year,aes(x=year,y=medianGdpperCap))+geom_line()+expand_limits(y=0)

##Barplot
by_species <- iris %>% filter(Sepal.Length>6) %>% group_by(Species) %>% summarize(medianPL=median(Petal.Length))
by_species
ggplot(by_species, aes(x=Species, y=medianPL)) + geom_col()

ggplot(iris_small, aes(x=Petal.Length))+ geom_histogram()

ggplot(iris_small, aes(x=Species, y=Sepal.Width))+ geom_boxplot()

view(iris_small)
install.packages("rmarkdown")
