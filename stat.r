#Read in the data
setwd("/home/dave/sa")
datain = read.csv("data.csv")

mean_g = with(datain, tapply(gvector_absolute_fitness, param_set, mean))
sort(mean_g)
max(mean_g)

#psets <- sample(levels(datain$param_set), 20, replace = FALSE)
few = subset(datain, param_set %in% c("set78", "set276", "set348", "set417"))
few$param_set <- droplevels(few$param_set)
#boxplot(gvector_absolute_fitness ~ param_set, data = s)
plotmeans(gvector_absolute_fitness ~ param_set, data = few, ylim = c(0, 11))
m = glm(gvector_absolute_fitness ~ num_epochs + bumps + c2 + num_organisms + sa_timesteps + activation_types + crossover_freq, data = datain, family = gaussian)
anove(m)
summarize(m)