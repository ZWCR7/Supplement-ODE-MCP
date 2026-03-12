library(coronavirus)

data('coronavirus')

ItalyIndex = which(coronavirus$country == 'Italy')
covid19Italy = coronavirus[ItalyIndex, ]

ItalyTime = unique(covid19Italy$date)
ItalyConfirmed = covid19Italy$cases[which(covid19Italy$type == 'confirmed')]
ItalyDeath = covid19Italy$cases[which(covid19Italy$type == 'death')]
ItalyRecovry = covid19Italy$cases[which(covid19Italy$type == 'recovery')]

ItalyTime = ItalyTime[31:561]
ItalyIndex = 1:531
ItalyConfirmed = ItalyConfirmed[31:561]
ItalyDeath = ItalyDeath[31:561]
ItalyRecovry = ItalyRecovry[31:561]


negative_index = unique(c(which(ItalyConfirmed < 0), which(ItalyDeath < 0), which(ItalyRecovry < 0)))
ItalyTime = ItalyTime[-negative_index]
ItalyIndex = ItalyIndex[-negative_index]
ItalyConfirmed = ItalyConfirmed[-negative_index]
ItalyDeath = ItalyDeath[-negative_index]
ItalyRecovry = ItalyRecovry[-negative_index]

Population = unique(covid19Italy$population)

ItalyTotalConfirm = cumsum(ItalyConfirmed)
ItalyRemove = cumsum(ItalyDeath) + cumsum(ItalyRecovry)
ItalyInfect = ItalyTotalConfirm - ItalyRemove
ItalySuscep = Population - ItalyTotalConfirm

save(list = c('ItalyInfect', 'ItalySuscep', 'Population', 'ItalyIndex', 'ItalyTime'), file = 'covid19Italy.RData')
rm(list = ls())
##########################################################################################################################