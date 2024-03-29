# install.packages("cranlogs")
library(cranlogs)
library(ggplot2)
library(dplyr)
library(tibble)

# Last week's downloads
cran_downloads(packages="nimbleAPT", when="last-week")


# How many overall downloads
mls <- cran_downloads(packages="nimbleAPT", from = "2021-11-11", to = Sys.Date()-1)
mls
sum(mls[,2])

# Plot
gr0 <- ggplot(mls, aes(date, count)) + geom_line(colour = "red",size=1)
gr0 + xlab("Time") + ylab("Nr. of downloads") + labs(title = paste0("nimbleAPT daily downloads ", Sys.Date()-1))

# Cumulate downloads
cumulative <- cumsum(mls[,2])
mls2 <- cbind(mls,cumulative)
# Date range
tibble(mls2) %>% filter(count>0) %>% select(date) %>% head(n=1) %>% pull() -> tStart
tibble(mls2) %>% filter(count>0) %>% select(date) %>% tail(n=1) %>% pull() -> tStop
# Plot
gr1 <- ggplot(mls2, aes(date, cumulative)) + geom_line(colour = "blue",size=1)
gr1 + xlab("Time") + ylab("Cumulative downloads") +
  labs(title = paste0("nimbleAPT downloads: ", tStart, " to ", tStop))
