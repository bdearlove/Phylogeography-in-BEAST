rootstate <- read.table('../output/location.states.log',skip=2,header=TRUE)
burnin <- 200
numstates <- dim(rootstate)[[1]]
rootstate.trim <- rootstate[burnin:numstates,]
as.data.frame.table(table(rootstate.trim$location)/(numstates-burnin))