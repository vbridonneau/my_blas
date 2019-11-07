times <- read.csv("times")
library(ggplot2)
mean <- aggregate(.~size, times[, 1-4], "mean")
min <- aggregate(.~size, times[, 1-4], "min")
max <- aggregate(.~size, times[, 1-4], "max")
stda <- aggregate(.~size, times[, 1-4], "sd")*sqrt(5/6)
names(stda)[names(stda) == "time"] = "std"
stat <- cbind(mean, stda[, 2])
names(stat)[names(stat) == "stda[, 2]"] = "std"
print(stat)
ggplot(stat, aes(x = size, y = time)) + 
	     geom_line() + 
	     geom_ribbon(aes(ymin = time - std,
                       	     ymax = time + std), alpha = 0.2)

#cover_min <- apply(mean, 2, "
#plot(mean, sub="Problem scaling", xlab="Vector size", ylab="GFlop/s", type="l", col="red")
#lines(min , sub="Problem scaling", xlab="Vector size", ylab="GFlop/s", type="l", col="blue")
#lines(max , sub="Problem scaling", xlab="Vector size", ylab="GFlop/s", type="l", col="green")
#lines(std , sub="Problem scaling", xlab="Vector size", ylab="GFlop/s", type="l", col="green")