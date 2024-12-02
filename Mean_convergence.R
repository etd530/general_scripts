data <- read.csv("test_results.csv")

measures_num = vector()
means_list = vector()
sd_list = vector()
subset_mean_lens = vector()
mean_of_means = vector()

for (i in 1:nrow(data)){
  for (n in 1:100){
    rownums = sample(nrow(data), i, replace = T)
    subset = data$Length[rownums]
    subset_mean_lens[length(subset_mean_lens)+1] = mean(subset)
    measures_num[length(measures_num)+1] <- i
  }
  mean_of_means[length(mean_of_means)+1] <- mean(subset_mean_lens[(length(subset_mean_lens)-99):length(subset_mean_lens)])
  means_list[length(means_list)+1] <- mean(subset_mean_lens)
  sd_list[length(sd_list)+1] <- sd(subset_mean_lens)
}

means_list

plot(measures_num, subset_mean_lens)
points(mean_of_means, pch = 17, col = "blue")
abline(a = mean_of_means[length(mean_of_means)], b = 0)
