data<-read.csv("13.DatasetMCI5AD2_TRAIN.csv")
wilcox.test(data[which(data$Classification=="MCI-p"),]$AGE,
            data[which(data$Classification=="MCI-s"),]$AGE,)

table(data$PTGENDER)
# 1 male, 2 female 
# Create a contingency table
data$PTGENDER<-ifelse(data$PTGENDER==1,"Male","Female")
contingency_table <- table(data$PTGENDER, data$Classification)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)


library(ggplot2)

ggplot(data, aes(x = PTGENDER, fill = Classification)) +
  geom_bar(position = "fill") +
  labs(title = "Classification Distribution by PTGENDER",
       y = "Proportion",
       fill = "Classification") +
  theme_minimal()

ggplot(data, aes(x = Classification, fill = PTGENDER)) +
  geom_bar(position = "fill") +
  labs(title = "Classification Distribution by PTGENDER",
       y = "Proportion",
       fill = "PTGENDER") +
  theme_minimal()
