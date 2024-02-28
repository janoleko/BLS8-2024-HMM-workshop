load("./data/elephant_data_prepped.RData")

# deleting moveHMM columns
newdata = data[,-c(1,2,3,8)]
# reordering columns
newdata = newdata[,c(3,4,5,1,2)]
# renaming coordinate columns
colnames(newdata)[4:5] = c("location.long", "location.lat")

head(newdata)

newdata = newdata[which(data$animalID=="2" & data$ID == "5"),]


write.csv(newdata, "./data/elephant_data.csv", row.names = FALSE)



