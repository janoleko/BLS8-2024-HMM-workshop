# Code for preparing the movebank data for the HMM analysis
# download at:
# https://www.datarepository.movebank.org/handle/10255/move.670

## setup----
library(lubridate) # version 1.8.0
library(moveHMM) # version 1.8
library(dplyr) # version 1.1.0

## read data----
data <- read.csv2("Movement syndromes across vertebrate taxa (data from Abrahms et al. 2017)-gps.csv", header = T, sep = ",")
# only use the elephants:
data <- data[data$individual.taxon.canonical.name == "Loxodonta africana", ]



## variable cleanup----
# remove some unnecessary variables
data <- subset(data, select = -c(visible, modelled, sensor.type, individual.taxon.canonical.name, study.name, tag.local.identifier))

# let event.id count from 1
data$event.id <- data$event.id - min(data$event.id) + 1

# make locations numeric
data$location.long <- as.numeric(data$location.long)
data$location.lat <- as.numeric(data$location.lat)

# timestamp variable in a more useful format:
data$timestamp <- strptime(data$timestamp, "%Y-%m-%d %H:%M:%S", tz = "GMT")
# separate date and time of day
data$date <- as.Date(data$timestamp)
data$tod <- as.numeric(format(data$timestamp, "%H"))
data$tod[data$tod == 0] <- 24 # doesn't actually make any difference, it's just easier for me this way


## irregular time series! separate into tracks----
# create a variable that measures the time between two observations
data <- data %>%
  arrange(individual.local.identifier, timestamp()) %>%
  group_by(individual.local.identifier) %>%
  mutate(diff = timestamp - lag(timestamp))
# set this variable to some high value when a time series of a new elephant begins
data$diff[which(data$individual.local.identifier != lag(data$individual.local.identifier))] <- 100

# I could impute missing observations for gaps of less than, say, 5 hours, but these don't happen.
# Instead, I'll separate the data into 152 tracks.

# variable with observations at which a new track should begin (when the lag is more than an hour)
new_tracks <- which(data$diff != 1)
# calculate lengths of the tracks
length_tracks <- new_tracks - lag(new_tracks)
length_tracks[1] <- new_tracks[1] - 1
length_tracks[length(length_tracks) + 1] <- (dim(data)[1] + 1) - new_tracks[length(new_tracks)]
# define track variable
data$track <- rep(c(1:length(length_tracks)), length_tracks)

# remove helper variables
rm(length_tracks)
rm(new_tracks)
data <- subset(data, select = -diff)


## moveHMM prep----
# rename elephant ID and track ID
data <- rename(data, animalID = individual.local.identifier)
data$animalID <- as.numeric(sub("elephant ", "", data$animalID))
data <- rename(data, ID = track)

# moveHMM doesn't seem to like variables of POSIX format, so I just use a subset with only necessary variables
data_step <- subset(data, select = c(ID, location.long, location.lat, event.id))
# get step lengths and turning angles
data_step <- prepData(data_step, coordNames = c("location.long", "location.lat")) # takes 1-2 minutes
# plot(data_step[data_step$ID==3,]) # could e.g. plot a track with moveHMM

# merge all variables together again
data_step$animalID <- data$animalID
data_step$tod <- data$tod
data_step$timestamp <- data$timestamp
data <- data_step
rm(data_step)

# calculate number of zero step lengths
no_zero <- sum(data$step == 0, na.rm = TRUE) # only 0.1% of the values
# substitute zero step lengths with small values, because they overcomplicate analysis otherwise
data$step[which(data$step == 0)] <- runif(no_zero, 0.0001, 0.00011)



## unrealistically long steps----
tail(sort(data$step)) # the two longest steps don't look right
# inspect and plot them
data[c(66876, 66877, 66878, 66879, 93273, 93274, 93275, 93276, 93277), ]
plot(data$x[data$ID == 148], data$y[data$ID == 148], type = "l")
plot(data$x[data$ID == 152], data$y[data$ID == 152], type = "l")

# I have no idea why these happen, but decide to just throw away the rest of the tracks after these long steps
data <- data[-c(66877:max(which(data$ID == 148)), 93275:max(which(data$ID == 152))), ]

# remove short tracks (<24 hours), because I don't know whether these might have a lot of GPS errors in them?
rl_trackID <- rle(data$ID) # run length encoding of track lengths
short_tracks <- rl_trackID$values[rl_trackID$lengths < 24] # identify short tracks
data <- data[-which(data$ID %in% short_tracks), ] # remove short tracks


## save data----
save(data, file = "elephant_data_prepped.RData")
