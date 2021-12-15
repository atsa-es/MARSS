###################################################
### code chunk number 2: Cs00_required_libraries
###################################################
library(MARSS)
library(maps)


###################################################
### code chunk number 3: Cs01_data2
###################################################
loggerheadNoisy[1:6, ]


###################################################
### code chunk number 4: Cs02_turtles2
###################################################
turtles <- levels(loggerheadNoisy$turtle)
turtles


###################################################
### code chunk number 5: Cs03_data3
###################################################
turtlename <- "BigMama"
theTurtle <- which(loggerheadNoisy$turtle == turtlename)
dat <- loggerheadNoisy[theTurtle, 5:6]
dat <- t(dat) # transpose


###################################################
### code chunk number 6: Cs04_fig-code
###################################################
# load the map package; you have to install it first
library(maps)
# Read in our noisy data (no missing values)
pdat <- loggerheadNoisy # for plotting
turtlename <- "BigMama"
theTurtle <- which(loggerheadNoisy$turtle == turtlename)
par(mai = c(0, 0, 0, 0), mfrow = c(1, 1))
map("state",
  region = c(
    "florida", "georgia", "south carolina",
    "north carolina", "virginia", "delaware", "new jersey", "maryland"
  ),
  xlim = c(-85, -70)
)
points(pdat$lon[theTurtle], pdat$lat[theTurtle],
  col = "blue", pch = 21, cex = 0.7
)


###################################################
### code chunk number 7: Cs05_badtagfig
###################################################
# load the map package; you have to install it first
library(maps)
# Read in our noisy data (no missing values)
pdat <- loggerheadNoisy # for plotting
turtlename <- "BigMama"
par(mai = c(0, 0, 0, 0), mfrow = c(1, 1))
map("state", region = c(
  "florida", "georgia", "south carolina", "north carolina",
  "virginia", "delaware", "new jersey", "maryland"
), xlim = c(-85, -70))
points(pdat$lon[which(pdat$turtle == turtlename)], pdat$lat[which(pdat$turtle == turtlename)],
  col = "blue", pch = 21, cex = 0.7
)


###################################################
### code chunk number 8: Cs06_setup2
###################################################
Z.model <- "identity"
U.model <- "unequal"
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"


###################################################
### code chunk number 9: Cs07_fitmodel
###################################################
kem <- MARSS(dat, model = list(
  Z = Z.model,
  Q = Q.model, R = R.model, U = U.model
))


###################################################
### code chunk number 10: Cs08_code-to-plot
###################################################
# Code to plot estimated turtle track against observations
# The estimates
pred.lon <- kem$states[1, ]
pred.lat <- kem$states[2, ]

par(mai = c(0, 0, 0, 0), mfrow = c(1, 1))
library(maps)
pdat <- loggerheadNoisy
turtlename <- "BigMama"
map("state",
  region = c(
    "florida", "georgia", "south carolina",
    "north carolina", "virginia", "delaware", "new jersey", "maryland"
  ),
  xlim = c(-85, -70)
)
points(pdat$lon[theTurtle], pdat$lat[theTurtle],
  col = "blue", pch = 21, cex = 0.7
)
lines(pred.lon, pred.lat, col = "red", lwd = 2)

goodturtles <- loggerhead
gooddat <- goodturtles[which(goodturtles$turtle == turtlename), 5:6]
points(gooddat[, 1], gooddat[, 2], col = "black", lwd = 2, pch = 3, cex = 1.1)
legend("bottomright", c(
  "bad locations", "estimated true location",
  "good location data"
),
pch = c(1, -1, 3), lty = c(-1, 1, -1),
col = c("blue", "red", "black"), bty = FALSE
)


###################################################
### code chunk number 11: Cs09_figbigmama
###################################################
# Code to plot estimated turtle track against observations
# The estimates
pred.lon <- kem$states[1, ]
pred.lat <- kem$states[2, ]

par(mai = c(0, 0, 0, 0), mfrow = c(1, 1))
library(maps)
pdat <- loggerheadNoisy
turtlename <- "BigMama"
map("state", region = c(
  "florida", "georgia", "south carolina", "north carolina",
  "virginia", "delaware", "new jersey", "maryland"
), xlim = c(-85, -70))
points(pdat$lon[which(pdat$turtle == turtlename)],
  pdat$lat[which(pdat$turtle == turtlename)],
  col = "blue", pch = 21, cex = 0.7
)
lines(pred.lon, pred.lat, col = "red", lwd = 2)

goodturtles <- loggerhead
gooddat <- goodturtles[which(goodturtles$turtle == turtlename), 5:6]
points(gooddat[, 1], gooddat[, 2],
  col = "black", lwd = 2, pch = 3, cex = 1.1
)
legend("bottomright", c(
  "bad locations", "estimated true location",
  "good location data"
),
pch = c(1, -1, 3), lty = c(-1, 1, -1),
col = c("blue", "red", "black"), bty = FALSE
)


###################################################
### code chunk number 12: Cs09b_GCDF
###################################################
GCDF <- function(lon1, lon2, lat1, lat2, degrees = TRUE, units = "miles") {
  temp <- ifelse(degrees == FALSE,
    acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1)),
    acos(sin(lat1 / 57.2958) * sin(lat2 / 57.2958) +
      cos(lat1 / 57.2958) * cos(lat2 / 57.2958) *
        cos(lon2 / 57.2958 - lon1 / 57.2958))
  )
  r <- 3963.0 # (statute miles) , default
  if ("units" == "nm") r <- 3437.74677 # (nautical miles)
  if ("units" == "km") r <- 6378.7 # (kilometers)
  return(r * temp)
}


###################################################
### code chunk number 13: Cs10_distance (eval = FALSE)
###################################################
## distance[i - 1] <- GCDF(
##   pred.lon[i - 1], pred.lon[i],
##   pred.lat[i - 1], pred.lat[i]
## )


###################################################
### code chunk number 14: Cs11_distance2
###################################################
distance <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance[i - 1] <- GCDF(
    pred.lon[i - 1], pred.lon[i],
    pred.lat[i - 1], pred.lat[i]
  )
}


###################################################
### code chunk number 15: Cs13_disfig
###################################################
par(mfrow = c(1, 1))
hist(distance) # make a histogram of distance traveled per day


###################################################
### code chunk number 16: Cs14_disfig2
###################################################
# Compare to the distance traveled per day if you used the raw data
distance.noerr <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance.noerr[i - 1] <- GCDF(dat[1, i - 1], dat[1, i], dat[2, i - 1], dat[2, i])
}
hist(distance.noerr) # make a histogram of distance traveled per day


###################################################
### code chunk number 17: Cs15_mean.miles
###################################################
# accounting for observation error
mean(distance)
# assuming the data have no observation error
mean(distance.noerr)


###################################################
### code chunk number 18: Cs16_turt-names
###################################################
levels(loggerheadNoisy$turtle)


###################################################
### code chunk number 20: Cs18_code
###################################################
###############################################################
#    GCDF FUNCTION
#    This function converts units of degrees lat/lon to miles,
#    nautical miles, or kilometers
###############################################################
GCDF <- function(lon1, lon2, lat1, lat2, degrees = TRUE, units = "miles") {
  # This is the function for the Great Circle Distance Formula
  # using decimal degrees or radians
  # Calculations at: http://www.nhc.noaa.gov/gccalc.shtml
  # This first component is only dependent on degrees or radians
  temp <- ifelse(degrees == FALSE, acos(sin(lat1) *
    sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1)),
  acos(sin(lat1 / 57.2958) * sin(lat2 / 57.2958) + cos(lat1 / 57.2958)
  * cos(lat2 / 57.2958) * cos(lon2 / 57.2958 - lon1 / 57.2958))
  )
  r <- 3963.0 # (statute miles) , default
  if ("units" == "nm") r <- 3437.74677 # (nautical miles)
  if ("units" == "km") r <- 6378.7 # (kilometers)
  return(r * temp)
}

library(maps) # must be installed from CRAN

# our noisy data (with no missing values)
turtles <- levels(loggerheadNoisy$turtle)
levels(loggerheadNoisy$turtle)

############## SPECIFY THE TURTLE NAME HERE
turtlename <- "BigMama"

dat <- loggerheadNoisy[which(loggerheadNoisy$turtle == turtlename), 5:6]
dat <- t(dat)

# Set up the MSSM model
# We are going to make a big approximation:
# We are going to pretend that 1 deg longitude
# change is equal to about the same distance (miles)
# over the range of latitudes that the turtles are moving
# That's not true.  There is about a 10% difference across
# their range of latitude movement.
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"
U.model <- "unequal"
Z.model <- "identity"

# Fit a random walk model for lon/lat to the lon/lat data
kem <- MARSS(dat, model = list(
  Z = Z.model,
  Q = Q.model, R = R.model, U = U.model
))
pred.lon <- kem$states[1, ]
pred.lat <- kem$states[2, ]

##########################################
# Plot the results
##########################################
op <- par(mai = c(0, 0, 0, 0))
map("state", region = c(
  "florida", "georgia",
  "south carolina", "north carolina", "virginia",
  "delaware", "new jersey", "maryland"
), xlim = c(-85, -65))
points(loggerheadNoisy$lon[which(loggerheadNoisy$turtle == turtlename)],
  loggerheadNoisy$lat[which(loggerheadNoisy$turtle == turtlename)],
  col = "blue", pch = 21, cex = 0.7
)
lines(pred.lon, pred.lat, col = "red", lwd = 2)

# add the good location data
goodturtles <- loggerhead
gooddat <- goodturtles[which(goodturtles$turtle == turtlename), 5:6]
points(gooddat[, 1], gooddat[, 2], col = "black", lwd = 2, pch = 3, cex = 1.1)
legend("bottomright", c("bad locations", "estimated true location", "good location data"), pch = c(1, -1, 3), lty = c(-1, 1, -1), col = c("blue", "red", "black"), bty = F)

# add the proposed fishery areas
lines(c(-77, -78, -78, -77, -77), c(33.5, 33.5, 32.5, 32.5, 33.5),
  col = "red", lwd = 2
)
lines(c(-75, -76, -76, -75, -75), c(38, 38, 37, 37, 38),
  col = "red", lwd = 2
)

###########################################
# Calculate the average miles traveled each day using
# the function GCDF defined above
# You must select and run the GCDF code first
###########################################
distance <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance[i - 1] <- GCDF(
    pred.lon[i - 1], pred.lon[i],
    pred.lat[i - 1], pred.lat[i]
  )
}

par(op) # reset plotting pars back to normal
par(mfrow = c(1, 2)) # make a 2 panel graph
hist(distance) # make a histogram of distance traveled per day
print(paste(
  "KalmanEM estimate of ave. mile per day for ",
  turtlename, " = ", mean(distance)
))

# Compare to the distance traveled per day if you used the raw data
distance <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance[i - 1] <- GCDF(dat[1, i - 1], dat[1, i], dat[2, i - 1], dat[2, i])
}
hist(distance) # make a histogram of distance traveled per day
mean(distance)
print(paste(
  "Raw estimate of ave. mile per day for ",
  turtlename, " = ", mean(distance)
))


