# This is a file to calulate sighting distance in best-case scenario for Sergestid shrimps. 

install.packages("data.table")
install.packages("gsl")

library(data.table)
library(gsl)


#Parameters: 
#   c = beam attenuation in 1/m 
#   dt = integration time in seconds 
#   A = eye size in meters
#   d = diameter of photoreceptors in um
#   E = emission intensity of bioluminescent source in photons/sec
#   N0 = number of photons absorbed by 1m pupil (from background light in horizontal direction)
#   k = diffuse attenuation (not needed when horizontal)
#   N0d = photons/s/m^2 in downward direction (used to calculate light blocked by cephalothorax)

#Create function to calculate sighting distance in meters

s.dist <- function(E,A,c,d,dt,N0){
  if(is.na(A)){return(NA)}else{
    b <- (E*dt)/(3*(1 + sqrt(1 + (2.8*(d^2)*N0)*dt)))
    r <- ((2/c)*lambert_W0((c*A/8)*sqrt(b)))
    return(r)}
}

#Set global parameter values: 

c <- 0.0468                         #Taken from Ruxton and Johnsen (2016)
dtlow  <- 0.0417                    #Taken from Frank (2000)
dthigh <- 0.0588                    #Taken from Frank (2000)
dt <- (dtlow+dthigh)/2              # best estimation of actual integration time
N0 <- 0                             #hypothetical situation of complete background darkness
d = 3*(10^-6)                       #Photoreceptor diameter was set to 3 μm (Land & Nilsson, 2012); everything expressed as meters

### Calculate Distances: 
## Assuming that actual aperture is equal to eye diameter and that all of the given environmental and bioluminescent conditions are the same. 

# Lindsay et al. 1999 estimate max bioluminescence emission from S. similis as 19 x 10^12/m2/s1. Needs to be photons/sec for the model below. To makes this conversion, I must estimate the surface area of the light organs and scale the intensity for the reduction in surface area. 
# From Foxton et al, estimated surface area of S. similis organs of pesta as 10 sq mm of light organ surface area. 
#1000mm in a meter
#so 100x reduction in brightness needed

#Brightness 19 x 10^12/m2/s1 divided by 100 = 19 x 10^10/mm2/s1 = 19 x 10^10 (190000000000) or roughly 10^10 photos/sec which is in line with biolum recordings in the deep sea (see) EE for ref.

#Best case scenario (full eye diameter; bright biolum; complete darkness in background)


## D. corniculum
# Full eye diameter as aperture
corniculum  <- s.dist(10^10, 0.0011, c, d, dt, 0)
print(corniculum)
# 2.380319 m

## D. henseni
# Full eye diameter as aperture
henseni  <- s.dist(10^10, 0.0010, c, d, dt, 0)
print(henseni)
# 2.174379 m
                           
## R. robusta
# Full eye diameter as aperture
robusta <- s.dist(10^10, 0.0018, c, d, dt, 0)
print(robusta)
# 3.770407 m


## C. hansjacobi
hansjacobi  <- s.dist(10^10, 0.0012, c, d, dt, 0)
print(hansjacobi)
# 2.584344 m


## S. atlanticus
atlanticus  <- s.dist(10^10, 0.0008, c, d, dt, 0)
print(atlanticus)
# 1.756593 m


## P. grandis
grandis  <- s.dist(10^10, 0.0015, c, d, dt, 0)
print(grandis)
# 3.185319 m


## R. regalis
regalis  <- s.dist(10^10, 0.0014, c, d, dt, 0)
print(regalis)
# 2.986806 m


## S. tenuiremis 
tenuiremis  <- s.dist(10^10, 0.0013, c, d, dt, 0)
print(tenuiremis)
# 2.786494 m


## P. armatus
armatus  <- s.dist(10^10, 0.0007, c, d, dt, 0)
print(armatus)
# 1.54466 m


## G. splendens
splendens  <- s.dist(10^10, 0.0012, c, d, dt, 0)
print(splendens)
# 2.584344 m


## G. edwardsii
edwardsii  <- s.dist(10^10, 0.0006, c, d, dt, 0)
print(edwardsii)
# 1.330641 m


## P. vigilax
vigilax  <- s.dist(10^10, 0.0006, c, d, dt, 0)
print(vigilax)
# 1.330641 m


## A. sargassi 
sargassi  <- s.dist(10^10, 0.0007, c, d, dt, 0)
print(sargassi)
# 1.54466 m


##A. pectinatus
pectinatus  <- s.dist(10^10, 0.0005, c, d, dt, 0)
print(pectinatus)
# 1.114491 m


## C. talismani 
talismani  <- s.dist(10^10, 0.0011, c, d, dt, 0)
print(talismani)
# 2.380319 m


## E. arcticus
arcticus  <- s.dist(10^10, 0.0014, c, d, dt, 0)
print(arcticus)
# 2.986806 m



# Next, we'll fit a regression to sighting distance in body lengths as a function of body lengths (to see if relative sighting distance changes with differences in species seizes) 
#Read in the data - sighting distance in BL ~ BL
data<-read.csv(file.choose(), header=TRUE)
head(data)

test<-lm(data$SD_BL ~ data$BL)
summary(test)


