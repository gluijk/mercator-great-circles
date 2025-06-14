# Mercator maps and great circles
# www.overfitting.net
# https://www.overfitting.net/2025/06/el-astuto-mercator-mapas-y-grandes.html

library(ggmap)  # map_data() provides (long, lat) pairs forming all countries borders
library(data.table)
library(Cairo)  # output antialiasing


# GIS functions (first two coded by ChatGPT)

# Conversion from (long, lat) to Mercator (x, y) coordinates in km
longlat_to_mercatorkm <- function(df, R = 6371.23) {
    # R=6371.23 is Earth's average radius (km)
    
    # About the Mercator projection:
    # The Mercator coordinates in km are not true Euclidean distances.
    # While the units are in km, they do not correspond directly
    # to real-world distances between points, especially as you move away
    # from the equator.
    # 
    # The Mercator projection distorts scale with latitude:
    #   It keeps angles and shapes locally correct (conformal).
    #   But distance and area are distorted, especially near the poles.
    # 
    # The x-axis (longitude) is linear, but the y-axis (latitude) is nonlinear,
    # stretched using log(tan(π/4 + φ/2)).
    # 
    # So while:
    #   x_km ≈ real horizontal distances near the equator
    #   y_km does not correspond to vertical distance traveled in km 
    
    # Check that required columns exist
    if (!all(c("long", "lat") %in% names(df))) {
        stop("Input dataframe must contain 'long' and 'lat' columns")
    }
    
    # Convert degrees to radians
    lon_rad <- df$long * pi / 180
    lat_rad <- df$lat * pi / 180
    
    # Apply Mercator projection
    x_km <- R * lon_rad
    y_km <- R * log(tan(pi / 4 + lat_rad / 2))
    
    # Return dataframe with original and projected coordinates
    return( cbind(df, x_km = x_km, y_km = y_km) )
}

# Equispaced interpolation of (long, lat) points along a great circle
# using spherical linear interpolation (slerp)
great_circle_points <- function(long1, lat1, long2, lat2, n = 100) {
    # Convert degrees to radians
    deg2rad <- function(deg) deg * pi / 180
    rad2deg <- function(rad) rad * 180 / pi
    
    # Convert lon/lat to radians
    lon1 <- deg2rad(long1)
    lat1 <- deg2rad(lat1)
    lon2 <- deg2rad(long2)
    lat2 <- deg2rad(lat2)
    
    # Convert to 3D Cartesian coordinates
    xyz1 <- c(cos(lat1) * cos(lon1), cos(lat1) * sin(lon1), sin(lat1))
    xyz2 <- c(cos(lat2) * cos(lon2), cos(lat2) * sin(lon2), sin(lat2))
    
    # Compute angle between the two vectors
    omega <- acos(sum(xyz1 * xyz2))
    
    # Interpolate using spherical linear interpolation (slerp)
    if (omega == 0) {
        # Points are the same; return a repeated point
        lon <- rep(long1, n)
        lat <- rep(lat1, n)
    } else {
        f <- seq(0, 1, length.out = n)
        sin_omega <- sin(omega)
        interpolated <- t(sapply(f, function(t) {
            a <- sin((1 - t) * omega) / sin_omega
            b <- sin(t * omega) / sin_omega
            xyz <- a * xyz1 + b * xyz2
            # Convert back to lon/lat
            lon <- atan2(xyz[2], xyz[1])
            lat <- asin(xyz[3] / sqrt(sum(xyz^2)))
            c(rad2deg(lon), rad2deg(lat))
        }))
        lon <- interpolated[,1]
        lat <- interpolated[,2]
    }
    
    return(data.frame(long = lon, lat = lat))
}

# Great circle distance in km between two given (long, lat) points
great_circle_distance=function(long1, lat1, long2, lat2, R=6371.23) {
    # R=6371.23 is Earth's average radius (km)
    
    # General great circle distance formula (spherical law of cosines)
    # calculates the shortest path between two points
    # on the surface of a sphere:
    # d = R * acos( sin(theta1) * sin(theta2) +
    #     cos(theta1) * cos(theta2) * cos(phi2 - phi1) )
    
    # Convert degrees to radians
    deg2rad=function(deg) deg * pi / 180
    
    d=R*acos(
        sin(deg2rad(lat1))*sin(deg2rad(lat2)) +
            cos(deg2rad(lat1))*cos(deg2rad(lat2))*cos(deg2rad(long2)-deg2rad(long1))
    )
    return(d)
}

# Function that provides the longitude span at a given latitude
# that covers a distance d in km along a great circle
longspan_latdistance=function(lat, d=1000, R=6371.23) {
    # R=6371.23 is Earth's average radius (km)
    
    # General great circle distance formula (spherical law of cosines)
    # calculates the shortest path between two points
    # on the surface of a sphere:
    # d = R * acos( sin(theta1) * sin(theta2) +
    #     cos(theta1) * cos(theta2) * cos(phi2 - phi1) )
    
    # We are solving for longspan_latdistance = phi2 - phi1
    # In our case (constant latitude): theta1 = theta2 = lat
    
    # Convert degrees to radians
    deg2rad=function(deg) deg * pi / 180
    rad2deg=function(rad) rad * 180 / pi

    delta=rad2deg(
        acos( (cos(d/R)-sin(deg2rad(lat))^2) / cos(deg2rad(lat))^2 )
    )
    
    return(delta)
}

# Initial compass bearing (also called azimuth) between two geographic
# coordinates (long, lat) when following the great circle
great_circle_initial_bearing <- function(long1, lat1, long2, lat2) {
    # Convert degrees to radians
    deg2rad <- function(deg) deg * pi / 180
    rad2deg <- function(rad) rad * 180 / pi
    
    long1 <- deg2rad(long1)
    lat1 <- deg2rad(lat1)
    long2 <- deg2rad(long2)
    lat2 <- deg2rad(lat2)
    
    delta_long <- long2 - long1
    
    x <- sin(delta_long) * cos(lat2)
    y <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(delta_long)
    
    initial_bearing <- atan2(x, y)
    compass_bearing <- (rad2deg(initial_bearing) + 360) %% 360
    
    return(compass_bearing)
}

# Loxodromic (rhumb line) bearing, the constant compass heading that maintains
# a fixed angle with meridians (i.e., straight lines on a Mercator projection)
loxodromic_rhumb_bearing <- function(long1, lat1, long2, lat2) {
    # Convert degrees to radians
    deg2rad <- function(deg) deg * pi / 180
    rad2deg <- function(rad) rad * 180 / pi
    
    long1 <- deg2rad(long1)
    lat1 <- deg2rad(lat1)
    long2 <- deg2rad(long2)
    lat2 <- deg2rad(lat2)
    
    delta_long <- long2 - long1
    
    # Handle crossing the antimeridian
    if (abs(delta_long) > pi) {
        delta_long <- ifelse(delta_long > 0, delta_long - 2 * pi, delta_long + 2 * pi)
    }
    
    delta_psi <- log(tan(pi/4 + lat2/2) / tan(pi/4 + lat1/2))
    
    bearing_rad <- atan2(delta_long, delta_psi)
    bearing_deg <- (rad2deg(bearing_rad) + 360) %% 360
    
    return(bearing_deg)
}


################################################################################

# IMAGE DIMENSIONS
DIMX=1920  # Full HD
DIMY=1080  # 1920 x 1080 pixels
NPOINTS=2000  # number of points of each great circle (flight route)
GAP=5

# READ WORLD COORDINATES
DT=data.table(map_data("world"))  # (long, lat) pairs for all countries
DT=unique(DT[, .(long, lat)])  # summarize to deduplicate points
DT$long[DT$long>180]=DT$long[DT$long>180]-360  # offset out of range points
DT=longlat_to_mercatorkm(DT)  # add Mercator (x,y) columns


############################
# 1. DRAW SOME FLIGHTS ALONG GREAT CIRCLES

# All flight destinations in (long, lat) provided by ChatGPT
cities=data.frame(
    city = c("Madrid", "New York", "Los Angeles", "Mexico City", "Bogota",
             "Buenos Aires", "Moscow", "Tokyo", "Johannesburg", "Sydney"),
    long = c(-3.7038, -74.0060, -118.2437, -99.1332, -74.0721,
             -58.3816, 37.6173, 139.6917, 28.0473, 151.2093),
    lat  = c(40.4168, 40.7128, 34.0522, 19.4326, 4.7110,
             -34.6037, 55.7558, 35.6895, -26.2041, -33.8688)
)
cities=longlat_to_mercatorkm(cities)
NCITIES=nrow(cities)

# Calculate and display great circle distances from Madrid
cities$km2madrid=great_circle_distance(cities$long[cities$city=="Madrid"],
                                       cities$lat[cities$city=="Madrid"],
                                       cities$long,
                                       cities$lat)
cities=cities[order(cities$km2madrid), ]  # order by distance asc
for (i in 2:NCITIES) print(paste0("Flight distance ",
                cities$city[cities$city=="Madrid"], "-",
                cities$city[i], ": ", round(cities$km2madrid[i], 1), " km"))
# Bar plot with distances from Madrid
MAXIMO=max(cities$km2madrid)
rango=2:NCITIES
barheights=barplot(cities$km2madrid[rango], names.arg=cities$city[rango],
        main="Flight distance in km from Madrid to...",
        ylim=c(0, MAXIMO*1.2), axes=FALSE, cex.names=0.6)
text(barheights, cities$km2madrid[rango]+MAXIMO/20,
     labels=round(cities$km2madrid[rango],1),
     col="red", cex=0.7)


# Plot Maps:
# Longitude/Latitude map
CairoPNG("Map_LongLat_FLIGHTS.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$long, DT$lat, main='Longitude/Latitude map', cex.main=4,
        xlim=c(-170,170), ylim=c(-75, 75),  # crop Antarctica a bit
        pch=20, cex=0.02, asp=1,
        xlab="Longitude (°)", ylab="Latitude (°)")
    abline(h=0, v=0, lty='dotted')  # draw equator and Greenwich meridian
    
    # Draw flight routes from Madrid
    for (i in 2:NCITIES) {
        flight=great_circle_points(cities$long[1], cities$lat[1],
                                   cities$long[i], cities$lat[i],
                                   n=NPOINTS)    
        points(flight$long, flight$lat, pch=20, cex=0.02, col='red')
    }
    # NY-Tokyo
    flight=great_circle_points(cities$long[cities$city=="New York"],
                               cities$lat[cities$city=="New York"],
                               cities$long[cities$city=="Tokyo"],
                               cities$lat[cities$city=="Tokyo"],
                               n=NPOINTS)
    points(flight$long, flight$lat, pch=20, cex=0.02, col='blue')
    # Cities names
    text(cities$long, cities$lat, cities$city, col='darkgreen', cex=3)
dev.off()

# Mercator map
CairoPNG("Map_Mercator_FLIGHTS.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$x_km, DT$y_km, main='Mercator map', cex.main=4,
        ylim=c(-12000, 16000),  # crop Antarctica a bit
        pch=20, cex=0.02, asp=1,
        xlab="Mercator km", ylab="Mercator km")
    abline(h=0, v=0, lty='dotted')  # draw equator and Greenwich meridian
    
    # Draw flight routes from Madrid
    for (i in 2:NCITIES) {
        flight=great_circle_points(cities$long[1], cities$lat[1],
                                   cities$long[i], cities$lat[i],
                                   n=NPOINTS)
        flightmerc=longlat_to_mercatorkm(flight)
        points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='red')
    }
    # NY-Tokyo
    flight=great_circle_points(cities$long[cities$city=="New York"],
                               cities$lat[cities$city=="New York"],
                               cities$long[cities$city=="Tokyo"],
                               cities$lat[cities$city=="Tokyo"],
                               n=NPOINTS)
    flightmerc=longlat_to_mercatorkm(flight)
    points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='blue')
    # Cities names
    text(cities$x_km, cities$y_km, cities$city, col='darkgreen', cex=3)
dev.off()


############################
# 2. DRAW SERIES OF 1000km GREAT CIRCLES AT DIFFERENT LATITUDES

NLATSL=1001  # number of latitudes to compute for side limits
NLATSC=21  # number of great circles to draw (odd number to include equator)
OFFSETLO=10  # longitude offset in degrees (as in original map)
dfos=data.frame(long=OFFSETLO, lat=0)  # lat value is irrelevant here
dfos=longlat_to_mercatorkm(dfos)
OFFSETKM=dfos$x_km  # offset conversion from long degrees to Mercator x_km
rm(dfos)

# Calculate 1000km side limits (m1) (only positive lat around Greenwich meridian)
m1=as.data.frame(matrix(nrow=NLATSL, ncol=2))
colnames(m1)=c('long', 'lat')
i=1
for (lat in seq(from=-85, to=85, length.out=NLATSL)) {
    m1$long[i]=longspan_latdistance(lat)/2  # /2 to replicate it L/R
    m1$lat[i]=lat
    i=i+1
}
m1=longlat_to_mercatorkm(m1)

# Calculate 1000km great circles (m2)
m2=as.data.frame(matrix(nrow=NLATSC, ncol=2))
colnames(m2)=c('long', 'lat')
i=1
for (lat in seq(from=-85, to=85, length.out=NLATSC)) {
    m2$long[i]=longspan_latdistance(lat)/2  # /2 to replicate it L/R
    m2$lat[i]=lat
    i=i+1
}
m2=longlat_to_mercatorkm(m2)
# Check all distances are 1000km
print(great_circle_distance(m2$long, m2$lat, -m2$long, m2$lat))


# Plot Maps:
# Longitude/Latitude map
CairoPNG("Map_LongLat_GC1000km.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$long, DT$lat, main='Longitude/Latitude map', cex.main=4,
         xlim=c(-170,170), ylim=c(-75, 75),  # crop Antarctica a bit
         pch=20, cex=0.02, asp=1,
         xlab="Longitude (°)", ylab="Latitude (°)")
    abline(h=0, v=0, lty='dotted')  # draw equator and Greenwich meridian
    
    # Draw side limits
    lines( m1$long+OFFSETLO, m1$lat, type='l', col='red', lty='longdash')
    lines(-m1$long+OFFSETLO, m1$lat, type='l', col='red', lty='longdash')
    
    # Draw great circles
    for (i in 2:(nrow(m2)-1)) {
        flight=great_circle_points(m2$long[i]+OFFSETLO, m2$lat[i],
                                  -m2$long[i]+OFFSETLO, m2$lat[i],
                                   n=NPOINTS)    
        points(flight$long, flight$lat, pch=20, cex=0.02, col='red')
    }
dev.off()

# Mercator map
CairoPNG("Map_Mercator_GC1000km.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$x_km, DT$y_km, main='Mercator map', cex.main=4,
         ylim=c(-13000, 16000),  # crop Antarctica a bit
         pch=20, cex=0.02, asp=1,
         xlab="Mercator km", ylab="Mercator km")
    abline(h=0, v=0, lty='dotted')  # draw equator and Greenwich meridian
    
    # Draw side limits
    lines( m1$x_km+OFFSETKM, m1$y_km, type='l', col='red', lty='longdash')
    lines(-m1$x_km+OFFSETKM, m1$y_km, type='l', col='red', lty='longdash')

    # Draw great circles    
    for (i in 2:(nrow(m2)-1)) {
        flight=great_circle_points(m2$long[i]+OFFSETLO, m2$lat[i],
                                  -m2$long[i]+OFFSETLO, m2$lat[i],
                                   n=NPOINTS)
        flightmerc=longlat_to_mercatorkm(flight)
        points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='red')
    }
dev.off()

