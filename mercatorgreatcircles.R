# Mercator maps and great circles
# www.overfitting.net
# https://www.overfitting.net/

library(ggmap)  # map_data() provides (x,y) pairs forming all countries borders
library(data.table)
library(Cairo)  # output antialiasing


# GIS functions (first two coded by ChatGPT)

# Conversion from (long, lat) to Mercator (x,y) coordinates in km
lonlat_df_to_mercator_km <- function(df, R=6371.23) {
    # R=6371.23 is Earth's average radius (km)
    
    # NOTE:
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
    cbind(df, x_km = x_km, y_km = y_km)
}

# Equispaced interpolation of (long,lat) points along a great circle
# using spherical linear interpolation (slerp)
great_circle_points_manual <- function(long1, lat1, long2, lat2, n = 100) {
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

deltalong=function(lat, d=1000, R=6371.23) {
    # R=6371.23 is Earth's average radius (km)
    # Function that provides the longitude span at a given latitude
    # that covers a distance d along a great circle
    
    # General great-circle distance formula (spherical law of cosines)
    # calculates the shortest path between two points
    # on the surface of a sphere:
    # d = R * acos( sin(theta1) * sin(theta2) +
    #     cos(theta1) * cos(theta2) * cos(phi2 - phi1) )
    
    # We are looking for deltalong = phi2 - phi1
    # theta1 = theta2 = lat
    
    # Convert degrees to radians
    deg2rad <- function(deg) deg * pi / 180
    rad2deg <- function(rad) rad * 180 / pi

    return(rad2deg(
        acos( (cos(d/R)-sin(deg2rad(lat))^2) / cos(deg2rad(lat))^2 )
        )
    )
}


################################################################################

# IMAGE DIMENSIONS
DIMX=1920  # Full HD animation
DIMY=1080  # 1920 x 1080 pixels
NPOINTS=2000  # number of points of each flight route
GAP=5

# READ WORLD COORDINATES
DT=data.table(map_data("world"))  # long/lat pairs for all countries
DT=unique(DT[, .(long, lat)])  # summarize to deduplicate points
DT$long[DT$long>180]=DT$long[DT$long>180]-360  # offset out of range points
DT=lonlat_df_to_mercator_km(DT)  # add Mercator (x,y) columns


############################
# 1. DRAW SOME FLIGHTS ALONG GREAT CIRCLES

# All flight destinations in  (long,lat)
cities=data.frame(
    city = c("Madrid", "New York", "Los Angeles", "Mexico City", "Bogotá",
             "Buenos Aires", "Moscow", "Tokyo", "Johannesburg", "Sydney"),
    long = c(-3.7038, -74.0060, -118.2437, -99.1332, -74.0721,
             -58.3816, 37.6173, 139.6917, 28.0473, 151.2093),
    lat  = c(40.4168, 40.7128, 34.0522, 19.4326, 4.7110,
             -34.6037, 55.7558, 35.6895, -26.2041, -33.8688)
)


# Plot Maps
CairoPNG("Map_LongLat_FLIGHTS.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$long, DT$lat, main='Longitude/Latitude map', cex.main=4,
        xlim=c(-170,170), ylim=c(-75, 75),  # crop Antarctica a bit
        pch=20, cex=0.02, asp=1,
        xlab="Longitude (°)", ylab="Latitude (°)")
    
    # Draw flight routes
    for (i in 2:nrow(cities)) {
        flight=great_circle_points_manual(cities$long[1], cities$lat[1],
                                          cities$long[i], cities$lat[i],
                                          n=NPOINTS)    
        points(flight$long, flight$lat, pch=20, cex=0.02, col='red')
    }
    # NY-Tokyo
    flight=great_circle_points_manual(cities$long[2], cities$lat[2],
                                      cities$long[8], cities$lat[8],
                                      n=NPOINTS)
    points(flight$long, flight$lat, pch=20, cex=0.02, col='blue')
    
    abline(h=0, v=0, lty='dotted')
dev.off()

CairoPNG("Map_Mercator_FLIGHTS.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$x_km, DT$y_km, main='Mercator map', cex.main=4,
        ylim=c(-12000, 16000),  # crop Antarctica a bit
        pch=20, cex=0.02, asp=1,
        xlab="Mercator km", ylab="Mercator km")
    
    # Draw flight routes
    for (i in 2:nrow(cities)) {
        flight=great_circle_points_manual(cities$long[1], cities$lat[1],
                                          cities$long[i], cities$lat[i],
                                          n=NPOINTS)
        flightmerc=lonlat_df_to_mercator_km(flight)
        points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='red')
    }
    # NY-Tokyo
    flight=great_circle_points_manual(cities$long[2], cities$lat[2],
                                      cities$long[8], cities$lat[8],
                                      n=NPOINTS)
    flightmerc=lonlat_df_to_mercator_km(flight)
    points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='blue')
    
    abline(h=0, v=0, lty='dotted')
dev.off()


############################
# 2. DRAW SERIES OF 1000km GREAT CIRCLES AT DIFFERENT LATITUDES

NLATSL=1001  # number of latitudes to compute for limits
NLATSC=21  # number of great circles to draw (odd number to include equator)
OFFSETLO=10  # longitude offset in degrees
dfos=data.frame(long=OFFSETLO, lat=0)  # lat value is irrelevant here
dfos=lonlat_df_to_mercator_km(dfos)
OFFSETKM=dfos$x_km  # offset conversion from long degrees to Mercator x_km
rm(dfos)

# Calculate 1000km limits (only positive side around Greenwich meridian)
m1=as.data.frame(matrix(nrow=NLATSL, ncol=2))
colnames(m1)=c('long', 'lat')
i=1
for (lat in seq(from=-85, to=85, length.out=NLATSL)) {
    m1$long[i]=deltalong(lat)/2  # divided by 2 because we'll replicate it L/R
    m1$lat[i]=lat
    i=i+1
}
m1=lonlat_df_to_mercator_km(m1)

# Calculate 1000km great circles
m2=as.data.frame(matrix(nrow=NLATSC, ncol=2))
colnames(m2)=c('long', 'lat')
i=1
for (lat in seq(from=-85, to=85, length.out=NLATSC)) {
    m2$long[i]=deltalong(lat)/2  # divided by 2 because we'll replicate it L/R
    m2$lat[i]=lat
    i=i+1
}
m2=lonlat_df_to_mercator_km(m2)


# Plot Maps
CairoPNG("Map_LongLat_GC1000km.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$long, DT$lat, main='Longitude/Latitude map', cex.main=4,
         xlim=c(-170,170), ylim=c(-75, 75),  # crop Antarctica a bit
         pch=20, cex=0.02, asp=1,
         xlab="Longitude (°)", ylab="Latitude (°)")
    
    # Draw side limits
    lines( m1$long+OFFSETLO, m1$lat, type='l', col='red', lty='longdash')
    lines(-m1$long+OFFSETLO, m1$lat, type='l', col='red', lty='longdash')
    
    # Draw great circles
    for (i in 2:(nrow(m2)-1)) {
        flight=great_circle_points_manual(m2$long[i]+OFFSETLO, m2$lat[i],
                                         -m2$long[i]+OFFSETLO, m2$lat[i],
                                          n=NPOINTS)    
        points(flight$long, flight$lat, pch=20, cex=0.02, col='red')
    }
    
    abline(h=0, v=0, lty='dotted')
dev.off()

CairoPNG("Map_Mercator_GC1000km.png", width=DIMX*2, height=DIMY*2, antialias="subpixel")
    # Draw map
    plot(DT$x_km, DT$y_km, main='Mercator map', cex.main=4,
         ylim=c(-13000, 16000),  # crop Antarctica a bit
         pch=20, cex=0.02, asp=1,
         xlab="Mercator km", ylab="Mercator km")
    
    # Draw side limits
    lines( m1$x_km+OFFSETKM, m1$y_km, type='l', col='red', lty='longdash')
    lines(-m1$x_km+OFFSETKM, m1$y_km, type='l', col='red', lty='longdash')

    # Draw great circles    
    for (i in 2:(nrow(m2)-1)) {
        flight=great_circle_points_manual(m2$long[i]+OFFSETLO, m2$lat[i],
                                         -m2$long[i]+OFFSETLO, m2$lat[i],
                                          n=NPOINTS)
        flightmerc=lonlat_df_to_mercator_km(flight)
        points(flightmerc$x_km, flightmerc$y_km, pch=20, cex=0.02, col='red')
    }
    
    abline(h=0, v=0, lty='dotted')
dev.off()

