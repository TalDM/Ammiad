#---
#title: "calculating plant density from survey done in 2020 on transect A"
#author: "Tal Dahan-Meir"
#date: "25/08/2021"
#---

require(geosphere)

csv = read.csv("data/2020_plant_densities.csv",colClasses="character")

#get dists between plants to next plants and add to data.frame
nPlants_a = nrow(csv)
coordinates_next_plants = data.frame("altitude_nextPlant" = rep(NA,nPlants_a),
                                     "latitude_nextPlant" = rep(NA,nPlants_a),
                                     "longitude_nextPlant" = rep(NA,nPlants_a))
for (i in 1:(nPlants_a-1)){
  coordinates_next_plants[i,] = csv[i+1,c("altitude","latitude","longitude")]
}
csv_distances = cbind (csv,coordinates_next_plants)
csv_distances$dist_meters = c( apply(csv_distances[1 : (nPlants_a-1), ], 1, function(tt) { distm(x = c(as.numeric(tt[6]), as.numeric(tt[5]) ), y = c(as.numeric(tt[9]), as.numeric(tt[8]) ) )} ), "NA")
csv_distances$position_meters = as.numeric(c( 0 ,cumsum (csv_distances[1:(nPlants_a-1), "dist_meters"] )))
csv_distances$plantsPerMeter = as.numeric(csv_distances$plants_until_next_point) / as.numeric(csv_distances$dist_meters)
csv_distances$plantsinMM_divided_by_pi = as.numeric(csv_distances$plants_in_1MM2) / pi
csv_distances$plantsinM_combined = as.numeric(csv_distances$plantsinMM_divided_by_pi)+as.numeric(csv_distances$plantsPerMeter)

head(csv_distances)


pdf("transect_a_density_2020.pdf", height=8, width=22)
plot (as.numeric(csv_distances$position_meters),
      as.numeric(csv_distances$plantsinMM_divided_by_pi),
      main = "Density transect A",
      cex.main=2,
      pch = 19,
      cex = 0.2,
      cex.lab = 1.5,
      cex.axis=1.4,
      ylim=c(0,25),
      xlim=c(0,320),
      col="black",
      xlab = "position (M)",
      ylab = "spikes per M^2")
#apply(csv_distances,1, function(xx){
 #     segments(x0 = as.numeric(xx[11]),
  #             y0 = as.numeric(xx[12]),
   #            x1 = as.numeric(xx[11])+as.numeric(xx[10]),
    #           y1 = as.numeric(xx[12]),
     #          col="grey") } )
apply(csv_distances,1, function(xx){
  rect(xleft = as.numeric(xx[11]),
      ybottom = 0,
      xright = as.numeric(xx[11])+as.numeric(xx[10]),
      ytop = as.numeric(xx[12])/2,
      col="#00ABA920",
      border= NA)} )
names_height_bonus = 0.5
apply(csv_distances,1, function(xx){
  text(x = as.numeric(xx[11]),
       y = as.numeric(xx[13])+names_height_bonus,
       as.character(xx[1]),
       col="black",
       cex=1)} )
dev.off()

###calculate mean density###
(sum(as.numeric(csv_distances$plantsinMM_divided_by_pi))+(sum(as.numeric(csv_distances$plants_until_next_point),na.rm=T)/2))/max(csv_distances$position_meters)
