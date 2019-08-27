

### convert GCM data to means and SDs ###

library(ecoclim)
library(tidyverse)
library(raster)


select <- dplyr::select



###### restoration parameters #####

# list of focal species
spp <- read.csv("focal_species.csv", stringsAsFactors = F)
spp <- unique(spp$Species)




### climate ensemble summaries ####

clim <- data.frame(path=list.files("F:/seeds/chelsa_derived", recursive=T, full.names=T),
                   stringsAsFactors=F) %>%
      mutate(set=sub("\\.tif", "", basename(path))) %>%
      separate(set, c("year", "model", "scenario"), remove=F, sep=" ")

ext <- extent(-125, -114, 32.5, 42)

# for each variable-year, mean and sd of model ensemble (summarize over models and scenarios)
vars <- c("PPT", "PET", "AET", "CWD", "RAR", "DJF", "JJA")
for(y in unique(clim$year)){
      for(v in 1:7){
            s <- stackBands(clim$path[clim$year==y], v)
            var <- vars[v]
            if(var=="PPT") s <- log10(s)

            avg <- mean(s)
            stdev <- calc(s, sd)
            s <- stack(avg, stdev)

            s <- crop(s, ext)

            writeRaster(s,
                        paste0("assets/seedsource_data/climate/ensemble ", y, " ", var, ".tif"),
                        overwrite=T)
      }
}






### species ranges #########

f <- list.files("E:/phycon/001_distributions/pred_s25km", full.names=T)
#f <- f[sub("\\.rds", "", basename(f)) %in% spp]

## transfer sdm to climate grid ##
template <- raster(clim$path[1])
for(x in spp){
      message(x)
      r <- readRDS(f[grepl(x, f)])
      m <- cellStats(r, max) * .2
      r <- r %>%
            reclassify(c(-1, m, NA,
                         m, 1, 1))
      d <- r %>%
            rasterToPoints() %>%
            as.data.frame()

      coordinates(d) <- c("x", "y")
      
      crs(d) <- crs(r)
      d <- spTransform(d, crs(template))

      r <- rasterize(d, template, field="layer") %>%
            crop(ext)
      saveRDS(r, paste0("assets/seedsource_data/ranges/", x, ".rds"))
}




#### soils ####

# soils data
soil <- list.files("F:/SoilGrids", full.names=T, pattern=".tif") %>%
      stack() %>%
      crop(ext)

# PCA to reduce dimensionality
sv <- values(soil)
na <- apply(sv, 1, function(x) any(is.na(x)))
sva <- sv[!na,]
pca <- prcomp(sva, center=T, scale=T)
npcs <- 5
svpc <- sv[,1:npcs]
svpc[!na,] <- pca$x[,1:npcs]
soil_pc <- soil[[1:npcs]]
soil_pc[] <- svpc

# aggregate to climate grid
template <- crop(template, ext)
template[] <- 1:ncell(template)

agg <- as.data.frame(svpc) %>%
      cbind(id=extract(template, coordinates(soil_pc)))
agg <- agg %>%
      group_by(id) %>%
      summarize_all(mean, na.rm=T)
agg <- select(agg, -id) %>% as.matrix()

soil <- stack(template, template, template,
              template, template)
for(i in 1:nlayers(soil)) soil[[i]][] <- agg[,i]

writeRaster(soil, "assets/seedsource_data/climate/soil_800m.tif",
            overwrite=T)




##### smoothed environmental rasters per species #####

vars <- sort(c("PPT", "AET", "CWD", "DJF", "JJA"))
clim_files <- data.frame(path=list.files("assets/seedsource_data/climate",
                                         full.names=T, pattern="ensemble"),
                         stringsAsFactors=F) %>%
      mutate(set=gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "var"), remove=F, sep=" ") %>%
      select(-junk) %>%
      filter(var %in% vars) %>% arrange(var)

soil_all <- stack("assets/seedsource_data/climate/soil_800m.tif")


# historic environments across species range
species_envt <- function(x){
      range <- readRDS(paste0("assets/seedsource_data/ranges/", x, ".rds"))
      
      clim <- stackBands(clim_files$path[clim_files$year=="1979-2013"], 1)
      names(clim) <- vars
      
      range <- crop(range, clim)
      clim <- clim %>% crop(range) %>% mask(range) %>% trim()
      soil <- soil_all %>% mask(range) %>% trim()
      
      list(soil = crop(soil, clim),
           clim = crop(clim, soil))
}

# smoothed historic environment representing gene flow
smoothed_envt <- function(e, radius){
      w <- matrix(1, 1+radius*2, 1+radius*2)
      center <- length(w) / 2 + .5
      
      m <- function(x, ...){
            if(is.na(x[center])) return(NA)
            mean(na.omit(x))
      }
      
      clim <- e$clim
      soil <- e$soil
      
      if(radius>0){
            for(i in 1:nlayers(clim)){
                  clim[[i]] <- focal(clim[[i]], w=w, fun=m, pad=T) %>%
                        mask(clim[[1]])
            }
            for(i in 1:nlayers(soil)){
                  soil[[i]] <- focal(soil[[i]], w=w, fun=m, pad=T) %>%
                        mask(soil[[1]])
            }
      }
      
      names(clim) <- vars
      names(soil) <- paste0("soil", 1:nlayers(soil))
      
      list(soil = soil,
           clim = clim)
}



for(x in spp){
      message(x)
      e <- species_envt(x)
      
      for(r in rev(c(0, 1, 2, 5, 10, 20, 50, 100))){
            message(paste("       r =", r))
            out <- paste0("assets/seedsource_data/smooth/", x, " ", r, ".rds")
            if(file.exists(out)) next()
            s <-  smoothed_envt(e, r)
            saveRDS(s, out)
      }
}
