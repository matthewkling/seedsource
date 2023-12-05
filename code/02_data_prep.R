

### convert GCM data to means and SDs ###

library(ecoclim)
library(tidyverse)
library(raster)


select <- dplyr::select



###### restoration parameters #####

# list of focal species
spp <- read.csv("assets/focal_species.csv", stringsAsFactors = F)
key_spp <- spp$Species


### climate ensemble summaries ####

clim <- data.frame(path=list.files("~/data/seeds/chelsa_derived", recursive=T, full.names=T),
                   stringsAsFactors=F) %>%
      mutate(set=sub("\\.tif", "", basename(path))) %>%
      separate(set, c("year", "model", "scenario"), remove=F, sep=" ")

ext <- extent(-125, -114, 32.5, 42)

# for each variable-year-scenario, mean and sd of model ensemble
vars <- c("PPT", "PET", "AET", "CWD", "RAR", "DJF", "JJA")
for(y in unique(clim$year)){
      for(v in 1:7){
            for(p in unique(clim$scenario)){
                  
                  paths <- clim$path[clim$year == y & clim$scenario == p]
                  if(length(paths) == 0) next()
                  s <- stackBands(paths, v)
                  var <- vars[v]
                  if(var=="PPT") s <- log10(s)
                  
                  avg <- mean(s)
                  stdev <- calc(s, sd)
                  s <- stack(avg, stdev)
                  
                  s <- crop(s, ext)
                  
                  terra::writeRaster(terra::rast(s),
                                     paste0("assets/climate/ensemble ", y, " ", p, " ", var, ".tif"),
                                     overwrite=T)
            }
      }
}

# collate individual models, plus mean (alternative to above)
msk <- rast("assets/all_species/smooth/NONE 0.tif")[[1]]
vars <- c("PPT", "PET", "AET", "CWD", "RAR", "DJF", "JJA")
for(y in unique(clim$year)){
      for(v in 1:7){
            for(p in unique(clim$scenario)){
                  
                  paths <- clim$path[clim$year == y & clim$scenario == p]
                  if(length(paths) == 0) next()
                  s <- stackBands(paths, v)
                  var <- vars[v]
                  s <- s %>% crop(ext)
                  
                  if(var == "PPT") s <- log10(s)
                  
                  avg <- mean(s)
                  # stdev <- calc(s, sd)
                  s <- stack(avg, s)
                  
                  terra::writeRaster(terra::rast(s) %>% terra::crop(msk) %>% terra::mask(msk),
                                     paste0("assets/models/ensemble ", y, " ", p, " ", var, ".tif"),
                                     overwrite = T)
            }
      }
}







#### soils ####

# soils data
soil <- list.files("F:/SoilGrids", full.names = T, pattern = ".tif") %>%
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

writeRaster(soil, "assets/climate/soil_800m.tif",
            overwrite=T)




### species ranges #########

#f <- f[sub("\\.rds", "", basename(f)) %in% spp]

## transfer sdm to climate grid ##
f <- list.files("~/data/ca_plant_ranges_philtransb/", full.names=T)
template <- raster(clim$path[1])
project_range <- function(x = "NONE", p = .1, outdir = "assets/all_species/ranges/"){
      outfile <- paste0(outdir, x, ".tif")
      if(file.exists(outfile)) return("skipped")
      message(x)
      fx <- f[grepl(x, f)]
      if(x == "NONE") fx <- f[[1]]
      r <- raster(fx)
      r <- projectRaster(r, template)
      m <- cellStats(r, max) * p
      if(x == "NONE") m <- -1
      r <- r %>%
            reclassify(c(-1, m, NA,
                         m, 1, 1)) %>%
            trim()
      writeRaster(r, outfile, overwrite = T)
      return("projected")
}




##### smoothed environmental rasters per species #####

vars <- sort(c("PPT", "AET", "CWD", "DJF", "JJA"))
clim_files <- data.frame(path=list.files("assets/climate",
                                         full.names=T, pattern="ensemble"),
                         stringsAsFactors=F) %>%
      mutate(set = gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "scenario", "var"), remove = F, sep = " ") %>%
      select(-junk) %>%
      filter(var %in% vars) %>% arrange(var)

soil_all <- stack("assets/climate/soil_800m.tif")


# historic environments across species range
species_envt <- function(x, dir = "assets/all_species/ranges/"){
      range <- raster(paste0(dir, x, ".tif"))
      
      clim <- stackBands(clim_files$path[clim_files$year=="1981-2010"], 1)
      names(clim) <- vars
      
      range <- crop(range, clim)
      clim <- clim %>% crop(range) %>% mask(range) %>% trim()
      soil <- soil_all %>% crop(range) %>% mask(range) %>% trim()
      
      list(soil = crop(soil, clim),
           clim = crop(clim, soil))
}

# smoothed historic environment representing gene flow
smoothed_envt <- function(e, radius){
      
      require(terra)
      e <- c(rast(e[[1]]), rast(e[[2]]))
      dime <- dim(e)
      
      rz <- res(e[[1]])[1]
      w <- focalWeight(e[[1]], rz/2 + rz*radius, "circle")
      w[w>0] <- 1
      
      if(ncol(w) > 2 * ncol(e) |
         nrow(w) > 2 * nrow(e)){
            e <- extend(e, max(dim(w)))
      }
      
      if(radius > 0) e <- terra::focal(e, w=w, fun=mean, na.rm=T, na.policy="omit")
      
      # names(e)[6:10] <- paste0("soil", 1:5)
      return(e)
}


all_spp <- list.files("~/data/ca_plant_ranges_philtransb/") %>%
      str_remove("\\.tif") %>%
      sample(length(.)) %>%
      c("NONE")

sp_smooth <- function(x, indir = "assets/all_species/range/", outdir = "assets/all_species/smooth/"){
      require(raster)
      if(file.exists(paste0(dir, x, " ", 100, ".tif"))) return("skipped")
      message(x)
      z <- species_envt(x, dir = indir)
      for(r in c(0, 1, 2, 5, 10, 20, 50, 100)){
            message(r)
            out <- paste0(outdir, x, " ", r, ".tif")
            if(file.exists(out)) next()
            terra::writeRaster(smoothed_envt(z, r), out, overwrite = T)
      }
      return("done")
}


# species niche breadths
niche_breadth <- function(x, dir = "assets/seedsource_data/smooth/"){
      require(terra)
      r <- rast(paste0(dir, x, " ", 0, ".tif"))
      m <- global(r, "mean", na.rm = T) %>% 
            as.data.frame() %>%
            rownames_to_column("var") %>%
            mutate(var = str_replace(var, "800m_", "PC"))
      s <- global(r, "sd", na.rm = T) %>% 
            as.data.frame()
      m$sd <- s$sd
      m$species <- x
      m
}


library(furrr)
plan(multisession, workers = 6)

future_map(all_spp, project_range)
sm <- future_map(all_spp, sp_smooth)
niches <- future_map_dfr(all_spp, niche_breadth)
saveRDS(niches %>% as_tibble(), "assets/range_stats.rds")

# would be better to just copy select species data to the appropriate folder instead of regenerating. but:
future_map(key_spp, project_range, outdir = "assets/select_species/ranges/")
sm <- future_map(key_spp, sp_smooth, indir = "assets/select_species/ranges/", outdir = "assets/select_species/smooth/")



# mask nevada from all climate layers (for computational efficiency)
f <- list.files("assets/climate/", full.names = T)
msk <- rast("assets/all_species/smooth/NONE 0.tif")[[1]]
map(f, function(x) x %>% rast() %>% crop(msk) %>% mask(msk) %>% writeRaster(x, overwrite = T) )



ff <- f[grepl("ssp", f)]
hh <- f[!grepl("ssp|soil", f)]
for(fi in ff){
      var <- str_remove(basename(fi), "\\.tif") %>%
            str_sub(-3, -1)
      hi <- hh[grepl(var, hh)]
      delta <- rast(fi)[[1]] - rast(hi)[[1]]
      writeRaster(delta, fi %>% str_replace("climate", "deltas") %>% str_replace("ensemble", "delta"), overwrite = T)
}
