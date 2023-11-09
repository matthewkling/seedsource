


library(raster)
library(tidyverse)



#### crop global rasters to study extent ####

# f <- c(list.files("F:/chelsa/cmip5", full.names=T),
#        list.files("F:/chelsa/monthly48", full.names=T))

cropNsave <- function(x){
      outfile <- paste0("~/data/seeds/chelsa_cropped/", basename(x))
      if(file.exists(outfile)) return("skipped")
      r <- raster(x)
      e <- extent(-125, -114, 28, 45)
      r <- crop(r, e)
      writeRaster(r, outfile, overwrite=T)
}

f <- c(list.files("~/data/CHELSA/v2/cmip/", full.names = T))
sapply(f, cropNsave)

f <- c(list.files("~/data/CHELSA/v2/raw/", full.names = T))
f <- f[grepl("_tas|_pr_", f)]
sapply(f, cropNsave)


#### parse metadata for raster files ####

md <- tibble(path = list.files("~/data/seeds/chelsa_cropped/", full.names = T),
             file = basename(path))
mdh <- filter(md, grepl("1981-2010", file)) %>%
      separate(file, c("junk1", "variable", "month", "junk2", "junk3"), sep = "_") %>%
      dplyr::select(-contains("junk")) %>%
      mutate(timeframe = "1981-2010")
mdf <- filter(md, !grepl("1981-2010|_bio", file)) %>%
      separate(file, c("junk1", "model", "junk2", "junk3", "scenario",  "variable",  
                       "month", "year", "year2", "junk4"), sep="_") %>%
      mutate(timeframe = paste0(year, "-", year2)) %>%
      dplyr::select(-contains("junk"), -year, -year2)
md <- full_join(mdf, mdh) %>%
      mutate(variable=case_when(variable=="temp10" ~ "tas",
                                variable=="tmax10" ~ "tasmax",
                                variable=="tmin10" ~ "tasmin",
                                variable=="prec" ~ "pr",
                                TRUE ~ variable),
             var_ord=case_when(variable=="pr" ~ 1, # necessary for water balance function
                               variable=="tas" ~ 2,
                               variable=="tasmax" ~ 3,
                               variable=="tasmin" ~ 4),
             month=str_pad(month, 2, "left", 0),
             set=paste(timeframe, model, scenario)) %>%
      arrange(set, var_ord, month)




#### calculate derived variables ####

derive <- function(data){
      
      message(data$set[1])
      outfile <- paste0("~/data/seeds/chelsa_derived/", data$set[1], ".tif")
      if(file.exists(outfile)) return("skipped")
      
      require(hydro)
      
      r <- stack(data$path)
      for(i in 1:12) r[[i]] <- r[[i]] / 10
      for(i in 13:48) r[[i]] <- r[[i]] / 10 - 273.15
      water <- hydro(r, already_latlong = T, ncores = 6)
      
      # seasonal temperature extremes
      djf <- mean(r[[c(37, 38, 48)]])
      jja <- mean(r[[30:32]])
      
      # export
      terra::writeRaster(terra::rast(stack(water, djf, jja)), 
                         outfile, 
                         overwrite = T)
}

mds <- split(md, md$set)
lapply(mds, derive)




