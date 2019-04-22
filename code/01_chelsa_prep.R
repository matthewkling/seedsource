


library(raster)
library(tidyverse)



#### crop global rasters to study extent ####

f <- c(list.files("F:/chelsa/cmip5", full.names=T),
       list.files("F:/chelsa/monthly48", full.names=T))

cropNsave <- function(x){
      outfile <- paste0("F:/seeds/chelsa_cropped/", basename(x))
      if(file.exists(outfile)) return("skipped")
      r <- raster(x)
      e <- extent(-125, -114, 28, 45)
      r <- crop(r, e)
      writeRaster(r, outfile, overwrite=T)
}

sapply(f, cropNsave)



#### parse metadata for raster files ####

md <- data.frame(path=list.files("F:/seeds/chelsa_cropped/", full.names=T), stringsAsFactors=F) %>%
      mutate(file=basename(path))
mdh <- filter(md, grepl("land", file)) %>%
      separate(file, c("junk1", "variable", "month", "junk2", "junk3")) %>%
      dplyr::select(-contains("junk")) %>%
      mutate(timeframe="1979-2013")
mdf <- filter(md, !grepl("land", file)) %>%
      separate(file, c("junk1", "variable", "junk2", "model", "scenario", 
                       "junk3", "junk4", "month", "timeframe"), sep="_") %>%
      dplyr::select(-contains("junk")) %>%
      mutate(timeframe=gsub("\\.tif", "", timeframe))
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
      
      outfile <- paste0("F:/seeds/chelsa_derived/", data$set[1], ".tif")
      if(file.exists(outfile)) return("skipped")
      
      message(data$set[1])
      
      # water balance variables
      source("e:/chilefornia/chilefornia/water_balance.R")
      lat <- latitude(data$path[1], already_latlong=T)
      r <- stack(c(data$path, lat)) %>%
            reclassify(c(-Inf, -1000, NA))
      names(r) <- c(paste0(data$variable, data$month), "latitude")
      water <- water_balance(r, temp_scalar=0.1, ncores=6)
      
      # seasonal temperature extremes
      djf <- mean(r[[c("tasmin12", "tasmin01", "tasmin02")]]/10)
      jja <- mean(r[[c("tasmax06", "tasmax07", "tasmax08")]]/10)
      
      # export
      writeRaster(stack(water, djf, jja), 
                  outfile, 
                  overwrite=T)
}

mds <- split(md, md$set)
lapply(mds, derive)




