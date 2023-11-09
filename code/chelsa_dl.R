
# "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/GFDL-ESM4/ssp585/pr/CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp585_pr_01_2071_2100_norm.tif"

library(tidyverse)

dl <- function(era = "2071-2100", 
               model = "GFDL-ESM4", 
               ssp = "ssp585", 
               variable = "pr", 
               month = "01"){
      
      root <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"
      filename <- paste0("CHELSA_", tolower(model), "_r1i1p1f1_w5e5_", ssp, "_", variable, "_", 
                         month, "_", str_replace(era, "-", "_"), "_norm.tif")
      url <- paste0(root, "/", era, "/", model, "/", ssp, "/", variable, "/", filename)
      
      dest <- paste0("~/data/climate/CHELSA/v2/cmip/", "/", filename)
      if(file.exists(dest)) return("skipped")
      
      options(timeout = 6000)
      download.file(url, dest)
}

expand_grid(era = c("2011-2040", "2041-2070", "2071-2100"), 
            model = c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"), 
            ssp = c("ssp126", "ssp370", "ssp585"), 
            variable = c("tasmax", "tasmin", "tas", "pr"), 
            month = str_pad(1:12, 2, "left", 0)) %>%
      pmap(dl)


