

# to do:
# preprocess raster stacks to avoid having to use stackBands



library(shiny)
library(rintrojs)
library(shinyalert)
library(dplyr)
library(tidyr)
library(ggplot2)
library(leaflet)
library(raster)
library(grid)
library(ecoclim)
select <- dplyr::select
extract <- raster::extract


# climate metadata
vars <- sort(c("PPT", "AET", "CWD", "DJF", "JJA"))
clim_files <- data.frame(path=list.files("../assets/seedsource_data",
                                         full.names=T, pattern="ensemble"),
                         stringsAsFactors=F) %>%
      mutate(set=gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "var"), remove=F, sep=" ") %>%
      select(-junk) %>%
      filter(var %in% vars) %>% arrange(var)

spps <- list.files("../assets/seedsource_data/ranges") %>% sub("\\.rds", "", .)

soil_all <- stack("../assets/seedsource_data/soil_800m.tif")



ll <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all_vars <- c(vars, paste0("soil_PC", 1:5), "clim_prob", "soil_prob", "prob")



ui <- fluidPage(
      
      useShinyalert(),
      
      column(3,
             br(),
             tags$img(src="logo2.png", width="50%", align="left"),
             br(),
             br(),
             br(),
             br(),
             br(),
             
             #fluidRow(column(10, selectizeInput("sp", "Focal species", spps)), 
             #         column(1, h1(" "), actionButton("i_species", "?"))),
             selectizeInput("sp", "Focal species", spps),
             
             shinyWidgets::sliderTextInput("radius", "Homogenization radius (km)",
                                           choices=c(0, 1, 2, 5, 10, 20, 50, 100), selected=5,
                                           grid = T),
             
             selectInput("time", "Target era", c("2041-2060", "2061-2080"), "2061-2080"),
             
             selectInput("climstat", "Climate similarity basis", c("GCM ensemble", "species niche"), "GCM ensemble"),
             
             sliderInput("pclim", "Climate importance (vs. soil)",
                         min=0, max=1, value=.5, width="100%"),
             
             selectizeInput("color", "Color variable", all_vars, "prob"),
             br(),
             plotOutput("scatter", height=300),
             selectizeInput("xvar", "x variable", vars, "PPT"),
             selectizeInput("yvar", "y variable", vars, "JJA")
      ),
      column(9,
             leafletOutput("map", height=1000)
      )
)


server <- function(input, output, session) {
      
      observeEvent(input$i_species, {
            shinyalert("Select a species", "This will load the geographic range map", type="info")
      })
      
      
      
      
      # default restoration site location: presidio
      s <- data.frame(species=NA, lat=37.801064, lon=-122.478557)
      coordinates(s) <- c("lon", "lat")
      crs(s) <- ll
      site <- reactiveValues(point = s)
      
      
      
      # update location when map is clicked
      observeEvent(input$map_click, {
            s <- data.frame(species=NA, lat=input$map_click$lat, lon=input$map_click$lng)
            coordinates(s) <- c("lon", "lat")
            crs(s) <- ll
            site$point <- s
      })
      
      
      # future environment at planting site
      site_future <- reactive({
            time <- input$time
            
            clim_mean <- stackBands(clim_files$path[clim_files$year==time], 1) %>% extract(site$point) %>% as.vector()
            clim_sd <- stackBands(clim_files$path[clim_files$year==time], 2) %>% extract(site$point) %>% as.vector()
            names(clim_mean) <- names(clim_sd) <- vars
            nclimdim <- length(clim_mean)
            
            soil_mean <- extract(soil_all, site$point)
            colnames(soil_mean) <- paste0("soil_PC", 1:ncol(soil_mean))
            nsoildim <- length(soil_mean)
            
            list(clim_mean = clim_mean,
                 clim_sd = clim_sd,
                 clim_ndim = nclimdim,
                 soil_mean = soil_mean,
                 soil_ndim = nsoildim)
      })
      
      # historic environments across species range
      species_envt <- reactive({
            range <- readRDS(paste0("../assets/seedsource_data/ranges/", input$sp, ".rds"))
            
            clim <- stackBands(clim_files$path[clim_files$year=="1979-2013"], 1)
            names(clim) <- vars
            #clim_hist <- extract(clim, site$point) %>% as.vector() ########## problem for siloing ######
            range <- crop(range, clim)
            clim <- clim %>% crop(range) %>% mask(range) %>% trim()
            
            soil <- soil_all %>% mask(range) %>% trim()
            
            list(soil = crop(soil, clim),
                 clim = crop(clim, soil))
      })
      
      # smoothed historic environment representing gene flow
      smoothed_envt <- reactive({
            radius <- input$radius
            w <- matrix(1, 1+radius*2, 1+radius*2)
            center <- length(w) / 2 + .5
            
            m <- function(x, ...){
                  if(is.na(x[center])) return(NA)
                  mean(na.omit(x))
            }
            
            clim <- species_envt()$clim
            soil <- species_envt()$soil
            
            if(radius>0){
                  for(i in 1:nlayers(clim)){
                        clim[[i]] <- focal(clim[[i]], w=w, fun=m) %>%
                              mask(clim[[1]])
                  }
                  for(i in 1:nlayers(soil)){
                        soil[[i]] <- focal(soil[[i]], w=w, fun=m) %>%
                              mask(soil[[1]])
                  }
            }
            
            list(soil = soil,
                 clim = clim)
      })
      
      # climate and soil differences between source populations and future planting site 
      sigmas <- reactive({
            
            # for every cell in species range, stdevs between historic clim and planting site future
            # (convert to probability and then sigma, correcting for dimensionality using chi sq distribution)
            sed <- function(x, ...){
                  if(is.na(x[1])) return(NA)
                  sum((x - site_future()$clim_mean)^2 / site_future()$clim_sd^2)# SED ^ 2
            }
            sigma <- function(x, ...){
                  if(is.na(x[1])) return(NA)
                  pchisq(x, site_future()$clim_ndim) %>% qchisq(1) %>% sqrt()
            }
            clim_prob <- calc(smoothed_envt()$clim, sed) %>% mask(smoothed_envt()$clim[[1]])
            clim_prob <- calc(clim_prob, sigma)
            clim_prob <- reclassify(clim_prob, c(10, Inf, 10))
            
            # for every cell in range, difference in soil, corrected for dimensionality
            sed <- function(x, ...){
                  if(is.na(x[1])) return(NA)
                  sum((x - site_future()$soil_mean)^2)# SED ^ 2
            }
            sigma <- function(x, ...){
                  if(is.na(x[1])) return(NA)
                  pchisq(x, site_future()$soil_ndim) %>% qchisq(1) %>% sqrt()
            }
            soil_prob <- calc(smoothed_envt()$soil, sed) %>% mask(smoothed_envt()$soil[[1]])
            soil_prob <- calc(soil_prob, sigma)
            soil_prob <- reclassify(soil_prob, c(10, Inf, 10))
            
            
            list(soil = soil_prob,
                 clim = clim_prob)
      })
      
      # overall dissimilarity combining soil and climate
      final <- reactive({
            
            prob <- (sigmas()$clim * input$pclim) + (sigmas()$soil * (1-input$pclim))
            
            d <- stack(smoothed_envt()$clim, smoothed_envt()$soil, 
                       sigmas()$clim, sigmas()$soil, prob)
            names(d) <- all_vars
            
            #cd <- rbind(clim_hist, clim_mean, clim_sd)
            
            list(rasters=d)#,
            #site_clim=cd,
            #site_soil=soil_mean)
      })
      
      
      
      
      colors <- c("red", "yellow", "seagreen1", "blue", "darkblue", "black")
      
      icon <- makeAwesomeIcon("leaf", markerColor="red")
      
      output$map <- renderLeaflet({
            
            r <- final()$rasters[[input$color]]
            
            pal <- colorNumeric(colors,
                                c(0, values(r)),
                                na.color = "transparent")
            
            latlon <- coordinates(site$point)
            
            #Stamen.TonerBackground
            #Esri.WorldTerrain
            #Esri.WorldGrayCanvas
            leaflet() %>%
                  setView(lng=latlon[1], lat=latlon[2], zoom=10) %>%
                  addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
                  addRasterImage(r, colors=pal, opacity=0.8) %>%
                  addAwesomeMarkers(lng=latlon[1], lat=latlon[2], icon=icon) %>%
                  addLegend(pal=pal, values=c(0, values(r)),
                            opacity=0.8, title="Sigma")
            
            
      })
      
      output$scatter <- renderPlot({
            
            vx <- input$xvar
            vy <- input$yvar
            
            d <- final()$rasters[[c(vx, vy, "prob")]] %>%
                  rasterToPoints() %>% as.data.frame() %>%
                  #arrange(desc(prob)) %>%
                  na.omit()
            names(d)[3:4] <- c("xvar", "yvar")
            
            #clim_sd <- final()$site_clim["clim_sd",]
            #clim_mean <- final()$site_clim["clim_mean",]
            #clim_h <- final()$site_clim["clim_hist",]
            #soil_mean <- final()$site_soil
            
            fill_col <- NA
            
            ggplot() +
                  geom_point(data=d, aes(xvar, yvar, color=prob), size=.5) +
                  scale_color_gradientn(colours=colors, limits=c(0, NA)) +
                  #annotate("polygon", color="red", fill=fill_col, alpha=.1,
                  #         x=clim_mean[vx] + clim_sd[vx] * cos(seq(0,2*pi,length.out=100)),
                  #         y=clim_mean[vy] + clim_sd[vy] * sin(seq(0,2*pi,length.out=100))) +
                  #annotate("polygon", color="red", fill=fill_col, alpha=.1,
                  #         x=clim_mean[vx] + clim_sd[vx]*2 * cos(seq(0,2*pi,length.out=100)),
                  #         y=clim_mean[vy] + clim_sd[vy]*2 * sin(seq(0,2*pi,length.out=100))) +
                  #annotate("polygon", color="red", fill=fill_col, alpha=.1,
                  #         x=clim_mean[vx] + clim_sd[vx]*3 * cos(seq(0,2*pi,length.out=100)),
                  #         y=clim_mean[vy] + clim_sd[vy]*3 * sin(seq(0,2*pi,length.out=100))) +
                  #annotate("segment", color="red",
                  #         x=clim_h[vx], y=clim_h[vy], xend=clim_mean[vx], yend=clim_mean[vy],
            #         arrow=grid::arrow(type="closed", angle=15, length=unit(.15, "in"))) +
            theme_minimal() +
                  theme(legend.position="none") +
                  labs(x=vx, y=vy)
      })
      
}

# Run the application
shinyApp(ui = ui, server = server)

