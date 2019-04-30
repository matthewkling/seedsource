

# to do:
# pre-load an initial set of results
# make smoothing window circular
# fix climate holes problem




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
clim_files <- data.frame(path=list.files("assets/seedsource_data/climate",
                                         full.names=T, pattern="ensemble"),
                         stringsAsFactors=F) %>%
      mutate(set=gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "var"), remove=F, sep=" ") %>%
      select(-junk) %>%
      filter(var %in% vars) %>% arrange(var)

spps <- list.files("assets/seedsource_data/ranges") %>% sub("\\.rds", "", .)

soil_all <- stack("assets/seedsource_data/climate/soil_800m.tif")

smoothed <- data.frame(path=list.files("assets/seedsource_data/smooth",
                                       full.names=T),
                       stringsAsFactors=F) %>%
      mutate(set=gsub("\\.rds", "", basename(path))) %>%
      separate(set, c("genus", "species", "radius"), remove=F, sep=" ") %>%
      mutate(gs=paste(genus, species),
             radius=as.integer(radius)) %>%
      select(path, gs, radius)



ll <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all_vars <- c(vars, paste0("soil_PC", 1:5), "clim_prob", "soil_prob", "prob")



ui <- fluidPage(
      
      useShinyalert(),
      
      column(4,
             fluidRow(
                   br(),
                   tags$img(src="logo2.png", width="50%", align="left"),
                   br()
             ),
             fluidRow(
                   column(6,
                          selectizeInput("sp", span("Focal species ", actionLink("i_species", "[?]")), 
                                         spps, selected="Quercus agrifolia"),
                          
                          shinyWidgets::sliderTextInput("radius", span("Neighborhood radius (km) ", actionLink("i_radius", "[?]")),
                                                        choices=sort(unique(smoothed$radius)), selected=10,
                                                        grid = T),
                          
                          selectInput("time", span("Target era", actionLink("i_time", "[?]")), 
                                      c("2041-2060", "2061-2080"), "2061-2080")
                   ),
                   column(6,
                          selectInput("climstat", span("Climatic similarity basis", actionLink("i_basis", "[?]")), 
                                      c("GCM ensemble", "species niche"), "species niche"),
                          
                          sliderInput("pclim", span("Climate importance (vs. soil)", actionLink("i_weight", "[?]")),
                                      min=0, max=1, value=.5, width="100%")
                   )
             ),
             hr(),
             fluidRow(
                   column(4,
                          selectizeInput("xvar", "X variable", vars, "PPT")
                   ),
                   column(4,
                          selectizeInput("yvar", "Y variable", vars, "JJA")
                   ),
                   column(4,
                          selectizeInput("color", "Color variable", all_vars, "prob")
                   )
             ),
             fluidRow(
                   plotOutput("scatter")
             )
      ),
      column(8,
             leafletOutput("map", height=1000)
      )
)


server <- function(input, output, session) {
      
      showModal(modalDialog(
            title="Seeds of Change",
            HTML("Welcome. This tool is aimed at helping to identify populations that may be suitable as seed sources for ecological restoration projects,",
                 "incorporating data on species ranges, soil variation, and future climate change.",
                 "<br><br>This version is an early prototype and still under development. Please contact mattkling@berkeley.edu with questions or bug reports."),
            easyClose = TRUE, footer = modalButton("Dismiss")
      ))
      
      observeEvent(input$i_species, 
                   {showModal(modalDialog(
                         title="Species range maps",
                         "The tool comes pre-loaded with geographic range maps for key species of California native plants.",
                         "Select a species to load its estimated California geographic range, which is modeled based on climate, distance to known observations, and landscape intactness.",
                         "We'll use the map to estimate the species' environmental tolerance, model gene flow among nearby populations, and hypothesize which populations will have suitable genetics for the planting site.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_radius, 
                   {showModal(modalDialog(
                         title="Incorporating gene flow",
                         "A source population's suitability as a genetic match for a planting site with the same environment depends on the population's historic balance between selection and gene flow.",
                         "The relative importance of gene swamping versus local adaptation is known to vary among species.",
                         "This parameter lets you set the size of the local neighborhood around each source population within which gene flow is expected to homogenize local adaptation.",
                         "Choose a small smoothing radius to model highly localized adaptation, or a large radius to model strong effects of widespread gene flow relative to local selection.",
                         "(Note that large values take a lot more computation time, so maybe walk down and grab a sandwich or something while you wait.)",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_time, 
                   {showModal(modalDialog(
                         title="Climate change is a moving target",
                         "This tool estimates how well the historic adaptive environment at each provenance site matches the projected future environment at the planting site.",
                         "For which future planting site time period would you like to estimate historic climatic similarity?",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_basis, 
                   {showModal(modalDialog(
                         title="Climate distance metric",
                         "The climatic difference 'sigma' between the planting site and each source population is expressed as a number of standard deviations.", 
                         "You can select whether this is based on variation among climate models, in which case sigma reflects the likelihood of a future climate match,",
                         "or on variation across the species California range, in which case sigma reflects difference in terms of realized niche breadth.",
                         "(Note that soil similarity is always calculated using the range method.)",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_weight, 
                   {showModal(modalDialog(
                         title="Climate versus soil",
                         "By default soil and climate are given equal weight when calculating similarity to the planting site environment.", 
                         "Adjust this slder to incorporate species-specific knowledge about their relative importance in shaping local adaptation,",
                         "or to explore how the results change based on assumptions about their importance.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      
      
      
      
      # default planting site location: presidio
      s <- data.frame(species=NA, lat=37.801064, lon=-122.478557)
      coordinates(s) <- c("lon", "lat")
      crs(s) <- ll
      site <- reactiveValues(point = s)
      
      # update planting site when map is clicked
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
      
      # historic climate at provenance sites
      smoothed_envt <- reactive({
            smoothed %>% 
                  filter(gs==input$sp,
                         radius==input$radius) %>%
                  pull(path) %>%
                  readRDS()
      })
      
      range_stats <- reactive({
            list(clim_mean = cellStats(smoothed_envt()$clim, "mean"),
                 clim_sd = cellStats(smoothed_envt()$clim, "sd", asSample=F),
                 soil_mean = cellStats(smoothed_envt()$soil, "mean"),
                 soil_sd = cellStats(smoothed_envt()$soil, "sd", asSample=F))
      })
      
      # climate differences between source populations and future planting site 
      clim_sigmas <- reactive({
            if(input$climstat=="GCM ensemble"){
                  site_mean <- site_future()$clim_mean
                  site_sd <- site_future()$clim_sd
                  ndim <- site_future()$clim_ndim
                  
                  sigma <- function(x, ...){
                        if(is.na(x[1])) return(NA)
                        sum((x - site_mean)^2 / site_sd^2) %>% 
                              pchisq(ndim) %>% 
                              qchisq(1) %>% 
                              sqrt()
                  }
                  clim_prob <- calc(smoothed_envt()$clim, sigma) %>%
                        reclassify(c(10, Inf, 10))
            }
            if(input$climstat=="species niche"){
                  clim <- smoothed_envt()$clim
                  range_mean <- range_stats()$clim_mean
                  range_sd <- range_stats()$clim_sd
                  clim <- (clim - range_mean) / range_sd
                  
                  site_mean <- (site_future()$clim_mean - range_mean) / range_sd
                  ndim <- site_future()$clim_ndim
                  
                  sigma <- function(x, ...){
                        if(is.na(x[1])) return(NA)
                        sum((x - site_mean)^2) %>% 
                              pchisq(ndim) %>% 
                              qchisq(1) %>% 
                              sqrt()
                  }
                  clim_prob <- calc(clim, sigma) %>%
                        reclassify(c(10, Inf, 10))
            }
            clim_prob
      })
      
      # soil differences
      soil_sigmas <- reactive({
            
            soil <- smoothed_envt()$soil
            range_mean <- range_stats()$soil_mean
            range_sd <- range_stats()$soil_sd
            soil <- (soil - range_mean) / range_sd
            
            site_mean <- (site_future()$soil_mean - range_mean) / range_sd
            ndim <- site_future()$soil_ndim
            
            sigma <- function(x, ...){
                  if(is.na(x[1])) return(NA)
                  sum((x - site_mean)^2) %>%
                        pchisq(ndim) %>% 
                        qchisq(1) %>% 
                        sqrt()
            }
            soil_prob <- calc(soil, sigma) %>%
                  reclassify(c(10, Inf, 10))
            
            soil_prob
      })
      
      # overall dissimilarity combining soil and climate
      final <- reactive({
            prob <- (clim_sigmas() * input$pclim) + (soil_sigmas() * (1-input$pclim))
            
            d <- stack(smoothed_envt()$clim, smoothed_envt()$soil, 
                       clim_sigmas(), soil_sigmas(), prob)
            names(d) <- all_vars
            list(rasters=d)
      })
      
      
      
      
      colors <- c("red", "yellow", "green", "blue", "darkblue", "black")
      colors <- c("red", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black")
      
      icon <- makeAwesomeIcon("leaf", markerColor="red")
      
      output$map <- renderLeaflet({
            
            withProgress(message = "Stand by.", 
                         value = 0, {
                               incProgress(.25, detail = "Computing climate similarities.")
                               x <- clim_sigmas()
                               incProgress(.5, detail = "Computing soil similarities.")
                               x <- soil_sigmas()
                               
                               incProgress(.75, detail = "Merging dimensions.")
                               r <- final()$rasters[[input$color]]
                               
                               incProgress(1, detail = "Generating plots.")
                               #browser()
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
      })
      
      output$scatter <- renderPlot({
            
            vx <- input$xvar
            vy <- input$yvar
            
            d <- final()$rasters[[c(vx, vy, "prob")]] %>%
                  rasterToPoints() %>% as.data.frame() %>%
                  #arrange(desc(prob)) %>%
                  na.omit()
            names(d)[3:4] <- c("xvar", "yvar")
            
            # future climate mean, and sd of either ensemble or niche 
            avg <- c(site_future()$clim_mean, site_future()$soil_mean)
            if(input$climstat == "GCM ensemble"){
                  std <- c(site_future()$clim_sd, range_stats()$soil_sd)
            }else{
                  std <- c(range_stats()$clim_sd, range_stats()$soil_sd)
            }
            
            fill_col <- NA
            
            ggplot() +
                  geom_point(data=d, aes(xvar, yvar, color=prob), size=.5) +
                  scale_color_gradientn(colours=colors, limits=c(0, NA)) +
                  annotate("polygon", color="red", fill=fill_col, alpha=.1,
                           x=avg[vx] + std[vx] * cos(seq(0,2*pi,length.out=100)),
                           y=avg[vy] + std[vy] * sin(seq(0,2*pi,length.out=100))) +
                  annotate("polygon", color="red", fill=fill_col, alpha=.1,
                           x=avg[vx] + std[vx]*2 * cos(seq(0,2*pi,length.out=100)),
                           y=avg[vy] + std[vy]*2 * sin(seq(0,2*pi,length.out=100))) +
                  annotate("polygon", color="red", fill=fill_col, alpha=.1,
                           x=avg[vx] + std[vx]*3 * cos(seq(0,2*pi,length.out=100)),
                           y=avg[vy] + std[vy]*3 * sin(seq(0,2*pi,length.out=100))) +
                  #annotate("segment", color="red",
                  #         x=hst[vx], y=hst[vy], xend=fut[vx], yend=fut[vy],
                  #         arrow=grid::arrow(type="closed", angle=15, length=unit(.15, "in"))) +
                  theme_minimal() +
                  theme(legend.position="none") +
                  labs(x=vx, y=vy)
      })
      
}

# Run the application
shinyApp(ui = ui, server = server)

