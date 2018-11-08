

library(shiny)
library(htmlwidgets)
library(dplyr)
library(raster)
library(geosphere)
library(caret)
library(grid)
select <- dplyr::select
extract <- raster::extract


# load data
d <- readRDS("../assets/seedsource_data/species_data.rds")
rsoil <- stack("../assets/seedsource_data/soil.tif")
clim_files <- list.files("../assets/seedsource_data/", pattern="climate", full.names=T)
clim <- lapply(clim_files, stack)


vars <- data.frame(abbv=c("geo_dist", "env_dist_t1", 
                          "climPC1", "climPC2", "climPC3", 
                          "soilPC1", "soilPC2", "soilPC3"),
                   description=c("Geographic distance", "Environmental dissimilarity",
                                 "Climate PC1", "Climate PC2", "Climate PC3",
                                 "Soil PC1", "Soil PC2", "Soil PC3"),
                   stringsAsFactors = F)


ui <- fluidPage(
  
  titlePanel("Predictive provenancing"),
  
  fluidRow(
    column(2, strong(p("Double-click map to select planting site"))),
    column(2, selectizeInput("sp", "Focal\nspecies", d$spp)),
    column(3, sliderInput("pclim", "Climate importance (vs. soil)", min=0, max=100, value=50, post="%", width="100%")),
    column(3, sliderInput("radius", "Homogenization radius", min=0, max=200, value=25, post="km", width="100%")),
    column(2, selectInput("color_var", "Point color", vars$description, "Environmental dissimilarity"))
  ),
  
  fluidRow(
    tags$head(tags$style(".shiny-plot-output{height:80vh !important;}")),
    column(4, 
           plotOutput("map", dblclick="map_dblclick", brush="map_brush")),
    column(4, 
           plotOutput("env_scatter", brush="env_brush"),
           fluidRow(
             column(6, selectInput("env_xvar", "x-axis variable", vars$description[grepl("PC", vars$abbv)], "Climate PC1")),
             column(6, selectInput("env_yvar", "y-axis variable", vars$description[grepl("PC", vars$abbv)], "Soil PC1")))),
    column(4, 
           plotOutput("dist_scatter", brush="dist_brush")    )#,
    #fluidRow(
    #      column(6, selectInput("dist_xvar", "x-axis variable", vars, "geo_dist")),
    #     column(6, selectInput("dist_yvar", "y-axis variable", vars, "env_dist1")))),
    
  )
  
)


server <- function(input, output) {
  
  # default restoration site location: presidio
  s <- data.frame(species=NA, lat=37.801064, lon=-122.478557)
  coordinates(s) <- c("lon", "lat")
  crs(s) <- crs(d$jeps)
  site <- reactiveValues(point = s)
  
  # update location when map is clicked
  observeEvent(input$map_dblclick, {
    s <- data.frame(species=NA, lat=input$map_dblclick$y, lon=input$map_dblclick$x)
    coordinates(s) <- c("lon", "lat")
    crs(s) <- crs(d$jeps)
    site$point <- s
  })
  
  site_clim <- reactive({ lapply(clim, extract, y=site$point) }) # slow
  
  data <- reactive({
    
    # geographic distances from restoration site
    j <- d$jeps[d$jeps$species==input$sp,]
    gdist <- as.vector(distm(site$point, j))
    
    # climate and soil data across species range
    jc <- d$jeps_clim[d$jeps$species==input$sp,]
    js <- d$jeps_soil[d$jeps$species==input$sp,]
    
    # climate and soils data for restoration site
    sc <- lapply(site_clim(), function(x){colnames(x) <- colnames(jc); return(x)})
    ss <- extract(rsoil, site$point)
    colnames(ss) <- colnames(js)
    
    # 3 climate PCs
    ctrans <- preProcess(jc, c("center", "scale", "pca"), pcaComp=3)
    tag <- function(d, tag){colnames(d) <- paste0(tag, colnames(d)); return(d)}
    sc <- lapply(sc, function(x) predict(ctrans, x) %>% tag("clim"))
    jc <- predict(ctrans, jc) %>% tag("clim")
    
    # 3 soil PCs
    strans <- preProcess(js, c("center", "scale", "pca"), pcaComp=3)
    ss <- predict(strans, ss) %>% tag("soil")
    js <- predict(strans, js) %>% tag("soil")
    
    # combine soil and climate data
    je <- cbind(jc, js)
    se <- lapply(sc, function(x) cbind(x, ss))
    
    # spatial smoothing to represent homogenizing gene flow
    if(input$radius > 0){
      dists <- distm(j) %>% as.matrix()
      je <- lapply(1:nrow(dists), function(i){
        jx <- je[which(dists[i,] <= input$radius*1000),]
        if(is.null(nrow(jx))) return(jx)
        as.data.frame(jx) %>% summarize_all(mean, na.rm=T)}) %>%
        do.call("rbind", .)
    }
    
    # apply climate:soil weights
    jeb <- je # keep an unweighted copy
    if(input$pclim != 50){
      je[,1:3] <- je[,1:3] * input$pclim / 50
      je[,4:6] <- je[,4:6] * (100 - input$pclim) / 50
    }
    
    # environmental distances
    edist <- lapply(se, function(x){
      rbind(x, je) %>% dist() %>% as.matrix() %>% "["((2:nrow(.)), 1) }) %>%
      do.call("cbind", .)
    colnames(edist) <- paste0("env_dist_t", 1:ncol(edist))
    
    # convert to data frames for plotting
    jd <- as.data.frame(j) %>%
      mutate(geo_dist = gdist/1000) %>%
      cbind(edist) %>%
      cbind(jeb) %>% # jeb
      mutate(type="restoration",
             env_dist_t1_rank=rank(env_dist_t1))
    
    return(list(jd=jd,
                je=je,
                se=se))
  })
  
  
  # convert variable descriptions to abbreviations
  col_var <- reactive({vars$abbv[vars$description==input$color_var]})
  env_xvar <- reactive({vars$abbv[vars$description==input$env_xvar]})
  env_yvar <- reactive({vars$abbv[vars$description==input$env_yvar]})
  
  
  # brushed data selection
  brushed <- reactiveValues(data = data.frame(latitude=0, longitude=0, env_dist_t1=0, geo_dist=0,
                                              climPC1=0, climPC2=0, climPC3=0, 
                                              soilPC1=0, soilPC2=0, soilPC3=0)[0,])
  observeEvent(input$map_brush, {
    brushed$data <- brushedPoints(data()$jd, input$map_brush,
                                  xvar="longitude", yvar="latitude")})
  observeEvent(input$dist_brush, {
    brushed$data <- brushedPoints(data()$jd, input$dist_brush,
                                  xvar="geo_dist", yvar="env_dist_t1")})
  observeEvent(input$env_brush, {
    brushed$data <- brushedPoints(data()$jd, input$env_brush,
                                  xvar=env_xvar(), yvar=env_yvar())})
  
  data_clr <- reactive({
    data <- data()
    data$jd <- data$jd %>% 
      mutate_(clr = col_var()) %>%
      mutate(clr=ecdf(clr)(clr))
    return(data)
  })
  
  output$map <- renderPlot({
    sd <- as.data.frame(site$point)
    ggplot() + 
      geom_polygon(data=d$md, aes(long, lat, group=group), fill="gray90") +
      geom_point(data=data_clr()$jd, aes_string("longitude", "latitude", color="clr")) +
      geom_point(data=brushed$data, aes(longitude, latitude), color="orange") +
      geom_point(data=sd, aes(lon, lat), color="magenta", shape=10, size=5) +
      scale_color_viridis_c() +
      theme_void() +
      theme(legend.position="none") +
      coord_fixed(ratio=1.2)
  })
  
  output$dist_scatter <- renderPlot({
    ggplot() + 
      geom_point(data=data_clr()$jd, aes_string("geo_dist", "env_dist_t1", color="clr")) +
      geom_point(data=brushed$data, aes_string("geo_dist", "env_dist_t1"), color="orange") +
      annotate(geom="point", x=0, y=0, color="magenta", shape=10, size=5) +
      scale_x_log10() +
      scale_color_viridis_c() +
      theme_minimal() +
      theme(legend.position="none") +
      labs(x="geographic distance (km)",
           y="environmental dissimilarity (climate and soil)")
  })
  
  output$env_scatter <- renderPlot({
    
    se <- data()$se %>%
      do.call("rbind", .) %>%
      as.data.frame() %>%
      mutate(metadata=basename(clim_files),
             metadata=sub("climate ", "", metadata),
             metadata=sub("\\.tif", "", metadata),
             timeframe=substr(metadata, 1, 9),
             model=substr(metadata, 11, nchar(metadata)-6),
             scenario=substr(metadata, nchar(metadata)-4, nchar(metadata)),
             group=paste(model, scenario))
    sh <- filter(se, timeframe=="1979-2013") %>%
      select(-metadata, -model, -scenario, -group)
    sf <- filter(se, timeframe=="2041-2060") %>%
      mutate(timeframe="1979-2013") %>%
      select(metadata, timeframe, model, scenario, group) %>%
      full_join(sh, .)
    se <- rbind(se, sf) %>%
      arrange(group, timeframe)
    
    
    ggplot() + 
      geom_point(data=data_clr()$jd, aes_string(env_xvar(), env_yvar(), color="clr")) +
      geom_point(data=brushed$data, aes_string(env_xvar(), env_yvar()), color="orange") +
      geom_path(data=se, aes_string(env_xvar(), env_yvar(), group="group", order="timeframe"), color="magenta") +
      geom_point(data=se, aes_string(env_xvar(), env_yvar()), color="magenta") +
      scale_color_viridis_c() +
      theme_minimal() +
      theme(legend.position="none") +
      labs(x=input$env_xvar,
           y=input$env_yvar)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

