


library(shiny)
library(rintrojs)
library(shinyalert)
library(dplyr)
library(tidyr)
library(ggplot2)
library(leaflet)
library(raster)
library(terra)
library(grid)
library(ecoclim)
library(purrr)
library(stringr)
select <- dplyr::select
extract <- terra::extract


assets <- ifelse(dir.exists("assets/all_species"), "assets/all_species/", "assets/select_species/")

# climate metadata
vars <- sort(c("PPT", "AET", "CWD", "DJF", "JJA"))
clim_files <- data.frame(path=list.files("assets/climate",
                                         full.names=T, pattern="ensemble"),
                         stringsAsFactors=F) %>%
      mutate(set=gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "ssp", "var"), remove=F, sep=" ") %>%
      select(-junk) %>%
      mutate(ssp = ifelse(ssp == "NA", "historic", toupper(ssp))) %>%
      filter(var %in% vars) %>% arrange(var)

times <- unique(clim_files$year)
ssps <- unique(clim_files$ssp)
ssps <- data.frame(ssp = ssps,
                   text = ifelse(ssps == "historic", ssps,
                                 paste0(str_sub(ssps, 1, 4), "-", str_sub(ssps, 5, 5), ".",
                                        str_sub(ssps, 6, 6))))

spps <- list.files(paste0(assets, "ranges")) %>% sub("\\.tif", "", .)
allspps <- readRDS("assets/all_species.rds")

soil_all <- rast("assets/climate/soil_800m.tif")
names(soil_all) <- paste0("soil_PC", 1:nlyr(soil_all))

smoothed <- data.frame(path=list.files(paste0(assets, "smooth"), full.names=T),
                       stringsAsFactors=F) %>%
      mutate(set = gsub("\\.tif", "", basename(path)),
             set = str_replace(set, "NONE ", "NONE NONE ")) %>%
      separate(set, c("genus", "species", "radius"), remove=F, sep=" ") %>%
      mutate(gs = paste(genus, species),
             gs = str_replace(gs, "NONE NONE", "NONE"),
             radius = as.integer(radius)) %>%
      select(path, gs, radius)

range_summaries <- readRDS("assets/range_stats.rds")

ll <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all_vars <- data.frame(abbv = c(vars, paste0("soil_PC", 1:5), "clim_prob", "soil_prob", "prob", "rgb", paste0("d", vars)))
all_vars$desc <- c("Actual evapotranspiration", "Climatic water deficit",
                   "Winter minimum temperature", "Summer maximum temperature",
                   "Annual precipitation",
                   "Soil PC1", "Soil PC2", "Soil PC3", "Soil PC4", "Soil PC5",
                   "Climate sigma", "Soil sigma", "Combined sigma",
                   "Ordination",
                   paste0("Change in ", c("Actual evapotranspiration", "Climatic water deficit",
                                          "Winter minimum temperature", "Summer maximum temperature", "Annual precipitation")))
all_vars$units <- c("(mm)", "(mm)",
                    "(deg. C)", "(deg. C)",
                    "(log10 mm)",
                    rep("(SD)", 5),
                    rep("(similarity to target, SD)", 3),
                    "(similar colors have\nsimilar enviroments)",
                    "(mm)", "(mm)", "(deg. C)", "(deg. C)", "(log10 mm)")
all_vars$min <- c(rep(NA, 10), rep(0, 3), NA, rep(NA, 5))
all_vars$palette <- c(rep("viridis", 10), rep("sigma", 3), NA, rep("delta", 5))

sigma <- function(x, site_mean){
      z <- sqrt(sum((x - site_mean)^2))
      z[] <- sqrt(qchisq(pchisq(z[], length(site_mean)), 1))
      classify(z, matrix(c(10, Inf, 10), nrow = 1))
}


# default focal site location: presidio
lonlat <- "-122.48, 37.80"
s <- data.frame(species=NA, lat=37.80, lon=-122.48)
coordinates(s) <- c("lon", "lat")
crs(s) <- ll


# define a keyPressed input that activates when RETURN is pressed
js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'

ui <- navbarPage("Seeds of Change [BETA]",
                 fluid = TRUE,
                 selected = "Home",
                 tags$script(js),
                 
                 tabPanel("Home",
                          column(4,
                                 fluidRow(
                                       column(4,
                                              selectizeInput("sp", 
                                                             span("Focal species ", actionLink("i_species", "[?]")), 
                                                             choices = spps, selected = "Quercus agrifolia"),
                                              textInput("lonlat",
                                                        span("Location ", actionLink("i_location", "[?]")), 
                                                        paste(s$lon, s$lat, sep = ", ")
                                              ),
                                              selectInput("mode", 
                                                          span("Focal site activity ", actionLink("i_mode", "[?]")), 
                                                          c("planting", "collection"), "planting")
                                       ),
                                       column(4,
                                              selectInput("time", 
                                                          span("Time period ", actionLink("i_time", "[?]")), 
                                                          times, "2071-2100"),
                                              selectInput("ssp", 
                                                          span("Scenario ", actionLink("i_ssp", "[?]")), 
                                                          ssps$text, "SSP5-8.5"),
                                              uiOutput("constraintControls")
                                       ),
                                       column(4,
                                              shinyWidgets::sliderTextInput("radius", 
                                                                            span("Smoothing radius", actionLink("i_radius", "[?]")),
                                                                            choices=sort(unique(smoothed$radius)), selected = 2,
                                                                            post = " km", grid = T),
                                              sliderInput("pclim", 
                                                          span("Soil versus Climate", actionLink("i_weight", "[?]")),
                                                          min=0, max=100, value=50, step=10, post = "% clim",
                                                          width="100%")
                                       )
                                       
                                 ),
                                 hr(),
                                 fluidRow(
                                       column(4, span(textOutput("target_label"), style="color:red")),
                                       column(4, textOutput("points_label")),
                                       column(4, span(textOutput("color_label"), style="color:darkblue")),
                                 ),
                                 fluidRow(
                                       plotOutput("scatter")
                                 ),
                                 fluidRow(
                                       column(4,
                                              selectizeInput("xvar", "X variable", setdiff(all_vars$desc, "ordination"), all_vars$desc[2])
                                       ),
                                       column(4,
                                              selectizeInput("yvar", "Y variable", setdiff(all_vars$desc, "ordination"), all_vars$desc[1])
                                       ),
                                       column(4,
                                              selectizeInput("color", "Color variable", 
                                                          all_vars$desc[c(13:11, 14, 1:10, 15:19)], 
                                                          all_vars$desc[13])
                                       )
                                 ),
                                 hr(),
                                 fluidRow( column(4, h5("Download ", actionLink("i_download", "[?]"))) ),
                                 fluidRow( column(4, downloadButton("download", "")) )
                                 
                          ),
                          column(8,
                                 leafletOutput("map", height=1000)
                          )
                 ),
                 
                 tabPanel("Instructions",
                          fluidRow(
                                column(12,
                                       htmltools::includeMarkdown("assets/instructions.md")
                                       # uiOutput('instructions')
                                )
                          )
                 ),
                 
                 tabPanel("About",
                          fluidRow(
                                column(12,
                                       htmltools::includeMarkdown("assets/about.md"),
                                       # tags$img(src="logo5.png", width="200px", align="center"),
                                       br(),
                                       tags$img(src="BIPPB_logo.jpg", width="300px", align="center"),
                                       
                                )
                          )
                 )
)



server <- function(input, output, session) {
      
      showModal(modalDialog(
            title="Seeds of Change",
            HTML("Welcome. This tool is aimed at comparing environmental patterns across California plant species ranges to support decisions around ecological restoration projects,",
                 "incorporating data on species ranges, soil variation, and future climate change. Details can be found on the 'Instructions' tab, and by clicking the '[?]' icons.",
                 "<br><br>This is a beta version still under active development. Documentation is incomplete, the tool has not been thoroughly tested, and bugs are possible. Please contact mattkling@berkeley.edu with questions or bug reports."),
            easyClose = TRUE, footer = modalButton("Dismiss")
      ))
      
      observeEvent(input$i_species, 
                   {showModal(modalDialog(
                         title="Species",
                         HTML("The tool comes pre-loaded with estimated geographic range maps for most California native plant species (only 25 species important for Bay Area restoration projects are included in this beta version, but all 5221 species listed below will ultimately be added).",
                         "Select a species to load its estimated California geographic range, which is modeled based on climate, distance to known observations, and landscape intactness.",
                         "We'll use this distribution to estimate the species' environmental tolerance, model gene flow among nearby populations, and hypothesize which populations may be adapted to the environment of the selected focal site.",
                         "Or for a generic species-agnostic analysis of environmental similarity between sites, enter 'NONE' in the species box.<br><br>"),
                         selectizeInput("allsp", "Full species list", choices = allspps), # not used server-side
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_radius, 
                   {showModal(modalDialog(
                         title="Smoothing radius (km)",
                         "A source population's suitability as a genetic match for a planting site with the same environment depends on the population's historic balance between selection and gene flow.",
                         "The relative importance of gene swamping versus local adaptation is known to vary among species.",
                         "This parameter lets you set the size of the local neighborhood around each source population within which gene flow is expected to homogenize local adaptation.",
                         "Choose a small smoothing radius to model highly localized adaptation, or a large radius to model strong effects of widespread gene flow relative to local selection.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_time, 
                   {showModal(modalDialog(
                         title="Time period",
                         "This tool estimates how well the historic adaptive environment at each provenance site matches the projected future environment at the planting site, or the inverse, if site activity is set to `collection.`",
                         "For which future time period would you like to estimate climatic similarity?",
                         "You can also select the historic period, in which case baseline climate data is used for both the focal site and species range",
                         "(Note that soil data does not change across time periods.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_ssp, 
                   {showModal(modalDialog(
                         title = "Intergovernmental Panel on Climate Change (IPCC) emissions scenarios (Shared Socio-economic Pathways (SSP))",
                         "The IPCC has defined a range of potential future climate trajectories based on different assumptions about climate change mitigation. Select which of these SSPs you would like to compare to historic climate.",
                         br(), br(), "SSP1-2.6 - effective sustainability policies, very low greenhouse gas emissions, net zero by 2050, global average temperature increase in year 2100 ~1.8 deg. C (range 1.3-2.4 deg. C) above pre-industrial",
                         br(), "SSP3-7.0 - regional imbalance, no additional policies, high greenhouse gas emissions, doubling 2020 to 2100, global average temperature increase in year 2100 ~3.6 deg. C (range 2.8-4.6 deg. C) above pre-industrial",
                         br(), "SSP5-8.5 - fossil fuel-intensive development, very high greenhouse gas emissions, doubling 2020 to 2050, global average temperature increase in year 2100 ~4.4 deg. C (range 3.3-5.7 deg. C) above pre-industrial",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      # observeEvent(input$i_basis, 
      #              {showModal(modalDialog(
      #                    title="Climate distance metric",
      #                    "The climatic difference 'sigma' between the planting site and each source population is expressed as a number of standard deviations.", 
      #                    "You can select whether this is based on variation among climate models, in which case sigma reflects the likelihood of a future climate match,",
      #                    "or on variation across the species California range, in which case sigma reflects difference in terms of realized niche breadth.",
      #                    "(Note that soil similarity is always calculated using the range method.)",
      #                    easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_weight, 
                   {showModal(modalDialog(
                         title="Soils versus climate",
                         "By default, soil and climate are given equal weight when calculating similarity to the focal site environment.", 
                         "Adjust this slider to incorporate species-specific knowledge about their relative importance in shaping local adaptation,",
                         "or to explore how the results change based on assumptions about their importance.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_constrain, 
                   {showModal(modalDialog(
                         title="Limit to species range",
                         "Check this box to limit prospective planting sites to locations within the species' current modeled range.",
                         "Uncheck it to consider all of California (which may be useful, e.g., if the modeled range omits areas where the species occurs).", 
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_mode, 
                   {showModal(modalDialog(
                         title = "Target site activity",
                         "Select whether the selected focal location is a 'planting' or 'collection' site.",
                         "If it's a planting site, its environment in the selected 'time period' will be compared to historic smoothed environments across the species range.", 
                         "If it's a collection site, its historic smoothed environment will be compared to range-wide environments from the selected time period.", 
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_download, 
                   {showModal(modalDialog(
                         title = "Download results",
                         "Click the 'download' button to save a raster file with the 'combined sigma' data for the current settings.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_location, 
                   {showModal(modalDialog(
                         title = "Select a focal location",
                         "To choose a focal site, click the map or enter 'Lon, Lat' in the box and press ENTER.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      
      
      
      # txt2html <- function(x){
      #       rawText <- readLines(x)
      #       splitText <- stringi::stri_split(str = rawText, regex = '\\n')
      #       replacedText <- lapply(splitText, shiny::p)
      #       return(replacedText)
      # }
      # output$instructions <- renderUI({ txt2html("assets/instructions") })
      # output$about <- renderUI({ txt2html("assets/about") })
      
      
      
      site <- reactiveValues(point = s)
      
      # update focal site when map is clicked
      observeEvent(input$map_click, {
            s <- data.frame(species=NA, lat=input$map_click$lat, lon=input$map_click$lng)
            updateTextInput(session, "lonlat", value = paste(round(s$lon, 2), round(s$lat, 2), sep = ", "))
            coordinates(s) <- c("lon", "lat")
            crs(s) <- ll
            site$point <- s
      })
      
      # update focal site via input text
      observeEvent(input$keyPressed, {
            crds <- as.numeric(stringr::str_split(input$lonlat, ", ")[[1]])
            s <- data.frame(species=NA, lat=crds[2], lon=crds[1])
            coordinates(s) <- c("lon", "lat")
            crs(s) <- ll
            site$point <- s
            # site$ll <- coordinates(s)
      })
      
      output$constraintControls <- renderUI({
            if(input$mode == "collection") checkboxInput("constrain", span("Limit results to species range", actionLink("i_constrain", "[?]")), TRUE)
      })
      
      ssp_sel <- reactiveValues(sel = "SSP585")
      observeEvent(input$ssp, { ssp_sel$sel <- c(ssp_sel$sel, ssps$ssp[ssps$text == input$ssp]) })
      
      observe({
            baseline <- as.character(input$time == "1981-2010")
            updateSelectInput(session, "ssp",
                              choices = switch(baseline, 
                                               "TRUE" = "historic",
                                               "FALSE" = setdiff(ssps$text, "historic")),
                              selected = ifelse(baseline == "TRUE", 
                                                "historic", ssps$text[ssps$ssp == tail(ssp_sel$sel[ssp_sel$sel != "historic"], 1)] ))
      })
      
      # for color variables, only show deltas if future timeframe is selected
      # observe({
      #       updateSelectInput(session, "color",
      #                         choices = switch(as.character(input$time == "1981-2010"),
      #                                          "TRUE" = all_vars$desc[c(13:11, 14, 1:10)],
      #                                          "FALSE" = all_vars$desc[c(13:11, 14, 1:10, 15:19)]))
      # })
      
      
      # historic
      smoothed_envt <- reactive({
            req(input$sp)
            y <- smoothed %>% 
                  filter(gs==input$sp,
                         radius==input$radius) %>%
                  pull(path) %>%
                  terra::rast()
            names(y) <- sub("800m_", "PC", names(y))
            return(list(clim = y[[6:10]],
                        soil = y[[1:5]]))
      })
      
      ref_envt <- list(clim = stackBands(clim_files$path[clim_files$year=="1981-2010"], 1) %>%
                             rast() %>%
                             setNames(vars))
      
      # future 
      future_envt <- reactive({
            time <- input$time
            ssp <- ssps$ssp[ssps$text == input$ssp]
            if(time == times[1]) ssp <- ssps$ssp[1]
            clim_mean <- stackBands(clim_files$path[clim_files$year==time & clim_files$ssp==ssp], 1) %>% rast()
            clim_sd <- stackBands(clim_files$path[clim_files$year==time & clim_files$ssp==ssp], 2) %>% rast()
            names(clim_mean) <- vars
            names(clim_sd) <- vars
            list(clim = clim_mean,
                 clim_sd = clim_sd,
                 soil = soil_all)
      })
      
      deltas <- reactive({
            time <- input$time
            ssp <- ssps$ssp[ssps$text == input$ssp]
            var <- all_vars$abbv[all_vars$desc == input$color]
            if(! str_detect(input$color, "Change in")) return(NULL)
            files <- list.files("assets/deltas", full.names = T)
            files <- files[grepl(time, files) & 
                                 grepl(tolower(ssp), files) &
                                 grepl(str_remove(var, "d"), files)]
            msk <- smoothed_envt()[[1]][[1]]
            rast(files) %>% crop(msk) %>% mask(msk) %>% setNames(var)
      })
      
      # reference climate across range (historic if planting mode)
      range_envt <- reactive({
            y <- switch(input$mode,
                        "planting" = smoothed_envt(),
                        "collection" = future_envt() )
            if(input$mode == "collection"){
                  if(!is.null(input$constrain)){
                        if(input$constrain) y <- y %>%
                                    map(crop, y = smoothed_envt()[[1]][[1]]) %>%
                                    map(terra::mask, mask = smoothed_envt()[[1]][[1]])
                  }
            }
            return(y)
      })
      
      # gcm_var <- reactive({
      #       x <- smoothed_envt()$clim[[1]]
      #       ss <- sample_n(as.data.frame(rasterToPoints(x)), 10)[,1:2]
      #       coordinates(ss) <- c("x", "y")
      #       crs(ss) <- ll
      #       ss <- rbind(ss, site$point)
      #       future_envt() %>% map(extract, y = ss)
      # })
      
      # target environment at selected site (future if planting mode)
      target_envt <- reactive({
            future <- future_envt() %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            historic <- smoothed_envt() %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            reference <- ref_envt %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            te <- list(focal = switch(input$mode, "planting" = future, "collection" = historic),
                       reference = switch(input$mode, "planting" = reference, "collection" = future))
            if(all(is.na(te$focal$clim))){
                  showModal(modalDialog(
                        title="Oops", "It seems your selected location is outside the allowed area. Please choose a site in California, and if focal site is set to 'collection', please choose a site within the species range.",
                        easyClose = TRUE, footer = modalButton("Dismiss") ))
                  updateSelectInput(session, "mode", selected = "planting")
            }
            return(te)
      })
      
      range_stats <- reactive({
            y <- filter(range_summaries,
                        species == input$sp)
            list(clim_mean = y$mean[6:10],
                 clim_sd = y$sd[6:10],
                 soil_mean = y$mean[1:5],
                 soil_sd = y$sd[1:5])
      })
      
      # soil differences
      soil_sigmas <- reactive({
            x <- range_envt()$soil
            range_mean <- range_stats()$soil_mean
            range_sd <- range_stats()$soil_sd
            x <- (x - range_mean) / range_sd
            site_mean <- (unlist(target_envt()$focal$soil) - range_mean) / range_sd
            sigma(x, site_mean)
      })
      
      # climate differences
      clim_sigmas <- reactive({
            x <- range_envt()$clim
            range_mean <- range_stats()$clim_mean
            range_sd <- range_stats()$clim_sd
            x <- (x - range_mean) / range_sd
            site_mean <- (unlist(target_envt()$focal$clim) - range_mean) / range_sd
            sigma(x, site_mean)
      })
      
      # overall dissimilarity combining soil and climate
      final <- reactive({
            pc <- input$pclim / 100
            req(clim_sigmas(), soil_sigmas())
            prob <- (clim_sigmas() * pc) + (soil_sigmas() * (1-pc))
            d <- c(range_envt()$clim, range_envt()$soil, clim_sigmas(), soil_sigmas(), prob)
            names(d) <- all_vars$abbv[1:13]
            if(!is.null(deltas())) d <- c(d, deltas())
            return(d)
      })
      
      
      
      sigma_pal <- c("red", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black")
      viridis_pal <- rev(c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black"))
      delta_pal <- rev(c("darkred", "orange", "gray", "dodgerblue", "darkblue"))
      
      icon <- makeAwesomeIcon("leaf", markerColor="red")
      
      rgb_raster <- reactive({
            req(final)
            f <- values(final()[[1:10]])
            a <- which(is.finite(rowSums(f)))
            v <- scale(na.omit(f))
            v[,1:5] <- v[,1:5] * input$pclim/100 * 2
            v[,6:10] <- v[,6:10] * (1-input$pclim/100) * 2
            v <- apply(prcomp(v)$x[,1:3], 2,
                       function(x) (rank(x)-1)/(length(x)-1))
            v <- v[, sample(1:3, 3)] # permute columns
            clr <- final()[[1:3]] %>% setValues(NA)
            for(i in 1:3){
                  if(sample(c(T, F), 1)){
                        clr[[i]][a] <- v[,i]
                  }else{clr[[i]][a] <- 1-v[,i]}
            }
            clr <- clr * 255 * .8 # scale to avoid white
            RGB(clr) <- 1:3
            names(clr) <- c("r", "g", "b")
            clr
      })
      
      output$map <- renderLeaflet({
            leaflet() %>%
                  setView(lng=coordinates(s)[1], lat=coordinates(s)[2], zoom=9) %>%
                  addProviderTiles(providers$Esri.WorldGrayCanvas)
      })
      
      observe({
            withProgress(
                  message = "Stand by.", 
                  value = 0, 
                  {
                        incProgress(.25, detail = "Computing climate similarities.")
                        x <- clim_sigmas()
                        incProgress(.5, detail = "Computing soil similarities.")
                        x <- soil_sigmas()
                        
                        incProgress(.75, detail = "Merging dimensions.")
                        
                        v <- all_vars[all_vars$desc == input$color, ]
                        f <- final()
                        incProgress(1, detail = "Generating plots.")
                        
                        latlon <- coordinates(site$point)
                        req(final())
                        
                        if(input$color == "Ordination"){
                              leafletProxy("map") %>% clearMarkers() %>% clearImages() %>%
                                    addRasterImage(rgb_raster(), opacity=0.8) %>%
                                    addAwesomeMarkers(lng=latlon[1], lat=latlon[2], icon=icon)
                        }else{
                              r <- f[[v$abbv]]
                              limits <- range(c(v$min, values(r)), na.rm = T)
                              if(str_detect(input$color, "Change in")) limits <- max(abs(values(r)), na.rm = T) * c(-1, 1)
                              pal <- colorNumeric(switch(v$palette, "sigma" = sigma_pal, "viridis" = viridis_pal, "delta" = delta_pal),
                                                  limits,
                                                  na.color = "transparent")
                              leafletProxy("map") %>% clearMarkers() %>% clearImages() %>%
                                    addRasterImage(r, colors=pal, opacity=0.8) %>%
                                    addAwesomeMarkers(lng=latlon[1], lat=latlon[2], icon=icon)
                        }
                        
                  })
      })
      
      output$scatter <- renderPlot({
            
            # if(str_detect(input$color, "Change in")) browser()
            
            
            # future climate mean, and sd of either ensemble or niche
            avg <- c(target_envt()$focal$clim, target_envt()$focal$soil) %>% unlist()
            ref <- c(target_envt()$reference$clim, target_envt()$focal$soil) %>% unlist()
            std <- c(range_stats()$clim_sd, range_stats()$soil_sd) %>% unlist()
            names(avg) <- names(std) <- names(ref) <- all_vars$abbv[1:length(avg)]
            
            vx <- all_vars$abbv[all_vars$desc == input$xvar]
            vy <- all_vars$abbv[all_vars$desc == input$yvar]
            vc <- all_vars$abbv[all_vars$desc == input$color]
            
            req(final())
            f <- final()
            if(vc == "rgb"){
                  d <- c(final(), rgb_raster())[[c(vx, vy, c("r", "g", "b"))]] %>% as.data.frame() %>% na.omit()
                  d$rgb <- rgb(d$r, d$g, d$b, maxColorValue = 255)
                  d <- select(d, -r, -g, -b)
            }else{
                  d <- f[[c(vx, vy, vc)]] %>% as.data.frame() %>% na.omit()
            }
            
            # compute range on full data set, and then subsample data
            minmax <- range(d[[vc]], na.rm = T)
            maxpoints <- 10000
            names(d) <- c("xvar", "yvar", "cvar")
            d <- d %>% sample_n(min(maxpoints, nrow(.)))
            
            
            # vc <- all_vars[all_vars$desc == input$color, ]
            fill_col <- NA
            
            
            
            
            x_label <- paste0(input$xvar, " ", all_vars$units[all_vars$desc == input$xvar])
            y_label <- paste(input$yvar, all_vars$units[all_vars$desc == input$yvar])
            
            
            
            # gcm uncertainty
            # gcmv <- gcm_var()
            # bubbles <- tibble(x = gcmv$clim[,vx],
            #                   y = gcmv$clim[,vy],
            #                   xsd = gcmv$clim_sd[,vx],
            #                   ysd = gcmv$clim_sd[,vy],
            #                   group = 1:length(x)) %>%
            #       expand_grid(id = 1:100,
            #                   sd = 1:2) %>%
            #       mutate(dx = cos(seq(0,2*pi,length.out=100))[id],
            #              dy = sin(seq(0,2*pi,length.out=100))[id]) %>%
            #       mutate(x = x + xsd * dx * sd,
            #              y = y + ysd * dy * sd)
            
            # plot
            req(d, vc)
            
            vci <- all_vars[all_vars$desc == input$color, ]
            scatter <- ggplot() +
                  geom_point(data = d, aes(xvar, yvar, color = cvar), size = .5) +
                  theme_minimal() +
                  theme(legend.position = "right") +
                  labs(x = x_label, y = y_label, color = NULL)
            
            if(vc == "rgb"){
                  scatter <- scatter + 
                        guides(color = guide_colorsteps()) +
                        scale_color_identity()
            }else{
                  limits <- c(ifelse(is.na(vci$min), minmax[1], vci$min), minmax[2])
                  if(str_detect(input$color, "Change in")) limits <- max(abs(minmax)) * c(-1, 1)
                  scatter <- scatter + 
                        guides(color = guide_colorbar(barheight = 12)) +
                        scale_color_gradientn(colours = switch(vci$palette, 
                                                               "sigma" = sigma_pal, "viridis" = viridis_pal, "delta" = delta_pal), 
                                              limits = limits)
            }
            
            
            
            # geom_polygon(data = bubbles, aes(x, y, group = paste(group, sd)),
            #              color = "black", fill = NA) + 
            
            # annotate("linerange", color="red", size=.5,
            #          x=avg[vx], 
            #          ymin=avg[vy] + std[vy]*c(-.5, .5), 
            #          ymax=avg[vy] + std[vy]*c(-2.5, 2.5)) +
            # annotate("segment", color="red", size=.5,
            #          y=avg[vy], yend=avg[vy], 
            #          x=avg[vx] + std[vx]*c(-.5, .5), 
            #          xend=avg[vx] + std[vx]*c(-2.5, 2.5)) +
            
            # annotate("point", color="red", shape=3, size=8,
            #          x=avg[vx], y=avg[vy]) +
            # annotate("polygon", color="red", fill=fill_col, size=.5,
            #          x=avg[vx] + std[vx] * cos(seq(0,2*pi,length.out=100)),
            #          y=avg[vy] + std[vy] * sin(seq(0,2*pi,length.out=100))) +
            # annotate("polygon", color="red", fill=fill_col, size=.5,
            #          x=avg[vx] + std[vx]*2 * cos(seq(0,2*pi,length.out=100)),
            #          y=avg[vy] + std[vy]*2 * sin(seq(0,2*pi,length.out=100))) +
            # annotate("segment", color="black",
            #          x=ref[vx], y=ref[vy], xend=avg[vx], yend=avg[vy],
            #          arrow=grid::arrow(type="closed", angle=15, length=unit(.15, "in"),
            #                            ends = ifelse(input$time == "1981-2010", "first", "last"))) +
            
            # coord_fixed(ratio = std[vx]/std[vy]) +
            
            
            if(!grepl("plot", paste(vx, vy))){
                  req(scatter)
                  scatter <- scatter +
                        annotate("segment", color="black", linewidth = 1.5,
                                 x=ref[vx], y=ref[vy], xend=avg[vx], yend=avg[vy],
                                 arrow=grid::arrow(type="closed", angle=15, length=unit(.15, "in"),
                                                   ends = ifelse(input$mode == "collection", "first", "last"))) +
                        annotate("segment", color="red", linewidth = 1,
                                 x=ref[vx], y=ref[vy], xend=avg[vx], yend=avg[vy],
                                 arrow=grid::arrow(type="closed", angle=15, length=unit(.15, "in"),
                                                   ends = ifelse(input$mode == "collection", "first", "last")))
            }
            
            return(scatter)
      })
      
      
      output$target_label <- renderText({ paste0("Target: ", input$mode, " site (",
                                                 ifelse(input$mode == "planting", 
                                                        paste0(input$time, ", ", input$ssp), 
                                                        paste0(times[1], ifelse(input$radius == 0, "", ", smoothed"))), 
                                                 # ", +/-2 range SD",
                                                 ")") })
      output$points_label <- renderText({ paste0("Points: potential ", setdiff(c("planting", "collection"), input$mode), " sites (",
                                                 ifelse(input$mode == "planting", 
                                                        paste0(times[1], ifelse(input$radius == 0, "", ", smoothed")), 
                                                        paste0(input$time, ", ", input$ssp)), ")") })
      output$color_label <- renderText({ paste("Color of map & points:", input$color,
                                               all_vars$units[all_vars$desc == input$color]) })
      
      
      
      # data download
      output$download <- downloadHandler(
            filename = function() {
                  ll <- coordinates(site$point)
                  paste0("sigma", 
                         "_", str_replace(input$sp, " ", "-"), 
                         "_", input$mode, 
                         "_pclim", input$pclim, 
                         "_r", input$radius, "km",
                         "_", input$time, 
                         "_", ssps$ssp[ssps$text == input$ssp], 
                         "_", round(ll[1], 2), "E",
                         "_", round(ll[2], 2), "N",
                         ".tif")
            },
            content = function(file) {
                  writeRaster(final()$prob, file)
            }
      )
      
}

# Run the application
shinyApp(ui = ui, server = server)
