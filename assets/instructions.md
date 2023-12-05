#### BACKGROUND:

The core function of the Seeds of Change tool is to measure the environmental similarity between a focal site (the location of a planting site, or the location where plant material was collected) and other sites across the range of a given plant species, typically spanning two points in time. The tool integrates data on plant species distributions, past and projected future climates, and soil characteristics.

This information can be used to help inform management activities including matching the environments of provenancing and planting sites for ecological restoration projects, collecting material from diverse environments for seed banking or restoration efforts, and assisting gene flow in the face of climate change.

The basic rationale for localized environmental information being relevant to management questions like these is the assumption that plant populations are genetically adapted to the historic environment of the locations where they live. Research suggests that populations often, but by no means always, exhibit true "local adaptation" meaning they exhibit better performance in their local environments than in other environments across the species range. Plants do of course at least have a demonstrated ability to survive on some timescale in the locations where they are found, which suggests that even in the absence of proper local adaptation, understanding patterns of climate departure from historic conditions can help to inform vulnerability assessments.

The goal of estimating the historic environment that a population inhabited and may be adapted to necessarily involves defining a population. This will depend on the rates of seed and pollen movement among nearby sites, relative to the rate of natural selection within a site: populations with strong selection and limited movement may be highly adapted to their local microenvironment, whereas populations with substantial seed and pollen movement relative to local selection may instead be adapted to the average environment across a larger swath of the species range. This tool uses spatial smoothing filters to approximate different scales of local adaptation.

#### SETTING MODEL INPUTS:

1.  To get started, select a ***species***. This will load a map of the modeled range of the species.

2.  Next, select a focal site ***location*** by clicking the map or entering a longitude and latitude in the input box.

3.  Also use the dropdown to select the ***focal site activity***, i.e. whether the site is a planting or collection location.

4.  Additional options include the future ***time period*** and emission ***scenario***, the ***smoothing radius*** that defines the scale of environmental smoothing applied to the historic data as discussed above, and the relative weight assigned to ***soil versus climate*** variables. These options can be set based on the particulars of your study species and use case, or they can be iteratively changed to explore the sensitivity of outputs to these choices.

Each time a setting is changed, results are recomputed and displayed on the map and scatterplot.

Click the **[?]** symbol next to a setting for additional information.

#### INTERPRETING RESULTS:

Results for a given set of model inputs are displayed on the map and scatterplot (which has a point for every 1 km grid cell on the map). The dropdown menus below the scatterplot can be used to change which output variables are displayed on the ***x*** and ***y*** dimensions of the scatterplot, and which variable is shown in ***color*** on the map and scatterplot. See below for descriptions of the variables.

The scatterplot includes an arrow indicating projected changes in the focal site's climate over time.

-   ***Sigma***: This is a measure of how different a location's environment is from the target site. This is the default variable represented by color. Sigma is a multidimensional z-score, representing the equivalent to how many one-dimensional standard deviations away from the focal site's environment a location is. The standard deviation is defined based on historic spatial variation across the species range. Users can view a combined version of sigma representing the mix of climate and soil defined on the input slider, or can view versions based solely on climate or soil.\
    ![illustration of sigma statistic](images/sigma-01.jpg){alt="illustration of sigma statistic" width="600"}

-   ***Climate variables***: The tool is based on 5 climate variables that are broadly relevant to the distribution and function of plants in California. These can be used as scatterplot axes as well as color variables, helping users to explore patterns of environmental variation across species ranges. The climate data used here are from [CHELSA](https://chelsa-climate.org/).

-   ***Soil principal components***: Soil similarity is based on a large number of chemical and physical soil variables relevant to plants, reduced to 5 uncorrelated dimensions via a principal component analysis. These principal components (PC) can be used as scatterplot dimensions and as color variables. Documentation will be added providing further detail about the particular soil properties associated with each PC. The soil data used here are from [SoilGrids](https://www.isric.org/explore/soilgrids).

-   ***Ordination***: This feature helps to visualize patterns of environmental variation across a species' range by using similar colors for sites with similar environments. This variable is not specific to the focal site and is based only on historic climate and soil (weighted based on the soil-climate slider). [NOTE: at last check this feature was not working properly on the web server currently hosting the application.]

-   ***Change in climate***: For each climate variable, users can also select the projected amount of climate change as a color variable. This represents the mean predicted change for an ensemble of 5 GCMs for the selected time period and SSP, relative to the baseline time period. No spatial smoothing is applied to these variables.

#### EXAMPLES:

-   ***Identifying potential collection sites***: The default settings that load when the app is initialized give an example geared toward identifying potential seed collection sites for coast live oak (Quercus agrifolia) to be planted at a focal site in San Francisco. The range model for this species estimates that it is relatively likely to be present in sites colored on the map. Sites colored in red (e.g. inland north of Santa Cruz and on the Big Sur coast) have (moderately-smoothed) historic environments similar to the planting site's projected late-21st century (2071-2100) environment under a high-emissions scenario (SSP5-8.5), potentially making them suitable as seed collection sites. The arrow on the scatterplot indicates that both climatic water deficit (CWD) and actual evapotranspiration (AET) are projected to increase at the focal site.

-   Identifying potential planting sites: [Forthcoming.]

#### CAVEATS:

All outputs are model estimates with inherent uncertainty, and should be used with the knowledge that they may not reflect true plant performance patterns.

-   ***Environmental data are uncertain.*** Climate and soil data for the historic and future time periods will not perfectly reflect conditions on the ground, due to a combination of factors including uncertainty in the data interpolation process, scale differences between the1 km resolution data and the microenvironment experienced by a plant, and uncertainties about future climate change. In addition, while the particular variables used here were selected for their general relevance to California plant ecology, they may not be the optimal variable set for any individual species.

-   ***Species ranges are only estimates.*** The species range maps reflect best efforts to model species distributions by associating known species occurrence localities with factors like climate and landscape intactness, true species distributions are not known.

-   ***Adaptive patterns are unknown.*** While the 'smoothing' setting in this tool is designed to let users explore potential patterns of local adaptation based on different assumed spatial scales of homogenizing gene flow, true adaptive genetic patterns are far more complex and are unknown for most species.
