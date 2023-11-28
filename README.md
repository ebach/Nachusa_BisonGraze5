# Nachusa_BisonGraze5
Analytical code for Chakravorty et al. (submitted) comparing plant community diversity with and without bison grazing at Nachusa Grasslands (IL).

Authors: Jennifer Chakravorty, John Herrington, Elizabeth Bach (elizabeth.bach@tnc.org)

The R code in this repository analyzes plant community composition with and without grazing in the first 5 years since bison (Bison bison) were reintroduced at The Nature Conservancy's Nachusa Grasslands. Nachusa Grasslands is located near Franklin Grove, IL (41.883706°N, -89.342410°W) and is owned and managed by The Nature Conservancy (TNC). Nachusa consists of more than 1,600 ha of native prairie remnants, planted restored of prairies, savannas and woodlands, and wetlands. In 2014, American bison (Bison bison) were reintroduced to Nachusa Grasslands in a phased approach. An initial 30 animals were released onto ~200 ha (the north unit) in fall 2014 and grazed exclusively in the north unit for the 2015 growing season. In addition to 20 calves born in spring 2015, an additional 20 adult animals joined the herd in fall 2015 and the grazing area was expanded to include the 400-ha south unit. Since winter 2015/16, the herd has had full access to the entire 600 ha grazing area. After 2018, the herd has been managed at approximately 100 animals overwinter, with new calves increasing its size each growing season.

To examine the impact of bison grazing on plant communities, twenty-two 10 m x 18 m bison fenced in exclosures were constructed and installed prior to bison reintroduction. The plots are stratified to represent the habitat types at Nachusa Grasslands. All treatment areas, but one, were included in areas that had been burned within 2 years of sampling. 
Where possible, locations of exclosures were chosen randomly within each community type. In cases where contiguous hectares of a community were sparse, as was the case for some of the exclosures located in remnant communities, exclosures had to be placed where adequate space was available. For each exclosure, three 12.5-meter permanent transects inside the fence (ungrazed) are paired with three transects of the same length outside the fence (grazed). Transect lines are located at least 2 m from fence lines and each transect is 3 m from the adjacent transect. Quarter-square-meter quadrats are placed on the east side of each transect line at five regular intervals – 0-0.5 m, 3-3.5 m, 6-6.5 m, 9-9.5 m, and 12-12.5 m – for a total of fifteen quadrats in each grazed and ungrazed plot. 

Please see Chakravorty et al. (submitted) for full experimental details, including maps of the preserve and experimental design schema.

Code in this repository was built and finalized in RStudio(version: 2023.09.1) running R 4.3.1.

NachusaBisonPlots_Analysis_final.R: This file applies linear mixed-effects models to compare plant richness, Shannon's diversity, grass:forb ratio, non-native:native ratio, and visual obstruction readings (VOR). It also produces the graphs showing these results (Figures 2, 3, and S2 in Chakravorty et al.). It uses packages lme4 (v. 1.1-34) and tidyverse (v. 2.0.0).

NachusaBisonExclosure_NMDS_final.R: This file applies multivariate statistics to visualize and analyze plant community data (Figures 4 and S3 in Chakravorty). It uses packages tidyverse (v. 2.0.0), vegan (v. 2.6-4), and ggplot2 (v. 3.4.3).
