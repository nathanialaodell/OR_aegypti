# OR_aegypti
This repository contains all code in order to reproduce the figures for the network analyses portion of "Synthesizing genetic and network methodologies to characterize invasion dynamics of Aedes aegypti" (In preparation) and its supplementary materials.

Code in this repository has been tested up to R v4.6.0 (2026-04-24).

## Directories setup
The code in this repository has been generalized to search for all data in a folder 
titled `~/Objects` for spatial network data and `~/Objects/kinship` in the parent directory. 

## Script descriptions

### Function sources

`From_to_functions.R` - Source script to turn a trap collection dataset into a dataframe with "from" and "to" columns to facilitate graph creation using `igraph`.

| Function name | Purpose | 
| ------------- | ------- |
| `nn.fun(x,y,z)` | Split data frame into sites that are origin (from) sites and destination sites and get their nearest neighbor to calculate Euclidean distance | 
| `new.fun(from, to, x, y)` | Merge from and to datasets into one with connections | 
| `gene_dist(data)` | Get the Euclidean distance for connections in a kinship dataset | 

`Creating_networks.R` - Source script to create and plot networks. 

| Function name | Purpose | 
| ------------- | ------- |
| `graph_wk_fun(data, wk)` | Create a directed graph from a data frame, filtered to a specific week (this is a roundabout way to create a temporal network) | 
| `get_nearest_neighbors(data, wk)` | Get the nearest neighbor for each site with a collection in a given week (used to create red, dashed lines in **S. figure 5**|
| `graph_wk_with_nn_fun(data, wk)` | Combination of the previous two functions | 
| `pn_nn_fun(data, graph, full_graph, title)` | Plot a nearest neighbor graph. The `data` argument is used to create a layout so that nodes are plotted according to their geocoordinates | 
| `pn_fun(data, graph, title, kin_legend, spat_legend, lab, temp)` | Plot a kinship or spatial graph | 
| `graph_wk_it_fun(data, wk)` | Create a single, static, weighted (by Euclidean distance) and directed graph |

### Analysis scripts

> [!NOTE]
> Scripts should be ran in the order signified by their placement in this markdown (alternatively, the number suffix preceding the script name). In general, scripts are ran according to the placement of their result(s) in the manuscript.

> [!NOTE]
> There are instances where a mismatch between genetic data and spatial data address representation occurs (i.e., there is a typo in one of the datasets that is being merged). In an effort to retain the datasets' original form(s), typos were corrected as needed within function calls, as opposed to manually changing the actual spreadsheet. Many traps were placed within or near residences, so to protect the privacy of homeowners, lines within R scripts where this occurred have been removed with an applicable note of "PII... obfuscated". Additionally, scripts where this process occurred are indicated by "... (obfuscated).R". 

`0_Study area map.R` - Create a map of the study areas.

`1_Sweep.R` - Preprocessing of collection data for Ashland and Talent, Oregon.

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
| Talent connectivity NA omitted.shp | (Esri) shapefile with "from" and "to" connections between traps (see manuscript) |
| Ashland connectivity NA omitted.shp | Same as above. | 


`2_Time ordered networks.R` - Creates and calculates centrality measures for 
static and time-ordered networks for Talent and Ashland, Oregon, 
\textit{Aedes aegypti} collection data from 2024 and 2025, respectively.

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
| Ashland sum table part 1 (spatial centrality).RDS | Chunk 1 to create **S. table 3**|

`3_Var_Extract.R` - Extract variables that could have been extracted at the time of analysis. Because the PRISM R package did not have the 
ability to pull data at 800 m resolution in early 2025, a separate perl script must be run (see the PRISM folder in this repo for instructions).
Note that this is no longer the case, so, alternatively, you could just extract the same weather variables using `prism::get_prism_dailys()`.

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
|Talent with climate.rds | Temperature max, mean; mean dew point; precipitation for dates with \textit{Aedes aegypti} collections|
| Talent with landcover.rds | Landcover data for Talent (ultimately not used in analysis) | 
| Talent with roads.rds | Road density for Talent (ultimately not used in analysis) |
| Talent NA with env.rds | All three above datasets merged | 
| Talent with rel humidity.rds | Relative humidity for dates with \textit{Aedes aegypti} collections |
| Talent with wind bearing.rds | Wind bearing for collection dates (ultimately not used in analysis) | 
| Ashland Clim attributes.rds | See above | 


`3.1_Merging land and climate.R`- Merging environmental variables to Talent and Ashland, Oregon, trap collection datasets for the purposes of homophily.

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
| Objects/Talent attributes.rds | Merge all datasets from `3_Var_Extract.R` | 

`4_Homophily.R` - Calculate environmental variable homophily for Talent and Ashland
network nodes. Homophily is calculated using a pearson coefficient for climate variables 
($\alpha = 0.05$).

`5_Spread analysis plot.R` - NOT RUN. Used to help create Figure 3.

`6_Hotspots.R` - Calculate outgoing contact chains and associated hotspot index
for Talent, Oregon. Importantly, the plot calls will not show properly in the R 
plot window. Save them to `.svg` files prior to viewing. 

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
| Talent sum table part 1 (spatial centrality).rds | Chunk 1 to create \textbf{S. table 2} |

`7_Dispersal kernels.R` - Calculate dispersal kernels for close-kin dyads. 

\textbf{Written objects}:

| File name | Description? | 
| --------- | ------------ | 
| Talent sum table part 3 (dispersal).RDS | Chunk 3 to create \textbf{S. table 2}|
| Ashland sum table part 3 (dispersal).RDS | Chunk 3 to create \textbf{S. table 3}|

Note that although the chunks to create \textbf{S.table 2 & 3} are created out of order, 
the scripts used to create each chunk are independent of one another (i.e., creating
part 3 prior to part 2 is a naming artifact and nothing more).

`8_Minimum spanning (obfuscated).R` - Calculate $F_{ST}$ between groups to acquire the null and alternative genetic graphs. There is one quirk to this script that needs mentioning: there are 3 address typos in 
the `~/Objects/kinship/Oregon_ALL_2026_GP.gen` file that required a copy of `~/Objects/Talent connectivity NA omitted.shp` to be made that matched the typos.

| File name | Description? | 
| --------- | ------------ | 
| Talent sum table part 2 (genomic centrality).RDS | Chunk 2 to create \textbf{S. table 2}|
| Ashland sum table part 2 (genomic centrality).RDS | Chunk 2 to create \texxtbf{S. table 3}|

`9_Relatedness.R (obfuscated)` - Get a quick distribution of relatedness between individuals.

`10_descriptive table big big (obfuscated).R` - Combine all three parts for \textbf{S. table 2 & 3} creation and analyze (pearson) correlations across and within graph-type centrality measures. 

`11_Triads.R` - Analyze the presence of triads within a gonotrophic cycle (~14 days).

*Network figures and table outputs:*

| Origin file | Figure | 
| ----------- | ----------- |
| `0_Study area map.R` | **Figure 1** |
| `5_Hotspots.R` | **Figure 3** |
| `8_Minimum spanning.R` | **Figure 4** |
| `11_Triads.R` | **Figure 5** |
| `1_Time ordered networks.R` | **S. figure 5** |
| `6_Hotspots.R` | **S. figure 6** |
| `9_Relatedness.R` | **S. figure 7** |
| `10_descriptive table big big (obfuscated).R` | **S. figure 8** |
| NA (conceptual framework) | **Table 1** |
| `7_Dispersal kernels.R` | **Table 2** |
| `6_Hotspots.R`; `7_Dispersal kernels.R`; `8_Minimum spanning (obfuscated).R` | **S. table 2** |
| `1_Time ordered networks.R`; `7_Dispersal kernels.R`; `8_Minimum spanning (obfuscated).R` | **S. table 3** |
| `10_descriptive table big big (obfuscated).R` | **S. table 4** |
| `3_Homophily.R` | **S.table 5** |
| `8_Relatedness.R` | **S. figure 6** |