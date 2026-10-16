# OR_aegypti
Code repository for network analyses in to-be-named OR aegypti invasion manuscript

## Summary
Code in this repository has been tested up to R v4.6.0 (2026-04-24).



## Directories setup
The code in this repository has been generalized to search for all data in a folder 
titled "~/Objects"--for spatial network data--and "~/Objects/kinship"--in the parent directory. 

## Script descriptions

### Function sources

### Analysis scripts

> [!NOTE]
> Scripts should be ran in the order signified by their placement in this markdown (alternatively, the number suffix preceding the script name). In general, scripts are ran according to the placement of their result(s) in the manuscript.

> [!NOTE]
> There are instances where a mismatch between genetic data and spatial data address representation occurs (i.e., there is a typo in one of the datasets that is being merged). In an effort to retain the datasets' original form(s), typos were corrected as needed within function calls, as opposed to manually changing the actual spreadsheet. Many traps were placed within or near residences, so to protect the privacy of homeowners, lines within R scripts where this occurred have been removed with an applicable note of "PII... obfuscated". Additionally, scripts where this process occurred are indicated by "... (obfuscated).R". 

`0_Study area map.R` - Create a map of the study areas.

`1_Sweep.R` - Preprocessing of collection data for Ashland and Talent, Oregon.

*Written objects:*

| File name | Description? | 
| --------- | ------------ | 
| Talent connectivity NA omitted.shp | (Esri) shapefile with "from" and "to" connections between traps (see manuscript) |
| Ashland connectivity NA omitted.shp | Same as above. | 


`2_Time ordered networks.R` - Creates and calculates centrality measures for 
static and time-ordered networks for Talent and Ashland, Oregon, 
**Aedes aegypti** collection data from 2024 and 2025, respectively.

*Written objects:*

| File name | Description? | 
| --------- | ------------ | 
| Talent sum table part 1 (spatial centrality).RDS | Chunk 1 to create **S. table 2**|
| Ashland sum table part 1 (spatial centrality).RDS | Chunk 1 to create **S. table 3**|

`3_Var_Extract.R` - Extract variables that could have been extracted at the time of analysis. Because the PRISM R package did not have the 
ability to pull data at 800 m resolution in early 2025, a separate perl script must be run (see the PRISM folder in this repo for instructions).
Note that this is no longer the case, so, alternatively, you could just extract the same weather variables using `prism::get_prism_dailys()`.

`3.1_Merging land and climate.R`- Merging environmental variables to Talent and Ashland, Oregon, trap collection datasets for the purposes of homophily.

`4_Homophily.R` - Calculate environmental variable homophily for Talent and Ashland
network nodes.

`5_Spread analysis plot.R` - NOT RUN. Used to help create Figure 3.

`6_Hotspots.R` - Calculate outgoing contact chains and associated hotspot index
for Talent, Oregon. Importantly, the plot calls will not show properly in the R 
plot window. Save them to `.svg` files prior to viewing. 

`7_Dispersal kernels.R` - Calculate dispersal kernels for close-kin dyads. 

*Written objects:*

| File name | Description? | 
| --------- | ------------ | 
| Talent sum table part 3 (dispersal).RDS | Chunk 3 to create **S. table 2**|
| Ashland sum table part 3 (dispersal).RDS | Chunk 3 to create **S. table 3**|

Note that although the chunks to create *S.table 2 & 3* are created out of order, 
the scripts used to create each chunk are independent of one another (i.e., creating
part 3 prior to part 2 is a naming artifact and nothing more).

`8_Minimum spanning (obfuscated).R` - There is one quirk to this script that needs mentioning: there are 3 address typos in 
the `~/Objects/kinship/Oregon_ALL_2026_GP.gen` file that required a copy of `~/Objects/Talent connectivity NA omitted.shp` to be made that matched the typos.

`9_Relatedness.R (obfuscated)` - Get a quick distribution of relatedness between individuals.

`10_descriptive table big big (obfuscated).R` - Combine all three parts for *S. table 2 & 3* creation and analyze (pearson) correlations across and within graph-type centrality measures. 

`11_Triads.R` - Analyze the presence of triads within a gonotrophic cycle (~14 days).

*Network figures and table outputs:*

| Origin file | Figure name | 
| ----------- | ----------- |
| `0_Study area map.R` | **Figure 1** |
| `5_Hotspots.R` | **Figure 3** |
| `8_Minimum spanning.R` | **Figure 4** |
| `11_Triads.R` | **Figure 5** |
| `1_Time ordered networks.R` | **S. figure 5** |
| `9_Relatedness.R` | **S. figure 6** |
| `6_Hotspots.R` | **S. figure 7** |
| `10_descriptive table big big (obfuscated).R | **S. figure 8** |
| NA (conceptual framework) | **Table 1** |
| `7_Dispersal kernels.R` | **Table 2** |
| `1_Time ordered networks.R`; `7_Dispersal kernels.R`; `8_Minimum spanning (obfuscated).R` | **S. table 2** |
| `1_Time ordered networks.R`; `7_Dispersal kernels.R`; `8_Minimum spanning (obfuscated).R` | **S. table 3** |
| `10_descriptive table big big (obfuscated).R` | **S. table 4** |
| `3_Homophily.R` | **S.table 5** |
| `8_Relatedness.R` | **S. figure 6** |