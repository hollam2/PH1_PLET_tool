library(dplyr)
library(rlang)
library(arrow)

source("search_data_lake/_fetch_occurrence_data.R")
source("search_data_lake/_open_parquet.R")
source("search_data_lake/_filter_parquet.R")

occ = fetch_occurrence_data()
print(occ)

# alternatives for occ parquet
# https://s3.waw3-1.cloudferro.com/emodnet/biology/eurobis_occurrence_data/eurobis_occurrences_geoparquet_2024-10-01.parquet
# https://s3.waw3-1.cloudferro.com/emodnet/emodnet_biology/12639/eurobis_parquet_2025-03-14.parquet


# my_parquet = occ[[1]]
my_parquet = "https://s3.waw3-1.cloudferro.com/emodnet/emodnet_biology/12639/eurobis_parquet_2025-03-14.parquet"
#my_parquet = occ[[1]]
print(my_parquet)


dataset = open_my_parquet(my_parquet)
print(dataset$schema)
column_names <- dataset$schema$names
print(column_names)



filter_params <- list(
  #datasetid = 8437,
  longitude = c(0, 1),
  latitude = c(50, 51),
  observationdate = c("2013-01-01", "2024-12-31"),
  #identifiedby = "Jonas Mortelmans"
  datasetid = 4687
)



# Apply filtering
my_selection <- filter_parquet(dataset, filter_params)
print(my_selection)


# select columns
my_subset = subset(my_selection, select=c(parameter,
                                          parameter_value,
                                          datasetid,
                                          observationdate
                                           
                                        
                                          )
                   )
















