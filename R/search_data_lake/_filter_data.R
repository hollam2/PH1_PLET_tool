source("search_data_lake/_fetch_occurrence_data.R")
source("search_data_lake/_open_parquet.R")


occ = fetch_occurrence_data()

print(occ)
my_parquet <- occ[[1]]

dataset = open_my_parquet(my_parquet)



# Filtering parameters
aphia_ID = 126417
sel_longitude = c(0, 1)
sel_latitude = c(50, 51)
start_date = "2019-01-01"
end_date = "2020-12-31"


# Apply filtering
my_selection <- eurobis |> 
  filter(aphiaidaccepted == aphia_ID,
         longitude > sel_longitude[1], 
         longitude < sel_longitude[2],
         latitude > sel_latitude[1], 
         latitude < sel_latitude[2],
         observationdate >= as.POSIXct(start_date),
         observationdate <= as.POSIXct(end_date)) |> 
  collect()

print(my_selection)