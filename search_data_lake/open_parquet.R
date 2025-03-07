library(arrow)
library(dplyr)
library(httr)

occurrence_data <- "https://s3.waw3-1.cloudferro.com/emodnet/emodnet_biology/12639/eurobis_gslayer_obisenv_19022025.parquet"


host <- "https://s3.waw3-1.cloudferro.com"


s3 <- S3FileSystem$create(endpoint_override = host)
s3_path <- "emodnet/emodnet_biology/12639/eurobis_gslayer_obisenv_19022025.parquet"

cat(s3_path)

cat("s3 opening\n")

# Open dataset
dataset <- open_dataset(s3_path, filesystem = s3, format = "parquet")

# Filtering the dataset
filtered_table <- dataset %>%
  filter(aphiaid == 1080, 
         latitude >= 51, latitude <= 51.5, 
         longitude >= 2.5, longitude <= 3.3)

# Convert to DataFrame
df <- as.data.frame(filtered_table)

print(head(df))
