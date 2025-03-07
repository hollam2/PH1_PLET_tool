import pyarrow.parquet as pq
import pyarrow.fs
import contextily
import pyarrow.dataset as ds
import pyarrow.compute as pc
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
from shapely.geometry import Point
from urllib.parse import urlparse

occurrence_data = ("https://s3.waw3-1.cloudferro.com/emodnet/emodnet_biology/"
                   "12639/eurobis_gslayer_obisenv_19022025.parquet")

parsed_url = urlparse(occurrence_data)
host = parsed_url.hostname
bucket_name = parsed_url.path.split('/')[1]
key = '/'.join(parsed_url.path.split('/')[2:])

print('s3 forming')
s3 = pyarrow.fs.S3FileSystem(endpoint_override=host)
s3_path = f"{bucket_name}/{key}"

print('s3 opening')
dataset = ds.dataset(s3_path, filesystem=s3, format="parquet")

filtered_table = dataset.to_table(
    filter=(
        (pc.field("aphiaid") == 1080) &
        # (pc.field("identifiedby") == "Jonas Mortelmans") &
        (pc.field("latitude") >= 51) &
        (pc.field("latitude") <= 51.5) &
        (pc.field("longitude") >= 2.5) &
        (pc.field("longitude") <= 3.3)
    )
)

df = filtered_table.to_pandas()
# df.to_csv("export_occurrences_Jonas_Mortelmans.csv")

print(f"{len(df)=}")

gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df.longitude, df.latitude),
    crs="EPSG:4326"  # WGS84 coordinate system
)
fig, ax = plt.subplots(figsize=(10, 10))
gdf = gdf.to_crs(epsg=3857)  # Reproject to Web Mercator for compatibility with contextily
gdf.plot(ax=ax, color="blue", markersize=10, alpha=0.6, label="Occurrences")
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, zoom=10)
ax.set_title("Map of Filtered AphiaID Occurrences with Background Map")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.legend()
plt.show()

