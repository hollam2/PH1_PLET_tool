import pystac_client

url = 'https://catalog.dive.edito.eu'
client = pystac_client.Client.open(url)
print(client)


collections = list(client.get_collections())
print(f"number of collections: {len(collections)}")

variable = "occurrence"

items = []
for collection in collections:
    if variable in collection.id:
        print(collection.id)
        for i, item in enumerate(collection.get_items()):
            items.append(item)

print(items)

for item in items:
    for key, value in item.assets.items():
        print(f"{key}: {value}")
        print("-"*25)
        if key == "parquet":
            occurrence_data = value.href

print("-_"*25)
print(occurrence_data)