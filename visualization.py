import pandas as pd
import geopandas as gpd
import folium
from pyproj import Transformer
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def read_urban_dataset(root_dir, dataset_name):
    # Load original dataset file
    # TODO: We need to unify the file_name, the type of objects saved in each file,
    #  and the column name of the geometry in the original data file and the gdf objects.
    # region (some of the regions could be 3D)
    if dataset_name == "NewYork":
        regions_file = f"{root_dir}/{dataset_name}/regions.pkl"
        regions_pd = pd.DataFrame(pd.read_pickle(regions_file))
        regions_pd = regions_pd.rename(columns={"shape": "geometry"})
        # TODO: Even if we claim geometry is this region_geometry_marker (e.g., "shape" for NewYork,
        #  the created gdf object still doesn't have "geometry" column but only "shape" column)
        regions_gdf = gpd.GeoDataFrame(regions_pd, geometry=r"geometry")
        regions_gdf["region_id"] = regions_gdf.index
    elif dataset_name == "Singapore":
        regions_file = f"{root_dir}/{dataset_name}/regions_updated.pkl"
        regions_gdf = pd.read_pickle(regions_file)
        regions_gdf["region_id"] = regions_gdf.index
    # buildings
    if dataset_name == "NewYork":
        buildings_file = f"{root_dir}/{dataset_name}/buildings.pkl"
        buildings_gdf = pd.read_pickle(buildings_file)
        buildings_gdf["building_id"] = buildings_gdf.index
    elif dataset_name == "Singapore":
        buildings_file = f"{root_dir}/{dataset_name}/building.pkl"
        buildings_pd = pd.DataFrame(pd.read_pickle(buildings_file))
        buildings_pd = buildings_pd.rename(columns={"shape": "geometry"})
        buildings_gdf = gpd.GeoDataFrame(buildings_pd, geometry="geometry")
        buildings_gdf["building_id"] = buildings_gdf.index
    # poi
    poi_file = f"{root_dir}/{dataset_name}/poi_updated.pkl"
    poi_gdf = pd.read_pickle(poi_file)
    poi_gdf["poi_id"] = poi_gdf.index
    # roads
    roads_file = f"{root_dir}/{dataset_name}/roads.pkl"
    if dataset_name == "NewYork":
        roads_gdf = pd.read_pickle(roads_file)
    elif dataset_name == "Singapore":
        roads_pd = pd.read_pickle(roads_file)
        #roads_pd = roads_pd.rename(columns={"shape": "geometry"})
        #TODO: For Singapore, its road.pkl contains two types of geometries, "geometry" is the GPS coordinates system, while "geometry" is its own coordinates system
        roads_gdf = gpd.GeoDataFrame(roads_pd, geometry="shape")
    roads_gdf["road_id"] = roads_gdf.index

    return regions_gdf, buildings_gdf, poi_gdf, roads_gdf


def get_transformer(city_name):
    if city_name == "Singapore":
        transformer = Transformer.from_crs("EPSG:6927", "EPSG:4326", always_xy=True)
    elif city_name == "NewYork":
        transformer = Transformer.from_crs("EPSG:2263", "EPSG:4326", always_xy=True)

    return transformer


def get_boundary(geometries, geometry_type):
    if geometry_type == "point":
        geometries = np.array(geometries)
        min_lat, max_lat = geometries[:, 0].min(), geometries[:, 0].max()
        min_lon, max_lon = geometries[:, 1].min(), geometries[:, 1].max()
    elif geometry_type == "polyline" or geometry_type == "polygon":
        points = np.array([point for polyline in geometries for point in polyline])
        min_lat, max_lat = points[:, 0].min(), points[:, 0].max()
        min_lon, max_lon = points[:, 1].min(), points[:, 1].max()

    return min_lat, max_lat, min_lon, max_lon


def create_point_map(points, city_name, color):
    min_lat, max_lat, min_lon, max_lon = get_boundary(points, geometry_type="point")
    city_center = [(min_lat + max_lat) / 2, (min_lon + max_lon) / 2]
    print(city_center)
    m = folium.Map(location=city_center, zoom_start=13, tiles='cartodbpositron')
    #m.fit_bounds([min_lat, min_lon], [max_lat, max_lon])
    # Add points
    for geometry in points:
        folium.CircleMarker(
            location=[geometry[0], geometry[1]],
            radius=5,
            color=color,
            fill=True,
            fill_color=color,
            fill_opacity=0.7,
            opacity=0.8  # Semi-transparent for visibility
        ).add_to(m)
    m.save(f"{city_name}_points.html")
    return m


def create_polyline_map(polylines, city_name, color):
    min_lat, max_lat, min_lon, max_lon = get_boundary(polylines, geometry_type="polyline")
    city_center = [(min_lat + max_lat) / 2, (min_lon + max_lon) / 2]
    m = folium.Map(location=city_center, zoom_start=13, tiles='cartodbpositron')  # Better contrast
    #m.fit_bounds([min_lat, min_lon], [max_lat, max_lon])
    for geometry in polylines:
        folium.PolyLine(
            geometry,
            color=color,
            weight=1,  # Thicker line
            opacity=0.8  # Semi-transparent for visibility
        ).add_to(m)

    m.save(f"{city_name}_polylines.html")
    return m


def create_polygon_map(polygons, city_name, color):
    min_lat, max_lat, min_lon, max_lon = get_boundary(polygons, geometry_type="polygon")
    city_center = [(min_lat + max_lat) / 2, (min_lon + max_lon) / 2]
    m = folium.Map(location=city_center, zoom_start=13, tiles='cartodbpositron')  # Better contrast
    #m.fit_bounds([min_lat, min_lon], [max_lat, max_lon])

    for geometry in polygons:
        folium.Polygon(
            geometry,
            color=color,
            weight=6,  # Thicker line
            opacity=0.8,  # Semi-transparent for visibility
            fill=True,
            fill_color=color,
            fill_opacity=0.6,
        ).add_to(m)
    m.save(f"{city_name}_polygons.html")
    return m


def create_region_map(city_name, polygons, colors):
    # Plot
    fig, ax = plt.subplots(figsize=(6, 6))
    for polygon, color in zip(polygons, colors):
        polygon_patch = patches.Polygon(polygon, closed=True, edgecolor=color, facecolor=color, alpha=0.6)
        ax.add_patch(polygon_patch)

    # Adjust Limits
    all_x = [x for poly in polygons for x, _ in poly]
    all_y = [y for poly in polygons for _, y in poly]
    ax.set_xlim(min(all_x) - 0.1, max(all_x) + 0.1)
    ax.set_ylim(min(all_y) - 0.1, max(all_y) + 0.1)
    ax.set_aspect('equal', adjustable='box')

    unique_colors = list(set(colors))

    # Legend
    legend_patches = [patches.Patch(color=color, label=f'Class {label}') for label, color in enumerate(unique_colors)]
    plt.legend(handles=legend_patches, loc='upper right')

    plt.title(f"{city_name}_RegionDCL")
    plt.show()


def dataset_visualization(gdf, geometry_type, city_name, color="red"):
    transformer = get_transformer(city_name)
    # Generate maps with clear colors
    if geometry_type == "point":
        geometries = gdf['geometry'].apply(lambda geom: list(geom.coords)).tolist()
        swapped_geometries = []
        for geometry in geometries:
            projected_geometry = transformer.transform(geometry[0][0], geometry[0][1])
            swapped_geometries.append([projected_geometry[1], projected_geometry[0]])
        map = create_point_map(swapped_geometries, city_name, color)
    elif geometry_type == "polyline":
        if city_name == "NewYork":
            geometries = gdf['geometry'].apply(lambda geom: list(geom.coords)).tolist()
            swapped_geometries = []
            for geometry in geometries:
                projected_polyline = [transformer.transform(coords[0], coords[1]) for coords in geometry]
                swapped_polyline = [[coords[1], coords[0]] for coords in projected_polyline]
                swapped_geometries.append(swapped_polyline)
        elif city_name == "Singapore":
            geometries = gdf['geometry'].apply(lambda geom: list(geom.coords)).tolist()
            swapped_geometries = []
            for geometry in geometries:
                projected_polyline = geometry
                swapped_polyline = [[coords[1], coords[0]] for coords in projected_polyline]
                swapped_geometries.append(swapped_polyline)
        map = create_polyline_map(swapped_geometries, city_name, color)
    elif geometry_type == "polygon":
        geometries = gdf['geometry'].apply(lambda geom: list(geom.exterior.coords)).tolist()
        geometries = random.sample(geometries, 100000)
        swapped_geometries = []
        for geometry in geometries:
            projected_polygon = [transformer.transform(coords[0], coords[1]) for coords in geometry]
            swapped_polygon = [[coords[1], coords[0]] for coords in projected_polygon]
            swapped_geometries.append(swapped_polygon)
        map = create_polygon_map(swapped_geometries, city_name, color)


def visualize_region(city_name, regions_gdf, buildings_gdf, labels, colors):
    assert regions_gdf.shape[0] == len(labels)
    assert len(set(labels)) == len(colors)
    label2color = {}
    label_unique = list(set(labels))
    for idx, label in enumerate(label_unique):
        label2color[label] = colors[idx]

    transformer = get_transformer(city_name)
    building_indices_list = regions_gdf['buildings'].tolist()
    geometries = []
    color_list = []

    for region_idx, building_indidces in enumerate(building_indices_list):
        for building_index in building_indidces:
            geometries.append(buildings_gdf.loc[building_index, "geometry"].exterior.coords)
            color_list.append(label2color[labels[region_idx]])
    swapped_geometries = []
    for geometry in geometries:
        projected_polygon = [transformer.transform(coords[0], coords[1]) for coords in geometry]
        swapped_polygon = [[coords[0], coords[1]] for coords in projected_polygon]
        swapped_geometries.append(swapped_polygon)
    create_region_map(city_name, swapped_geometries, color_list)


if __name__ == "__main__":
    city_color = {"NewYork": "green", "Singapore": "blue"}
    for dataset_name in ["NewYork", "Singapore"]:
        print(dataset_name)
        root_dir = "/Users/jiali/Desktop/Poly2Vec/submission/submission_ICML2025/data"
        regions_gdf, buildings_gdf, poi_gdf, roads_gdf = read_urban_dataset(root_dir, dataset_name)
        #dataset_visualization(poi_gdf, "point", city_name=dataset_name, color=city_color[dataset_name])
        #dataset_visualization(roads_gdf, "polyline", city_name=dataset_name, color=city_color[dataset_name])
        #dataset_visualization(buildings_gdf, "polygon", city_name=dataset_name, color=city_color[dataset_name])
    for dataset_name in ["Singapore"]:
        labels = random.choices([0, 1, 2, 3, 4], k=regions_gdf.shape[0])
        colors = ["green",  "blue", "red", "pink", "grey"]
        full_regions_gdf = pd.read_pickle(f"{root_dir}/{dataset_name}/regions_with_pois_buildings_roads.pkl")
        visualize_region(dataset_name, full_regions_gdf, buildings_gdf, labels, colors)





