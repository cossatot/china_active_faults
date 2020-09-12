import pandas as pd
import numpy as np
import json
from copy import deepcopy

pdf = pd.read_csv("../block_data/block_poles_eur_rel.csv")
pdf = pdf.set_index(pdf['mov'].values)

pdf['time'] = 0

gdf_10 = pdf[["mov", "time", "lat", "lon", "rotrate", "fix"]]

gdf_10["rotrate"] *= -10
gdf_10["time"] = 10.

#gdf_10.to_csv("../gplates/rots_10.rot", index=False, header=False, sep="\t")

gdf_0 = gdf_10.copy(deep=True)
gdf_0["time"] = 0.
gdf_0["rotrate"] = 0.
gdf_0["lat"] = 0.
gdf_0["lon"] = 0.

gdf_11 = gdf_0.copy(deep=True)
gdf_11["time"] = 11.


gdf_all = pd.concat((
    gdf_0, 
    gdf_10, 
    #gdf_11
))


gdf_all = gdf_all.sort_values(["mov", "time"])

gdf_all["comment"] = "!"
gdf_10["comment"] = "!"

gdf_all.to_csv("../gplates/rots_all.rot", index=False, header=False, sep="\t")
gdf_10.to_csv("../gplates/rots_10.rot", index=False, header=False, sep="\t")



def sind(ang):
    return np.sin(np.radians(ang))

def cosd(ang):
    return np.cos(np.radians(ang))

def atand2(y,x):
    return np.degrees(np.arctan2(y,x))

def atand(x):
    return np.degrees(np.arctan(x))


def sphere_to_cart(lon, lat):
    
    r = 1.
    
    x = r * cosd(lat) * cosd(lon)
    y = r * cosd(lat) * sind(lon)
    z = r * sind(lat) 
    
    return np.array([x,y,z])


def cart_to_sphere(x, y, z):
    lon = atand2(y, x)
    lat = atand(z / np.sqrt(x**2 + y**2))
    
    return [lon, lat]

def make_rot_matrix(cpole, ang):
    
    Ex, Ey, Ez = cpole[0], cpole[1], cpole[2]
    
    r11 = Ex * Ex * (1 - cosd(ang)) + cosd(ang)
    r12 = Ex * Ey * (1 - cosd(ang)) - Ez * sind(ang)
    r13 = Ex * Ez * (1 - cosd(ang)) + Ey * sind(ang)
    
    r21 = Ey * Ex * (1 - cosd(ang)) + Ez * sind(ang)
    r22 = Ey * Ey * (1 - cosd(ang)) + cosd(ang)
    r23 = Ey * Ez * (1 - cosd(ang)) - Ex * sind(ang)
    
    r31 = Ez * Ex * (1 - cosd(ang)) - Ey * sind(ang)
    r32 = Ez * Ey * (1 - cosd(ang)) + Ex * sind(ang)
    r33 = Ez * Ez * (1 - cosd(ang)) + cosd(ang)
    
    return [[r11, r12, r13],
            [r21, r22, r23],
            [r31, r32, r33]]


def rotate_block(block, time=1):
    trace = block["geometry"]["coordinates"][0]
    pole = pdf.loc[block['properties']['fid']]
    
    cpole = sphere_to_cart(pole.lon, pole.lat)
    R = make_rot_matrix(cpole, pole.rotrate * time)
    
    trace_cart = [sphere_to_cart(*tr) for tr in trace]
    
    rot_pts = [(tr @ R) for tr in trace_cart]
    
    new_trace_sphere = [cart_to_sphere(*rp) for rp in rot_pts]
    
    return new_trace_sphere




with open('../block_data/chn_blocks.geojson') as f:
    gj = json.load(f)


gj_plus_5 = deepcopy(gj)
gj_minus_5 = deepcopy(gj)

for feat in gj_plus_5["features"]:
    feat["geometry"]["coordinates"] = [rotate_block(feat, 5)]

for feat in gj_minus_5["features"]:
    feat["geometry"]["coordinates"] = [rotate_block(feat, -5)]

for feat in gj["features"]:
    feat["properties"]["end_time"] = -9999.
    feat["properties"]["start_time"] = 20.

for feat in gj_plus_5["features"]:
    feat["properties"]["end_time"] = -9999.
    feat["properties"]["start_time"] = 20.
    #feat["properties"]["start_time"] = 5.

for feat in gj_minus_5["features"]:
    feat["properties"]["TOAGE"] = -9999.
    feat["properties"]["FROMAGE"] = 20.
    feat["properties"]["PLATEID"] = feat["properties"]["fid"]
    feat["properties"]["IMPORT_AGE"] = 0.



with open('../gplates/gp_chn_blocks.geojson', 'w') as gf:
    json.dump(gj, gf)

with open('../gplates/gp_chn_blocks_plus_5.geojson', 'w') as gf:
    json.dump(gj_plus_5, gf)

with open('../gplates/gp_chn_blocks_minus_5.geojson', 'w') as gf:
    json.dump(gj_minus_5, gf)
