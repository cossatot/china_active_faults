using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta
using ArchGDAL
const AG = ArchGDAL

using PyPlot

fault_file = "../block_data/chn_faults.geojson"
#gnss_vels_file = "./block_vels.csv"
# load data

# faults
fault_json = JSON.parsefile(fault_file);

data_gpkg = "../block_data/chn_faults_blocks.gpkg"
geol_slip_rates_file = "../block_data/geol_slip_rate_pts.geojson"

fault_df = Oiler.IO.gis_vec_file_to_df(data_gpkg; layername="chn_faults")
block_df = Oiler.IO.gis_vec_file_to_df(data_gpkg; layername="blocks")
gnss_df_all = Oiler.IO.gis_vec_file_to_df(data_gpkg;layername="gnss_vels")
geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(geol_slip_rates_file)

println("n blocks: ", size(block_df, 1))

fault_weight = 1.


function vel_nothing_fix(vel)
    if vel == ""
        return 0.
    elseif isnothing(vel) | ismissing(vel)
        return 0.
    else
        if typeof(vel) == String
            return parse(Float64, vel)
        else
            return vel
        end
    end
end

function err_nothing_fix(err)
    if err == ""
        return 20.
    elseif isnothing(err) | ismissing(err)
        return 20.
    else
        if typeof(err) == String
            return parse(Float64, err)
        else
            return err
        end
    end
end


function row_to_fault(row)
    trace = Oiler.IO.get_coords_from_geom(row[:geometry])

    Oiler.Fault(trace = trace,
        dip_dir = row[:dip_dir],
        extension_rate = vel_nothing_fix(row[:v_ex]),
        extension_err = err_nothing_fix(row[:e_ex]) * fault_weight,
        dextral_rate =vel_nothing_fix(row[:v_rl]),
        dextral_err = err_nothing_fix(row[:e_rl]) * fault_weight,
        dip = row[:dip],
        name = row[:name],
        hw = row[:hw],
        fw = row[:fw],
        usd=row[:usd],
        lsd=row[:lsd],
        )
end

faults = [row_to_fault(fault_df[i,:]) for i in 1:size(fault_df, 1)]

#faults = map(feat_to_Fault, fault_json["features"]);
fault_vels_ = map(Oiler.fault_to_vels, faults);
fault_vels = reduce(vcat, fault_vels_)
println("n faults: ", length(faults))

#fault_vels = [f for f in fault_vels if (f.vn != 0.) & (f.ve != 0.)]
println("n faults vels: ", length(fault_vels))


# geol slip rates
function make_vel_from_slip_rate(slip_rate_row, fault_df)
    fault_idx = parse(Int, slip_rate_row[:fault_seg])
    fault_row = @where(fault_df, findall(x->(x==fault_idx), fault_df.fid))[1,:]
    fault = row_to_fault(fault_row)

    extension_rate = vel_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err])
    dextral_rate =vel_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:extension_err])

    ve, vn = Oiler.Faults.fault_slip_rate_to_ve_vn(dextral_rate, 
                                                   extension_rate,
                                                   fault.strike)

    ee, en = Oiler.Faults.fault_slip_rate_to_ve_vn(dextral_err, 
                                                   extension_err,
                                                   fault.strike)

    pt = Oiler.IO.get_coords_from_geom(slip_rate_row[:geometry])
    lon = pt[1]
    lat = pt[2]
    
    VelocityVectorSphere(lon=lon, lat=lat, ve=ve, vn=vn, fix=fault.hw,
                         mov=fault.fw, vel_type="fault")

end


geol_slip_rate_vels = []
for i in 1:size(geol_slip_rate_df, 1)
    slip_rate_row = geol_slip_rate_df[i,:]
    if slip_rate_row[:include] == true
        push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row, 
                                                           fault_df))
    end
end

geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))

# gnss


gnss_block_idx = Oiler.IO.get_block_idx_for_points(gnss_df_all, block_df)

function gnss_vel_from_row(row, block)
    pt = Oiler.IO.get_coords_from_geom(row[:geometry])
    lon = pt[1]
    lat = pt[2]
    Oiler.VelocityVectorSphere(lon = lon, lat = lat, ve = row.e_vel,
        vn = row.n_vel, ee = row.e_err, en = row.n_err, name = row.station,
        fix = "1111", mov = string(block), vel_type="GNSS")
end

gnss_vels = []

for (i, block) in enumerate(gnss_block_idx)
    if !ismissing(block)
        gv = gnss_vel_from_row(gnss_df_all[i,:], block) 
        push!(gnss_vels, gv)
    end
end

gnss_vels = convert(Array{VelocityVectorSphere}, gnss_vels)

vels = vcat(fault_vels, gnss_vels, geol_slip_rate_vels);
# vels = vcat(fault_vels, gnss_vels);

println("n gnss vels: ", length(gnss_vels))

vel_groups = Oiler.group_vels_by_fix_mov(vels);

# solve
poles = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults = faults,
                                               sparse_lhs=true,
                                               weighted = true)

# look at outputs

rates  = Oiler.Utils.get_fault_slip_rates_from_poles(faults, poles)

new_fault_json = deepcopy(fault_json)

for (i, rate) in enumerate(rates)
    new_fault_json["features"][i]["properties"]["v_rl"] = rate[1]
    new_fault_json["features"][i]["properties"]["e_rl"] = 0.
    new_fault_json["features"][i]["properties"]["v_ex"] = rate[2]
    new_fault_json["features"][i]["properties"]["e_ex"] = 0.
end

open("../block_data/chn_faults_out.geojson", "w") do ff
    JSON.print(ff, new_fault_json)
end


#block_centroids = CSV.read("./block_centroids.csv")
block_centroids = [AG.centroid(block_df[i, :geometry]) for i in 1:size(block_df, 1)]

pole_arr = collect(values(poles))
pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]

centroids_lon = []
centroids_lat = []
centroids_ve = []
centroids_vn = []
eur_rel_poles = Array{Oiler.PoleCart}(undef, size(block_centroids,1))
for (i, b_cent) in enumerate(block_centroids)
    #try
        #row = block_centroids[i,:]
        bc_lon = AG.getx(b_cent, 0)
        bc_lat =  AG.gety(b_cent, 0)
        push!(centroids_lon, bc_lon)
        push!(centroids_lat, bc_lat)
        PvGb = Oiler.BlockRotations.build_PvGb_deg(bc_lon, bc_lat)
        
        pole = Oiler.Utils.get_path_euler_pole(pole_arr, "1111", string(block_df[i, :fid]))
        
        ve, vn, vu = PvGb * [pole.x, pole.y, pole.z]
        push!(centroids_ve, ve)
        push!(centroids_vn, vn)
        eur_rel_poles[i] = pole
    #catch
    #    row = block_centroids[i,:]
    #    fid = row.fid
    #    warn_msg = "no pole found for block $fid" 
    #    @warn warn_msg
    #end
end

centroids = DataFrame()
centroids.lon = centroids_lon
centroids.lat = centroids_lat
#centroids.fid =

CSV.write("../block_data/block_poles_eur_rel.csv", 
          Oiler.IO.poles_to_df(eur_rel_poles, convert_to_sphere=true))

vlon = [v.lon for v in gnss_vels]
        
vlat = [v.lat for v in gnss_vels]
vve = [v.ve for v in gnss_vels]
vvn = [v.vn for v in gnss_vels]

vg_keys = sort(collect(Tuple(keys(vel_groups))))
fault_groups = Oiler.Utils.group_faults(faults, vg_keys)

locking_partial_groups = Dict()
for vg in vg_keys
    if haskey(fault_groups, vg)
        locking_partial_groups[vg] = sum([
            Oiler.Elastic.calc_locking_effects_segmented_fault(fault, vlon, vlat)
            for fault in fault_groups[vg]
        ])
    #else
    #    println("Aint no fault here")
    end
end

pred_lock_vels_dict = Dict()

for vg in keys(locking_partial_groups)
    pv = [poles[vg].x, poles[vg].y, poles[vg].z]
    pred_lock_vels_dict[vg] = [part * pv for part in
    locking_partial_groups[vg]]
end

pred_lock_vels = sum(values(pred_lock_vels_dict))


ple, pln = [v[1] for v in pred_lock_vels], [v[2] for v in pred_lock_vels]


gnss_pred_vels = [
    Oiler.Utils.get_vel_vec_from_pole(vel, poles[(vel.fix, vel.mov)])
    for vel in gnss_vels
]

pred_block_ve = [v[1] for v in gnss_pred_vels]
pred_block_vn = [v[2] for v in gnss_pred_vels]

pred_gnss_ve = ple + pred_block_ve
pred_gnss_vn = pln + pred_block_vn


figure(figsize=(14,10))
for fault in faults
    plot(fault.trace[:,1], fault.trace[:,2], "k-", lw=0.3)
end
quiver(centroids_lon, centroids_lat, centroids_ve, centroids_vn, color="C1", scale=300)
#quiver(vlon, vlat, vve, vvn, color="k", scale=300)
#quiver(vlon, vlat, ple, pln, color="r", scale=300)
quiver(vlon, vlat, pred_gnss_ve, pred_gnss_vn, color="r", scale=300)
quiver(vlon, vlat, vve-pred_gnss_ve, vvn-pred_gnss_vn, color="b", scale=300)
axis("equal")

figure(figsize=(10,4))
subplot(1,2,1)
errorbar([f.dextral_rate for f in faults], [r[1] for r in rates],
         xerr = [f.dextral_err for f in faults], fmt="o")
title("dextral")
subplot(1,2,2)
errorbar([f.extension_rate for f in faults], [r[2] for r in rates],
         xerr = [f.extension_err for f in faults], fmt="o")
title("extension")

show()
