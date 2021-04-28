using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta
using ArchGDAL
const AG = ArchGDAL

using PyPlot

fault_file = "../block_data/chn_faults.geojson"
#gnss_vels_file = "../block_data/gnss_vels.csv"
gnss_vels_file = "../block_data/gnss_vels.geojson"
# load data

# faults
fault_json = JSON.parsefile(fault_file);

data_gpkg = "../block_data/chn_faults_blocks.gpkg"
geol_slip_rates_file = "../block_data/geol_slip_rate_pts.geojson"
block_file = "../block_data/chn_blocks.geojson"

#fault_df = Oiler.IO.gis_vec_file_to_df(data_gpkg; layername="chn_faults")
fault_df = Oiler.IO.gis_vec_file_to_df(fault_file)
#block_df = Oiler.IO.gis_vec_file_to_df(data_gpkg; layername="blocks")
block_df = Oiler.IO.gis_vec_file_to_df(block_file)

gnss_df_all = Oiler.IO.gis_vec_file_to_df(gnss_vels_file)
#gnss_df_all = DataFrame(CSV.File(gnss_vels_file))
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
        return 10.
    elseif isnothing(err) | ismissing(err)
        return 10.
    else
        if typeof(err) == String
            err = parse(Float64, err)
            if iszero(err)
                return 10.
            else
                return err
            end
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
#fault_vels_ = map(Oiler.fault_to_vel, faults);
fault_vels = reduce(vcat, fault_vels_)
println("n faults: ", length(faults))

#fault_vels = [f for f in fault_vels if (f.vn != 0.) & (f.ve != 0.)]
println("n faults vels: ", length(fault_vels))


# geol slip rates
function make_vel_from_slip_rate(slip_rate_row, fault_df)
    fault_seg = slip_rate_row[:fault_seg]
    fault_idx = parse(Int, fault_seg)
    fault_row = @where(fault_df, findall(x->(x==fault_idx), fault_df.fid))[1,:]
    fault = row_to_fault(fault_row)

    extension_rate = vel_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err])
    dextral_rate =vel_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:dextral_err])

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
                         mov=fault.fw, vel_type="fault", name=fault_seg)

end


geol_slip_rate_vels = []
for i in 1:size(geol_slip_rate_df, 1)
    slip_rate_row = geol_slip_rate_df[i,:]
    if (slip_rate_row[:include] == true) | (slip_rate_row[:include] == "1")
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
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults = faults,
                                               sparse_lhs=true,
                                               weighted = true,
                                               predict_vels=true,
                                               pred_se=true)

# look at outputs

poles = results["poles"]

new_fault_json = deepcopy(fault_json)

for (i, rate) in enumerate(results["predicted_slip_rates"])
    new_fault_json["features"][i]["properties"]["v_rl"] = round.(rate.dextral_rate, digits=3)
    new_fault_json["features"][i]["properties"]["e_rl"] = round.(rate.dextral_err, digits=3)
    new_fault_json["features"][i]["properties"]["v_ex"] = round.(rate.extension_rate, digits=3)
    new_fault_json["features"][i]["properties"]["e_ex"] = round.(rate.extension_err, digits=3)
end

open("../block_data/chn_faults_out.geojson", "w") do ff
    JSON.print(ff, new_fault_json)
end


block_centroids = [AG.centroid(block_df[i, :geometry]) 
                   for i in 1:size(block_df, 1)]

pole_arr = collect(values(results["poles"]))
pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]

centroids_lon = []
centroids_lat = []
centroids_ve = []
centroids_vn = []
centroids_ee = []
centroids_en = []
centroids_cen = []
eur_rel_poles = Array{Oiler.PoleCart}(undef, size(block_centroids,1))
for (i, b_cent) in enumerate(block_centroids)
    bc_lon = AG.getx(b_cent, 0)
    bc_lat =  AG.gety(b_cent, 0)
    push!(centroids_lon, bc_lon)
    push!(centroids_lat, bc_lat)
    #PvGb = Oiler.BlockRotations.build_PvGb_deg(bc_lon, bc_lat)
    
    pole = Oiler.Utils.get_path_euler_pole(pole_arr, "1111", 
                                           string(block_df[i, :fid]))
    
    #ve, vn, vu = PvGb * [pole.x, pole.y, pole.z]
    block_vel = Oiler.BlockRotations.predict_block_vel(bc_lon, bc_lat, pole)
    push!(centroids_ve, block_vel.ve)
    push!(centroids_vn, block_vel.vn)
    push!(centroids_ee, block_vel.ee)
    push!(centroids_en, block_vel.en)
    push!(centroids_cen, block_vel.cen)
    eur_rel_poles[i] = pole
end

centroids_pred_df = DataFrame()
centroids_pred_df.lon = round.(centroids_lon, digits=5)
centroids_pred_df.lat = round.(centroids_lat, digits=5)
centroids_pred_df.ve =  round.(centroids_ve, digits=3)
centroids_pred_df.vn =  round.(centroids_vn, digits=3)
centroids_pred_df.ee =  round.(centroids_ee, digits=3)
centroids_pred_df.en =  round.(centroids_en, digits=3)
centroids_pred_df.cen = round.(centroids_cen, digits=3)

CSV.write("../block_data/block_vels.csv", centroids_pred_df)

CSV.write("../block_data/block_poles_eur_rel.csv", 
          Oiler.IO.poles_to_df(eur_rel_poles, convert_to_sphere=true))

pred_geol_slip_rates = []
for (i, rate) in enumerate(geol_slip_rate_vels)
    fault_idx = parse(Int, rate.name)
    fault_row = @where(fault_df, findall(x->(x==fault_idx), fault_df.fid))[1,:]
    fault = row_to_fault(fault_row)
    
    if haskey(poles, (rate.fix, rate.mov))
        pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault, 
                        poles[(rate.fix, rate.mov)])
    else
        pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault,
                        poles[(rate.mov, rate.fix)])
    end
    
    push!(pred_geol_slip_rates, pred_rate)
end


obs_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(vel_groups)]
pred_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(results["predicted_vels"])]

vlon, vlat = Oiler.Utils.get_coords_from_vel_array(obs_vels)
ove = [v.ve for v in obs_vels]
ovn = [v.vn for v in obs_vels]

pve = [v.ve for v in pred_vels]
pvn = [v.vn for v in pred_vels]

rve, rvn = ove - pve, ovn - pvn

pred_gnss_df = DataFrame()
pred_gnss_df.lon = vlon
pred_gnss_df.lat = vlat
pred_gnss_df.ve =  round.(pve, digits=3)
pred_gnss_df.vn =  round.(pvn, digits=3) 
pred_gnss_df.ee =  round.([v.ee for v in pred_vels], digits=3)
pred_gnss_df.en =  round.([v.en for v in pred_vels], digits=3)
pred_gnss_df.cen = round.([v.cen for v in pred_vels], digits=3)
pred_gnss_df.re =  round.(rve, digits=3)
pred_gnss_df.rn =  round.(rvn, digits=3)
pred_gnss_df.ree = round.(sqrt.(pred_gnss_df.ee.^2 + [v.ee^2 for v in obs_vels]), digits=3)
pred_gnss_df.ren = round.(sqrt.(pred_gnss_df.en.^2 + [v.en^2 for v in obs_vels]), digits=3)

CSV.write("../block_data/pred_gnss.csv", pred_gnss_df)


inc = map(!, iszero.(geol_slip_rate_df[!,:include]))

function check_missing(val)
    if val == ""
        return false
    elseif val == 0.
        return false
    elseif ismissing(val) | isnothing(val)
        return false
    else
        return true
    end
end

dex_inc = [check_missing(d) for d in geol_slip_rate_df[!,:dextral_rate]]
ext_inc = [check_missing(d) for d in geol_slip_rate_df[!,:extension_rate]]

dex_geol_obs = geol_slip_rate_df[!,:dextral_rate][dex_inc .* inc]
dex_geol_err = parse.(Float64, geol_slip_rate_df[!,:dextral_err][dex_inc .* inc])
dex_geol_pred = [p[1] for p in pred_geol_slip_rates[dex_inc[inc]]]
dex_geol_pred_err = [p[3] for p in pred_geol_slip_rates[dex_inc[inc]]]

ext_geol_obs = parse.(Float64, geol_slip_rate_df[!,:extension_rate][ext_inc .* inc])
ext_geol_err = parse.(Float64, geol_slip_rate_df[!,:extension_err][ext_inc .* inc])
ext_geol_pred =  [p[2] for p in pred_geol_slip_rates[ext_inc[inc]]]
ext_geol_pred_err = [p[4] for p in pred_geol_slip_rates[ext_inc[inc]]]

all_obs = vcat(ove, ovn, dex_geol_obs, ext_geol_obs)
all_errs = vcat([v.ee for v in obs_vels], [v.en for v in obs_vels], 
                dex_geol_err, ext_geol_err)

all_pred = vcat(pve, pvn, dex_geol_pred, ext_geol_pred)

chi_sq = sum(((all_obs .- all_pred).^2 ) ./ all_errs.^2 )

n_param = length(keys(vel_groups)) / 3.

red_chi_sq = chi_sq / (length(all_obs) - n_param)


println("Reduced Chi Sq: $red_chi_sq")


figure(figsize=(14,10))
for fault in faults
    plot(fault.trace[:,1], fault.trace[:,2], "k-", lw=0.3)
end
quiver(centroids_lon, centroids_lat, centroids_ve, centroids_vn, color="C1", scale=300)
quiver(vlon, vlat, pve, pvn, color="r", scale=300)
quiver(vlon, vlat, rve, rvn, color="b", scale=300)
axis("equal")



figure(figsize=(5,9))
suptitle("Observed vs. Modeled Quaternary Slip Rates")

subplot(2,1,1)
data_max = maximum([maximum(dex_geol_obs), maximum(dex_geol_pred)])
data_min = minimum([minimum(dex_geol_obs), minimum(dex_geol_pred)])

plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5)

axis("equal")
errorbar(dex_geol_obs, dex_geol_pred, xerr =  dex_geol_err, yerr = dex_geol_pred_err,
         fmt=",", elinewidth=0.3)

autoscale(false)

fill_between([data_min, 0., -data_min], 
             [data_min, 0., data_min], 
             [data_min, data_min, data_min],
             color="grey",
             lw=0.,
             alpha=0.1)

fill_between([data_min, 0., -data_min], 
             [-data_min, -data_min, -data_min],
             [-data_min, 0., -data_min], 
             color="grey",
             lw=0.,
             alpha=0.1)

xlabel("observed")
ylabel("modeled")
title("dextral")

subplot(2,1,2)

data_max = maximum([maximum(ext_geol_obs), maximum(ext_geol_pred)])
data_min = minimum([minimum(ext_geol_obs), minimum(ext_geol_pred)])
plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5)

axis("equal")
errorbar(ext_geol_obs, ext_geol_pred, xerr=ext_geol_err, yerr=ext_geol_pred_err, 
         fmt=",", elinewidth=0.3)

autoscale(false)

fill_between([data_min, 0., -data_min], 
             [data_min, 0., data_min], 
             [data_min, data_min, data_min],
             color="grey",
             lw=0.,
             alpha=0.1)

fill_between([data_min, 0., -data_min], 
             [-data_min, -data_min, -data_min],
             [-data_min, 0., -data_min], 
             color="grey",
             lw=0.,
             alpha=0.1)


xlabel("observed")
ylabel("modeled")
title("extension")

tight_layout()

savefig("../../../china_faults_paper/figs/fault_rates.pdf")

show()
