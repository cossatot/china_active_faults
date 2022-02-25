using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta
#using ArchGDAL
#const AG = ArchGDAL

plots = true

if plots == true
    using PyPlot
end

save_results = false


# load data
fault_file = "../block_data/chn_faults.geojson"
gnss_vels_file = "../block_data/gnss_vels.geojson"

# faults
data_gpkg = "../block_data/chn_faults_blocks.gpkg"
geol_slip_rates_file = "../block_data/geol_slip_rate_pts.geojson"
block_file = "../block_data/chn_blocks.geojson"

block_df = Oiler.IO.gis_vec_file_to_df(block_file)
block_df[!,:fid] = string.(block_df[!,:fid])

gnss_df_all = Oiler.IO.gis_vec_file_to_df(gnss_vels_file)

println("n blocks: ", size(block_df, 1))

fault_weight = 2.

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        fault_file,
                                                        block_df=block_df)
fault_df[!,:fid] = string.(fault_df[!,:fid])
println("n faults: ", length(faults))
println("n faults vels: ", length(fault_vels))

geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(geol_slip_rates_file)
geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                    geol_slip_rate_df,
                                                    fault_df)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))

gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df_all, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

vels = vcat(fault_vels, gnss_vels, geol_slip_rate_vels);

println("n gnss vels: ", length(gnss_vels))

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                               sparse_lhs=true,
                                               weighted=true,
                                               predict_vels=true,
                                               pred_se=true)

# look at outputs

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df)

println(results["stats_info"])


poles = results["poles"]

if save_results == true
    Oiler.IO.write_fault_results_to_gj(results, 
    "../results/chn_faults_out.geojson",
    name="China fault slip rates")
end

#block_centroids = [AG.centroid(block_df[i, :geometry]) 
#                   for i in 1:size(block_df, 1)]
#
#pole_arr = collect(values(results["poles"]))
#pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]
#
#centroids_lon = []
#centroids_lat = []
#centroids_ve = []
#centroids_vn = []
#centroids_ee = []
#centroids_en = []
#centroids_cen = []
#eur_rel_poles = Array{Oiler.PoleCart}(undef, size(block_centroids, 1))
#for (i, b_cent) in enumerate(block_centroids)
#    bc_lon = AG.getx(b_cent, 0)
#    bc_lat =  AG.gety(b_cent, 0)
#    push!(centroids_lon, bc_lon)
#    push!(centroids_lat, bc_lat)
#    # PvGb = Oiler.BlockRotations.build_PvGb_deg(bc_lon, bc_lat)
#    
#    pole = Oiler.Utils.get_path_euler_pole(pole_arr, "1111", 
#                                           string(block_df[i, :fid]))
#    
#    # ve, vn, vu = PvGb * [pole.x, pole.y, pole.z]
#    block_vel = Oiler.BlockRotations.predict_block_vel(bc_lon, bc_lat, pole)
#    push!(centroids_ve, block_vel.ve)
#    push!(centroids_vn, block_vel.vn)
#    push!(centroids_ee, block_vel.ee)
#    push!(centroids_en, block_vel.en)
#    push!(centroids_cen, block_vel.cen)
#    eur_rel_poles[i] = pole
#end
#
#centroids_pred_df = DataFrame()
#centroids_pred_df.lon = round.(centroids_lon, digits=5)
#centroids_pred_df.lat = round.(centroids_lat, digits=5)
#centroids_pred_df.ve =  round.(centroids_ve, digits=3)
#centroids_pred_df.vn =  round.(centroids_vn, digits=3)
#centroids_pred_df.ee =  round.(centroids_ee, digits=3)
#centroids_pred_df.en =  round.(centroids_en, digits=3)
#centroids_pred_df.cen = round.(centroids_cen, digits=3)


if save_results == true
    CSV.write("../results/block_vels.csv", centroids_pred_df)

    CSV.write("../results/block_poles_eur_rel.csv", 
              Oiler.IO.poles_to_df(eur_rel_poles, convert_to_sphere=true))

    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups,
                                        name="../results/chn_gnss_results.csv")
                                           
end

if plots == true
    map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)
    rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                               geol_slip_rate_vels, 
                                               fault_df, results)
    
    if save_results == true
        savefig("../../../china_faults_paper/figs/fault_rates.pdf")
    end
    
    show()
end

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")
