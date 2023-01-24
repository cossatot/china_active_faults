using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta
using PyPlot

save_results = false

fault_weight = 2.

fault_file = "../block_data/chn_faults.geojson"
gnss_vels_file = "../block_data/gnss_vels.geojson"
tibet_vel_field_file = "../geod/tibet_vel_field_2021_11_05.geojson"
# load data

# faults
geol_slip_rates_file = "../block_data/geol_slip_rate_pts.geojson"
block_file = "../block_data/chn_blocks.geojson"

fault_df = Oiler.IO.gis_vec_file_to_df(fault_file)
block_df = Oiler.IO.gis_vec_file_to_df(block_file)


println("n blocks: ", size(block_df, 1))


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

tib_vel_df = Oiler.IO.gis_vec_file_to_df(tibet_vel_field_file)
tib_vel_df[!, :fid] = string.(tib_vel_df.fid)

@info " doing comet insar vels"
@time tib_vels = Oiler.IO.make_vels_from_gnss_and_blocks(tib_vel_df, block_df;
    fix="1111", epsg=102016, name=:fid)

vels = vcat(fault_vels, 
            gnss_vels, 
            #tib_vels, 
            geol_slip_rate_vels);

println("n gnss vels: ", length(gnss_vels))

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                               sparse_lhs=true,
                                               weighted=true,
                                               predict_vels=true,
                                               pred_se=true)

# look at outputs
Oiler.ResultsAnalysis.get_block_centroid_vels(results, block_df; fix="1111")
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

    Oiler.IO.write_block_centroid_vels_to_csv(results;
                                        outfile="../results/block_vels_velfield.csv",
                                        fix="1111")
    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups,
                                        name="../results/chn_gnss_results_velfield.csv")
                                           
end

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

if save_results == true
    savefig("../../../china_faults_paper/figs/fault_rates_velfield.pdf")
end

show()

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")
