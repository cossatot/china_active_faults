using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta

using PyPlot

save_results = true


fault_file = "../block_data/chn_faults.geojson"
gnss_vels_file = "../block_data/gnss_vels.geojson"
vel_field_file = "../geod/tibet_vel_field.geojson"

bound_file = "../block_data/tibet_vel_field_bounds.geojson"

cea_block_file = "../../c_asia_blocks/block_data/c_asia_blocks.geojson"
cea_fault_file = "../../c_asia_blocks/block_data/c_asia_faults.geojson"
cea_gnss_file = "../../c_asia_blocks/gnss_data/c_asia_vels_rollins.geojson"
cea_slip_rate_file = "../../c_asia_blocks/block_data/c_asia_geol_slip_rates.geojson"


# load data
geol_slip_rates_file = "../block_data/geol_slip_rate_pts.geojson"
block_file = "../block_data/chn_blocks.geojson"

chn_block_df = Oiler.IO.gis_vec_file_to_df(block_file)
cea_block_df = Oiler.IO.gis_vec_file_to_df(cea_block_file)

chn_block_df[!,:fid] = string.(chn_block_df[!,:fid])

bound_df = Oiler.IO.gis_vec_file_to_df(bound_file)
block_df = vcat(chn_block_df,
                cea_block_df;
                cols=:union)

block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)


println("n blocks: ", size(block_df, 1))

fault_weight = 2.


fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        fault_file,
                                                        cea_fault_file,
                                                        block_df=block_df,
                                                        lsd_default=10.,
                                                        dip_adj_remainder=0.1,
                                                        e_default=5.0,
                                                        subset_in_bounds=true,
                                                        check_blocks=true)
fault_df[!,:fid] = string.(fault_df[!,:fid])
println("n faults: ", length(faults))
println("n faults vels: ", length(fault_vels))

chn_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(geol_slip_rates_file)
cea_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)
geol_slip_rate_df = vcat(chn_geol_slip_rate_df, cea_geol_slip_rate_df)


geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                    geol_slip_rate_df,
                                                    fault_df)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))


non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(
    x->Oiler.Boundaries.boundary_to_vels(x; ee=5.0, en=5.0),
    non_fault_bounds)...)




println("n non-fault-bound vels: ", length(bound_vels))

gnss_df_all = Oiler.IO.gis_vec_file_to_df(gnss_vels_file)
cea_gnss_df = Oiler.IO.gis_vec_file_to_df(cea_gnss_file)

gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df_all, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

comet_vels = Oiler.IO.make_vels_from_gnss_and_blocks(cea_gnss_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:site,
    fix="1111"
)

vel_field_df = Oiler.IO.gis_vec_file_to_df(vel_field_file)
vel_field_df[!,"station"] = string.(vel_field_df[!,:fid])


vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(vel_field_df, block_df;
    fix="1112")

vels = vcat(fault_vels, 
            gnss_vels, 
            comet_vels,
            vel_field_vels,
            geol_slip_rate_vels, 
            bound_vels,
            );

println("n gnss vels: ", length(gnss_vels))

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                               sparse_lhs=true,
                                               weighted=true,
                                               predict_vels=true,
                                               pred_se=false)

# look at outputs
Oiler.ResultsAnalysis.get_block_centroid_vels(results, block_df; fix="1111")
Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df)
Oiler.ResultsAnalysis.calculate_resid_block_strain_rates(results)

println(results["stats_info"])


poles = results["poles"]

if save_results == true
    Oiler.IO.write_fault_results_to_gj(results, 
    "../results/chn_faults_out_insar.geojson",
    name="China fault slip rates")

    Oiler.IO.write_geol_slip_rate_results_to_csv(results;
                                outfile="../results/geol_slip_rates_insar.csv")

    Oiler.IO.write_block_centroid_vels_to_csv(results;
                                        outfile="../results/block_vels_insar.csv",
                                        fix="1111")
    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups,
                                    name="../results/chn_gnss_results_insar.csv")
                                           
end

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

show()


Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")
