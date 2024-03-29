using Revise

using Oiler

using JSON
using CSV
using DataFrames, DataFramesMeta

using PyPlot

save_results = false



fault_file = "../block_data/chn_faults.geojson"
gnss_vels_file = "../geod_data/gnss_vels.geojson"

# load data
geol_slip_rates_file = "../geol_data/geol_slip_rate_pts.geojson"
block_file = "../block_data/chn_blocks_mod.geojson"

block_df = Oiler.IO.gis_vec_file_to_df(block_file)
block_df[!,:fid] = string.(block_df[!,:fid])


println("n blocks: ", size(block_df, 1))

fault_weight = 2.


fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        fault_file,
                                                        block_df=block_df,
                                                        lsd_default=10.,
                                                        adjust_err_by_dip=false,
                                                        e_default=5.0,
                                                        check_blocks=true)
fault_df[!,:fid] = string.(fault_df[!,:fid])
println("n faults: ", length(faults))
println("n faults vels: ", length(fault_vels))

geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(geol_slip_rates_file)
geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                    geol_slip_rate_df,
                                                    fault_df)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))


non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(
    x->Oiler.Boundaries.boundary_to_vels(x; ee=5.0, en=5.0),
    non_fault_bounds)...)


println("n non-fault-bound vels: ", length(bound_vels))


gnss_df = Oiler.IO.gis_vec_file_to_df(gnss_vels_file)

gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)
println("n gnss vels: ", length(gnss_vels) )

vels = vcat(fault_vels, 
            gnss_vels, 
            bound_vels,
            );

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                               sparse_lhs=true,
                                               weighted=true,
                                               predict_vels=true,
                                               factorization="lu",
                                               pred_se=true)

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
    "../results/chn_faults_out_mod.geojson",
    name="China fault slip rates")

    #Oiler.IO.write_geol_slip_rate_results_to_csv(results;
    #                outfile="../results/geol_slip_rates.csv")

    #Oiler.IO.write_block_centroid_vels_to_csv(results;
    #                outfile="../results/block_vels.csv",
    #                                    fix="1111")
    #Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups,
    #                name="../results/chn_gnss_results.csv")
                                           
end

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

#if save_results == true
#    if strain_part
#        fig_dir = "/Users/itchy/Desktop/strain_part/"
#    else
#        fig_dir = "/Users/itchy/Desktop/no_strain_part/"
#    end
#    savefig(fig_dir * "/fault_rates.pdf")
#end

show()


#Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
#                                 ref_pole="1111", directory="../web_viewer")

#run(`ipython /Users/itchy/research/gem/china_faults_paper/ss_profiles.py`)
#run(`ipython misfit_plots.py`)
