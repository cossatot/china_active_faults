using JSON

block_json = JSON.parsefile("/Users/itchy/research/gem/fault_data/china/block_data/chn_blocks.geojson");
fault_json = JSON.parsefile("/Users/itchy/research/gem/fault_data/china/block_data/chn_faults.geojson");

blocks = block_json["features"];
faults = fault_json["features"];

function check_winding_order(coords::Array{Float64,2})

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i,:], coords[i - 1,:]) for i in 2:size(coords, 1)])))
end


function check_winding_order(coords)

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i], coords[i - 1]) for i in 2:size(coords, 1)])))
end

#check_winding_order(blocks[7]["geometry"]["coordinates"][1])

for block in blocks
    if check_winding_order(block["geometry"]["coordinates"][1]) == -1
        block["geometry"]["coordinates"][1] = reverse(block["geometry"]["coordinates"][1])
    end
end

function line_to_segs(line)
    [round.([line[i][1] line[i][2];line[i+1][1] line[i+1][2]]; digits=5) for i in 1:length(line)-1]
end

#line_to_segs(blocks[1]["geometry"]["coordinates"][1][1:6])

block_segs = Dict("nea001" => [], "nea002"=>[])
for block in blocks
    block_segs[string(block["properties"]["fid"])] = line_to_segs(block["geometry"]["coordinates"][1])
end

fault_segs = Dict()
for fault in faults
    fault_segs[fault["properties"]["fid"]] = line_to_segs(fault["geometry"]["coordinates"])
end

function check_hw_fw(fault, block_segs=block_segs, fault_segs=fault_segs)
    hw = fault["properties"]["hw"]
    fw = fault["properties"]["fw"]
    hw_adj = false
    fw_adj = false
    trace_segs = fault_segs[fault["properties"]["fid"]]
    
    for trace_seg in trace_segs
        if haskey(block_segs, fw)
            if (trace_seg in block_segs[fw])
                fw_adj = true
            end
        end
        if haskey(block_segs, hw)
            if (reverse(trace_seg, dims=1) in block_segs[hw])
                hw_adj = true
            end
        end
    end

    results = Dict{String, Any}("hw" => hw_adj, "fw" => fw_adj)
    
    if !fw_adj
        for (block, segs) in block_segs
            if (trace_segs[1] in segs)
                results["good_fw"] = block
                #println(block)
            end
        end
    end
    if !hw_adj
        for (block, segs) in block_segs
            if (reverse(trace_segs[1], dims=1) in segs)
                #println(block)
                results["good_hw"] = block
            end
        end
    end
            
    results
end

bad_fault_counter = 0
for fault in faults
    hw_fw_check = check_hw_fw(fault)

    if hw_fw_check["hw"] + hw_fw_check["fw"] < 2
        global bad_fault_counter = bad_fault_counter + 1
    end

    if hw_fw_check["hw"] == false
        println(fault["properties"]["fid"], ": hw: ", fault["properties"]["hw"])
        if haskey(hw_fw_check, "good_hw") 
            println("  good hw: ", hw_fw_check["good_hw"])
        end
    end
    if hw_fw_check["fw"] == false
        println(fault["properties"]["fid"], ": fw: ", fault["properties"]["fw"])
        if haskey(hw_fw_check, "good_fw") 
            println("  good fw: ", hw_fw_check["good_fw"])
        end
    end
end
    
println("n bad faults: ", bad_fault_counter)
