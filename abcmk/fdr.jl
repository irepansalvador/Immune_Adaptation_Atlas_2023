using CSV, DataFrames, MKtest, Suppressor, ProgressMeter, JLD2,RCall, Random,StatsBase

function control_pvalue(param::MKtest.parameters,path::String,genes::Matrix,data_tgp::String=labstorage * "/annotations/MKdata_may2023.txt",rates::String=labstorage*"/rates_hpc.jld2")


    tmp_folder,tmp_file = splitdir(path)

    mkpath(tmp_folder* "/revision/"*tmp_file)
    tmp_folder = tmp_folder * "/revision/" * tmp_file

    tmp_genes = tmp_folder *"/" * splitdir(tempname())[end]

    CSV.write(tmp_genes,DataFrame(genes,:auto),header=false)

    alpha,sfs,divergence = MKtest.parse_sfs(param,data=data_tgp,gene_list=tmp_genes)

    out = @suppress begin
        try

            @time summstat   = MKtest.summary_statistics(param,sfs,divergence,h5_file=rates,output_folder = tmp_folder,summstat_size=10^5);
            @time posteriors = MKtest.ABCreg(output_folder=tmp_folder,S=length(param.dac),tol=2500,rm_summaries=true);

            out        = MKtest.summary_abc(posteriors,stat="mode");

            insertcols!(out[:inference],1,:type=>"case")
            insertcols!(out[:inference],1,:cell=>tmp_file)
            insertcols!(out[:inference],1,:line=>splitpath(tmp_folder)[end-3])
            out[:inference]
        catch
            DataFrame()
        end
    end

    rm.(filter(x -> occursin("gz",x),readdir(tmp_folder,join=true)))

    return(out)
end

labstorage = "/labstorage/jmurgamoreno/immune_adaptation_atlas/raw_data/"

param = adap = MKtest.parameters(n=661,dac=[2,4,5,10,20,50,200,661],cutoff=[0.0,0.7]);


df_abcmk = vcat(map(x->CSV.read(x,DataFrame)[:,1:6],filter(x-> occursin("inference",x),readdir(labstorage * "/fdr",join=true)))...)
df_abcmk = df_abcmk[df_abcmk.type .== "case",:]
df_abcmk.cell .= df_abcmk.line .* "_" .* df_abcmk.cell
df_abcmk.cell .= replace.(df_abcmk.cell,"exdef_folder_"=>"")

files_c = filter(x-> !isdir(x),readdir(labstorage * "/fdr/cell_lines_outputs/",join=true))[1:end-1]
files_a = filter(x-> !isdir(x),readdir(labstorage * "/fdr/adult_outputs/",join=true))[1:end-1]
files_m = filter(x-> !isdir(x),readdir(labstorage * "/fdr/macrophages_outputs/",join=true))[1:end-1]
files_l = filter(x-> !isdir(x),readdir(labstorage * "/fdr/lung_outputs/",join=true))[1:end-1]
files = Dict("developmental"=>files_c,"adult"=>files_a,"activated_macrophages"=>files_m,"lung"=>files_l)



out_fdr = DataFrame[]
out_error = String[]
@showprogress for dataset in keys(files)

    dataset = replace.(files[key],"_control.txt"=>"")
    dataset .= replace.(dataset,"_case.txt"=>"")
    dataset = unique(dataset)
    @show dataset
    for i in dataset
        @show i

        control = CSV.read(i * "_control.txt",Tables.matrix,header=false)[:,2:end];
        case = CSV.read(i * "_case.txt",Tables.matrix,header=false);

        tmp = shuffle(unique(vcat(case,vec(control))));

        fdr_genes = permutedims(hcat(sample.(fill(tmp,10),length(case),replace=true)...));

        try
            df_tmp = control_pvalue(param,i,fdr_genes);

            df_case = df_abcmk[df_abcmk.cell .== splitdir(i)[end],:];
            df_out = hcat(unique(df_case[:,1:2]),DataFrame(:type1=>sum( df_tmp.α .>= df_case.α)/size(df_tmp,1)));
            insertcols!(df_out,1,:dataset=>key)

            push!(out_fdr,df_out)
        catch e
            push!(out_error,i)
        end
    end
end

df_fdr = vcat(out_fdr...)
df_fdr.cell .= replace.(df_fdr.cell,"Myeloid_ALL_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"Myeloid_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"Lymphoid_ALL_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"Lymphoid_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"HSC_progenitors_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"Tcells_"=>"")
df_fdr.cell .= replace.(df_fdr.cell,"Bcells_"=>"")
df_fdr.line .= replace.(df_fdr.line,"exdef_folder"=>"Activated macrophages")

CSV.write(labstorage * "/fdr/fdr.txt",df_fdr)

df_fdr = CSV.read(labstorage * "/fdr/fdr_all.txt",DataFrame)


@rput df_fdr
@rput labstorage
R"""
df = as.data.table(df_fdr)
df$fdr = 0

target = c("developmental","adult","activated_macrophages","lung")
out=list();for(i in target){
    tmp = df[dataset==i]
    tmp$fdr = p.adjust(tmp$type1,'fdr')
    out[[i]] = tmp[order(type1)]
}

df_dataset = rbindlist(out)
fwrite(df_dataset,paste0(labstorage,"/fdr/fdr_dataset.txt"))
"""
