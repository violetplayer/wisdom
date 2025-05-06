##==================================================================================================================
## Main file for parallel job execution
##==================================================================================================================

#### minimal imports
using Distributed
using SharedArrays

#### initialize workers
n_cores = parse(Int64, ARGS[2])
if length(workers()) > 1
    rmprocs(workers())
end
addprocs(n_cores)
proc = workers()
println("Workers added: $n_cores")

#### run set up script
@everywhere input_deck = $(ARGS[1])
@everywhere include("setup.jl")

#### kludge fix for now
dir = pwd()
file_id = ""

##==================================================================================================================
## Parallel engines for jobs
##==================================================================================================================

##--------------------------
## Neutral Beam MC
##--------------------------

#### axes: beam number, energy component [1, 1/2, 1/3], x axis, y axis, z axis
fastn = SharedArray{Float64}(length(beams), 3, length(RSgrid.x), length(RSgrid.y), length(RSgrid.z))
nb_src = SharedArray{Float64}(OSgrid.N)

function MCbeams!(n::Integer)
    @sync @distributed for i=1:n
        for j=eachindex(beams)
            ni,nj,nk,w = throw_neutral(beams[j], state, RSgrid)
            fastn[j,:,ni,nj,nk] += w / n
        end
    end
end

function MCsource!(n::Integer)
    @sync @distributed for _=1:n
        for beam in beams
            i,w = throw_neutral(beam, state, OSgrid)
            nb_src[i] += w / n
        end
    end
end

##--------------------------
## Weight functions
##--------------------------

results = SharedArray{Float64}(N + 1,3) # +1 to ensure dimensionality

function engine(n::Integer)

    #### CX method
    function MCweights!(view::CXView, view_index::Integer, n::Integer)
        start_index = (view_index - 1) * n
        if view.beam_index == -1
            @sync @distributed for j=1:n
                I,J,V = throw_particle(view, state, OSgrid, RSgrid) # COO format for sparse matrices
                results[start_index + j, :] .= [I, J, V * OSgrid.V / N]
                #put!(results, (i,j,w*OSgrid.V/N))
            end
        else
            @sync @distributed for j=1:n
                I,J,V = throw_particle(view, state, beams, fn_itp, OSgrid, RSgrid)
                results[start_index + j, :] .= [I, J, V * OSgrid.V / N]
            end
        end
    end

    #### fusion product method
    function MCweights!(view::FusionProductView, view_index::Integer, n::Integer)
        start_index = (view_index - 1) * n
        @sync @distributed for j=1:n
            I,J,V = throw_particle(view, state, RSgrid, OSgrid)
            results[start_index + j, :] .= [I, J, V * OSgrid.V / N]
        end
    end

    for (i,view) in enumerate(views)
        MCweights!(view, i, n)
        pdone = round(i / length(views) * 100, digits=1)
        println("$(pdone)% completed")
    end

    results[end,:] .= [length(signals),OSgrid.N,0.0]
    W = sparse(Integer.(results[:,1]), Integer.(results[:,2]), results[:,3]) ## SparseMatrixCSC
    save_weights(dir, file_id, W)

    return W
end


##==================================================================================================================
## Do the computation
##==================================================================================================================

## timing
t = time()

## Initialize output file
initialize_output_files(dir, file_id, OSgrid, signals)

## if there are no states, exit nicely
if OSgrid.N == 0
    println("Exiting: provided orbit-space axes don't span a non-trivial space.")
    println("\nFinished!")
    exit()
end

##--------------------------
## Compute or load fast neutrals
##--------------------------

if have_fastn
    println("\nLoading fast neutral density from file...")
    load_fastn(run_params["files"]["fastn_file"])
else
    println("\nNo fast neutral file provided. Computing fast neutral distribution...")
    MCbeams!(Nfastn)
    println("\nSaving fast neutrals...")
    #save_fastn(dir, file_id, beams, fastn)
end
@everywhere _fn_itp = interpolate_fast_neutrals(fastn, beams, RSgrid)
@everywhere const fn_itp = _fn_itp
tdone = round((time() - t) / 60, digits=2)
stat = have_fastn ? "loaded" : "computed"
println("\nFast neutral density $(stat), t = $(tdone) minutes.")

if build_source
    println("\nBuilding source distribution...")
    MCsource!(Nfastn)
    F = sparse(nb_src)
end
tdone = round((time() - t) / 60, digits=2)
println("\nSource distribution computed, t = $(tdone) minutes.")

##--------------------------
## Compute + output weight functions + signals
##--------------------------

println("\nComputing weights...")
W = engine(Nion)
tdone = round((time() - t) / 60, digits=2)
println("\nWeight computation done, t = $(tdone) minutes.")

## signals!
S = Vector(W * F) # project to dense vector, doesn't need to be sparse anymore
save_signals(dir, file_id, S)
tdone = round((time() - t) / 60, digits=2)
println("\nSignals computed, t = $(tdone) minutes.")

### no more parallelism
rmprocs(workers())

##--------------------------
## Make plots of weight functions & compute signals
##--------------------------

if make_plots
    #println("\nCannot make plots, because I haven't built this yet!")
    include("plots.jl")
    
    # load
    println("\nMaking plots...")
    plot_view_weights!(W, dir, views, OSgrid)
    plot_fast_ions!(F, dir, OSgrid)
    plot_npa_spectra!(S, dir, views)
    tdone = round((time() - t) / 60, digits=2)
    println("\nPlots made, t = $(tdone) minutes.")
else
    println("\nNot making plots...")
end

tdone = round((time() - t) / 60, digits=2)
println("\nFinished! Total time : $(tdone) minutes.")