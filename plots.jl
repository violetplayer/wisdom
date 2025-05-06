using Plots
using LaTeXStrings
pythonplot()
default(guidefont=(20), tickfont=(18), titlefont=(28), legendfontsize=16, legend=:none)

function plot_view_weights!(W::SparseMatrixCSC, dir::String, views::Vector{AbstractView}, OSgrid::OrbitSpaceGrid)

    _, _, Pstar, _ = get_orbit_bounds(state, OSgrid)
    for view in views
        wv = zeros(OSgrid.N)
        for sig in view.signals
            wv .+= W[sig.index,:]
        end
        w3d = OSvector_to_OSarray(wv, OSgrid)
        w_PH = [sum(w3d[i,j,:]) for i=eachindex(OSgrid.Pθ), j=eachindex(OSgrid.H)]
        w_PZ = [sum(w3d[i,:,j]) for i=eachindex(OSgrid.Pθ), j=eachindex(OSgrid.ζ)]
        w_HZ = [sum(w3d[:,i,j]) for i=eachindex(OSgrid.H), j=eachindex(OSgrid.ζ)]

        p = plot(size=(600,600),dpi=200)
        contourf!(OSgrid.Pθ, OSgrid.H, transpose(w_PH) / maximum(w_PH), levels=20, cmap=:inferno)
        vline!([Pstar], lw=3, color=:white, style=:dash, label=:none)
        xlabel!(L"$P_\theta / P_0$")
        ylabel!("H (keV)")
        savefig("$(dir)/Wph_$(view.name).png")

        p = plot(size=(600,600),dpi=200)
        contourf!(OSgrid.Pθ, OSgrid.ζ, transpose(w_PZ) / maximum(w_PZ), levels=20, cmap=:inferno)
        vline!([Pstar], lw=3, color=:white, style=:dash, label=:none)
        xlabel!(L"$P_\theta / P_0$")
        ylabel!(L"$\vartheta$ (deg.)")
        savefig("$(dir)/Wpz_$(view.name).png")

        p = plot(size=(600,600),dpi=200)
        contourf!(OSgrid.H, OSgrid.ζ, transpose(w_HZ) / maximum(w_HZ), levels=20, cmap=:inferno)
        xlabel!("H (keV)")
        ylabel!(L"$\vartheta$ (deg.)")
        savefig("$(dir)/Whz_$(view.name).png")
    end
end

function plot_fast_ions!(F::SparseVector, dir::String, OSgrid::OrbitSpaceGrid)

    _, _, Pstar, _ = get_orbit_bounds(state, OSgrid)

    f3d = OSvector_to_OSarray(Vector(F), OSgrid)
    f_PH = [sum(f3d[i,j,:]) for i=eachindex(OSgrid.Pθ), j=eachindex(OSgrid.H)]
    f_PZ = [sum(f3d[i,:,j]) for i=eachindex(OSgrid.Pθ), j=eachindex(OSgrid.ζ)]
    f_HZ = [sum(f3d[:,i,j]) for i=eachindex(OSgrid.H), j=eachindex(OSgrid.ζ)]

    p = plot(size=(600,600),dpi=200)
    contourf!(OSgrid.Pθ, OSgrid.H, transpose(f_PH), levels=20, cmap=:inferno)
    vline!([Pstar], lw=3, color=:white, style=:dash, label=:none)
    xlabel!(L"$P_\theta / P_0$")
    ylabel!("H (keV)")
    savefig("$(dir)/Fph.png")

    p = plot(size=(600,600),dpi=200)
    contourf!(OSgrid.Pθ, OSgrid.ζ, transpose(f_PZ), levels=20, cmap=:inferno)
    vline!([Pstar], lw=3, color=:white, style=:dash, label=:none)
    xlabel!(L"$P_\theta / P_0$")
    ylabel!(L"$\vartheta$ (deg.)")
    savefig("$(dir)/Fpz.png")

    p = plot(size=(600,600),dpi=200)
    contourf!(OSgrid.H, OSgrid.ζ, transpose(f_HZ), levels=20, cmap=:inferno)
    xlabel!("H (keV)")
    ylabel!(L"$\vartheta$ (deg.)")
    savefig("$(dir)/Fhz.png")
end

function plot_npa_spectra!(S::Vector{Float64}, dir::String, views::Vector{AbstractView})
    for (i,view) in enumerate(views)
        name = view.name
        if view isa NPAChord
            ΔE = [sig.ΔE_ch for sig in view.signals]
            x = [sig.E_ch for sig in view.signals]
            y = S[[sig.index for sig in view.signals]]
            
            p = plot(size=(500,500),dpi=200)
            plot!(x, normalize!(y ./ ΔE), lw=3, color=:black)
            xlabel!("E (keV)")
            ylabel!("S (arb.)")
            savefig("$(dir)/spectrum_$(name).png")
        end
    end
end
# function get_KL_map(signals::Array{AbstractSignal})
#     for i=eachindex(signals)
#         for j=eachindex(signals)
#             if i != j

#             end
#         end
#     end
# end