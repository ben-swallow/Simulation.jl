function _compartment_idx(compartment, names)
    idx = findfirst(names .== compartment)
    isnothing(idx) && throw(ArgumentError("Compartment $compartment not in $names"))
    return idx
end

function _default_steps(abuns)
    N_steps = size(abuns, 3)
    return Int.(floor.(range(1, N_steps; length=4)))
end

"""
    plot_epiheatmaps(
        epi::AbstractEpiSystem,
        abuns::AbstractArray{<:Integer, 3};
        compartment="Exposed",
        steps=[],
    )

Plot heatmaps of `abuns` for `compartment` at `steps`.

## Arguments
- `epi`: The `AbstractEpiSystem` to plot.
- `abuns`: The array of abundances to plot, of size Ncompartments x Ncells x Nsteps

## Keyword arguments
- `compartment`: The compartment to plot
- `steps`: A list of steps to plot (one heatmap for each step). If empty, plots 4
    equally-spaced steps.
"""
plot_epiheatmaps
@userplot Plot_EpiHeatmaps
function _check_args(h)
    correct_args = (
        length(h.args) == 2 &&
        isa(h.args[1], AbstractEpiSystem) &&
        isa(h.args[2], AbstractArray{<:Integer, 3})
    )
    if !correct_args
        throw(ArgumentError(
            "$(typeof(h)) requires (AbstractEpiSystem, abuns); got: $(typeof(h.args))"
        ))
    end
end
@recipe function f(
    h::Plot_EpiHeatmaps;
    compartment="Exposed",
    steps=[],
)
    _check_args(h)
    epi, abuns = h.args
    idx = _compartment_idx(compartment, epi.epilist.names)
    if isempty(steps)
        steps = _default_steps(abuns)
    end

    layout := length(steps)

    subplot = 1
    gridsize = size(epi.epienv.habitat.matrix)
    for step in steps
        data = Float64.(reshape(abuns[idx, :, step], gridsize...))
        data[.!epi.epienv.active] .= NaN
        @series begin
            seriestype := :heatmap
            title := "Day $step ($compartment)"
            subplot := subplot
            background_color --> :lightblue
            background_color_outside --> :white
            aspect_ratio --> 1
            grid --> false
            data
        end
        subplot += 1
    end
end

"""
    plot_epidynamics(
        epi::AbstractEpiSystem,
        abuns::AbstractArray{<:Integer, 3};
        category_map=nothing,
    )

Plot the dynamics of `abuns` summed over space, as a function of time.

## Arguments
- `epi`: The `AbstractEpiSystem` to plot.
- `abuns`: The array of abundances to plot, of size Ncompartments x Ncells x Nsteps

## Keyword arguments
- `category_map`: An iterable of key-value pairs where the keys are category names, and the
    values are a list of compartment indices associated with that category. These
    compartments will be summed in the plot. For example, the following will plot the sum of
    compartments 1 and 2 as the `Susceptible` category, and the sum of compartments 3 and 4
    as the `Infected` category.

        category_map = ("Susceptible" => [1, 2], "Infected" => [3, 4])

    If `category_map` is `nothing`, all compartments in `epi` will be plotted separately
    with their corresponding names.
"""
plot_epidynamics
@userplot Plot_EpiDynamics
@recipe function f(h::Plot_EpiDynamics; category_map=nothing)
    _check_args(h)
    epi, abuns = h.args

    if isnothing(category_map)
        # Make each compartment its own category
        category_map = (name => [idx] for (idx, name) in enumerate(epi.epilist.names))
    end

    for (name, idx) in category_map
        data = vec(mapslices(sum, abuns[idx, :, :], dims = (1, 2)))
        title --> "Infection dynamics"
        xguide --> "Step"
        yguide --> "Totals"
        @series begin
            label := name
            data
        end
    end
end
