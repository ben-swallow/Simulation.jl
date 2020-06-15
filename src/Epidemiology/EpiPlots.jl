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

!!! note
    Heatmaps are transposed by default. Pass in `transpose=false` to turn this off.
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
    class=:Exposed,
    steps=[],
)
    _check_args(h)
    epi, abuns = h.args
    if isempty(steps)
        steps = _default_steps(abuns)
    end

    layout := length(steps)
    # Transpose by default
    # `match_dimensions` is the same as `transpose`
    match_dimensions --> true
    seriescolor --> :heat

    Nsteps = size(abuns)[end]
    ax = (AxisArrays.axes(epi.epienv.initial_population)..., Axis{:step}(1:Nsteps))
    abuns_grid = AxisArray(reshape(abuns, length.(ax)), ax)

    subplot = 1
    gridsize = size(epi.epienv.habitat.matrix)
    for step in steps
        data = abuns_grid[class=class, step=step].data
        reducedims = Tuple(1:ndims(data)-2)
        data = sum(data, dims=reducedims)
        data = dropdims(data, dims=reducedims)
        data = Matrix{Any}(data)
        data[.!epi.epienv.active] .= NaN
        x = 1:size(data, 2)
        y = 1:size(data, 1)
        if plotattributes[:match_dimensions]
            x, y = y, x
        end
        @series begin
            seriestype := :heatmap
            title := "Step $step ($class)"
            subplot := subplot
            background_color --> :lightblue
            background_color_outside --> :white
            aspect_ratio --> 1
            grid --> false
            xlims --> extrema(x)
            ylims --> extrema(y)
            x, y, data
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
@recipe function f(h::Plot_EpiDynamics; classes=[:Susceptible])
    _check_args(h)
    epi, abuns = h.args

    Nsteps = size(abuns)[end]
    ax = (AxisArrays.axes(epi.epienv.initial_population)..., Axis{:step}(1:Nsteps))
    abuns_grid = AxisArray(reshape(abuns, length.(ax)), ax)

    for class in classes
        # Need .data here due to https://github.com/JuliaArrays/AxisArrays.jl/issues/113
        data = abuns_grid[class=class].data
        data = vec(mapslices(sum, data, dims = Tuple(1:ndims(data)-1)))
        title --> "Infection dynamics"
        xguide --> "Step"
        yguide --> "Totals"
        @series begin
            label := class
            data
        end
    end
end
