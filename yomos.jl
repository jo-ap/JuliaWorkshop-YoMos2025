## Multiple Dispatch: A key feature of julia 


## [1] Simulations and Plots

# loading packages
using DifferentialEquations 
using GLMakie

# choose ode solver (this is a Runge-Kutta solver that switches to Rodas4P if
# stiffness is detected)
ode_alg = AutoTsit5(Rodas4P());

# define RHS of ODE system (function must have this signature and mutates du)
function f!(du, u, p, t)
  x, y = u
  ρ, α, θ, δ, γ = p 
  du .= [
    ρ*x*(1-x) - α/(θ + x)*x*y
    -δ*y + γ*α/(θ + x)*x*y
  ]
  return nothing
end

# inital conditions, parameters and time span
u0 = [0.1,0.2]
p = (ρ=1.0, α=2.0, θ=1.0, δ=0.5, γ=1.0)
tspan = (0.0, 50.0)

# Solve ODE problem
problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem, ode_alg; saveat=LinRange(tspan..., 501))

# simple plot 
lines(sol)

# Alternatively: Plot with some more control 
# plot ODESolution type with Makie as series
Makie.convert_arguments(T::Type{<:Series}, sol::ODESolution) = Makie.convert_arguments(T, sol.t, hcat(sol.u...))

fig, ax, s = series(sol; labels=["x", "y"], color=[Makie.wong_colors()[1], Makie.wong_colors()[2]]);
axislegend(ax)
fig

## Integration with other packages is generally very nice in julia
using Measurements

p = (ρ=1.0 ± 0.04, α=2.0 ± 0.05, θ=1.0 ± 0.08, δ=0.5 ± 0.01, γ=1.0 ± 0.09)
p = (ρ=2.4 ± 0.05, α=5.4 ± 0.04, θ=0.8 ± 0.08, δ=0.4 ± 0.01, γ=1.0 ± 0.09)
p = [1.0 ± 0.05, 2.0 ± 0.1, 1.0 ± 0.08, 0.5 ± 0.01, 1.0 ± 0.09]
# p = [1.0 ± 0.05, 2.0 ± 0.1, 1.0 ± 0.08, 0.5 ± 0.01, 1.0 ± 0.09]
u0 = [0.1, 0.2] .± 0.0
tspan = (0.0, 30.0) .± 0.0
problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem; saveat=LinRange(tspan..., 501))

fig = Figure();
ax = Axis(fig[1,1])
band!(ax, sol.t, sol[1,:]; color=Makie.wong_colors()[1], alpha=0.3)
band!(ax, sol.t, sol[2,:]; color=Makie.wong_colors()[2], alpha=0.3)
lines!(ax, sol.t, sol[1,:]; color=Makie.wong_colors()[1])
lines!(ax, sol.t, sol[2,:]; color=Makie.wong_colors()[2])
fig

## [2] Fit model to data 

# load additional packages
using DifferentialEquations
using Optimization
using OptimizationPolyalgorithms
using SciMLSensitivity
using ForwardDiff

# Define dummy data
t_data = 0:1:50
# we can get data at every time point from the ODESolution object by interpolation 
x_data = sol(t_data)[1,:] .+ 0.03.*randn(length(t_data))
y_data = sol(t_data)[2,:] .+ 0.03.*randn(length(t_data))
xy_data = hcat(x_data, y_data)'

## Optimization

# initial guess for parameters
p0 = [values(p)...] .+ 0.5.*rand(length(p))

# Define a loss metric function to be minimized
function loss(p)
    newprob = remake(problem, p = p)
    sol = solve(newprob, saveat = t_data)
    loss = sum(abs2, sol .- xy_data)
    return loss
end

# keep track of parameters and loss during fit
p_evolve = [p0]
l_evolve = [loss(p0)]

# Define a callback function to monitor optimization progress
function callback(state, l)
    # save intermediate parameters and loss values
    push!(p_evolve, state.u)
    push!(l_evolve, l)
    return false
end

# Set up the optimization problem with our loss function and initial guess
adtype = AutoForwardDiff()
optf = OptimizationFunction((x, _) -> loss(x), adtype)
optprob = OptimizationProblem(optf, p0)

# Optimize the ODE parameters for best fit to our data
pfinal = solve(optprob, PolyOpt(),
    callback = callback,
    maxiters = 200)
p_fit = round.(pfinal, digits = 2)

## Visualize fit process

# Makie can create responsive plots using Observables
i = Observable(1)
title_str = @lift "iter = $($i-1), loss = $(round(l_evolve[$i]; digits = 4))"
p_current = @lift p_evolve[$i]

# data and true solution
fig = Figure(size=(1200,800));
Label(fig[1,1:2], title_str; fontsize=20)
ax1 = Axis(fig[2,1], xlabel="time", title="Simulation")
scatter!(ax1, t_data, x_data; label="x data")
scatter!(ax1, t_data, y_data; label="y data")
series!(ax1, sol; linestyle=:dot, labels=["x true soluton", "y true solution"], color=[Makie.wong_colors()[1], Makie.wong_colors()[2]])
problem_current = @lift ODEProblem(f!, u0, tspan, $p_current);
sol_current = @lift solve($problem_current; saveat=LinRange(t_data[1], t_data[end], 101));
series!(ax1, sol_current; labels=["x fit", "y fit"], color=[Makie.wong_colors()[1], Makie.wong_colors()[2]])
axislegend(ax1)

# parameters
ax2 = Axis(fig[2,2], title="Parameters")
barplot!(0.75:2:(2*length(p)-0.25), [values(p)...]; gap=0.5, label="true")
barplot!(1.25:2:(2*length(p)+0.25), p_current; gap=0.5, label="fit")
axislegend(ax2)
ylims!(ax2, high=1.05*maximum(maximum.(p_evolve)))
fig

for k in eachindex(p_evolve)
  i[] = k 
  sleep(0.1)
end


## [3] GUI 

fig = Figure(size=(1200,800));
ax = Axis(fig[1,1]);

# sliders for interactivity
slider_names = [string.([keys(p)...]); ["x(0)", "y(0)"]]
sliders = [(
  label = name,
  range = 0:0.1:10,
  format = x -> "$x"
) 
  for name in slider_names
]
lsgrid = SliderGrid(fig[1:2,2], sliders..., tellheight=false)
sliderobservables = [s.value for s in lsgrid.sliders]
slider_values = lift(sliderobservables...) do slvalues...
  [slvalues...]
end
# set to inital values
for i in 1:length(sliderobservables)
  set_close_to!(lsgrid.sliders[i], [values(p)..., u0...][i])
end

# add slider for tspan 
T_end_slider = SliderGrid(fig[2,1], (label = "T", range = [10, 100, 200, 500, 1000, 2000, 5000, 10000], format = x->"$x"))
T_end = T_end_slider.sliders[1].value

# solve ode problem depending on slider
_p = @lift $(slider_values)[1:length(p)]
_u0 = @lift $(slider_values)[length(p)+1:end]
_tspan = @lift (0.0, $T_end)
prob = @lift remake(problem; u0=$_u0, tspan=$_tspan, p=$_p);
sol = @lift solve($prob, ode_alg; saveat=LinRange(0.0, $T_end, 501));
series!(ax, sol, labels=["x", "y"], color=Makie.wong_colors()[1:2], linewidth=2)
axislegend(ax)

# update y axis range on change
on(sol) do sol
  y_range = [minimum(sol), maximum(sol)]
  y_pad = 0.05*(y_range[2] - y_range[1])
  x_pad = 0.05*(sol.t[end] - sol.t[1])
  xlims!(ax, 0.0 - x_pad, sol.t[end] + x_pad)
  ylims!(ax, y_range[1] - y_pad, y_range[2] + y_pad)
end
  
# layout
rowgap!(lsgrid.layout, 7)
colsize!(fig.layout, 1, Relative(2/3))

fig

