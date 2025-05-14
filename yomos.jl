## [1] Warm up 
## [1.1] Some basic collection types 

# constructor
x = [1, 2, 3] # Vector{Int64}
y = [1.0, 2, 3] # Vector{Float64}
nt = (a=1.9, b=2, c="hello world") # named tuple (immutable), has no fixed element type
M = [1 2; 3 4] # Matrix{Int64}
[i^2 for i in 1:3] # list comprehension a la python

# indexing
x[1]
M[1,2]
nt[1]
nt.c

# [1.2] composite types (objects) 
# with parametric type (real and imaginary part have the same numeric type)
mutable struct ComplexNumber{T<:Real}
  re::T 
  im::T
end

# a constructor function is automatically added
z = ComplexNumber(1, 0)

# not ok, because real and imaginary part do not have the same type
z = ComplexNumber(1.0, 0)

# [1.3] Multiple Dispatch: A key feature of julia 

# ...but we can add constructors to allow promotion to the closest super type
function ComplexNumber(r::Real, i::Real)
  r, i = promote(r, i)
  return ComplexNumber(r, i)
end
ComplexNumber(r) = ComplexNumber(r, zero(r))

# Lets add a function for pretty printing of our type to the REPL
Base.show(io::IO, ::MIME"text/plain", z::ComplexNumber) = print(io, "$(z.re) + $(z.im)i")

# now this works
ComplexNumber(2.0, 1)
ComplexNumber(1)

# we can also define the behaviour for existing functions 
import Base: +
+(x::ComplexNumber,y::ComplexNumber) = ComplexNumber(x.re + y.re, x.im + y.im)

z1, z2 = ComplexNumber(1, 2), ComplexNumber(3.0, -1.2)
z1 + z2 

# Another illustrating example 
abstract type Animal end
struct Dog  <: Animal
  name::String
  age::Int64
end
struct Cat <: Animal
  name::String
  age::Int64
end
struct Mouse <: Animal
  name::String
  age::Int64
end

# define interactions
meet(x::Dog, y::Cat) = println("Dog $(x.name) chases cat $(y.name) up a tree")
meet(x::Mouse, y::Dog) = println("Mouse $(x.name) scares dog $(y.name) away")
meet(x::Cat, y::Mouse) = println(x.age > 10 ? "Cat $(x.name) is to old to see mouse $(y.name)" : "Cat $(x.name) sneaks up on mouse $(y.name)")
# interactions are symmetric
meet(x::Animal, y::Animal) = meet(y, x)

meet(Cat("Bob", 10), Dog("Alice", 1))
meet(Dog("Alice", 1), Mouse("Jerry", 10))
meet(Cat("Tom", 15), Mouse("Jerry", 10))
meet(Cat("Bill", 8), Mouse("Jerry", 10))

# the compiler automatically dispatches to the best match
# we can check that with the macro `@which`
@which meet(Cat("Bob", 10), Dog("Alice", 1))
@which meet(Dog("Alice", 1), Cat("Bob", 10))

## [2] ODE Simulations and Plots

# loading packages (this makes all exported symbols of a package available in the Main namespace)
using DifferentialEquations 
using GLMakie

# define RHS of ODE system (function must have this signature)
# note that f! mutates du (in julia mutating functions are typically indicated with `!`)
function f!(du, u, p, t)
  x, y = u
  ρ, α, θ, δ, γ = p 
  du .= [
    ρ*x*(1-x) - α/(θ + x)*x*y
    -δ*y + γ*α/(θ + x)*x*y
  ]
  return nothing # return `nothing`, update `du` directly (in-place function)
end

# inital conditions, parameters and time span
u0 = [0.1,0.2]
p = (ρ=1.0, α=2.0, θ=1.0, δ=0.5, γ=1.0)
tspan = (0.0, 50.0)

# Solve ODE problem
problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem; saveat=LinRange(tspan..., 501))

# simple plot 
fig, ax, l = lines(sol; label=["x", "y"])
Legend(fig[1,2], l, ["x", "y"])

# Alternatively: Plot with some more control 
# tell Makie how to plot ODESolution as series
Makie.convert_arguments(T::Type{<:Series}, sol::ODESolution) = Makie.convert_arguments(T, sol.t, hcat(sol.u...))

fig, ax, s = series(sol; labels=["x", "y"], color=[Makie.wong_colors()[1], Makie.wong_colors()[2]]);
axislegend(ax)
fig

## Integration with other packages is generally very nice in julia (duo to multiple dispatch)
using Measurements

# p = (ρ=2.4 ± 0.05, α=5.4 ± 0.04, θ=0.8 ± 0.08, δ=0.4 ± 0.01, γ=1.0 ± 0.09)
p = (ρ=1.8 ± 0.05, α=0.9 ± 0.04, θ=0.8 ± 0.08, δ=0.2 ± 0.01, γ=1.0 ± 0.09)
u0 = [0.1, 0.2] .± 0.05
tspan = (0.0, 100.0) .± 0.0
problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem; saveat=LinRange(tspan..., 501))

fig = Figure();
ax = Axis(fig[1,1])
band!(ax, sol.t, sol[1,:]; color=Makie.wong_colors()[1], alpha=0.3)
band!(ax, sol.t, sol[2,:]; color=Makie.wong_colors()[2], alpha=0.3)
lines!(ax, sol.t, sol[1,:]; color=Makie.wong_colors()[1], label="x")
lines!(ax, sol.t, sol[2,:]; color=Makie.wong_colors()[2], label="y")
axislegend(ax)
fig

## [3] GUI 

p = (ρ=1.8, α=2.1, θ=0.4, δ=0.3, γ=0.3)
u0 = [0.1, 0.2]
tspan = (0.0, 100.0)

fig = Figure(size=(1200,800));
ax = Axis(fig[1,1]);

# create sliders for interactivity
slider_names = [string.([keys(p)...]); ["x(0)", "y(0)"]]
sliders = [(
  label = name,
  range = 0:0.1:10,
  format = x -> "$x"
) 
  for name in slider_names
]
# organize sliders in grid
lsgrid = SliderGrid(fig[1:2,2], sliders..., tellheight=false)
# sliders values as Vector of observables
sliderobservables = [s.value for s in lsgrid.sliders]
# slider values as observable Vector (listens to changes in any slider)
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

# solve ode problem with parameters defined by slider
_p = @lift $(slider_values)[1:length(p)]
_u0 = @lift $(slider_values)[length(p)+1:end]
_tspan = @lift (0.0, $T_end)
# remake and solve ODEProblem with new arguments
prob = @lift remake(problem; u0=$_u0, tspan=$_tspan, p=$_p);
sol = @lift solve($prob; saveat=LinRange(0.0, $T_end, 501));
# plot solution (Makie recognizes observables and updates the plot on change)
series!(ax, sol, labels=["x", "y"], color=Makie.wong_colors()[1:2], linewidth=2)
axislegend(ax)

# update axes range on change
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


