## 1 Fundamentals
## 1.1 Some Collection Types 
x = [1, 2, 3] # Vector{Int64}
y = [1.0, 2, 3] # Vector{Float64}
nt = (a=1.9, b=2, c="hello world") # named tuple (immutable), has no fixed element type
M = [1 2; 3 4] # Matrix{Int64}

[i^2 for i in 1:3]
x[1]
M[1,2]
nt.c

for i in eachindex(nt)
    println(nt[i])
end

## 1.2 Functions
function hello(x)
    println("hello $(x)!")
    return nothing
end

f(x) = x^2
g = x -> sqrt(x)
hello("world")
g(4)

## 1.3 Composite Types (i.e. Objects) 

# Type def
mutable struct ComplexNumber{T<:Real}
  re::T 
  im::T
end

# Constructor and field access
z = ComplexNumber(1, 2)
z.re
z.im

# pretty display
Base.show(io::IO, ::MIME"text/plain", z::ComplexNumber) = print(io, "$(z.re) + $(z.im)i")

# Constructors for mixed types
ComplexNumber(1.0, 0.0) # works
ComplexNumber(1.0, 0)   # does not work (type mismatch)

## 1.4 Multiple Dispatch: A Key Feature of Julia 

## 1.4.1 Adding Functionality to our `ComplexNumber` Type
# define constructor for mixed types
function ComplexNumber(r::Real, i::Real)
  r, i = promote(r, i)
  return ComplexNumber(r, i)
end;

ComplexNumber(0.3, true)

ComplexNumber(r::Real, i::Real) = ComplexNumber(promote(r,i)...)
ComplexNumber(r::Real) = ComplexNumber(r, zero(r))
ComplexNumber(1)

# add functionality
import Base: +
+(x::ComplexNumber, y::ComplexNumber) = ComplexNumber(x.re + y.re, x.im + y.im)
+(x::ComplexNumber, y::Real) = ComplexNumber(x.re + y, x.im)

ComplexNumber(1, 2) + ComplexNumber(3.0, -1.2)
ComplexNumber(1, 2) + 1

## 1.4.2 Another Illustrating Example 
# define objects/types
abstract type Animal end
struct Mouse <: Animal
    name::String
end
struct Cat <: Animal
    name::String
end
struct Dog  <: Animal
    name::String
    age::Int64
end

# define behaviour
meet(x::Mouse, y::Dog) = println("Mouse \"$(x.name)\" scares dog \"$(y.name)\" away.")
meet(x::Cat, y::Mouse) = println("Cat \"$(x.name)\" sneaks up on mouse \"$(y.name)\".")
meet(x::Dog, y::Cat) = println(x.age > 10 ? "Dog \"$(x.name)\" is to old to run fast: Cat \"$(y.name)\" strolls away." : "Dog \"$(x.name)\" chases cat \"$(y.name)\" up a tree.")

# symmetry
meet(x::Animal, y::Animal) = meet(y, x)

# create instances
alice = Cat("Alice")
bob = Dog("Bob", 1)
jerry = Mouse("Jerry")
bill = Dog("Bill", 10)

# what's happening?
meet(alice, bob)
meet(bob, jerry)
meet(alice, jerry)
meet(bill, alice)

# dispatch
@which meet(alice, bob) # meet(::Animal, ::Animal)
@which meet(bob, alice) # meet(::Dog, ::Cat)


## 2 ODE Simulations and Plots

## 2.1 First Steps

# load packages
using DifferentialEquations 
using GLMakie;

# def RHS
function f!(du, u, p, t)
    x, y = u # we can unpack vectors or tuples like this
    ρ, α, θ, δ, γ = p # julia has support for unicode 🥳
    du .= [
        ρ*x*(1-x) - α/(θ + x)*x*y
        -δ*y + γ*α/(θ + x)*x*y
    ]
    return nothing # we return `nothing` because we update `du` in-place 
end;

# initial conditions, parameters and time span
u0 = [0.1,0.2]
p = (ρ=1.0, α=2.0, θ=1.0, δ=0.5, γ=1.0)
tspan = (0.0, 50.0);

# def problem and solve with default algorithm
problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem; saveat=LinRange(tspan..., 501));

# simple plot
lines(sol)

# plot with more control 
# tell Makie what to do with ODESolution object
Makie.convert_arguments(T::Type{<:Series}, sol::ODESolution) = Makie.convert_arguments(T, sol.t, hcat(sol.u...))

fig = Figure();
ax = Axis(fig[1,1], xlabel="time", ylabel="population density", title="The MacArthur-Rosenzweig Model")
series!(ax, sol; 
        labels=["x", "y"],
        color=color=Makie.wong_colors()[1:2],
        linewidth=2);
axislegend(ax)
fig

# sol object has interpolation (this is a function call)
sol(10.12345)

## 2.2 Integration with `Measurements.jl`
# now parameters have uncertainty
using Measurements
p = (ρ=1.8 ± 0.2, α=1.1 ± 0.04, θ=0.8 ± 0.08, δ=0.2 ± 0.01, γ=1.0 ± 0.09)
u0 = [0.1, 0.2] .± 0.00
tspan = (0.0, 100.0) .± 0.0

problem = ODEProblem(f!, u0, tspan, p)
sol = solve(problem; saveat=LinRange(tspan..., 501))

# plot solution with uncertainty band
labels = ["prey", "predator"]
fig = Figure();
ax = Axis(fig[1,1])
for i in 1:2
    band!(ax, sol.t, sol[i,:]; color=Makie.wong_colors()[i], alpha=0.3)
    lines!(ax, sol.t, sol[i,:]; color=Makie.wong_colors()[i], label=labels[i])
end
axislegend(ax)
fig

## 3 Interactive GUI for ODE Simulations

function odegui(f!::Function, u0::NamedTuple, p::NamedTuple; 
                size::Tuple{Int64, Int64}=(800,500), # optional arguments
                range=0:0.1:10)

    # create an empty figure and axis
    fig = Figure(size=size);
    ax = Axis(fig[1,1]);

    # create sliders for interactivity
    # use keys of named tuples for state variables and parameters
    slider_names = [string.([keys(p)...]); string.([keys(u0)...]) .* "(0)"]
    sliders = [(
        label = name,
        range = 0:0.1:10,
        format = x -> "$x"
    ) 
        for name in slider_names
    ]
    # organize sliders in a grid
    lsgrid = SliderGrid(fig[1:2,2], sliders..., tellheight=false)
    # sliders values as Vector{Observable{Float64}}
    sliderobservables = [s.value for s in lsgrid.sliders]
    # slider values as Observable{Vector{Float64}} (listens to changes in any slider)
    slider_values = lift(sliderobservables...) do slvalues...
        [slvalues...]
    end
    # set sliders close to inital values
    for i in 1:length(sliderobservables)
        set_close_to!(lsgrid.sliders[i], [values(p)..., u0...][i])
    end

    # add slider for tspan with hard coded choices
    T_end_slider = SliderGrid(fig[2,1], (label = "T", range = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], format = x->"$x"))
    T_end = T_end_slider.sliders[1].value
    set_close_to!(T_end_slider.sliders[1], 50)

    # solve ode problem with parameters defined by slider
    # the @lift macro can be used to define variables that depend on other observable variables and need to be updated
    _p = @lift $(slider_values)[1:length(p)]
    _u0 = @lift $(slider_values)[length(p)+1:end]
    _tspan = @lift (0.0, $T_end)
    # remake and solve ODEProblem with new arguments
    prob = @lift remake(problem; u0=$_u0, tspan=$_tspan, p=$_p);
    sol = @lift solve($prob; saveat=LinRange(0.0, $T_end, 501));
    # plot solution (Makie recognizes observables and updates the plot on change)
    series!(ax, sol, labels=["x", "y"], color=Makie.wong_colors()[1:2], linewidth=2)
    axislegend(ax)

    # update plot ranges on change
    on(sol) do sol
        y_range = [minimum(sol), maximum(sol)]
        y_pad = 0.05*(y_range[2] - y_range[1])
        x_pad = 0.05*(sol.t[end] - sol.t[1])
        xlims!(ax, 0.0 - x_pad, sol.t[end] + x_pad)
        ylims!(ax, y_range[1] - y_pad, y_range[2] + y_pad)
    end

    # update layout
    rowgap!(lsgrid.layout, 7)
    colsize!(fig.layout, 1, Relative(2/3))

    return fig
end

u0 = (x=0.1, y=0.2)
p = (ρ=1.0, α=2.0, θ=1.0, δ=0.5, γ=1.0)
tspan = (0.0, 50.0)

odegui(f!, u0, p)

