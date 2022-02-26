
struct IpoptWrapper <: NLPSolverWrapper
    # Ipopt problem
    prob::IpoptProblem

    # Options dictionary
    options::Dict
end

# Function evaluations should use a different update such that the NLP is not reevaluated
# if the same decision variable is passed to geval! after feval!, etc...
"""
    IpoptWrapper(feval!, geval!, gradfeval!, jacgeval!, n, x_L, x_U,
        m, g_L, g_U, nele_jac)

    Constructs IpoptWrapper type. Instantiates Ipopt problem and sets 
    hessian_aproximation to limited-memory as second derivates are not supported (yet).

    # Arguments
    - feval::Function : Objective function. In the form of feval(x) -> f(x)
    - geval!::Function : Constraint function. In the form geval!(g,x), mutates g
    - gradfeval!::Function : Evaluates gradient of objective function. In the form gradfeval!(grad,x), mutates grad 
    - jacgeval!::Function : Evaluates Jacobian of constraints. In the form jacgeval(values,rows,cols,x), mutates 
                            values, rows, and cols. Values::Union{Nothing,Vector{Float64}}
    - n::Int : Number of decision variables
    - x_L::Vector{Float64} : Lower bounds on decision variables 
    - x_U::Vector{Float64} : Upper bounds on decision variables 
    - m::Int : Number of constraint functions 
    - g_L::Vector{Float64} : Lower bounds on constraints 
    - g_U::Vector{Float64} : Upper bounds on constriants
    - nele_jac::Int : Number of non-zero elements in the Jacobian
"""
function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, n::Int, x_L::Vector{Float64}, x_U::Vector{Float64},
    m::Int, g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int)

    # Create Ipopt problem
    prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Set hessian approximation (required unless we decide to support hessian computation)
    Ipopt.AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")

    # Create options Dict (For resetting user options durring problem recreation after mesh refinement)
    options = Dict("hessian_approximation" => "limited-memory")

    # Craete Ipopt Wrapper 
    return IpoptWrapper(prob, options)
end

function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, n::Int, x_L::Vector{Float64}, x_U::Vector{Float64},
    m::Int, g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, options::Dict)

    # Create Ipopt problem
    prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Set options 
    SetOptions!(prob, options)

    # Create IpoptWrapper 
    return IpoptWrapper(prob, options)
end

function SetInitialGuess!(wrapper::IpoptWrapper, x::Vector{Float64})
    if length(x) != wrapper.prob.n 
        error("Decision vector guess passed to SetInitialGuess! with incorrect length!")
    end
    @inbounds wrapper.prob.x .= x
    return nothing
end

function Optimize!(wrapper::IpoptWrapper)
    @warn "Optimize!(wrapper::IpoptWrapper does not do anything after solving!"
    solvestat = Ipopt.IpoptSolve(wrapper.prob)
end

function SetStringOption!(wrapper::IpoptWrapper, str1::String, str2::String)
    push!(wrapper.options, str1 => str2)
    Ipopt.AddIpoptStrOption(wrapper.prob, str1, str2)
    return nothing
end

function SetIntOption!(wrapper::IpoptWrapper, str::String, int::Int)
    push!(wrapper.options, str => int)
    Ipopt.AddIpoptIntOption(wrapper.prob, str, int)
    return nothing
end

function SetFloatOption!(wrapper::IpoptWrapper, str::String, float::AbstractFloat)
    push!(wrapper.options, str => float)
    Ipopt.AddIpoptNumOption(wrapper.prob, str, float)
    return nothing
end

# Functions to set options via dictionary
function SetOptions!(prob::IpoptProblem, options::Dict)
    # Itterate through keys in options
    for key in keys(options)

        # Get value of key
        value = options[key]

        # Call option setters depending on type of value
        if value isa String 
            Ipopt.AddIpoptStrOption(prob, key, value)
        elseif value isa Integer
            Ipopt.AddIpoptIntOption(prob, key, value)
        elseif value isa Float64
            Ipopt.AddIpoptNumOption(prob, key, value)
        else
            error("Option " * key * " passed with value of invalid type " * typeof(value))
        end
    end
end
SetOptions!(wrapper::IpoptWrapper, options::Dict) = SetOptions!(wrapper.prob, options)
