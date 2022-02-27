
mutable struct IpoptWrapper <: NLPSolverWrapper
    # Ipopt problem
    prob::IpoptProblem

    # Initial guess set flag 
    initGuessSet::Bool 

    # Has optimized flag 
    hasOptimized::Bool

    # User options
    strOpts::Dict{String,String}
    intOpts::Dict{String,Int64}
    numOpts::Dict{String,Float64}
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
    - x_L::Vector{Float64} : Lower bounds on decision variables 
    - x_U::Vector{Float64} : Upper bounds on decision variables 
    - g_L::Vector{Float64} : Lower bounds on constraints 
    - g_U::Vector{Float64} : Upper bounds on constriants
    - nele_jac::Int : Number of non-zero elements in the Jacobian
"""
function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int)

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Create Ipopt problem
    prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Set hessian approximation (required unless we decide to support hessian computation)
    Ipopt.AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")

    # Create string options (For resetting user options durring problem recreation after mesh refinement)
    strOpts = Dict("hessian_approximation" => "limited-memory")

    # Create integer and floating point options 
    intOpts = Dict{String,Int64}()
    numOpts = Dict{String,Float64}()

    # Craete Ipopt Wrapper 
    return IpoptWrapper(prob, false, false, strOpts, intOpts, numOpts)
end

function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, 
    options::Dict{String,T}) where {T}

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Create Ipopt problem
    prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Set options 
    SetOptions!(prob, options)

    # Get concrete dictionaries
    strDict, intDict, numDict = ConcreteDicts(options)

    # Create IpoptWrapper 
    return IpoptWrapper(prob, false, false, strDict, intDict, numDict)
end

function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, 
    strOpts::Dict{String,String}, intOpts::Dict{String,Int64}, numOpts::Dict{String,Float64})

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Create Ipopt problem
    prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Set options 
    SetOptions!(prob, strOpts)
    SetOptions!(prob, intOpts)
    SetOptions!(prob, numOpts)

    # Create IpoptWrapper 
    return IpoptWrapper(prob, false, false, strOpts, intOpts, numOpts)
end

function IpoptWrapper(feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, wrapper::IpoptWrapper)

    return IpoptWrapper(feval, geval!, gradfeval!, jacgeval!, x_L, x_U,
        g_L, g_U, nele_jac, wrapper.strOpts, wrapper.intOpts, wrapper.numOpts)
end

# Function to reset Ipopt optimizer after mesh refinement
function ResetIpoptWrapper!(wrapper::IpoptWrapper, feval::Function, geval!::Function, gradfeval!::Function, 
    jacgeval!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int)

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In ResetIpoptWrapper!, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In ResetIpoptWrapper!, g_L and g_U must be of the same length")

    # Create new Ipopt problem
    wrapper.prob = Ipopt.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 0, 
            feval, (x,g) -> geval!(g,x), (x,grad) -> gradfeval!(grad,x), 
            (x,rows,cols,values) -> jacgeval!(values,rows,cols,x), (args...) -> false)

    # Reset flags 
    wrapper.initGuessSet = false
    wrapper.hasOptimized = false

    # Set options
    SetOptions!(wrapper.prob, wrapper.strOpts)
    SetOptions!(wrapper.prob, wrapper.intOpts)
    SetOptions!(wrapper.prob, wrapper.numOpts)

    return nothing
end

function SetInitialGuess!(wrapper::IpoptWrapper, x::Vector{Float64})
    if length(x) != wrapper.prob.n 
        error("Decision vector guess passed to SetInitialGuess! with incorrect length!")
    end
    @inbounds wrapper.prob.x .= x
    wrapper.initGuessSet = true
    return nothing
end

function Optimize!(wrapper::IpoptWrapper)
    if wrapper.initGuessSet
        solvestat = Ipopt.IpoptSolve(wrapper.prob)
        wrapper.hasOptimized = true
    else
        error("In Optimize!(wrapper::IpoptWrapper), cannot optimize without setting initial guess!")
    end
end

function GetSolution(wrapper::IpoptWrapper)
    if wrapper.hasOptimized == true
        sol = wrapper.prob.x
    else
        error("In GetSolution(wrapper::IpoptWrapper), cannot get solution before optimization!")
    end
    return sol
end

function SetStringOption!(wrapper::IpoptWrapper, str1::String, str2::String)
    push!(wrapper.strOpts, str1 => str2)
    Ipopt.AddIpoptStrOption(wrapper.prob, str1, str2)
    return nothing
end

function SetIntOption!(wrapper::IpoptWrapper, str::String, int::Int64)
    push!(wrapper.intOpts, str => int)
    Ipopt.AddIpoptIntOption(wrapper.prob, str, int)
    return nothing
end

function SetFloatOption!(wrapper::IpoptWrapper, str::String, float::Float64)
    push!(wrapper.numOpts, str => float)
    Ipopt.AddIpoptNumOption(wrapper.prob, str, float)
    return nothing
end

# Functions to set options via dictionary
function SetOptions!(prob::IpoptProblem, options::Dict{String,Any})
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
    return nothing
end
function SetOptions!(prob::IpoptProblem, options::Dict{String,String})
    # Itterate through keys in options 
    for key in keys(options)
        Ipopt.AddIpoptStrOption(prob, key, options[key])
    end
end
function SetOptions!(prob::IpoptProblem, options::Dict{String,Int64})
    # Itterate through keys in options 
    for key in keys(options)
        Ipopt.AddIpoptIntOption(prob, key, options[key])
    end
end
function SetOptions!(prob::IpoptProblem, options::Dict{String,Float64})
    # Itterate through keys in options 
    for key in keys(options)
        Ipopt.AddIpoptNumOption(prob, key, options[key])
    end
end
SetOptions!(wrapper::IpoptWrapper, options::Dict) = SetOptions!(wrapper.prob, options)
