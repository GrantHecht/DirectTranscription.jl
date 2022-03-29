_pfts = Union{PathFunction, Nothing}
struct PathFunctionSet{DF <: _pfts,CF <: _pfts,AF <: _pfts}
    # Flags indicating contents of set 
    hasDynamics::Bool 
    hasCost::Bool 
    hasAlgebraic::Bool

    # Path functions
    dynFuncs::DF
    costFuncs::CF
    algFuncs::AF

    # Number of state, control, and static parameters
    nStates::Int
    nControls::Int
    nStatic::Int

    function PathFunctionSet(funcs...) 
        # Get number of arguments
        n = length(funcs)

        # Loop through funcs to get function types
        dynFuncsIdx     = 0
        costFuncsIdx    = 0
        algFuncsIdx     = 0
        nStates         = 0
        nControls       = 0
        nStatic         = 0
        numDynFuncs     = 0
        numAlgFuncs     = 0
        numCostFuncs    = 0
        for i = 1:n
            if funcs[i] isa PathFunction
                # Get function indicies
                if GetFunctionType(funcs[i]) <: Dynamics 
                    dynFuncsIdx = i
                    numDynFuncs += 1
                elseif GetFunctionType(funcs[i]) <: Algebraic 
                    algFuncsIdx = i 
                    numAlgFuncs += 1
                elseif GetFunctionType(funcs[i]) <: Cost 
                    costFuncsIdx = i
                    numCostFuncs += 1
                    # Check that integral cost function only has a single function
                    if GetNumberOfFunctions(funcs[i]) > 1
                        error("Integral cost function cannot be a vector function.")
                    end
                end

                # Check that there is only one of each function type
                if numDynFuncs > 1 || numAlgFuncs > 1 || numCostFuncs > 1
                    if numDynFuncs > 1
                        str = "dynamics"
                    elseif numAlgFuncs > 1
                        str = "algebraic"
                    else
                        str = "cost"
                    end
                    error("Path function set can only contain one of each function type. Extra " * 
                        str * " functions passed to path function set constructor.")
                end
                
                # Get and verify number of input parameters
                if i == 1
                    nStates     = GetNumberOfStates(funcs[i])
                    nControls   = GetNumberOfControls(funcs[i])
                    nStatic     = GetNumberOfStatics(funcs[i])
                else
                    if nStates != GetNumberOfStates(funcs[i])
                        error("All functions in path function set must accept the same number of state variables.")
                    end
                    if nControls != GetNumberOfControls(funcs[i])
                        error("All functions in path function set must accept the same number of control variables.")
                    end
                    if nStatic != GetNumberOfStatics(funcs[i])
                        error("All functions in path function set must accept the same number of static variables.")
                    end
                end
            else
                error("Arguments to PathFunctionSet must all be of type PathFunction.")
            end
        end

        # Set flags
        hasDynamics     = dynFuncsIdx != 0
        hasCost         = costFuncsIdx != 0
        hasAlgebraic    = algFuncsIdx != 0
        
        # Grab individual functions
        if hasDynamics 
            dynFuncs = funcs[dynFuncsIdx]
        else
            #dynFuncs = nothing
            error("Path function set must contain dynamics functions.")
        end
        if hasCost 
            costFuncs = funcs[costFuncsIdx]
        else
            costFuncs = nothing 
        end
        if hasAlgebraic
            algFuncs = funcs[algFuncsIdx]
        else
            algFuncs = nothing
        end

        # Instantiate object 
        dft = typeof(dynFuncs)
        cft = typeof(costFuncs)
        aft = typeof(algFuncs)
        new{dft,cft,aft}(hasDynamics, hasCost, hasAlgebraic, dynFuncs, costFuncs, algFuncs, nStates, nControls, nStatic)
    end 
end

# Methods to get number of variables
GetNumberOfStates(pfs::PathFunctionSet)     = pfs.nStates
GetNumberOfControls(pfs::PathFunctionSet)   = pfs.nControls
GetNumberOfStatics(pfs::PathFunctionSet)    = pfs.nStatic

# Methods to check if function set has dynamics, algebraic, or cost functions
HasDynamicsFunctions(pfs::PathFunctionSet)  = pfs.hasDynamics
HasAlgebraicFunctions(pfs::PathFunctionSet) = pfs.hasAlgebraic
HasCostFunctions(pfs::PathFunctionSet)      = pfs.hasCost

# Methods to get number of dynamics, algebraic, or cost functions
GetNumberOfDynamicsFunctions(pfs::PathFunctionSet)  = HasDynamicsFunctions(pfs) ? GetNumberOfFunctions(pfs.dynFuncs) : 0
GetNumberOfAlgebraicFunctions(pfs::PathFunctionSet) = HasAlgebraicFunctions(pfs) ? GetNumberOfFunctions(pfs.algFuncs) : 0
GetNumberOfCostFunctions(pfs::PathFunctionSet)      = HasCostFunctions(pfs) ? GetNumberOfFunctions(pfs.costFuncs) : 0

# Methods to get dynamics, algebraic, or cost functions
function GetDynamicsFunctions(pfs::PathFunctionSet) 
    HasDynamicsFunctions(pfs) || error("Path function set does not contain dynamics functions.")
    return pfs.dynFuncs
end
function GetAlgebraicFunctions(pfs::PathFunctionSet)
    HasAlgebraicFunctions(pfs) || error("Path function set does not contain dynamics functions.")
    return pfs.algFuncs
end
function GetCostFunctions(pfs::PathFunctionSet)
    HasCostFunctions(pfs) || error("Path function set does not contain cost functions.")
    return pfs.costFuncs
end