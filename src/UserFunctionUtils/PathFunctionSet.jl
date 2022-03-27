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
        for i = 1:n
            if funcs[i] isa PathFunction
                # Get function indicies
                if GetFunctionType(funcs[i]) <: Dynamics 
                    dynFuncsIdx = i
                elseif GetFunctionType(funcs[i]) <: Algebraic 
                    algFuncsIdx = i 
                elseif GetFunctionType(funcs[i]) <: Cost 
                    costFuncsIdx = i
                end
                
                # Get and verify number of input parameters
                if i == 1
                    nStates     = GetNumberOfStates(funcs[i])
                    nControls   = GetNumberOfControls(funcs[i])
                    nStatic     = GetNumberOfStatics(funcs[i])
                else
                    if nStates != GetNumberOfStates(funcs[i])
                        error("All functions in path function set but accept the same number of state variables.")
                    end
                    if nControls != GetNumberOfControls(funcs[i])
                        error("All functions in path function set but accept the same number of control variables.")
                    end
                    if nStatic != GetNumberOfStatics(funcs[i])
                        error("All functions in path function set but accept the same number of static variables.")
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
            dynFuncs = nothing
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
