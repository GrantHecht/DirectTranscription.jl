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

    function PathFunctionSet(funcs...) 
        # Get number of arguments
        n = length(funcs)

        # Loop through funcs to get function types
        dynFuncsIdx = 0
        costFuncsIdx = 0
        algFuncsIdx = 0
        for i = 1:n
            if funcs[i] isa PathFunction
                if GetFunctionType(funcs[i]) <: Dynamics 
                    dynFuncsIdx = i
                elseif GetFunctionType(funcs[i]) <: Algebraic 
                    algFuncsIdx = i 
                elseif GetFunctionType(funcs[i]) <: Cost 
                    costFuncsIdx = i
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
        new{dft,cft,aft}(hasDynamics, hasCost, hasAlgebraic, dynFuncs, costFuncs, algFuncs)
    end 
end

