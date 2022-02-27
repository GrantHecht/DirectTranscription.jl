# Function to get dictionaries of concrete types 
function ConcreteDicts(options::Dict{String,Any})
    # Create dictionaries
    strDict = Dict{String,String}()
    intDict = Dict{String, Int64}()
    numDict = Dict{String, Float64}()

    # Itterate through keys in options 
    for key in keys(options)
        # Get value of key 
        value = options[key]

        # Push to correct dictionary
        if value isa String 
            push!(strDict, key => value)
        elseif value isa Integer 
            push!(intDict, key => value)
        elseif value isa Float64 
            push!(numDict, key => value)
        else
            error("Option " * key * " passed with value of invalid type " * typeof(value))
        end
    end

    return (strDict, intDict, numDict)
end

function CombineDicts(strDict::Dict{String,String},intDict::Dict{String,Int},numDict::Dict{String,Float64})
    cDict = Dict{String,Any}()
    for key in keys(strDict)
        push!(cDict, key => strDict[key])
    end
    for key in keys(intDict)
        push!(cDict, key => intDict[key])
    end
    for key in keys(numDict)
        push!(cDict, key => numDict[key])
    end

    return cDict
end