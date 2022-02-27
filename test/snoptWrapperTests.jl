using DirectTranscription
using Snopt
#using Suppressor

function feval!(g, df, dg, x, deriv)
    # Compute f
    f = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]

    # Compute g 
    g[1] = x[1] * x[2] * x[3] * x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2

    if deriv 
        # Gradient of f 
        df[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
        df[2] = x[1] * x[4]
        df[3] = x[1] * x[4] + 1
        df[4] = x[1] * (x[1] + x[2] + x[3])

        # Jacobian of g
        # Constraint (row) 1
        dg[1] = x[2] * x[3] * x[4]  # 1,1
        dg[2] = x[1] * x[3] * x[4]  # 1,2
        dg[3] = x[1] * x[2] * x[4]  # 1,3
        dg[4] = x[1] * x[2] * x[3]  # 1,4
        # Constraint (row) 2
        dg[5] = 2 * x[1]  # 2,1
        dg[6] = 2 * x[2]  # 2,2
        dg[7] = 2 * x[3]  # 2,3
        dg[8] = 2 * x[4]  # 2,4
    end
    return f, false
end

x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]

g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

# Create SnoptWrapper 
wrapper = DirectTranscription.SnoptWrapper(feval!, x_L, x_U, g_L, g_U, 8, Matrix{Float64}(undef, 0, 0))

# Set initial guess 
DirectTranscription.SetInitialGuess!(wrapper, [1.0, 5.0, 5.0, 1.0])

# Set sparsity pattern 
DirectTranscription.SetSparsity!(wrapper, [1,1,1,1,2,2,2,2],[1,2,3,4,1,2,3,4])

# Set options 

# Optimize 
DirectTranscription.Optimize!(wrapper)

