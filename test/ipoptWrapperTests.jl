using DirectTranscription
using Suppressor

function eval_f(x)
    return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

function eval_g!(g,x)
    g[1] = x[1] * x[2] * x[3] * x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
    return nothing
end

function eval_grad_f!(grad_f,x)
    grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
    grad_f[2] = x[1] * x[4]
    grad_f[3] = x[1] * x[4] + 1
    grad_f[4] = x[1] * (x[1] + x[2] + x[3])
    return nothing
end

function eval_jac_g!(values,rows,cols,x)
    if values === nothing
        # Constraint (row) 1
        rows[1] = 1
        cols[1] = 1
        rows[2] = 1
        cols[2] = 2
        rows[3] = 1
        cols[3] = 3
        rows[4] = 1
        cols[4] = 4
        # Constraint (row) 2
        rows[5] = 2
        cols[5] = 1
        rows[6] = 2
        cols[6] = 2
        rows[7] = 2
        cols[7] = 3
        rows[8] = 2
        cols[8] = 4
    else
        # Constraint (row) 1
        values[1] = x[2] * x[3] * x[4]  # 1,1
        values[2] = x[1] * x[3] * x[4]  # 1,2
        values[3] = x[1] * x[2] * x[4]  # 1,3
        values[4] = x[1] * x[2] * x[3]  # 1,4
        # Constraint (row) 2
        values[5] = 2 * x[1]  # 2,1
        values[6] = 2 * x[2]  # 2,2
        values[7] = 2 * x[3]  # 2,3
        values[8] = 2 * x[4]  # 2,4
    end
    return
end

n   = 4
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]

m   = 2
g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

# Create IpoptWrapper
wrapper = DirectTranscription.IpoptWrapper(eval_f, eval_g!, 
    eval_grad_f!, eval_jac_g!, n, x_L, x_U, m, g_L, g_U, 8)

# Set initial guess
DirectTranscription.SetInitialGuess!(wrapper, [1.0, 5.0, 5.0, 1.0])

# Set options 
DirectTranscription.SetFloatOption!(wrapper, "tol", 1e-10)
DirectTranscription.SetIntOption!(wrapper, "print_level", 1)
DirectTranscription.SetStringOption!(wrapper, "check_derivatives_for_naninf", "yes")

# Optimize 
@suppress DirectTranscription.Optimize!(wrapper)

# Check optimization results 
@test DirectTranscription.GetSolution(wrapper) ≈ 
    [0.9999999900091949, 4.742999643601108, 3.8211499789170844, 1.3794082932197205]

# Reset wrapper 
DirectTranscription.ResetIpoptWrapper!(wrapper, eval_f, eval_g!, 
    eval_grad_f!, eval_jac_g!, n, x_L, x_U, m, g_L, g_U, 8)

# Set initial guess
DirectTranscription.SetInitialGuess!(wrapper, [1.0, 5.0, 5.0, 1.0])

# Optimize 
@suppress DirectTranscription.Optimize!(wrapper)

# Check optimization results 
@test DirectTranscription.GetSolution(wrapper) ≈ 
    [0.9999999900091949, 4.742999643601108, 3.8211499789170844, 1.3794082932197205]

