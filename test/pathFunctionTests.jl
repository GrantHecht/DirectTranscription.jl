using DirectTranscription
using LinearAlgebra
using SparseArrays

# Simple keplarian dynamics with thrust acceleration ODE function
function keplar!(out, x, u, p, t)
    # Earth gravitational parameter (likely would never be a static parameter,
    # but using here for testing purposes...)
    μ = p[1]

    # Requirements
    @views r    = norm(x[1:3])

    # Derivatives
    @views out[1:3] .= x[4:6]
    @views out[4:6] .= -(μ / r^3).*x[1:3] .+ u
    return nothing
end

# Analytical Jacobians
function stateJac!(out, x, u, p, t)
    μ = p[1]
    out .= 0.0

    # ddr/dv 
    for i in 1:3
        out[i, 3 + i] = 1.0
    end

    # ddv/dr
    @views r = norm(x[1:3])
    @views mul!(out[4:6, 1:3], x[1:3], transpose(x[1:3]))
    @views out[4:6, 1:3] .*= (3.0 / r^2)
    for i in 1:3; out[3 + i, i] -= 1.0; end
    @views out[4:6, 1:3] .*= (μ/r^3)
    return nothing
end

function controlJac!(out, x, u, p, t)
    out .= 0.0
    for i in 1:3
        out[3 + i, i] = 1.0
    end
    return nothing
end

function staticJac!(out, x, u, p, t)
    out .= 0.0
    @views r = norm(x[1:3])
    @views out[4:6] .= -(1/r^3).*x[1:3]
end

# ===== Construct Path Functions
nFuncs      = 6
nStates     = 6
nControls   = 3
nStatic     = 1
stateSP     = sparse([1, 2, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6],
                     [4, 5, 6, 1, 2, 3, 1, 2, 3, 1, 2, 3],
                     [true for _ in 1:12])
controlSP   = sparse([4,5,6],[1,2,3],[true for _ in 1:3])
staticSP    = sparse([4,5,6],[1,1,1],[true for _ in 1:3])
adpf        = PathFunction(Dynamics(), keplar!, nFuncs, nStates, nControls, nStatic)
anpf        = PathFunction(Dynamics(), keplar!, stateJac!, controlJac!, staticJac!, nothing,
                nFuncs, nStates, nControls, nStatic, stateSP, controlSP, staticSP, sparse([],[],[]))

# ===== Test function evaluation
# Create inputs 
adout       = Vector{Float64}(undef,DirectTranscription.GetNumberOfFunctions(adpf))
anout       = Vector{Float64}(undef,DirectTranscription.GetNumberOfFunctions(anpf))
states      = rand(DirectTranscription.GetNumberOfStates(adpf))
controls    = rand(DirectTranscription.GetNumberOfControls(adpf))
static      = rand(DirectTranscription.GetNumberOfStatics(adpf))
time        = 0.0

# Test that nuber of parameters are correct
@test length(adout)     == nFuncs
@test length(states)    == nStates
@test length(controls)  == nControls
@test length(static)    == nStatic

# Evaluate functions with method
DirectTranscription.EvaluateFunction(adpf, adout, states, controls, static, time)
DirectTranscription.EvaluateFunction(anpf, anout, states, controls, static, time)

# Evaluate wrapped function
outWrapped  = Vector{Float64}(undef, nFuncs)
keplar!(outWrapped, states, controls, static, time)

# Test for equality
for i in 1:nFuncs 
    @test adout[i] == outWrapped[i]; 
    @test anout[i] == outWrapped[i];
end

# ===== Test Sparsity Detection and Jacobian evaluation
# State Jacobian
(rs,cs,_)       = findnz(stateSP)
adstateJac      = zeros(6,6)
anstateJac      = zeros(6,6)
adstateJacSP    = sparse(rs,cs,ones(length(rs)))
anstateJacSP    = sparse(rs,cs,ones(length(rs)))
stateJacWrapped = zeros(6,6)
DirectTranscription.EvaluateJacobian(DirectTranscription.State(), 
    adpf, adstateJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.State(),
    anpf, anstateJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.State(),
    adpf, adstateJacSP, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.State(),
    anpf, anstateJacSP, states, controls, static, time)
stateJac!(stateJacWrapped, states, controls, static, time)
for col in 1:6
    for row in 1:6
        @test adpf.stateSP[row,col] == stateSP[row,col]
        @test adstateJac[row, col] ≈ stateJacWrapped[row, col]
        @test anstateJac[row, col] ≈ stateJacWrapped[row, col]
        @test adstateJacSP[row, col] ≈ stateJacWrapped[row, col]
        @test anstateJacSP[row, col] ≈ stateJacWrapped[row, col]
    end
end

# Control Jacobian
(rs,cs,_)       = findnz(controlSP)
adcontrolJac    = zeros(6,3)
ancontrolJac    = zeros(6,3)
adcontrolJacSP  = sparse(rs,cs,ones(length(rs)))
ancontrolJacSP  = sparse(rs,cs,ones(length(rs)))
controlJacWrap  = zeros(6,3)
DirectTranscription.EvaluateJacobian(DirectTranscription.Control(),
    adpf, adcontrolJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Control(),
    anpf, ancontrolJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Control(),
    adpf, adcontrolJacSP, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Control(),
    anpf, ancontrolJacSP, states, controls, static, time)
controlJac!(controlJacWrap, states,  controls, static, time)
for col in 1:3
    for row in 1:6
        @test adpf.controlSP[row,col] == controlSP[row,col]
        @test adcontrolJac[row, col] ≈ controlJacWrap[row, col]
        @test ancontrolJac[row, col] ≈ controlJacWrap[row, col]
        @test adcontrolJacSP[row, col] ≈ controlJacWrap[row, col]
        @test ancontrolJacSP[row, col] ≈ controlJacWrap[row, col]
    end
end

# Static Jacobian
(rs,cs,_)       = findnz(staticSP)
adstaticJac     = zeros(6,1)
anstaticJac     = zeros(6,1)
adstaticJacSP   = sparse(rs,cs,ones(length(rs)))
anstaticJacSP   = sparse(rs,cs,ones(length(rs)))
staticJacWrap   = zeros(6,1)
DirectTranscription.EvaluateJacobian(DirectTranscription.Static(),
    adpf, adstaticJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Static(),
    anpf, anstaticJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Static(),
    adpf, adstaticJacSP, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Static(),
    anpf, anstaticJacSP, states, controls, static, time)
staticJac!(staticJacWrap, states, controls, static, time)
for row in 1:6
    @test adpf.staticSP[row, 1] == staticSP[row, 1]
    @test adstaticJac[row, 1] ≈ staticJacWrap[row, 1]
    @test anstaticJac[row, 1] ≈ staticJacWrap[row, 1]
    @test adstaticJacSP[row, 1] ≈ staticJacWrap[row, 1]
    @test anstaticJacSP[row, 1] ≈ staticJacWrap[row, 1]
end

# Time Jacobian
adtimeJac       = zeros(6,1)
antimeJac       = zeros(6,1)
DirectTranscription.EvaluateJacobian(DirectTranscription.Time(),
    adpf, adtimeJac, states, controls, static, time)
DirectTranscription.EvaluateJacobian(DirectTranscription.Time(),
    anpf, antimeJac, states, controls, static, time)
for row in 1:6
    @test adtimeJac[row, 1] ≈ 0.0
    @test antimeJac[row, 1] ≈ 0.0
end