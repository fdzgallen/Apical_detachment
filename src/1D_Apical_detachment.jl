# Standard code to run a simulation using the mechanochemical model presented in article "INSERT DOI"
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# One can work with pure local inhibition by changing the following mechanical parameters:
# vCTE=0   α=0   β=0   σₐ₀=0
# Written by Andreu F Gallen working in Turlier lab and in collaboration with Orion Weiner's lab

using DelimitedFiles
using Gridap
using Plots
using Plots.PlotMeasures
include("Plots_RhoRacA.jl")

function plotting(ylab, po, pPNG, i)
    plot(po)
    xlabel!("ξ[μm]")
    ylabel!(ylab)
    savefig(pPNG*"rho"*i*".png")
end

function conservation(sρ0, sρ, Minitial)
    return 0.1*(Minitial-(sρ0+sρ))
end

function threshold(x, x₀, xth)
    return (0.5 * (1 .+ tanh.(x/x₀ .- xth/x₀)))
end

#function to run a single simulation with a few given parameters
function run(χ, λ⁻², η, T, Δt, part)

    # Time discretisation parameters
    t₀ = 0.0
    t = t₀
    nΔt = trunc(Int, T/Δt)

    domain = (0, L)
    partition = (part)
    model = CartesianDiscreteModel(domain, partition)
    Δx = L/partition

    # Lets write down some parameters in a txt just in case 
    pVTU="./output/"*simulation*"/VTU/"
    mkpath(pVTU)
    pPNG="./output/"*simulation*"/"
    mkpath(pPNG)
    mkpath(pPNG*"ezrin_time/")
    mkpath(pPNG*"ezrin_unbound_time/")
    mkpath(pPNG*"Rac_time/")
    mkpath(pPNG*"Rho_time/")

    # Lets copy the code in the output folder to be able to check code used for each simulation
    cp(@__FILE__, pPNG*split(@__FILE__, "/")[end], force = true)


    order = 1
    #Starting Boundary conditions
    x₀(x) = 0.0
    xₗ(x) = 0.0

    v₀(x) = 0.0
    vₗ(x) = 0.0

    # Creating FE space
    reffe = ReferenceFE(lagrangian, Float64, order)
    W0 = TestFESpace(model, reffe; dirichlet_tags = ["tag_1"])
    WD0 = TestFESpace(model, reffe; dirichlet_tags = ["tag_1", "tag_2"])
    X = TrialFESpace(WD0, [x₀, xₗ])
    X2 = TrialFESpace(W0, [x₀])
    V = TrialFESpace(WD0, [v₀, vₗ])

    degree = 2
    #triangulate FE space
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    neumanntags = ["tag_2"]
    Γ = BoundaryTriangulation(model, tags = neumanntags)
    dΓ = Measure(Γ, degree)

    #Neumann and defining trial spaces and test space for ρ_b 
    Q0 = TestFESpace(model, reffe)
    RHO = TrialFESpace(Q0)
    RHO0 = TrialFESpace(Q0)

    #defining Starting conditions for some variables and dummy variables to be able to build the equations
    uh_ρb = interpolate_everywhere(kon*M0/partition/(koff+kon), RHO)
    uh_ρb_old = uh_ρb
    uh_ρu = interpolate_everywhere(koff*M0/partition/(koff+kon), RHO0)
    uh_ρu_old = uh_ρu
    tension = interpolate_everywhere(0, RHO)
    uh_x = zero(X)
    uh_x_old = uh_x
    uh_v = zero(V)
    f1(x) = 0.01*x[1]
    f2(x) = α₀+0.5*exp(-x[1]^2/wrac^2)
    Minitial=sum(get_free_dof_values(uh_ρu))+sum(get_free_dof_values(uh_ρb))

    print(
        "Start1 sum(ρ_b) "*string(sum(get_free_dof_values(uh_ρb)))*", sum(ρ_u) "*string(
            sum(get_free_dof_values(uh_ρu)),
        )*", and sum(ρ_b+ρ_u)"*string(
            sum(get_free_dof_values(uh_ρu))+sum(get_free_dof_values(uh_ρb)),
        )*"\n",
    )

    #variables to store temporal information that we would like to plot later
    tensiont = zeros(trunc(Int, T/Δt)+1, partition+1)
    vt = zeros(trunc(Int, T/Δt)+1, partition-1)
    xt = zeros(trunc(Int, T/Δt)+1, partition)
    Mt = zeros(trunc(Int, T/Δt))
    #We start wriiting down the values of some parameters at time=0
    λ = 0.0
    ρbt = zeros(trunc(Int, T/Δt)+1, partition+1)
    ρbt[1, :] = get_free_dof_values(uh_ρb)[:]
    ρut = zeros(trunc(Int, T/Δt)+1, partition+1)
    ρut[1, :] = get_free_dof_values(uh_ρu)[:]

    #we write down the weak form of the membrane equation
    # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a 
    a(k, x, w) = ∫(k*(∇(x)⋅∇(w)))dΩ
    m(Δt, x, w) = ∫((χ)*(x*w)/Δt)dΩ
    t=0
    h(x) = ten_max * threshold(t, 10, 100) * (sin((t-100)*ω) + 1) / 2 #tension oscillations

    aₓ(x, w) = ∫(k*(∇(x)⋅∇(w)))dΩ + ∫(ηₘ*((∇(x)/Δt) ⋅ ∇(w)))dΩ + m(Δt, x, w)
    bₓ(uh_v, x_old, w) =
        ∫(ηₘ * ((∇(x_old)/Δt) ⋅ ∇(w)))dΩ + ∫(w*h)*dΓ + m(Δt, x_old, w) + m(Δt, uh_v, w)

    # Now for MAC bound
    # mass term for the temporal evolution ρ_b
    mρ(Δt, ρ_b, w) = ∫((ρ_b*w)/Δt)dΩ
    aρb(ρ_b, v, w) = ∫((ρ_b*w)/Δt)dΩ + ∫(k0 * exp(-eb/kBT) * (ρ_b*w))dΩ
    bρb(w, ρ_u, ρ_b_old, ten) =
        ∫((ρ_b_old*w)/Δt)dΩ + ∫(
            k0*exp(-ten/kBT)*(ρ_u*w) +
            λ/L*(k0*exp(-ten/kBT)/(k0*exp(-ten/kBT) + k0*exp(-eb/kBT))),
        )dΩ

    #Now for MCA unbound
    # mass term for the temporal evolution ρ_b0
    aρu(ρ_u, x, x_old, w) = mρ(Δt, ρ_u, w) + ∫(k0*exp(-ten/kBT)*(ρ_u*w))dΩ
    bρu(w, ρ_b, ten) =
        ∫(
            k0 *
            exp(-eb/kBT) *
            (ρ_b*w) *
            (ρ_b*w)+λ/L*(
                k0 * exp(-eb/kBT) * (ρ_b*w)/(k0*exp(-ten/kBT)+k0 * exp(-eb/kBT) * (ρ_b*w))
            ),
        )dΩ + mρ(Δt, uh_ρu_old, w)

    #Now for v
    aᵥ(ρ_b, v, w) = ∫(η*(∇(v)⋅∇(w)))dΩ + ∫(χ*(v*w))dΩ
    bᵥ(w) = m(Δt, uh_x, w) - m(Δt, uh_x_old, w) # +  ∫( w*σₐ₀*∇(uh_rho)⋅VectorValue(1.0) )dΩ #no feedback is ∫( w*∇σₐ )dΩ
    Aᵥ(v, w) = aᵥ(uh_ρb, v, w)
    Bᵥ(w) = bᵥ(w)
    #We can now use ρ_b and x to solve v
    op_v = AffineFEOperator(Aᵥ, Bᵥ, V, WD0)

    #SOLVING X AT t=0
    Aₓ(x, w) = aₓ(x, w)
    Bₓ(w) = bₓ(uh_v, uh_x_old, w)
    op_x = AffineFEOperator(Aₓ, Bₓ, X2, W0)
    uh_x = solve(op_x)
    uh_x_old = uh_x

    ten = get_free_dof_values(tension)[1]
    #SOLVE ρ_b AT t=0 vien initial velocity zero
    Aρ(ρ_b, w) = aρb(ρ_b, uh_v, w)
    Bρ(w) = bρb(w, uh_ρu, uh_ρb_old, ten)
    op_ρ = AffineFEOperator(Aρ, Bρ, RHO, Q0)
    uh_ρb = solve(op_ρ)
    uh_ρb_old = uh_ρb

    #SOLVE ρ_u AT t=0
    Aρ0(ρ_u, w) = aρu(ρ_u, uh_x, uh_x_old, w)
    Bρ0(w) = bρu(w, uh_ρb_old, ten)
    op_ρ0 = AffineFEOperator(Aρ0, Bρ0, RHO0, Q0)
    uh_ρu = solve(op_ρ0)
    uh_ρu_old = uh_ρu

    #SOLVING V AT t=0 using x, xold and rho
    uh_v=zero(V)

    #we define the tension for a spring
    mten(u, v) = ∫((u*v))dΩ
    mten2(u, v) = ∫((k * (∇(u)⋅VectorValue(1.0))) * v)dΩ
    #sten(u,v) = ∫( 10*γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

    Aten(u, v) = mten(u, v) #+ sten(u,v)
    bten(v) = mten2(uh_x, v)
    op_ten = AffineFEOperator(Aten, bten, RHO, Q0)
    tension = solve(op_ten)

    vt[1, :] = get_free_dof_values(uh_v)[:]
    xt[1, :] = get_free_dof_values(uh_x)[:]
    tensiont[1, :] = get_free_dof_values(tension)[:]

    m_api(u, v) = ∫(kon*(M-u)*exp(-F*β)*v - koff*u*exp(-u*β)*v)dΩ
    #s(u,v) = ∫( 10*γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩ

    A_api(u, v) = m_api(u, v) #+ sten(u,v) 
    B_api(v) = m_api(xₕ, v)
    #op_api = AffineFEOperator(A_api , B_api , RHO0 , Q0)
    #c_api = solve(op_api)

    i = 0
    t=0
    dummyx0=0

    writevtk(
        Ω,
        pVTU*"VTU$i",
        cellfields = ["x"=>uh_x, "v"=>uh_v, "rho"=>uh_ρb, "rho0"=>uh_ρu],
    )
    for ti = t₀:Δt:(T-Δt)
        i1 = sum(get_free_dof_values(uh_ρb))
        i2 = sum(get_free_dof_values(uh_ρu))
        λ = conservation(i2, i1, Minitial)
        i3 = i1+i2
        i4 = trunc(t)
        i = i + 1
        t = t + Δt

        #boundary conditions for v and x 
        if i > 1
            vₗ = (xt[i, end] - xt[i-1, end])/Δt
        end

        #we introduce a slight relaxation for the membrane, decreases 2% x at the Boundary condition only
        # dummyx0=dummyx0*0.98+Δt*vCTErac/(1+ten[end]*ten[end]/ten0)
        # xₗ = dummyx0
        # X = TrialFESpace(WD0,[x₀,xₗ])

        @info "Time step $i/$nΔt, time $i4, sum(ρ_b+ρ_u) $i3"
        Mt[i]=i3

        ρbt[i+1, :] = get_free_dof_values(uh_ρb)[:]
        ρut[i+1, :] = get_free_dof_values(uh_ρu)[:]

        # Updating v to solve ρ and x
        op_x = AffineFEOperator(Aₓ, Bₓ, X2, W0)
        uh_x = solve(op_x)

        po=get_free_dof_values(uh_ρb)
        poo = get_free_dof_values(uh_ρu)

        op_ten = AffineFEOperator(Aten, bten, RHO, Q0)
        tension = solve(op_ten)

        op_ρ = AffineFEOperator(Aρ, Bρ, RHO, Q0)
        uh_ρb = solve(op_ρ)
        uh_ρb_old = uh_ρb

        op_ρ0 = AffineFEOperator(Aρ0, Bρ0, RHO0, Q0)
        uh_ρu = solve(op_ρ0)
        uh_ρu_old = uh_ρu

        V = TrialFESpace(WD0, [v₀, vₗ])
        op_v = AffineFEOperator(Aᵥ, Bᵥ, V, WD0)
        uh_v=solve(op_v)

        #updating x_old per time derivarive
        uh_x_old = uh_x
        writevtk(
            Ω,
            pVTU*"VTU$i",
            cellfields = ["x"=>uh_x, "v"=>uh_v, "rho"=>uh_ρb, "rho0"=>uh_ρu],
        )

        plotting("ρ_b", po, pPNG*"ezrin_time/", "$i")
        plotting("ρ_u", poo, pPNG*"ezrin_unbound_time/", "$i")


        vt[i+1, :] = get_free_dof_values(uh_v)[:]
        xt[i+1, :] = get_free_dof_values(uh_x)[:]
        tensiont[i+1, :] = get_free_dof_values(tension)[:]
        ten = get_free_dof_values(tension)[1]

    end
    plots_apical(nΔt, Δx, vt, xt, tensiont, ρbt, ρut, pPNG, ten_max)

    return tensiont, ρbt, vt
end


###############CONSTANT MECHANICAL PARAMETERS ##############
partition=100       #size of the FE system
χ = 10.0             # [pN s/um^3] friction coefficient
η = 100000.0               # [pN s/ um] 2D viscosity of the cortex
L = 20.0           # [um] System size length
k = 1000.0           # [pN /um] elastic constant for the membrane
σₐ₀ = 100.0        # [pN /um^2] maximum active "pressure"
ηₘ = 0.0 # [pN s/um] viscosity of the membrane 
ten_max = 100.0
ω = 0.5 #frequency of the tension oscillations

koff = 1.4  # ezrin koff [1/s] toff=0.7s Fritzsche et al
kon = 5.0  # ezrin kon [1/s] ton =0.2s Fritzsche et al
D = 0.3   # diffusion of the membrane [um^2/s] 0.003  Fritzsche et al
M0 = 10.0 #total amount of ezrin
wrac = L/5 #Opto signal use this for the Gaussian distribution width
#in case one wants nonlinear deactication if the system would be unstable 

k0 = 1.0 #unbinding rate at zero force
eb = 5.0 #energy scale for the force dependence of the unbinding rate
kBT = 1.0 #thermal energy scale for the force dependence of the unbinding rate

#adimensional numbers that define the system
τₐ = η/σₐ₀ #for σₐ₀=100
τₑ = η/k #almost always?
Pe⁻¹ = D*η/(σₐ₀*L^2) #D*η/(σₐ₀*L^2) 0.00075 for σₐ₀=100
λ = sqrt(η*L / (χ*L*L*M0)) #\bar λ in the suppl material. Adimensional lenthscale

rac0 = 0.01# ADIMENSINAlizes rac switch for protrusion, changes the slope of the threshold function
vCTE = 0.15 #Scale for the polimerizatio velocity
ten0 = 10.0 #denominator for protrusion dependence on tension
tenth = 5.0 #threshold for tension to activate rho
sig0 = 1.00  #adimensionalizes the tension term for rho
ρth = 0.068 #MCA/ezrin value below which protrusion starts
rho0 = 0.01 #adimensionalizes the MCA concentration term for rac

#time variables
topto=50 #time to start opto signal
Δt = 1.0 #timestep
T = 500 #time to finish simulation
#if we wanna do multiple simulations in a row with different variables
setsv=1
setsx=2

#output folder name
len=trunc(λ, digits = 2)
simulation = "T=$T deltat=$Δt/v4/ten_max=$ten_max k0=$k0 eb=$eb kBT=$kBT freq=$ω friction=$χ visco_memb=$ηₘ η=$η k=$k M0=$M0/"
print("Simulation name: "*simulation*"\n")

#IF we run multiple simulations with changing parameters, we use this arrays to store the results of each simulation
xχ = zeros(setsx)
xη = zeros(setsv)
tension = zeros(setsx, setsv, Int32(T/Δt)+1, partition+1)
v = zeros(setsx, setsv, Int32(T/Δt)+1, partition-1)
ρ_b = zeros(setsx, setsv, Int32(T/Δt)+1, partition+1)
a = zeros(setsx, setsv, Int32(T/Δt)+1, partition+1)
b = zeros(setsx, setsv, Int32(T/Δt)+1, partition+1)
tension2 = zeros(setsx, setsv, Int32(T/Δt)+1, partition-2)
v2 = zeros(setsx, setsv, Int32(T/Δt)+1, partition-1)
ρ2 = zeros(setsx, setsv, Int32(T/Δt)+1, partition+1)
αt = zeros(setsx, setsv, Int32(T/Δt)+1, partition)
βt = zeros(setsx, setsv, Int32(T/Δt)+1, partition)

λ⁻² = 1/(λ*λ)
run(χ, λ⁻², η, T, Δt, partition)

#FOR for the possible multiple simulations
# for i = 2:1:setsx
#     for j = 1:1:setsv
#         χ = 2*10^((i))    # [pN s/um^3] friction coefficient
#         χ = trunc(χ,digits=2)
#         λ⁻² = 1/(λ*λ)
#         xχ[i]=χ
#         η= 0+0*j+10^(1*(j+1))              # [pN s/ um] 2D viscosity of the cortex
#         lη = η
#         xη[j]=lη
#         tension[i, j, :, :], ρ_b[i, j, :, :], v[i, j, :, :] =
#             run(χ, λ⁻², lη, T, Δt, partition)
#     end
# end
