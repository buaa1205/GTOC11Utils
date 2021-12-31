function solve_with(vars, p; alg=DEFAULT_ALG, kwargs...)
	@unpack λ, t = vars
	u0 = ComponentArray([p.chaser; p.target; λ], getaxes(p.prob.u0))
	prob = remake(p.prob; u0, tspan=(zero(t), t), p)
	return solve(prob, alg; DEFAULT_SIM_ARGS..., kwargs...)
end

function loss(vars, p)
	sol = solve_with(vars, p)
	@unpack x, λ = sol[end]
	@unpack chaser, target = x
    dr = SVector{3}(chaser.r) - SVector{3}(target.r)
    dṙ = SVector{3}(chaser.ṙ) - SVector{3}(target.ṙ)
    return 5e6*1/2*dr'dr + 1e6*1/2*dṙ'dṙ
	# return sum(abs2, SVector{6}(chaser) - SVector{6}(target))
end

"""
    low_thrust_transfer(chaser, target; kwargs...)

Calculate low thrust transfer between a chaser with initial conditions `chaser` and a target with
initial conditions `target`. Initial conditions should be given as `[rx, ry, rz, vx, vy, vz]`

## Keyword arguments
| Name | Default | Description |
|:---- |:-------|:----------- |
| `λ` | `-1e12rand(6).-5e11` | Initial guess for costate initial conditions
| `λ_lb` | `fill(-Inf,6)` | Lower bound for costate initial condition
| `λ_ub` | `fill(Inf,6)` | Upper bound for costate initial condition
| `t` | `8.0` | Initial guess for stop time
| `t_lb` | `0.0` | Lower bound for stop time
| `t_ub` | `15.0` | Upper bound for stop time
| `Γ` | `ustrip(Γ)` | Acceleration magnitude (in AU/yr)
| `μ` | `ustrip(μ)` | Sun's gravitational constant (in AU³/yr²)
| `alg`| `Fminbox(NelderMead())` | Optimization algorithm (from Optim.jl or similar)
| `autodiff` | `GalacticOptim.AutoFiniteDiff()` | Autodiff method (from GalacticOptim.jl)
"""
function low_thrust_transfer(chaser, target;
                            λ=-1e12rand(6).-5e11, λ_lb=fill(-Inf,6), λ_ub=fill(Inf,6),
                            t=8.0, t_lb=0.0, t_ub=15.0,
                            Γ=ustrip(Γ), μ=ustrip(μ),
                            opt_alg=Fminbox(NelderMead()), autodiff=GalacticOptim.AutoFiniteDiff(),
                            opt_kwargs...
                            )

	x = ComponentArray(OptInput(; λ, t))
	lb = ComponentArray(OptInput(; λ=λ_lb, t=t_lb))
	ub = ComponentArray(OptInput(; λ=λ_ub, t=t_ub))

    p = (; Γ, μ, chaser, target, prob=opt_rel_prob)
    f = OptimizationFunction(loss, autodiff)
    prob = OptimizationProblem(f, x, p; lb, ub)

    return solve(prob, opt_alg; opt_kwargs...)
end




## Everything above this is garbage

const nl_var_ax = getaxes(ComponentArray(λ=state_vec(zeros(6)), t=0.0))

function scale_to_requirements(x; r_req=10km, v_req=0.01m/s)
    r, v = x[1:3], x[4:6]
    r_scaled = unit(r_req).(r*DEFAULT_DISTANCE_UNIT)./r_req
    v_scaled = unit(v_req).(v*DEFAULT_DISTANCE_UNIT/DEFAULT_TIME_UNIT)./v_req
    return [r_scaled; v_scaled]
end

distance_metric(a, b) = sqrt(sum(abs2, scale_to_requirements(SVector{6}(a) - SVector{6}(b))))

internal_norm(u, t) = norm(scale_to_requirements(SVector{12}(u)[1:6]))

nl_fun(out, u, p) = nl_fun(ComponentArray(out, nl_var_ax), ComponentArray(u, nl_var_ax), p)
function nl_fun(out::ComponentArray, u::ComponentArray, p)
	@unpack λ, t = u
	@unpack station_state_final, asteroid_initial, min_t0, prob, α, alg = p

    tf = typeof(t)(prob.tspan[2])
    t0 = typeof(t)(clamp(t, min_t0, tf-ustrip(yr, 1d)))
    # t0 = typeof(t)(max(t, min_t0))

	# u0 = copy(prob.u0)
    # u0.x .= propagate(t, asteroid_initial)
	# u0.λ .= λ
    u0 = ComponentArray([propagate(t0, asteroid_initial); λ], getaxes(prob.u0))
	prob = remake(prob; u0=u0, tspan=(t0, tf))
	sol = solve(prob, alg; saveat=[tf], DEFAULT_SIM_ARGS...)
	# sol = solve(prob, Tsit5(); saveat=[tf])
	uf = sol[end]

    out.λ .= scale_to_requirements(station_state_final - SVector{6}(uf.x))
    # @. out.λ.r = 1e-10(1e10station_state_final[1:3] - 1e10uf.x.r)*DEFAULT_DISTANCE_UNIT |> km |> ustrip
    # @. out.λ.ṙ = (station_state_final[4:6] - uf.x.ṙ)*(DEFAULT_DISTANCE_UNIT/DEFAULT_TIME_UNIT) |> m/s |> ustrip
    out.t = (1 + uf.λ'uf.x) * α
    return out
end

get_candidate_solutions(station, asteroids::AbstractMatrix, args...; kwargs...) = get_candidate_solutions(station, collect(eachrow(asteroids)), args...; kwargs...)
function get_candidate_solutions(station, asteroids, tf, t0_guess;
                                 n_candidates=1,
                                 trans_scale=1,
                                 alg=DEFAULT_ALG,
                                 nl_alg=NLSolveJL(autoscale=false),
                                #  nl_alg=NewtonRaphson(),
                                 saveat=ustrip(DEFAULT_TIME_UNIT(1d)),
                                 kwargs...)
    ## Solve reverse problem
    Δt = tf - t0_guess
    min_t0 = tf-1.5Δt

    # Propagate the station to the final time
    station_state_final = propagate(tf, station)

    # Choose the first five final costates at random and calculate last from Hamiltonian
    λf = @SVector(rand(6)) .- 0.5
    max_i = argmax(station_state_final)
    # max_i = 6
    @set! λf[max_i] = -(1 + station_state_final[Not(max_i)]'*λf[Not(max_i)]) / station_state_final[max_i]

    # Set up problem
    uf = ComponentArray(OneVehicleSimState(x=collect(station_state_final), λ=λf))
    back_prob = remake(opt_prob; u0=uf, tspan=(tf, t0_guess))

    # Solve
    back_sol = solve(back_prob, Tsit5())

    # Get initial state and costate
    u0 = back_sol[end]
    back_station = u0.x
    λ0 = u0.λ

    forward_prob = remake(opt_prob; tspan=(t0_guess, tf))


    ## Get closest asteroids at that point
    # sorted = sort(asteroids; by=asteroid->sum(abs2, propagate(t0_guess, asteroid) - back_station))
    sorted = sort(asteroids; by=asteroid->distance_metric(propagate(t0_guess, asteroid),  back_station))
    besties = sorted[1:n_candidates]


    ## Loop over chosen asteroids
    nl_sols = map(besties) do bestie
        # u0.x = propagate(0, bestie)

        # ## Solve for the optimal trajectory
        # # Set up forward ODE problem
        # tspan = (t0, tf)
        # forward_prob = remake(back_prob; u0=u0, tspan=tspan)

        # Set up nonlinear problem
        nl_u = [λ0; t0_guess]# ComponentArray(; λ=λ0, t=t0_guess)
        nl_p = (
            station_state_final = station_state_final,
            asteroid_initial = bestie,
            min_t0 = min_t0,
            prob = forward_prob,
            α = trans_scale,
            alg = alg,
        )
        nl_prob = NonlinearProblem{true}(nl_fun, nl_u, nl_p)

        # Solve
        nl_sol = (; sol=solve(nl_prob, nl_alg; kwargs...), bestie=bestie)
        # # u0.λ = nl_sol.u.λ
        # u0.λ = nl_sol.u[1:6]
        # u0.x = propagate(nl_sol.u[end], bestie)
        # # u0
        # println((nl_sol.u[end], tf))
        # solve(remake(forward_prob; u0=u0, tspan=(nl_sol.u[end], tf)), alg; saveat=saveat)
    end

    filter!(x->x.sol[end]<tf, nl_sols)

    return map(nl_sols) do (sol, bestie)
        u0.λ = sol.u[1:6]
        u0.x = propagate(sol.u[end], bestie)
        new_prob = remake(forward_prob; u0=u0, tspan=(max(sol.u[end], min_t0), tf))
        solve(new_prob, DEFAULT_ALG; saveat=saveat, DEFAULT_SIM_ARGS...)
    end
end
