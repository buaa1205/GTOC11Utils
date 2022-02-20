function spacecraft!(D, vars, params, t; a=0)
	@unpack r, ṙ = vars
	@unpack μ = params

	r = SVector{3}(r)
	R = sqrt(sum(abs2, r))

	D.r .= ṙ
	D.ṙ .= -μ/R^3*r .+ a
	
	return D
end

function opt_spacecraft!(D, vars, p, t)
	@unpack x, λ = vars
	@unpack Γ, μ = p
	
	θ = atan(λ[5], λ[4])
	ϕ = atan(sqrt(λ[5]^2 + λ[4]^2), λ[6])
	@log θ rad2deg(θ)
	@log ϕ rad2deg(ϕ)
	
	sθ, cθ = sincos(θ)
	sϕ, cϕ = sincos(ϕ)
	@log a = Γ*SVector(sϕ*cθ, sϕ*sθ, cϕ)
	
	spacecraft!(D.x, x, p, t; a)
	
	r = x.r
	R = sqrt(sum(abs2, r))
	n = μ/R^3
	m = 3μ/R^5
	
	D.λ.r.x = n*λ.ṙ.x + m*r.x*(-λ.ṙ.x + λ.ṙ.y + λ.ṙ.z)
	D.λ.r.y = n*λ.ṙ.y + m*r.y*( λ.ṙ.x - λ.ṙ.y + λ.ṙ.z)
	D.λ.r.z = n*λ.ṙ.z + m*r.z*( λ.ṙ.x + λ.ṙ.y - λ.ṙ.z)
	D.λ.ṙ .= -λ.r
	
	return D
end

function opt_rel_spacecraft!(D, vars, p, t)
	@unpack x, λ = vars
	@unpack Γ, μ = p
	
	θ = atan(λ[5], λ[4])
	ϕ = atan(sqrt(λ[5]^2 + λ[4]^2), λ[6])
	@log θ rad2deg(θ)
	@log ϕ rad2deg(ϕ)
	
	sθ, cθ = sincos(θ)
	sϕ, cϕ = sincos(ϕ)
	@log a = Γ*SVector(sϕ*cθ, sϕ*sθ, cϕ)
	
	spacecraft!(D.x.chaser, x.chaser, p, t; a)
	spacecraft!(D.x.target, x.target, p, t)
	
	r = x.chaser.r
	R = sqrt(sum(abs2, r))
	n = μ/R^3
	m = 3μ/R^5
	
	D.λ.r.x = n*λ.ṙ.x + 3*m*r.x*(-λ.ṙ.x*r.x - λ.ṙ.y*r.y - λ.ṙ.z*r.z)
	D.λ.r.y = n*λ.ṙ.y + 3*m*r.y*(-λ.ṙ.x*r.x - λ.ṙ.y*r.y - λ.ṙ.z*r.z)
	D.λ.r.z = n*λ.ṙ.z + 3*m*r.z*(-λ.ṙ.x*r.x - λ.ṙ.y*r.y - λ.ṙ.z*r.z)
	D.λ.ṙ .= -λ.r
	
	return D
end



const opt_rel_prob = let
    ic = sim_proto_two_vehicles
	p = (; Γ=ustrip(Γ), μ=ustrip(μ))
    ODEProblem(opt_rel_spacecraft!, ic, (0.0, 10.0), p)
end

const opt_prob = let
	ic = sim_proto_one_vehicle
	p = (; Γ=ustrip(Γ), μ=ustrip(μ))
	ODEProblem(opt_spacecraft!, ic, (0.0, 10.0), p)
end

const spacecraft_prob = let
	ic = state_vec(zeros(6))
	p = (; Γ=ustrip(Γ), μ=ustrip(μ))
	ODEProblem(spacecraft!, ic, (0.0, 10.0), p)
end