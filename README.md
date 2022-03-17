# GTOC11Utils
## Installation
```julia
(@v1.6)> add https://github.com/buaa1205/GTOC11Utils
```

## Use
Note: λ can be typed in the Julia REPL in VSCode as `\lambda<TAB>` and ṙ can be typed as `r\dot<TAB>`.
```julia
using GTOC11Utils


# Define state vectors
chaser = [-1.0, -0.0, -0.0, -0.0, -6.324154185, -0.0]

target = [1.524, -0.0, -0.0, -0.0, 5.109916581, -0.0]


# Solve the optimization problem 
out = low_thrust_transfer(chaser, target)


# Inspect returned optimal values in the `u` property
out.u.t     # 1.6238469777682614
out.u.λ.r.x   # -2.4348070681482498e11
out.u.λ.ṙ.y   # -2.2321843456764603e11


# Run the problem with the optimal values to get the full solution
sol = run_solution(out)


# Plot the solution
chase_plot(sol; show_accel=true)    # Plot the trajectories in 3D
acceleration_plot(sol)              # Plot the accelerations vs. time




```
