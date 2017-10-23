module OSQPSolver
using Formatting
export solveOSQP, qpResult, qpSettings, test


# -------------------------------------
# HELPER FUNCTIONS
# -------------------------------------
  # compute projection of x onto a box defined by l and u
  function project(x::Array{Float64},l::Array{Float64},u::Array{Float64})
    return min.( max.(x,l), u)
  end

# -------------------------------------
# TYPE DEFINITIONS
# -------------------------------------
  type qpResult
    x::Array{Float64}
    y::Array{Float64}
    cost::Float64
    iter::Int64
  end

  type qpSettings
    rho::Float64
    sigma::Float64
    alpha::Float64
    eps_abs::Float64
    eps_rel::Float64
    eps_prim_inf::Float64
    eps_dual_inf::Float64
    max_iter::Int64
    verbose::Bool

    #constructor
    function qpSettings(;rho=0.1,sigma=10e-6,alpha=1.5,eps_abs=1e-3,eps_rel=1e-3,eps_prim_inf=1e-5,eps_dual_inf=1e-5,max_iter=2500,verbose=false)
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose)
    end
  end



# SOLVER ROUTINE
# -------------------------------------
  function solveOSQP(P::Array{Float64},q::Array{Float64},A::Array{Float64},l::Array{Float64},u::Array{Float64},settings::qpSettings)

    #Load algorithm settings
    σ = settings.sigma
    α = settings.alpha
    ρ = settings.rho
    ϵ_abs = settings.eps_abs
    ϵ_rel = settings.eps_rel
    ϵ_prim_inf = settings.eps_prim_inf
    ϵ_dual_inf = settings.eps_dual_inf

    # instantiate variables
    iter = 0
    r_prim = 100
    r_dual = 100

    # determine size of decision variables
    m = size(A,1)
    n = size(A,2)
    x = zeros(n)
    y = zeros(m)
    z = zeros(m)

    xhat = copy(x)
    zhat = copy(z)
    cost = Inf

    # create initial guess vectors
    xPrev = zeros(n,1)
    yPrev = zeros(m,1)
    zPrev = zeros(m,1)
    nuNew = zeros(m,1)
    xNew = zeros(n,1)
    yNew = zeros(m,1)
    zNew = zeros(m,1)

    # print information about settings to the screen
    println("-"^50 * "\n" * " "^8 * "OSQP Implementation in pure Julia\n" * " "^18 * "Michael Garstka\n"  * " "^8 * "University of Oxford, October 2017\n" * "-"^50 * "\n")
    println("Problem:  variables n = $(n), constraints m = $(m)")
    println("Settings: ϵ_abs = $(ϵ_abs), ϵ_rel = $(ϵ_rel),\n" * " "^10 * "ϵ_prim_inf = $(ϵ_prim_inf), ϵ_dual_inf = $(ϵ_dual_inf),\n" * " "^10 * "ρ = $(ρ), σ = $(σ), α = $(α),\n" * " "^10 * "max_iter = $(settings.max_iter)\n\n")

    tic()

    # KKT matrix M
    M = [P+σ*eye(n) A';A -1/ρ*eye(m)]
    M = sparse(M)

    # Do LDLT Factorization: A = LDL^T
    F = ldltfact(M)

    for iter = 1:1:settings.max_iter

      # assign previous variables
      xPrev = xNew
      yPrev = yNew
      zPrev = zNew

      # construct right hand side
      b = [σ*xPrev-q; zPrev - 1/ρ*yPrev]

      #solve linear system M*k = b with help of factorization matrix
      k = F\b

      #deconstruct solution vector k = [xt_(k+1);nu_(k+1)]
      xt = k[1:n]
      nuNew = k[n+1:end]
      zt = zPrev + 1/ρ * (nuNew - yPrev)

      # Projection steps
      xNew = α*xt + (1-α)*xPrev
      zNew = project( α*zt + (1-α)*zPrev + 1/ρ*yPrev,l,u)

      # update dual variable
      yNew = yPrev + ρ* (α*zt + (1-α)*zPrev - zNew)

      # update cost
      cost = (1/2 * xNew'*P*xNew + q'*xNew)[1]

      # compute residuals to check for termination condition
      r_prim = norm(A*xNew - zNew)
      r_dual = norm(P*xNew + q + A'*yNew)

      # compute deltas
      δx = xNew - xPrev
      δy = yNew - yPrev

      # print iterations steps
      # TODO: Might be slow
      if settings.verbose
        if iter == 1
          println("Iter:\tObjective:\tPrimal Res:\tDual Res:")
        end
        if mod(iter,100) == 0 || iter == 1 || iter == settings.max_iter
          printfmt("{1:d}\t{2:.2f}\t\t{3:.4f}\t\t{4:.4f}", iter,cost,r_prim,r_dual)
       end
      end
      # check primal infeasibility (2-norm or inf-norm?)
    #   norm_δy = norm(δy)
    #   if norm_δy > ϵ_prim_inf^2
    #    δy = δy/norm_δy
    #    # second condition
    #    if (u'*max(δy,0) + l'*min(δy,0) ) <= - ϵ_prim_inf
    #     # first condition
    #     # FIXME: isnt there a *norm(δy) missing?
    #     if norm(A'*δy) <= ϵ_prim_inf
    #       cost = Inf
    #       xNew = NaN*ones(n,1)
    #       yNew = NaN*ones(m,1)
    #       break
    #     end
    #   end
    # end

    # #check dual infeasibility
    # norm_δx = norm(δx)
    # if norm_δx > ϵ_dual_inf^2
    #   δx = δx/norm_δx
    #   if q'*δx < - ϵ_dual_inf
    #     if norm(P*δx) < ϵ_dual_inf
    #       Aδx = A * δx

    #       for i = 1:m
    #           if (u[i] < 1e18) && (Aδx[i] > ϵ_dual_inf) || (l[i] > -1e18) && (Aδx[i] < - ϵ_dual_inf)
    #             break
    #           end
    #       end
    #       cost = -Inf
    #       xNew = NaN*ones(n,1)
    #       yNew = NaN*ones(m,1)
    #       break
    #     end
    #   end
    # end



      # check convergence with residuals
      ϵ_pri = ϵ_abs *sqrt(m) + ϵ_rel * norm(zNew)
      ϵ_dual = ϵ_abs *sqrt(n) + ϵ_rel * ρ * norm(A'*yNew)
      if ( r_prim < ϵ_pri && r_dual < ϵ_dual  )
        break
      end
    end

    # print solution to screen
    rt = toq()
    println("\n\n" * "-"^50 * "\nRESULT: \nTotal Iterations: $(iter), Cost: $(cost), x: $(xNew), y: $(yNew),\nPrimal Res: $(r_prim), Dual Res: $(r_dual)\nRuntime: $(rt)s ( $(rt*1000)ms)\n" * "-"^50 )

    # create result object
    result = qpResult(xNew,yNew,cost,iter);

    return result;

  end

end




