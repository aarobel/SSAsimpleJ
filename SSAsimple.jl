using ForwardDiff
using NLsolve
#using SparseArrays

#Prescribed parameters in a dictionary
params = Dict();
params["b0"] = -100;
params["bx"] = -1e-3;
params["A"] = 4.227e-25;
params["year"] = 3600 * 24 * 365;
params["n"] = 3;
params["C"] = 1e6;
params["rho_i"] = 917;
params["rho_w"] = 1028;
params["g"] = 9.81;
params["B"] = params["A"]^(-1 / params["n"]);
params["m"] = 1 / params["n"];
params["accum"] = 1 / params["year"];

#Scaling parameters
params["hscale"] = 1000;
params["ascale"] = 0.1 / params["year"];
params["uscale"] =(params["rho_i"] * params["g"] * params["hscale"] * params["ascale"] / params["C"])^(1 / (params["m"] + 1));
params["xscale"] = params["uscale"] * params["hscale"] / params["ascale"];
params["tscale"] = params["xscale"] / params["uscale"];
params["eps"] = params["B"] * ((params["uscale"] / params["xscale"])^(1 / params["n"])) / (2 * params["rho_i"] * params["g"] * params["hscale"]);
params["lambda"] = 1 - (params["rho_i"] / params["rho_w"]);
params["transient"] = 0;

#Grid parameters
params["NT"] = 1000;
params["TF"] = params["year"] * 1000;
params["dt"] = params["TF"] / params["NT"];

params["N1"] = 150;
params["N2"] = 2000;
params["sigGZ"] = 0.94;
params["NX"] = params["N1"] + params["N2"];

#sigma = vcat(LinRange(0, 1, convert(Integer, params["N1"])));
sigma1=LinRange(params["sigGZ"]/(params["N1"]+0.5), params["sigGZ"], convert(Integer, params["N1"]));
sigma2=LinRange(params["sigGZ"], 1, convert(Integer, params["N2"]+1));
sigma = vcat(sigma1,sigma2[2:params["N2"]+1]);
grid = Dict("sigma" => sigma);
grid["sigma_elem"] = vcat(0,(sigma[1:params["NX"]-1] + sigma[2:params["NX"]]) ./ 2)
grid["dsigma"] = diff(grid["sigma"]);

#Define bed function
function bed(x,params::Dict{Any,Any})
    return params["b0"].+params["bx"].*x;
end

#Define bed deriv
function dbdx(x,params::Dict{Any,Any})
    return params["bx"].*ones(Float64,size(x));
end

#Define function with schoof velocity
function u_schoof(hg,params::Dict{Any,Any})
    us = (((params["A"]*(params["rho_i"]*params["g"])^(params["n"]+1) * params["lambda"]^params["n"])/(4^params["n"] * params["C"]))^(1/(params["m"]+1)))*(hg^(((params["m"]+params["n"]+3)/(params["m"]+1))-1));
end

#calculate function with h, u, xg
function flowline!(F, varin, varin_old, params, grid, bedfun::Function)
    #Unpack grid
    NX = params["NX"];
    N1 = params["N1"];
    dt = params["dt"]/params["tscale"];
    ds = grid["dsigma"];
    sigma = grid["sigma"];
    sigma_elem = grid["sigma_elem"];

    #Unpack parameters
    xscale = params["xscale"];
    hscale = params["hscale"];
    lambda = params["lambda"];
    m      = params["m"];
    n      = params["n"];
    a      = params["accum"]/params["ascale"];
    eps    = params["eps"];
    transient = params["transient"];

    #Unpack variables
    h = varin[1:NX]
    u = varin[NX+1:2*NX]
    xg = varin[2*NX+1]

    h_old = varin_old[1:NX];
    xg_old = varin_old[2*NX+1];

    #Calculate bed
    hf = -bedfun(xg*xscale,params)/(hscale*(1-lambda));
    hfm = -bedfun(xg*last(sigma_elem)*xscale,params)/(hscale*(1-lambda));
    b =  -bedfun(xg.*sigma.*xscale,params)./hscale;

    #Calculate thickness functions
    F[1]      = transient*(h[1]-h_old[1])/dt + (2*h[1]*u[1])/(ds[1]*xg) - a;
    F[2]      = transient*(h[2]-h_old[2])/dt -
                    transient*sigma_elem[2]*(xg-xg_old)*(h[3]-h[1])/(2*dt.*ds[2]*xg) +
                        (h[2]*(u[2]+u[1]))/(2*xg*ds[2]) - a;
    F[3:NX-1] = transient.*(h[3:NX-1] .- h_old[3:NX-1])./dt .-
                    transient.*sigma_elem[3:NX-1].*(xg - xg_old).*(h[4:NX].-h[2:NX-2])./(2 .* dt .* ds[3:NX-1] .* xg) .+
                        (h[3:NX-1] .* (u[3:NX-1] .+ u[2:NX-2]) .- h[2:NX-2] .* (u[2:NX-2] .+ u[1:NX-3]))./(2 .* xg .* ds[3:NX-1]) .- a;
#    F[N1] = transient.*(h[N1] .- h_old[N1])./dt .-
#                    transient.*sigma[N1].*(xg - xg_old).*(h[N1].-h[N1-1])./(dt .* ds[N1] .* xg) .+
#                        (h[N1] .* (u[N1] .+ u[N1-1]) .- h[N1-1] .* (u[N1-1] .+ u[N1-2]))./(2 .* xg .* ds[N1]) .- a;
    F[N1] = (1+0.5*(1+(ds[N1]/ds[N1-1])))*h[N1] - 0.5*(1+(ds[N1]/ds[N1-1]))*h[N1-1] - h[N1+1];
    F[NX]     = transient*(h[NX]-h_old[NX])/dt -
                    transient*sigma[NX]*(xg-xg_old)*(h[NX]-h[NX-1])/(dt*ds[NX-1]*xg) +
                        (h[NX]*(u[NX]+u[NX-1]) - h[NX-1]*(u[NX-1]+u[NX-2]))/(2*xg*ds[NX-1]) - a;

    #Calculate velocity functions
    F[NX+1]      = (4*eps/(xg*ds[1])^((1/n)+1))*(h[2]*(u[2]-u[1])*abs(u[2]-u[1])^((1/n)-1) -
                  h[1]*(2*u[1])*abs(2*u[1])^((1/n)-1)) - u[1]*abs(u[1])^(m-1) -
                  0.5*(h[1]+h[2])*(h[2]-b[2]-h[1]+b[1])/(xg*ds[1]);
    F[NX+2:2*NX-1] = (4 .* eps ./(xg .* ds[2:NX-1]).^((1/n)+1)) .* (h[3:NX] .* (u[3:NX] .- u[2:NX-1]) .* abs.(u[3:NX].-u[2:NX-1]).^((1/n)-1) .-
                  h[2:NX-1] .* (u[2:NX-1] .- u[1:NX-2]) .* abs.(u[2:NX-1] .- u[1:NX-2]).^((1/n)-1)) .-
                  u[2:NX-1] .* abs.(u[2:NX-1]).^(m-1) .- 0.5 .* (h[2:NX-1] .+ h[3:NX]) .* (h[3:NX] .- b[3:NX] .- h[2:NX-1] .+ b[2:NX-1])./(xg .* ds[2:NX-1]);
    F[NX+N1] = (u[N1+1]-u[N1])/ds[N1] - (u[N1]-u[N1-1])/ds[N1-1];
    F[2*NX]     = (1/(xg*ds[NX-1])^(1/n))*(abs(u[NX]-u[NX-1])^((1/n)-1))*(u[NX]-u[NX-1]) - lambda*hf/(8*eps);

    #Calculate grounding line functions
    F[2*NX+1]        = 3*h[NX] - h[NX-1] - 2*hf;
end


#Initial steady-state calculation
xg = 200e3/params["xscale"];
hf = (-bed(xg*params["xscale"],params)/params["hscale"])/(1-params["lambda"]);
h  = 1 .- (1-hf).*grid["sigma"];
u  = 0.3.*(grid["sigma_elem"].^(1/3)) .+ 1e-3;

huxg_old = vcat(h,u,xg);
huxg     = huxg_old;

#Use automatic differentiation to calculate Jacobian for nonlinear solver (more accurate and faster than finite difference!)
f = varin -> (F = fill(zero(promote_type(eltype(varin), Float64)), 2*params["NX"]+1); flowline!(F, varin, huxg_old, params, grid, bed); return F)
J = fill(0.0, 2*params["NX"]+1, 2*params["NX"]+1);
Jf! = (J,varin) -> ForwardDiff.jacobian!(J, f, varin);
@time huxg_out=nlsolve((F,varin) ->flowline!(F, varin, huxg_old, params, grid, bed), Jf!, huxg);

#@time huxg_out=nlsolve((F,varin) ->flowline!(F, varin, huxg_old, params, grid, bed), huxg);

#flowline!(F, huxg_out.zero, huxg_old, params, grid, bed);
u_bl = u_schoof(-bed(last(huxg_out.zero)*params["xscale"],params)/(1-params["lambda"]),params);
u_num=huxg_out.zero[2*params["NX"]]*params["uscale"];
err = 100*(u_num-u_bl)/u_bl
