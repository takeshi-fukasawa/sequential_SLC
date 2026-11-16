using LinearAlgebra
using Distributions
using ForwardDiff
using BenchmarkTools
using CSV
using DataFrames
using Zygote
using ReverseDiff

## Parameters
n=100;
P0=Matrix{Int64}(I, n, n);
P1=Matrix{Int64}(I, n, n);

#P0=rand(Float64, n,n);
#P1=rand(Float64, n,n);

P0=P0./sum(dims=2,P0);
#P1=P1./sum(dims=2,P1);


u0=zeros(n,1);
u1=zeros(n,1);
beta_param=0.9;
V0=ones(n,1);
vec=ones(n,1);

## Function to compute Inclusive value
function Phi_func(V)
    EV0=P0*V;
    out=log.(exp.(u0.+beta_param.*EV0)+exp.(u1));
    return out
end

function Phi2_func(V)
    util0=u0.+beta_param.*(P0*V);
    out=log.(exp.(util0)+exp.(u1));
    return out
end

function jvp_func(func,primal,tangent)
    g(t)=func(primal+t*tangent);
    jvp_result=ForwardDiff.derivative(g,0.0);
    #jvp_result=ReverseDiff.gradient(g,0.0);
    
    #jvp_result=Zygote.jacobian(g,0.0)[1];
    return jvp_result
end


function g(t)
    return Phi_func(V0+t*vec);
end

function jvp_func2()
    jvp_result=ForwardDiff.derivative(g,0.0);
    #jvp_result=ReverseDiff.gradient(g,0.0);
    
    #jvp_result=Zygote.jacobian(g,0.0)[1];
    return jvp_result
end

function jvp_func3()
    #jvp_result=(ForwardDiff.jacobian(Phi_func,V0))*vec;
    jvp_result=(ReverseDiff.jacobian(Phi_func,V0))*vec;
    #jvp_result=(Zygote.jacobian(Phi_func,V0))*vec;# ten times slower??
    
    return jvp_result
end

function Jacobian_Phi_func(V)
    util_0=u0.+beta_param.*P0*V;
    denom=exp.(util_0)+exp.(u1);
    CCP0=exp.(util_0)./denom;
    CCP1=exp.(u1)./denom;

    Jacobian=CCP0.*(beta_param*P0);

    return Jacobian
end

function JVP_Phi_func(V)
    util_0=u0.+beta_param.*P0*V;
    denom=exp.(util_0)+exp.(u1);
    CCP0=exp.(util_0)./denom;
    CCP1=exp.(u1)./denom;

    Jacobian=CCP0.*(beta_param*P0);

    JVP=Jacobian*vec;
    return JVP
end

function JVP_Phi_func_no_Jacobian(V)
    util_0=u0.+beta_param.*P0*V;
    denom=exp.(util_0)+exp.(u1);
    CCP0=exp.(util_0)./denom;
    CCP1=exp.(u1)./denom;

    JVP=CCP0.*(beta_param*(P0*vec));

    return JVP
end

function jacobian_func()
    #jacobian=(ForwardDiff.jacobian(Phi_func,V0));
    jacobian=(ReverseDiff.jacobian(Phi_func,V0));
    
    return jacobian
end

function jacobian_func2()
    #jacobian=(ForwardDiff.jacobian(Phi2_func,V0));
    jacobian=(ReverseDiff.jacobian(Phi2_func,V0));

    return jacobian
end

ratio_jacobian_check=Jacobian_Phi_func(V0)[:]./jacobian_func()[:];


function jvp_numer_diff_func(func,primal,tangent,eps_val)
    jvp_result=(func(primal+eps_val*tangent)-func(primal-eps_val*tangent))./(2*eps_val);
    return jvp_result
end

## Numerical differentiation
eps_val=1e-5;

#Jac=ForwardDiff.jacobian(Phi_func,V0)
if 1==0
    @time for i=1:1000
        global jvp_numer_diff=jvp_numer_diff_func(Phi_func,V0,vec,eps_val);
    end
    #@time for i=1:1000
    #    global jvp_auto_diff=jvp_func(Phi_func,V0,vec);
    #end
    @time for i=1:1000
        global jvp_auto_diff2=jvp_func2();
    end
end
#difference=maximum(abs.(jvp_numer_diff-jvp_auto_diff2));

## Experiment 1: Compute Jacobians
println("Experiment 1")

if 1==1
    result_Jacobian_auto_diff=@benchmark Jac_auto_diff=jacobian_func();
    result_Jacobian_analy_diff=@benchmark Jac_analy_diff=Jacobian_Phi_func(V0);
end


## Experiment2: Compute JVP
println("Experiment 2")
if 1==1
    result_JVP_numer_diff=@benchmark jvp_numer_diff=jvp_numer_diff_func(Phi_func,V0,vec,eps_val);
    
    result_JVP_auto_diff=@benchmark jvp_auto_diff2=jvp_func2();
    result_JVP_auto_diff_compute_Jacobian=@benchmark jvp_auto_diff3=jvp_func3();

    result_JVP_analy_diff_compute_Jacobian=@benchmark JVP_Phi_func(V0);
    result_JVP_analy_diff=@benchmark JVP_Phi_func_no_Jacobian(V0);
    
end

function result_func(result)
    out1=[BenchmarkTools.mean(result), BenchmarkTools.std(result)]
    out1=map(x -> (x.time)/(1000^2), copy(out1));
    out=[out1;BenchmarkTools.memory(result)/(1024^2)];
    return out
end

result_Jacobian=cat(dims=2,
    result_func(result_Jacobian_auto_diff),
    result_func(result_Jacobian_analy_diff)
);
result_Jacobian=permutedims(result_Jacobian,(2,1));

result_JVP=cat(dims=2,
    result_func(result_JVP_numer_diff),
    result_func(result_JVP_auto_diff),
    result_func(result_JVP_auto_diff_compute_Jacobian),
    result_func(result_JVP_analy_diff),
    result_func(result_JVP_analy_diff_compute_Jacobian)
);
result_JVP=permutedims(result_JVP,(2,1));

if 1==1
    if n==100
        CSV.write("result_n_100.csv", DataFrame(round.(cat(dims=1,result_Jacobian,result_JVP),digits=2),:auto));
    elseif n==1000
        CSV.write("result_n_1000.csv", DataFrame(round.(cat(dims=1,result_Jacobian,result_JVP),digits=2),:auto));
    end
end
