# include(Cluster)
using CSV
using JuMP, Ipopt
cd(pwd())
mutable struct Generator
    id::Int
    p_min::Float64
    p_max::Float64
    q_min::Float64
    q_max::Float64
end

mutable struct Bus
    id::Int # i
    # a::Int # a
    # v0::Float64 # VM
    # theta0::Float64 # VA
    v_max::Float64 # NVHI
    v_min::Float64 # NVLO
    Pd::Float64
    Qd::Float64
    # vK_max::Float64 # EVHI
    # vK_min::Float64 # EVLO

    # Bus(Data) = bus_constr(new(), Data)
end

mutable struct Line
    id::Int64
    from::Int64
    to::Int64
    # r::Float64
    # x::Float64
    b::Float64
    rateA::Int64
    y::Complex

end
mutable struct Network
    Generators::Array{Generator}
    # Generators_id::Array{Int64}
    Buses::Array{Bus}
    # Buses_id::Array{Int64}
    Lines::Array{Line}
    Slacks::Array{Int64}
    G::Array{Float64, 2}
    B::Array{Float64, 2}

end

branches= copy(CSV.read("9_buses/branch_9.csv"))
bus = copy(CSV.read("9_buses/bus_9.csv"))
gen = copy(CSV.read("9_buses/gen_9.csv"))
gencost = copy(CSV.read("9_buses/gencost_9.csv"))
gencost[!, :ID] = 1:3

BusData = DataFrame(bus = 1:6,
                    Pd = [ 0; 0; 0; 70; 70; 70],
                    Qd = [ 0; 0; 0; 70; 70; 70],
                    Vmin = [1.05; 1.05; 1.07; 0.95; 0.95; 0.95],
                    Vmax = [1.05; 1.05; 1.07; 1.05; 1.05; 1.05])

# Setting BranchData
BranchData = DataFrame(fbus = [1; 1; 1; 2; 2; 2; 2; 3; 3; 4; 5],
                        tbus = [2; 4; 5; 3; 4; 5; 6; 5; 6; 5; 6],
                        r = [0.1; 0.05; 0.08; 0.05; 0.05; 0.1; 0.07; 0.12; 0.02; 0.2; 0.1],
                        x = [0.2; 0.2 ;0.3; 0.25; 0.1; 0.3; 0.2; 0.26; 0.1; 0.4; 0.3],
                        b = [0.04; 0.04; 0.06; 0.06; 0.02; 0.04; 0.05; 0.05; 0.02; 0.08; 0.06],
                        rate = [40; 60; 40; 40; 60; 30; 90; 70; 80; 20; 40])

# Setting Generation Data
GenData = DataFrame(bus = [1; 2; 3],
                    Pmin = [50; 37.5; 45],
                    Pmax = [200; 150; 180],
                    Qmin = [-100; -100; -100],
                    Qmax = [100; 100; 100],
                    c2 = [0.00533; 0.00889; 0.00741],
                    c1 = [11.669; 10.333; 10.833],
                    c0 = [213.1; 200; 240])
branches = BranchData
bus = BusData
gen = GenData
generators = Generator[]

Sbase = 100
for i in 1:size(gen, 1)
    generator = Generator(gen[i, :bus], gen[i, :Pmin] / Sbase,
     gen[i, :Pmax] / Sbase, gen[i, :Qmin] / Sbase, gen[i, :Qmax] / Sbase)
    # generator = Generator(gen[i, :bus], gen[i, :Pmin],
    #  gen[i, :Pmax], gen[i, :Qmin], gen[i, :Qmax])
    push!(generators, generator)
end
buses  = Bus[]
for j in 1:size(bus, 1)
    bus_j = Bus(bus[j, :bus], bus[j, :Vmax],
    bus[j, :Vmin], bus[j, :Pd] / Sbase, bus[j, :Qd]/Sbase)
    # bus_j = Bus(bus[j, :bus_i], bus[j, :Vmax],
    #  bus[j, :Vmin], bus[j, :Pd], bus[j, :Qd])
    push!(buses, bus_j)
end

lines = Line[]
for i in 1:size(branches, 1)
    line = Line(i, branches[i, :fbus],
    branches[i, :tbus], branches[i, :b], branches[i, :rate],
    1 ./ (branches[i, :r] + im*branches[i, :x]))
    push!(lines, line)
end

# Y = zeros(Complex, size(buses, 1), size(buses, 1))
function ybus(BranchData, buses)
    NBuses = size(buses, 1)
    NLines=size(BranchData,1);
    Buses=1:NBuses;
    Lines=1:NLines;
    # Y = 1/Z
    # Branch Impedance and Admittance Parameters
    Y = zeros(Complex,NLines, NLines)
    BranchData[:z] = BranchData[:r] + im*BranchData[:x]
    BranchData[:y] = 1 ./ BranchData[:z];
    Y = zeros(Complex,NBuses, NBuses)
    #Non-diagonal Y
    bus_idx = [bus.id for bus in buses]
    for l in Lines
        fb = findall(x->x==BranchData[l,:fbus], bus_idx)[1]
        tb = findall(x->x==BranchData[l, :tbus], bus_idx)[1]
        Y[fb,tb] = -BranchData[l,:y]
        Y[tb,fb] = -BranchData[l,:y]
    end
    # Diagonal elements
    for i = Buses
        for j = Lines
            fb = findall(x->x==BranchData[j,:fbus], bus_idx)[1]
            tb = findall(x->x==BranchData[j, :tbus], bus_idx)[1]
            if (i==fb)||(i==tb)
                Y[i,i] += BranchData[j, :y]+ 0.5*im*BranchData[j,:b]

            end
        end
    end
    return Y
end
#now we can make admitance matrix for AC or DC network with function ybus
Y=ybus(branches, buses)
network =  Network(generators, buses, lines, [1], real(Y), imag(Y))
model = nothing
model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0))
p_g = Array{JuMP.VariableRef, 1}()
q_g = Array{JuMP.VariableRef, 1}()
e_g = Array{JuMP.VariableRef, 1}()
f_g = Array{JuMP.VariableRef, 1}()
nodal_constr_p = Array{JuMP.ConstraintRef, 1}()
nodal_constr_q = Array{JuMP.ConstraintRef, 1}()
v_constr = Array{JuMP.ConstraintRef, 1}()
for i in 1:size(network.Buses, 1)
    id = network.Buses[i].id
    if id in [k.id for k in network.Generators]
        push!(p_g, @variable(model, lower_bound=network.Generators[i].p_min,
         upper_bound=network.Generators[i].p_max, base_name="p_$id"))
        push!(q_g, @variable(model, lower_bound=network.Generators[i].q_min,
         upper_bound=network.Generators[i].q_max, base_name="q_$id"))
    else
        push!(p_g, @variable(model, lower_bound=0, upper_bound=0, base_name="p_$id", start=0.0))
        push!(q_g, @variable(model, lower_bound=0, upper_bound=0, base_name="q_$id", start=0.0))
    end
    # push!(e_g, @variable(model, lower_bound=network.Buses[i].v_min,
    #  upper_bound=network.Buses[i].v_max, base_name="v_$id"))
    if id in network.Slacks
        push!(f_g, @variable(model, lower_bound=0, upper_bound=0, base_name="f_$id"))
    else
        push!(f_g, @variable(model, base_name="f_$id", start=0.0))
    end
    push!(e_g, @variable(model, base_name="e_$id", start=1.0))
    # push!(f_g, @variable(model, base_name="f_$id"))

    push!(v_constr, @NLconstraint(model, network.Buses[i].v_min^2 <=
    e_g[end]^2 + f_g[end]^2 <= network.Buses[i].v_max^2))# base_name="v_constr_$i"))
end


G = network.G
B = network.B
for i in 1:size(network.Buses, 1)
    push!(nodal_constr_p, @NLconstraint(model, p_g[i]  - network.Buses[i].Pd ==
     sum(e_g[i] * (G[i, j] * e_g[j] - B[i, j] * f_g[j] )
       + f_g[i] * (G[i, j] * f_g[j] + B[i, j] * e_g[j]) for j in 1:size(network.Buses, 1))))

    push!(nodal_constr_q, @NLconstraint(model, q_g[i]  - network.Buses[i].Qd ==
    sum(f_g[i] * (G[i, j] * e_g[j] - B[i, j] * f_g[j] )
      - e_g[i] * (G[i, j] * f_g[j] + B[i, j] * e_g[j]) for j in 1:size(network.Buses, 1))))
end


# pnodal = Array{JuMP.ConstraintRef, 1}()
# qnodal = Array{JuMP.ConstraintRef, 1}()
# for i in 1:size(network.Buses, 1)
#     push!(pnodal, @NLconstraint(model, p_g[i] - network.Buses[i].Pd = sum()))
# end

@objective(model, Min, sum(gencost[g,:c2]*(p_g[g]* Sbase)^2+gencost[g,:c1]*(p_g[g]*Sbase)+gencost[g,:c0] for g in 1:3))
result = optimize!(model)
print(termination_status(model))
print(JuMP.objective_value((model)))
print(value.(p_g))
print(value.(q_g))
print(value.(e_g))
print(value.(f_g))
# e = convert(DataFrame, rand(3,3))

# (branches[:, 2] in [4, 5]) .| (branches[:, 1] .== 2)
nodes = [1,4,9,5]
cl_branches = branches[([ x in nodes for x in branches[!, 1]]) .|
    ([ x in nodes for x in branches[!, 2]]), :]
cl_numb = collect(Set(append!([cl_branches[i, :tbus] for i in 1:size(cl_branches, 1)],
    [cl_branches[i, :tbus] for i in 1:size(cl_branches, 1)])))
cl_bus = bus[[x in cl_numb for x in bus[!, 1]], :]
cl_gen = gen[[x in nodes for x in gen[!, 1]], :]
cl_gencost = gencost[[x in nodes for x in gencost[!, :ID]], :]

cl_generators = Generator[]

Sbase = 100
for i in 1:size(cl_gen, 1)
    generator = Generator(cl_gen[i, :bus], cl_gen[i, :Pmin] / Sbase,
     cl_gen[i, :Pmax] / Sbase, cl_gen[i, :Qmin] / Sbase, cl_gen[i, :Qmax] / Sbase)
    # generator = Generator(gen[i, :bus], gen[i, :Pmin],
    #  gen[i, :Pmax], gen[i, :Qmin], gen[i, :Qmax])
    push!(cl_generators, generator)
end
cl_buses  = Bus[]
for j in 1:size(cl_bus, 1)
    bus_j = Bus(cl_bus[j, :bus_i], cl_bus[j, :Vmax],
    cl_bus[j, :Vmin], cl_bus[j, :Pd] / Sbase, cl_bus[j, :Qd]/Sbase)
    # bus_j = Bus(bus[j, :bus_i], bus[j, :Vmax],
    #  bus[j, :Vmin], bus[j, :Pd], bus[j, :Qd])
    push!(cl_buses, bus_j)
end

cl_lines = Line[]
for i in 1:size(cl_branches, 1)
    line = Line(i, cl_branches[i, :fbus],
    cl_branches[i, :tbus], cl_branches[i, :b], cl_branches[i, :rateA],
    1 ./ (cl_branches[i, :r] + im*cl_branches[i, :x]))
    push!(cl_lines, line)
end

NBuses = size(cl_buses, 1)
NLines=size(cl_branches,1);
Buses=1:NBuses;
Lines=1:NLines;
bus_idx = [bus.id for bus in cl_buses]
for l in Lines

    fb = findall(x->x==cl_branches[l,:fbus], bus_idx)[1]
    print(fb)
    tb = findall(x->x==cl_branches[l, :tbus], bus_idx)[1]
    print("elif")
    Y[fb,tb] = -BranchData[l,:y]
    Y[tb,fb] = -BranchData[l,:y]
end

cl_Y = ybus(cl_branches, cl_buses)
cl_1 =  Network(generators, buses, lines, [1], real(Y), imag(Y))

cl_1 = Network(nodes, cl_branches, 1)
