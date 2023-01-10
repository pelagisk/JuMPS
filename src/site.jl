struct SiteType
    name::String
    d::Int
    operators::Dict #Dict{Symbol, Union{UniformScaling{Bool}, Array{ComplexF64, 2}}}
end

GenericType(d) = SiteType("generic", d, Dict(:I => I))

# TODO define some common types

Base.getindex(st::SiteType, label::Symbol) = st.operators[label]

mutable struct Site
    Γ::Array{ComplexF64, 3}
    Λ::Diagonal{Float64, Array{Float64, 1}}
    type::SiteType
end

# TODO prefer to use Site(Γ, Λ, type) directly
Site(Γ, Λ) = Site(Γ, Λ, GenericType(size(Γ, 2)))

get_physical_dim(site::Site) = site.type.d

get_right_bond_dim(site::Site) = size(site.Λ, 2)

get_type(site::Site) = site.type

get_operator(site::Site, key) = site.type[key]

get_Γ(site::Site) = site.Γ

get_Λ(site::Site) = site.Λ

function set_Γ(site::Site, Γ)
    site.Γ = Γ
end

function set_Λ(site::Site, Λ)
    site.Λ = Λ
end

function contract_left_Γ(U, site::Site)
    @tensor begin
        T[α, a, β] := U[α, γ] * site.Γ[γ, a, β]
    end
    return T
end

function contract_right_Γ(site::Site, U)
    @tensor begin
        T[α, a, β] := site.Γ[α, a, γ] * U[γ, β]
    end
    return T
end

function contract_left_right_Γ(U, site::Site, V)
    @tensor begin
        T[α, a, β] := U[α, γ] * site.Γ[γ, a, δ] * V[δ, β]
    end
    return T
end

get_A(site::Site) = contract_right_Γ(site, site.Λ)

function get_entanglement(site::Site)
    Λ = filter(λ -> λ > 0, get_Λ(site))
    return -sum(Λ.^2 .* log.(Λ.^2))
end
