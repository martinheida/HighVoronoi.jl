
function cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator::II) where {II<:Geometry_Integrator}
    return 0.0 
end

function cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator::II) where {II<:Montecarlo_Integrator}
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
    Integral = Integrator.Integral
    cdw = cell_data_writable(Integral,_Cell,dfvb,dfvvb)
    old_neighbors = cdw.neighbors
    activate_data_cell(data,_Cell,old_neighbors)

    xs = data.extended_xs
    vector = xs[_Cell]
    dim = length(vector)
    V = 0.0
    for i in 1:length(old_neighbors)
        n = old_neighbors[i]
        #hascelldata(Integrator.Integral,n)
        if !(n in calculate) && hascelldata(Integrator.Integral,n)
            neigh_data = cell_data_writable(Integrator.Integral,n,dfvb,dfvvb)
            _Cell_index = findfirstassured(_Cell,neigh_data.neighbors)
            distance=0.5*norm(vector-xs[n])
            neigh_area = neigh_data.area[_Cell_index] 
            if neigh_area>0 # otherwise pointless
                factor = distance/dim
                vol2 = 0.5*( neigh_area-cdw.area[i])
                vol = vol2*factor
                cdw.volumes[1] += vol
                neigh_data.volumes[1] -= vol
                #V+=vol
                cdw.area[i] += vol2
                neigh_data.area[_Cell_index] -= vol2
                if !(typeof(Integrator.interface)==Nothing || length(neigh_data.interface_integral)<length(neigh_data.neighbors) || Integrator.heuristic) 
                    cdw.bulk_integral .-= cdw.interface_integral[i] .* factor
                    neigh_data.bulk_integral .-= neigh_data.interface_integral[_Cell_index] .* factor
                    cdw.interface_integral[i] .+= neigh_data.interface_integral[_Cell_index]
                        cdw.interface_integral[i] .*= 0.5
                        neigh_data.interface_integral[_Cell_index] .= cdw.interface_integral[i]
                        cdw.bulk_integral .+= cdw.interface_integral[i] .* factor
                        neigh_data.bulk_integral .+= neigh_data.interface_integral[_Cell_index] .* factor
                    end
            end
        end
    end
    return V
end

function cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator::II) where {II<:Union{Polygon_Integrator,Fast_Polygon_Integrator}}
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
    Integral = Integrator.Integral
    cdw = cell_data_writable(Integral,_Cell,dfvb,dfvvb)
    old_neighbors = cdw.neighbors
    activate_data_cell(data,_Cell,old_neighbors)

    xs = data.extended_xs
    vector = xs[_Cell]
    dim = length(vector)
    V = 0.0
    for i in 1:length(old_neighbors)
        n = old_neighbors[i]
        if !(n in calculate)
            distance=0.5*norm(vector-xs[n])
            has_area_data = cdw.area[i]>0
            vol2 = has_area_data ? cdw.area[i] : get_area(Integral,n,_Cell) 
            vol = vol2*distance/dim
            cdw.volumes[1] += vol
            V+=vol
            cdw.area[i] = vol2
            if inter 
                    AA = has_area_data ? cdw.interface_integral[i] : get_integral(Integral,n,_Cell)
                    cdw.interface_integral[i] .= AA
                    cdw.bulk_integral .+= AA.*(distance/dim)
            end
        end
    end
    return V
end
function cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator::II) where {II<:Union{Heuristic_Integrator}}
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
    Integral = Integrator.Integral
    cdw = cell_data_writable(Integral,_Cell,dfvb,dfvvb)
    old_neighbors = cdw.neighbors
    activate_data_cell(data,_Cell,old_neighbors)

    xs = data.extended_xs
    vector = xs[_Cell]
    dim = length(vector)
    V = 0.0
    for i in 1:length(old_neighbors)
        n = old_neighbors[i]
        if !(n in calculate) && hascelldata(Integrator.Integral,n)
            distance=0.5*norm(vector-xs[n])
            AA = get_integral(Integral,n,_Cell)
            cdw.interface_integral[i] .= AA
            cdw.bulk_integral .+= AA.*(distance/dim)
        end
    end
    return V
end
@inline function cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator::II) where {II<:Union{HeuristicMCIntegrator}}
    cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator.mc)
    cleanup_cell(vol,ar,bulk,inter,_Cell,iterate, calculate, data,Integrator.heuristic)
end

