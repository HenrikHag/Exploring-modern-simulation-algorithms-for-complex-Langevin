module UsersGuide

using Plots, Statistics

export Plots, LastRowFromFile, PlotTwoPointCorrelation, PlotExp, plot_x
export GetColumn, GetLastMean, GetTwoPointData, GetTP1data, Err

# Creating functions to work on results from User's Guide to M C Methods

"""Return elements separated by "," in last row of file as Vector{Float64}  
O(1) last row file lookup
"""
function LastRowFromFile(filename)
    # row = split(last(readlines(filename)),",")
    row = []
    open(filename) do f
        # a = first(Iterators.reverse(eachline(f)))
        # println(iterate(a))
        seekend(f)
        seek(f,position(f)-2)
        while Char(peek(f)) != '\n'
            seek(f,position(f)-1)
        end
        Base.read(f, Char)
        row = append!(row, split(Base.read(f, String),","))
    end
    rowData = Vector{Float64}(undef,length(row))
    for i = 1:length(row)
        # println(parse.(Float64,row[i]))
        rowData[i] = parse.(Float64,row[i])
    end
    return rowData
end
# lastRow = LastRowFromFile("results/expfull.csv")

"""Returns column(s) of file delimited by ","
"""
function GetColumn(col,filename)
    if length(col) == 1
        al = Vector{Float64}(undef,0)
        for c = col
            for r = readlines(filename)
                push!(al,parse.(Float64,split(r,",")[c]))
            end
        end
        return al
    end
    all1 = Matrix{Float64}(undef,countlines(filename),length(col))
    r1 = 1
    for r = readlines(filename)
        all1[r1,:] = parse.(Float64,split(r,",")[col])
        r1 += 1
    end
    return all1
end

"""Get last row of file with means  
Returns Float64 of elements in range 2:n_tau+1
"""
function GetLastMean(meanf, n_tau)
    return parse.(Float64,split(last(readlines(meanf)),","))[2:n_tau+1]
end

"""Gets the last n_tau columns from "filename" (measuredObs-file)
"""
function GetTwoPointData(filename)
    ind = div(length(LastRowFromFile(filename))-1,4)
    return GetColumn((3*ind+2:4*ind+1),filename)#,ind,(3*ind:4*ind+1)
end

"""Gets the second to last n_tau columns from "filename" (measuredObs-file)
"""
function GetTP1data(filename)
    ind = div(length(LastRowFromFile(filename))-1,4)
    return GetColumn((2*ind+2:3*ind+1),filename)
end

"""Calculate the error  
for observables in obsf  
with mean in meanf
"""
function Err(n_tau, meanf, obsf)
    mean = GetLastMean(meanf, n_tau)
    variance = zeros(Float64,n_tau,1)
    err = zeros(Float64,n_tau,1)
    numb = 0
    @time for r = readlines(obsf)
        for i = 1:n_tau
            variance[i] += (parse.(Float64,split(r,",")[2:n_tau+1])[i] - mean[i])^2
        end
        numb += 1
        if numb%20==0
            println(numb)
        end
    end
    for i=1:n_tau
        variance[i] /= numb-1
    end
    for i = 1:n_tau
        err[i] = √(variance[i]/numb)
    end
    return err
end



"""Plots two-point correlation from array exp_xjxi, normalized by first element
"""
function PlotTwoPointCorrelation(exp_xjxi)
    # println(exp_xjxi)
    plot([(i,exp_xjxi[i]/exp_xjxi[1]) for i=1:length(exp_xjxi)], title="Two-Point Correlation", label="⟨x₁xᵢ⟩", yaxis=:log10, size=[700,500])
end

"""Plots expectationvalues
"""
function PlotExp(exp,n)
    if n == 1
        title = "Expectationvalue x"
        label = "⟨xᵢ⟩"
    elseif n == 2
        title = "Expectationvalue x²"
        label = "⟨(xᵢ)²⟩"
    elseif n == 3
        title = "Expectationvalue x, x²"
        label = ["⟨xᵢ⟩" "⟨(xᵢ)²⟩"]
        println("Plot?")
        println("Plotted?")
    else
        hline([0,mean(exp)])
        println("n != {1, 2}")
        return 1
    end
    plot!(exp, title=title, label=label)
    # plot([1:length(exp),exp], title=title, label=label)
end

"""Plots expectationvalues of x₁ or those in array "number"
"""
function plot_x(n_tau, meanf,number)
    rowData=[]
    for r = readlines(meanf)
        push!(rowData,parse.(Float64,split(r,",")[2:n_tau+1]))
    end
    N_meas = length(rowData)
    matrixData = zeros(Float64, N_meas, n_tau)
    for i=1:n_tau
        for ii=1:N_meas
            matrixData[ii,i]=rowData[ii][i]
        end
    end
    # reshape(rowData,N_meas,1)
    if number == 0
        plot(matrixData[:,1],label="⟨x⟩ₘₑₐ")
    else
        plot(matrixData[:,number[1]], label="⟨x$(number[1])⟩ₘₑₐ")
        for i = number[2:length(number)]
            println(i)
            plot!(matrixData[:,i], label="⟨x$(i)⟩ₘₑₐ")
        end
    end
    # return matrixData
    return
end
# PlotTwoPointCorrelation(LastRowFromFile("results/expfull.csv")[33:48])
# plotattr(:Plot)

# """G(Δτ)=⟨x(τ)x(τ+Δτ)⟩
# """
# function AutoCorrelation(array, n)
#     autoCorrelation = 0
#     for i = 1:n
#         autoCorrelation += array[i]
#     end
#     return autoCorrelation
# end

# lastRow = LastRowFromFile("results/expfull.csv")
# AutoCorrelation(lastRow, 4)

end