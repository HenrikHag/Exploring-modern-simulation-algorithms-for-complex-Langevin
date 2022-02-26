module UsersGuide

using Plots, Statistics, StatsBase

export GetColumn, GetData, GetTwoPointData, GetTP1data, GetExpXData, GetLastMean, Err1
export Jackknife1
export LastRowFromFile, PlotExp, plot_x, PlotProbDD, PlotProbDDe, PlotTPCF, PlotAC

# Creating functions to work on results from User's Guide to M C Methods






#                                               #
#               File handling (reading)         #
#                                               #
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
        append!(row, split(Base.read(f, String),","))
    end
    rowData = Vector{Float64}(undef,length(row))
    for i = 1:length(row)
        rowData[i] = parse.(Float64,row[i])
    end
    return rowData
end
# lastRow = LastRowFromFile("results/expfull.csv")



"""Returns column(s) of file delimited by ","
"""
function GetColumn(col::Number,filename::String)
    al = Vector{Float64}(undef,0)
    for c = col
        for r = readlines(filename)
            push!(al,parse.(Float64,split(r,",")[c]))
        end
    end
    return al
end
function GetColumn(col,filename::String)
    all1 = Matrix{Float64}(undef,countlines(filename),length(col))
    r1 = 1
    for r = readlines(filename)
        all1[r1,:] = parse.(Float64,split(r,",")[col])
        r1 += 1
    end
    return all1
end


"""Get columns from "filename" throwing away the index at column 1,  
then taking a group of columns 1:n_tau corresponding to a dataset
"""
function GetData(filename,Nn_tau,n)
    ind = div(length(LastRowFromFile(filename))-1,Nn_tau)
    return GetColumn(((n-1)*ind+2:n*ind+1),filename)
end


"""Gets the last n_tau columns from "filename" (measuredObs-file)
"""
function GetTwoPointData(filename)
    return GetData(filename,4,4)
    #GetColumn((3*ind+2:4*ind+1),filename)#,ind,(3*ind:4*ind+1)
end

"""Gets the second to last n_tau columns from "filename" (measuredObs-file)
"""
function GetTP1data(filename)
    return GetData(filename,4,3)
end


"""Gets the n-th n_tau columns from "filename" (expfull/measuredObs-file)  
n = 1: ⟨x̂⟩  
n = 2: ⟨x̂²⟩  
n = 3: ⟨x₁xᵢ⟩  
Specify number (array or Int) for a specific column in the n-th n_tau column  
"""
function GetExpXData(filename)
    return GetData(filename,3,1)
end
function GetExpXData(filename, n)
    return GetData(filename,3,n)
end
function GetExpXData(filename, n, number)
    ind = div(length(LastRowFromFile(filename))-1,3)
    # if ∉(max(number),[1:ind])
    #     return ErrorException("Number not in range of 1:n_tau")
    # end
    return GetColumn(number.+((n-1)*ind+1),filename)
end


"""Get last row of file with means  
Returns Float64 of elements in range 2:n_tau+1
"""
function GetLastMean(meanf, n_tau)
    return parse.(Float64,split(last(readlines(meanf)),","))[2:n_tau+1]
end




#                                               #
#           Mean and Error Estimation           #
#                                               #
"""Calculate mean and error of elements in an array  
Pass a matrix to calculate mean and error for each column  
"""
function Err1(array1::AbstractArray)
    mean1 = mean(array1)
    err1 = sum((array1.-mean1).^2)
    # for i=1:length(array1)
    #     err1 += (array1[i]-mean1)^2
    # end
    err1 /= length(array1)*(length(array1)-1)
    return [mean(array1), √(err1)]#std(array1)/√length(array1)] #  √(var(array1)/length(array1)),
end
function Err1(matrix1::AbstractMatrix)
    err1 = Matrix{Float64}(undef,length(matrix1[1,:]),2)
    for i=1:length(matrix1[1,:])
        err1[i,:] = Err1(matrix1[:,i])
    end
    return err1
end



"""Does Jackknife analysis with `binsize = n-1`, where `n = length(array1)`  
If type is matrix it returnes a matrix of results of Jacknife analysis for each column
"""
function Jackknife1(array1::AbstractArray)
    jf = [mean(array1[2:length(array1)])]
    for i=2:length(array1)-1
        append!(jf,mean(append!(array1[1:i-1],array1[(i+1):length(array1)])))
    end
    append!(jf,mean(array1[1:length(array1)-1]))
    jf = jf.-mean(array1)
    jf = jf.^2
    jfvm = mean(jf)
    # jfvm = 0
    # for i=1:length(jf)
    #     jfvm += (jf[i]-mean(jf))^2
    # end
    jfvm *= (length(array1)-1)#/length(jf)
    # println("Done")
    return [mean(array1),√(jfvm)]
end

function Jackknife1(matrix1::AbstractMatrix)
    jf = Matrix{Float64}(undef,length(matrix1[1,:]),2)
    for i=1:length(matrix1[1,:])
        a = Jackknife1(matrix1[:,i])
        jf[i,:] = a
    end
    return jf
end






#                                               #
#               Plotting of data                #
#                                               #
"""Plots final expectationvalues  
n = 1: ⟨x̂ᵢ⟩  
n = 2: ⟨x̂ᵢ²⟩  
n = 3: ⟨x̂ᵢ⟩, ⟨x̂ᵢ²⟩  
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
        println("n != {1, 2, 3}")
        return 1
    end
    plot!(exp, title=title, label=label)
    # plot([1:length(exp),exp], title=title, label=label)
end

"""Plots running expectationvalues of x₁ or those in array "number"
n = 1: ⟨x̂ᵢ⟩
n = 2: ⟨x̂ᵢ²⟩
n = 3: ⟨x₁xᵢ⟩
"""
function plot_x(meanf)
    matrixData = GetColumn(2,meanf)
    plot(matrixData[:,1],label="⟨x⟩ₘₑₐ")
    return
end
function plot_x(meanf, n)
    matrixData = GetExpXData(meanf,n)
    plot(matrixData[:,1],label="⟨x⟩ₘₑₐ")
    return
end
function plot_x(meanf, n, number)
    matrixData = GetExpXData(meanf, n, number)
    labl = ""
    if n == 1
        labl = "⟩ₘₑₐ"
    elseif n == 2
        labl = "²⟩ₘₑₐ"
    else
        labl = "x₁⟩ₘₑₐ"
    end
    plot(matrixData[:,1], label=string("⟨x_$(number[1])",labl))
    for (i,val) = enumerate(number[2:length(number)])
        println(i+1)
        plot!(matrixData[:,i+1], label=string("⟨x_$(val)",labl))
    end
    return
end
function plot_x(meanf, n, number::Number)
    matrixData = GetExpXData(meanf, n, number)
    labl = ""
    if n == 1
        labl = "⟩ₘₑₐ"
    elseif n == 2
        labl = "²⟩ₘₑₐ"
    else
        labl = "x₁⟩ₘₑₐ"
    end
    plot(matrixData, label=string("⟨x_$(number)",labl))
    return
end


#                                       #
#      Probability density diagram      #
#                                       #
"""Plots Probability Density Diagram from data in column 2:n_tau+1  
"""
function PlotProbDD(file,incsize1)
    arr1 = GetColumn(2:Int((length(LastRowFromFile(file))-1)/4)+1,file)
    arr1 = reshape(arr1,:)
    histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:incsize1:(floor(maximum(arr1)*10)+1)/10],normed=true,xlabel="x",ylabel="|ψ_0|²")#,weights=repeat([1/length(arr1)],length(arr1))
end

"""Calculates the analytical Probability Density Diagram for the HO, and appends it to a plot  
"""
function PlotProbDDe(m,ω,ħ,range1)
    println("m = ",m,", ω = ",ω,", ħ = ",ħ)
    plot!([x for x=-range1:0.01:range1],[((m*ω/(π*ħ))^(1/4)*exp(-m*ω*x^2/(2*ħ)))^2 for x=-range1:0.01:range1],linewidth=2)
end



#                                       #
#         Two-Point Correlation         #
#                                       #
function PlotTPCF(filename)
    tpcd = Jackknife1(GetTwoPointData(filename))
    tpcr = Err1(GetTwoPointData(filename))
    println(tpcd[:,2])
    println(tpcr[:,2])
    # open("results/Twopointdata.csv","a") do file
    #     for i = 1:length(tpcd[:,1])
    #         write(file,string(i/2-1/2," ",tpcd[i,1]," ",tpcd[i,2],"\n"))
    #     end
    # end
    plot(tpcd[:,1],yerr=tpcd[:,2],yrange=[1.4*10^-3,10^2],yaxis=:log,title="Two-Point Correlation", label="⟨x₍ᵢ₊ₓ₎xᵢ⟩")
    plot!(tpcr[:,1],yerr=tpcr[:,2],yrange=[1.4*10^-3,10^2],yaxis=:log,title="Two-Point Correlation", label="⟨x₍ᵢ₊ₓ₎xᵢ⟩")
end



#                                       #
#           Auto Correlation            #
#                                       #
function PlotAC(filename,leng)
    data1 = GetData(filename,4,1)
    if leng > length(data1[:,1])
        leng = length(data1[:,1])
        println("PlotAC: Length specified to large, using length(data1[:,1]) = N_meas")
    end
    autocorrdata = transpose(StatsBase.autocor(data1,[i for i=0:leng-1]))
    jkf1 = Jackknife1(autocorrdata)
    jkf1[:,1]
    plot(jkf1[:,1],yerr=jkf1[:,2],title="AutoCorr by StatsBase package")
end


end