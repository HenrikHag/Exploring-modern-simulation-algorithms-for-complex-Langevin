
using Dates

export findDate

"""
Finds todays date and returns it in format: "yy.mm.dd"  
"""
function findDate()
    y = string(Dates.year(Dates.today()))[3:4]
    m = string(Dates.month(Dates.today()))
    d = string(Dates.day(Dates.today()))
    if length(m)==1 m = string("0",m) end
    if length(d)==1 d = string("0",d) end
    return string(y,".",m,".",d)
end
