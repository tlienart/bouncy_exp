f = open("$(ARGS[1])", "r")

rows  = open("rows.csv",  "w")
cols  = open("cols.csv",  "w")
rates = open("rates.csv", "w")

for line in eachline(f)
    ss = split(line,",")
    write(rows,  ss[1]*"\n")
    write(cols,  ss[2]*"\n")
    write(rates, ss[3]*"\n")
end
close(f); close(rows); close(cols); close(rates)
