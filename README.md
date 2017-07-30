# Experiments with the Bouncy Particle Sampler

Using the BPS for specific models: experiments repo (assumes you have a working `PDMP.jl`)

## PMF

Tech docs:

```
(cd pmf/docs; pdflatex main.tex; )
```

#### How data was obtained

The ratings in `/pmf/data` are obtained from the 1M movielens data. Assume the file is called `ratings.dat`. The delimiters are changed from `::` to `,` like so:

```bash
a=ratings.dat; sed s/::/,/g $a > _tmp; mv _tmp $a;
```

The format of every line of the file is now `i,j,rij,hash`.
To make readings more efficient / faster, we produce three files one with just
the rows, one with just the columns, one with just the ratings.
(The file `splitmat.jl` is available in the `pmf/data` folder).

```bash
julia splitmat.jl ratings.dat;
```

Then this can be read

```julia
rows  = vec(readdlm("rows.csv",  Int))
cols  = vec(readdlm("cols.csv",  Int))
rates = vec(readdlm("rates.csv", Int))
```
