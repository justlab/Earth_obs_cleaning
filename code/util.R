in.sf = function(x, y, crs, region)
# Return a logical vector indicating whether each (x, y) point is in the
# specified region (an `sf` object).
    drop(st_intersects(sparse = F,
        convert.crs(cbind(x, y), crs, st_crs(region), sf = T),
        region))

pretty.table.numbers = function(d)
   {d = copy(d)
    for (col in names(d)) d[, (col) :=
       (if (is.integer(get(col)))
            scales::comma(get(col), accuracy = 1)
        else if (str_detect(col, "\\ABias"))
            sprintf("%+.03f", get(col))
        else if (is.double(get(col)))
            sprintf("%.03f", get(col))
        else
            get(col))]
    d[]}
