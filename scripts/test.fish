#!/usr/bin/env fish
mkdir tmp
for i in (seq 10)
    ./two_poly_ring_simple_c/main < ./testcases/$i.in > ./tmp/out
    if not diff -q ./testcases/$i.out ./tmp/out
        rm -r tmp
        exit -1
    end
end

rm -r tmp
