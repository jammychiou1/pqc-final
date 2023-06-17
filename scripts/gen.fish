#!/usr/bin/env fish
for i in (seq 10)
    ./scripts/gen.py > ./testcases/$i.in
    ./sage_ref/sage_ref.sage < ./testcases/$i.in > ./testcases/$i.out
end
