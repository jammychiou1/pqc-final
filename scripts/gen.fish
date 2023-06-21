#!/usr/bin/env fish
for i in (seq 10)
    ./scripts/gen.py > ./scripts/testcases/$i.in
    ./scripts/sage_ref/sage_ref.sage < ./scripts/testcases/$i.in > ./scripts/testcases/$i.out
end
