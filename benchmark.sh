date

echo "nand35-5 old vs new:"
~/storm/build/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/nand/nand.prism  -prop 's=4 & z/N<0.1' -const N=35,K=5 -bisim | grep preprocessing
~/storm/build3/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/nand/nand.prism  -prop 's=4 & z/N<0.1' -const N=35,K=5 -bisim | grep preprocessing

echo "crowds10-10 old vs new:"
~/storm/build/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/crowds/crowds.prism  --build-all-labels  -const TotalRuns=10,CrowdSize=10 -bisim --prop 'observe0>1' | grep preprocessing
~/storm/build3/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/crowds/crowds.prism  --build-all-labels  -const TotalRuns=10,CrowdSize=10 -bisim --prop 'observe0>1' | grep preprocessing

echo "herman15 old vs new:"
~/storm/build/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/herman/herman.15.prism --build-all-labels -bisim | grep preprocessing
~/storm/build3/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/herman/herman.15.prism --build-all-labels -bisim | grep preprocessing

echo "cluster100 old vs new:"
~/storm/build/bin/storm --prism ~/git/qcomp/benchmarks/ctmc/cluster/cluster.prism -pc --build-all-labels -bisim -const N=100 | grep preprocessing
~/storm/build3/bin/storm --prism ~/git/qcomp/benchmarks/ctmc/cluster/cluster.prism -pc --build-all-labels -bisim -const N=100 | grep preprocessing

echo "philosophers16 old vs new"
~/storm/build/bin/storm --jani ~/git/qcomp/benchmarks/ctmc/philosophers/philosophers.16.jani -pc --build-all-labels -bisim -const TIME_BOUND=100 | grep preprocessing
~/storm/build3/bin/storm --jani ~/git/qcomp/benchmarks/ctmc/philosophers/philosophers.16.jani -pc --build-all-labels -bisim -const TIME_BOUND=100 | grep preprocessing

echo "haddasmonmege 5e6 old vs new"
~/storm/build/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/haddad-monmege/haddad-monmege.pm --build-all-labels -bisim -const N=5000000,p=0.5 | grep preprocessing
~/storm/build3/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/haddad-monmege/haddad-monmege.pm --build-all-labels -bisim -const N=5000000,p=0.5 | grep preprocessing

# xctrace record --template 'Time Profiler' --output /tmp/storm.trace --launch -- ~/storm/build2/bin/storm --prism ~/git/qcomp/benchmarks/dtmc/herman/herman.15.prism --build-all-labels -bisim
#
# Fri Jul 24 10:55:17 CEST 2026
# nand35-5 old vs new:
# Time for model preprocessing: 3.615s.
# Time for model preprocessing: 0.850s.
# crowds10-10 old vs new:
# Time for model preprocessing: 2.015s.
# Time for model preprocessing: 0.965s.
# herman15 old vs new:
# Time for model preprocessing: 0.423s.
# Time for model preprocessing: 0.300s.
# cluster100 old vs new:
# Time for model preprocessing: 13.762s.
# Time for model preprocessing: 0.180s.
# philosophers16 old vs new
# Time for model preprocessing: 11.441s.
# Time for model preprocessing: 0.965s.
# haddasmonmege 5e6 old vs new
# Time for model preprocessing: 4.053s.
# Time for model preprocessing: 2.980s.