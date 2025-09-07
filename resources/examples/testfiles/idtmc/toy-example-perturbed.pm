dtmc

module toy
    s : [0..3] init 0; // 0 = init, 1,2 = intermediates, 3 = final

    // initial to intermediates (deterministic choice split evenly)
    [] s=0 -> 0.5 : (s'=1) + 0.5 : (s'=2);

    // intermediates to final with interval probabilities
    [] s=1 -> [0.71,0.89] : (s'=3) + [0.12,0.28] : (s'=1);
    [] s=2 -> [0.69,0.91] : (s'=3) + [0.08,0.32] : (s'=2);

    // final absorbing
    [] s=3 -> 1 : (s'=3);
endmodule

label "intermediate" = (s=1) | (s=2);
label "final"        = (s=3);