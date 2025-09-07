dtmc

module toy
    s : [0..3] init 0; // 0 = init, 1,2 = intermediates, 3 = final

    // initial to intermediates (deterministic choice split evenly)
    [] s=0 -> 0.5 : (s'=1) + 0.5 : (s'=2);

    // intermediates to final with interval probabilities
    [] s=1 -> [0.7,0.9] : (s'=3) + [0.1,0.3] : (s'=1);
    [] s=2 -> [0.7,0.9] : (s'=3) + [0.1,0.3] : (s'=2);

    // final absorbing
    [] s=3 -> 1 : (s'=3);
endmodule

label "intermediate" = (s=1) | (s=2);
label "final"        = (s=3);