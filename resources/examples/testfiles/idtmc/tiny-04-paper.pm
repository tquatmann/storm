dtmc

module tiny04
    // state variable
    s : [0..3];

    // transitions
    [] s=0 -> [0.5, 0.69] : (s'=1) + [0.3, 0.51] : (s'=2) + [0.2, 0.79] : (s'=3);
    [] s=1 -> [0.5, 0.7] : (s'=1) + [0.3, 0.5] : (s'=2) + [0.2, 0.8] : (s'=3);
    [] s=2 -> [0.78, 1] : (s'=2) + [0.2, 0.8] : (s'=3);
    [] s=3 -> [1, 1] : (s'=3);
endmodule

init (s = 0) endinit

// atomic propositions
label "a" = (s=3);
