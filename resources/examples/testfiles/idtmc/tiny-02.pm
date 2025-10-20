dtmc

module tiny02
    // state variable
    s : [0..3];

    // transitions
    [] s=0 -> [1, 1] : (s'=2);                                  // from 0 to 2 with prob 1
    [] s=1 -> [1, 1] : (s'=3);                                  // from 1 to 3 with prob 1
    [] s=2 -> [0.3, 0.3] : (s'=0) + [0.7, 0.7] : (s'=2);          // [0.3,0.3] and [0.7,1] ⇒ 0.3, 0.7
    [] s=3 -> [0.2, 0.2] : (s'=1) + [0.8, 0.8] : (s'=3);          // [0.2,0.2] and [0.8,1] ⇒ 0.2, 0.8
endmodule

init (s = 0) | (s = 1) endinit

// atomic propositions
label "a" = (s=2) | (s=3);
label "b" = (s=0) | (s=1);
