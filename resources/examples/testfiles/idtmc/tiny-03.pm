// ./storm --prism ../../resources/examples/testfiles/idtmc/tiny-03.pm --prop "P=? [F \"a\"];P=? [F \"b\"]" --uncertainty-resolution min

dtmc

module tiny03
    // state variable
    s : [0..3];

    // transitions
    [] s=0 -> [0.2, 0.3] : (s'=1) + [0.2, 0.2] : (s'=2) + [0.5, 0.6] : (s'=3);
    [] s=1 -> [1, 1] : (s'=1);                            
    [] s=2 -> [0.2, 0.2] : (s'=0) + [0.3, 0.3] : (s'=1) + [0.1, 0.1] : (s'=2) + [0.4, 0.4] : (s'=3);
    [] s=3 -> [1, 1] : (s'=3);
endmodule

init (s = 0) endinit

// atomic propositions
label "a" = (s=1);
label "b" = (s=3);
