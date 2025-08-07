fresh

C = structor("bulk","first-fields-first");

C.str.a1.a2 = ["aa11","aa21"; "aa12","aa22"];
C.str.a1.b2 = "ab";
C.str.b1 = "b";
C.str.c1.a2.a3 = "caa";
C.str.c1.a2.b3 = "cab";
C.str.c1.b2 = "cb";
C.str.a1.c2 = "ac";
C.str.c1.a2.c3.a4 = ["caca11", "caca12"];

disp("====== str:")
C.str
disp("====== data:")
C.data
disp("====== vec:")
C.vec
disp("====== map:")
C.map
disp("====== node_tree:")
C.node_tree


%%

C.build_structure("shallow-values-first")
C.vec

C.build_structure("bredth-to-first","bulk",1);
C.vec

C.build_structure("bredth-to-first","bulk",0);
C.vec


%%

C.vec(4) = "new value!";

C.str.ne.b2.nonexistant = "another new value!!";

C.vec

C.str.a1.a2 = "setting the whole thing at once";

%%

C.str.a1.a2(1,2) = "indexed!";
C.vec

%%

C.build_structure("shallow-values-first","column");
C.vec


%%

C_copy = C.copy;
C_copy.vec

C_zeros = C.zeros;
C_zeros.str

C_ones = C.ones;
C_ones.vec

C_subcopy = C.str.a1();
C_subcopy.vec


%%

C = structor("column");
C.str.X = [ 1  2  3  4;
           10 20 30 40];
C.str.U = [200 400 600 800];

disp('Non-interpolated:')
C.vec

C_interp = C.interp([1:4],[1:4]+0.5);

disp('interpolated:')
C_interp.vec