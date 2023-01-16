
// Display warning for floating point exception
ieee(1)

clear;
//load oono
//Y=newsp108(1:3:108,:);
//Y=Y+1;
mtlb_load("Y");

// ! L.7: mtlb(Y) can be replaced by Y() or Y whether Y is an M-file or not
[m,n] = size(mtlb_double(mtlb(Y)));

C = rand(m,2);

for t = 1:10;
  // !! L.12: Unknown function pinv, original calling sequence used
  // ! L.12: mtlb(Y) can be replaced by Y() or Y whether Y is an M-file or not
  A = mtlb_double(pinv(C))*mtlb_double(mtlb(Y));
  // Y=CA A=pinv(C)*Y
  A = mtlb_max(0,A);
  // ! L.14: mtlb(Y) can be replaced by Y() or Y whether Y is an M-file or not
  // !! L.14: Unknown function pinv, original calling sequence used
  C = mtlb_double(mtlb(Y))*mtlb_double(pinv(A));
  C = mtlb_max(0,C);
  // ! L.16: mtlb(Y) can be replaced by Y() or Y whether Y is an M-file or not
  E(1,t) = norm(mtlb_s(mtlb_double(mtlb(Y)),C*A));
end;
