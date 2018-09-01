\ Model matlab
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2.5e-14 C0 + 2.5e-14 C1 + 2.5e-14 C2 + 2.5e-14 C3 + 2.5e-14 C4
   + 2.5e-14 C5 + [ 1.34449e-08 C0 * C1 + 2.17708e-08 C0 * C2
   + 1.54952e-09 C1 * C2 + 1.34449e-08 C3 * C4 + 2.17708e-08 C3 * C5
   + 1.54952e-09 C4 * C5 ] / 2 
Subject To
 R0: C0 + C3 = 1
 R1: C1 + C4 = 1
 R2: C2 + C5 = 1
Bounds
End
