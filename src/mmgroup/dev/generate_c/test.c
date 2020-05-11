/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%GEN h
// // comment for header
// %%GEN c

// %%EXPORT_TABLE 
uint32_t my_table[3] = {
  // %%TABLE y, uint{TABLE_SIZE}
0x00000001UL,0x00000002UL,0x00000003UL
};

// %%EXPORT
void int f1(int a, int b)
{
    // %%AddTo a,  {mult}*3 
    a = a + 51;
    a +=  my_table[b % 3] ;
    // %%FOR i in range(3)
       // %%IF {i} % 2
       // %%ELSE 
          // Add zero to variable a
       // %%END IF 
       a += 0;

       // %%FOR j in range(4)
          a += 0;
          a += 0;
          a += 0;
          a += 0;
       // %%END FOR # loop j
       // %%IF {i} % 2
          // Add an odd number to variable a
       // %%END IF 
       a += 1;

       // %%FOR j in range(4)
          a += 0;
          a += 1;
          a += 2;
          a += 3;
       // %%END FOR # loop j
       // %%IF {i} % 2
       // %%ELSE IF {i} > 0
          // Add an even number to variable a
       // %%END IF 
       a += 2;

       // %%FOR j in range(4)
          a += 0;
          a += 2;
          a += 4;
          a += 6;
       // %%END FOR # loop j
    // %%END FOR  # loop i
    // %%FOR t0, t1 in pair_list 
        a += my_table[0] * 7;
        a += my_table[1] * 4;
        a += my_table[2] * 9;
    // %%END FOR
    // %%FOR k in range(1,3)
        a += 
            0 * y_table[0] + 
            2 * y_table[1] + 
            4 * y_table[2];
        a += 
            0 * y_table[0] + 
            4 * y_table[1] + 
            8 * y_table[2];
    // %%END FOR
    return a;
}


   