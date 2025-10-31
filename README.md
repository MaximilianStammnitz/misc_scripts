                .-^-.
     /\/\      /( X )\      _ _._        ___
    (o  o)----<  \_/  >---/ /|/|\ \---- / _ \   ΔΔG
     \__/  . .  /| |\  . ./_/ |_| \_\  | / \ | fitness→
      ||  /_\  /_|_|_\ /_\   . . .     | |_| | mut eff.
   ___||_/___\_|__DMS__|/___\__|___     \___/
  /__/_/____/_\_|____|_/____\_/_\__\

        .-.
   .---(   )---.   volcano plot
  / .--`-'--.  \      ^     *
 | /  p value \ |     |   *
 | |  *    *  | |     |     *       effect size→
 |  \   *    /  |     +-----------→
  \  '.___.''  /
   '----------'

  [ mutation matrix ]
   AA AC AG AT | fractured grid
   CA CC CG CT | |\ |\ |\ |\
   GA GC GG GT | | \| \| \| \
   TA TC TG TT | |__|__|__|__

     /\   split transcript
  __/  \__  broken reads →  ||||||||||  ||| |||

  faces in fragments:   /\   /\   /\ 
                       (==) (<> ) (..)
                        \/   \/   \/
                 antibodies^ variants^
                ────────────┴────────

 dose–response (log10 dose on x, response on y)

   resp ↑
        |                     *
     1.0|                  ***     Hill~1.2
        |               ***  *
        |            ***      *    asymptote
        |         ***          *
     0.5| ** ** **             * ← EC50
        |*   *                   *
        +-----------------------------→ log10 dose
          -3   -2   -1    0    1

 variants & IC50 shifts

   resp ↑
     1.0|  ─*─***─*──────*──────  A (right shift)
        | ─* *   ***  *──*───     B (shallow)
     0.5|*      *  **   *         C (partial)
        +-----------------------------→ log10 dose
           EC50_A  EC50_B  EC50_C

  quick legend:
    * data points     ─ stitched fit
    EC50/IC50 mid-response; right shift = lower potency
