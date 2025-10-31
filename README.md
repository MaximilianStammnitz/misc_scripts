"Deep mutational scanning data analysis using R" by P. ChaGPicasso
──────────────────────────────────────────────────────────────────
                 
                 .-^-.
     /\/\      /( X )\        _ _._           ___
    (o  o)----<  \_/  >-----/ /|/|\ \------ / _  \   ΔΔG
     \__/  . .  /| |\  . .  /_/ |_| \_\    | / \ |  fitness
      ||  /_\  /_|_|_\  /_\    . . .       | |_| |   →
   ___||_/___\_|__DMS__|/___\____|____      \___/  mut effect
  /__/_/____/_\_|____|_/____\_/_\___\_\

         .-.
   .----(   )----.     volcano plot
  /  .---`-'---.  \      ^
 |  /  p  value  \ |     |   *
 | |   *      *   | |    |      *
 |  \     *     /  |     +----------→ effect size
  \  '.___*__.''  /
   '--------------'

  [ mutation matrix ]
   AA  AC  AG  AT  |   ..fractured grid..
   CA  CC  CG  CT  |   |\  |\  |\  |\
   GA  GC  GG  GT  |   | \ | \ | \ | \
   TA  TC  TG  TT  |   |__\|__\|__\|__\

      /\      split transcript
   __/  \__   and broken reads →  ||||||||||
  /  /\    \       ||||||||||     |||| ||||
  \_/  \____/       |||   |||     |||   ||

   faces in fragments:   /\     /\      /\ 
                        (==)   (<> )   (..)
                         \/     \/      \/
                    antibodies ^  variants ^
                   ────────────┴───────────


    fractured dose–response curves (log10 dose on x, response on y)

        response ↑
                 |                         *
             1.0 |                      ***      Hill ~ 1.2
                 |                   ***   *
                 |                ***       *
                 |             ***           *            asymptote
                 |          ***               *
             0.5 |  ***  ***                  *  ← EC50 (midpoint)
                 | **  **                      *
                 |*    *                        *
                 +--------------------------------------------→ log10 dose
                  -3     -2     -1      0      1      2

        replicate curves / variants (IC50 shifts & slope changes)

        response ↑
             1.0 |     ──*───***────*──────────*────────  variant A (right-shifted, IC50 ↑)
                 |   ─*─*  *     ***   *───*────          variant B (shallower slope)
                 | *─*    *         ***    *             variant C (partial efficacy)
             0.5 |*       *            **    *
                 +------------------------------------------------→ log10 dose
                    EC50_A →|      EC50_B →|   EC50_C →|

        quick legend:
          * data points            ─ stitched “cubist” fit
          EC50/IC50 = mid response / half-max inhibition
          Hill slope = steepness; right shift = lower potency
