#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"

//-----------------------------------------------------------------------------
void ProgressBar::loadbar(unsigned int x, unsigned int n, unsigned int w) {
   if ( (x != n) && (x % (n/100) != 0) ) return;

   float ratio  =  x/(float)n;
   int   c      =  ratio * w;

   cout << setw(3) << (int)(ratio*100) << "% [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "]\r" << flush;
}//loadbar

//-----------------------------------------------------------------------------
void ProgressBar::loadbar2(unsigned int x, unsigned int n, unsigned int w) {
   if ( (x != n) && (x % (n/100) != 0) ) return;
 
   float ratio  =  x/(float)n;
   int   c      =  ratio * w;
 
   cout << setw(3) << (int)(ratio*100) << "% [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "] (" << x << "/" << n << ")\r" << flush;
}//loadbar2

//-----------------------------------------------------------------------------
void ProgressBar::loadbar3(string prefix, unsigned int x, unsigned int n, unsigned int w) {
   if ( (x != n) && (x % (n/100) != 0) ) return;
 
   float ratio  =  x/(float)n;
   int   c      =  ratio * w;
 
   cout << prefix << " [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "] " << setw(3) << (int)(ratio*100) << "% (" << x << "/" << n << ")\r" << flush;
}//loadbar3

//-----------------------------------------------------------------------------
// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
void ProgressBar::loadBar(int x, int n, int r, int w) {
   // Only update r times.
   if ( x % (n/r) != 0 ) return;
 
   // Calculuate the ratio of complete-to-incomplete.
   float ratio = x/(float)n;
   int   c     = ratio * w;
 
   // Show the percentage complete.
   printf("%3d%% [", (int)(ratio*100) );
 
   // Show the load bar.
   for (int x=0; x<c; x++)
      printf("=");
 
   for (int x=c; x<w; x++)
      printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    //printf("]\n\033[F\033[J");
    printf("]\n33[1A33[2K");
}

//ClassImp(ProgressBar)
